/*
 *  GPU depth renderer: the MPR pipeline as three GL 4.3 compute
 *  passes (doc/TAPE-DESIGN.md Rounds 5-6).  The shaders are the
 *  proven ones from tests/gpu.cpp, generated per deck with the blob
 *  offsets, slot count, and C opcode values baked into the source.
 *
 *  libEGL is dlopen'd at first use, so SbFab carries no GL link
 *  dependency; without a usable stack every entry point returns
 *  false and callers keep the CPU path.
 *
 *  Correctness contract (referee'd in tests/gpu.cpp): pixels are a
 *  pure function of field signs, device min/max choices follow the
 *  CPU's exact conditions (operand-taint veto included), interval
 *  ops are either `precise` mirrors of math_i or conservative
 *  widenings - sound choices leave shortened-tape values exact, so
 *  the image is bit-identical to the CPU renderer's.
 */

#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <vector>

#ifdef STIBIUM_HAS_EGL

#include <dlfcn.h>
#include <EGL/egl.h>
#include <GL/glcorearb.h>

#endif

#include "fab/tree/gpu.h"
#include "fab/tree/render.h"
#include "fab/tree/node/opcodes.h"

#ifdef STIBIUM_HAS_EGL

namespace {

constexpr uint32_t GPU_TILE = 16;
constexpr uint32_t GPU_SLAB = 16;   // z voxels per classified slab

////////////////////////////////////////////////////////////////////////////
//  EGL via dlopen, GL via eglGetProcAddress

#define EGL_FUNCS(X) \
    X(PFNEGLGETDISPLAYPROC,      eglGetDisplay) \
    X(PFNEGLINITIALIZEPROC,      eglInitialize) \
    X(PFNEGLBINDAPIPROC,         eglBindAPI) \
    X(PFNEGLCHOOSECONFIGPROC,    eglChooseConfig) \
    X(PFNEGLCREATECONTEXTPROC,   eglCreateContext) \
    X(PFNEGLMAKECURRENTPROC,     eglMakeCurrent) \
    X(PFNEGLGETPROCADDRESSPROC,  eglGetProcAddress)

typedef EGLDisplay (*PFNEGLGETDISPLAYPROC)(EGLNativeDisplayType);
typedef EGLBoolean (*PFNEGLINITIALIZEPROC)(EGLDisplay, EGLint*, EGLint*);
typedef EGLBoolean (*PFNEGLBINDAPIPROC)(EGLenum);
typedef EGLBoolean (*PFNEGLCHOOSECONFIGPROC)(EGLDisplay, const EGLint*,
                                             EGLConfig*, EGLint, EGLint*);
typedef EGLContext (*PFNEGLCREATECONTEXTPROC)(EGLDisplay, EGLConfig,
                                              EGLContext, const EGLint*);
typedef EGLBoolean (*PFNEGLMAKECURRENTPROC)(EGLDisplay, EGLSurface,
                                            EGLSurface, EGLContext);
typedef void (*(*PFNEGLGETPROCADDRESSPROC)(const char*))(void);

#define DECLARE(type, name) type p_##name = nullptr;
EGL_FUNCS(DECLARE)
#undef DECLARE

#define GL_FUNCS(X) \
    X(PFNGLGETSTRINGPROC,             glGetString) \
    X(PFNGLCREATESHADERPROC,          glCreateShader) \
    X(PFNGLSHADERSOURCEPROC,          glShaderSource) \
    X(PFNGLCOMPILESHADERPROC,         glCompileShader) \
    X(PFNGLGETSHADERIVPROC,           glGetShaderiv) \
    X(PFNGLGETSHADERINFOLOGPROC,      glGetShaderInfoLog) \
    X(PFNGLCREATEPROGRAMPROC,         glCreateProgram) \
    X(PFNGLATTACHSHADERPROC,          glAttachShader) \
    X(PFNGLLINKPROGRAMPROC,           glLinkProgram) \
    X(PFNGLGETPROGRAMIVPROC,          glGetProgramiv) \
    X(PFNGLUSEPROGRAMPROC,            glUseProgram) \
    X(PFNGLDELETESHADERPROC,          glDeleteShader) \
    X(PFNGLDELETEPROGRAMPROC,         glDeleteProgram) \
    X(PFNGLGENBUFFERSPROC,            glGenBuffers) \
    X(PFNGLDELETEBUFFERSPROC,         glDeleteBuffers) \
    X(PFNGLBINDBUFFERPROC,            glBindBuffer) \
    X(PFNGLBUFFERDATAPROC,            glBufferData) \
    X(PFNGLBINDBUFFERBASEPROC,        glBindBufferBase) \
    X(PFNGLDISPATCHCOMPUTEPROC,       glDispatchCompute) \
    X(PFNGLMEMORYBARRIERPROC,         glMemoryBarrier) \
    X(PFNGLGETBUFFERSUBDATAPROC,      glGetBufferSubData) \
    X(PFNGLFINISHPROC,                glFinish) \
    X(PFNGLGETERRORPROC,              glGetError) \
    X(PFNGLGETUNIFORMLOCATIONPROC,    glGetUniformLocation) \
    X(PFNGLUNIFORM1UIPROC,            glUniform1ui)

#define DECLARE(type, name) type name = nullptr;
GL_FUNCS(DECLARE)
#undef DECLARE

struct Gpu
{
    EGLDisplay dpy = EGL_NO_DISPLAY;
    EGLContext ctx = EGL_NO_CONTEXT;
    bool ok = false;
    std::string renderer;

    Gpu()
    {
        void* egl = dlopen("libEGL.so.1", RTLD_NOW | RTLD_LOCAL);
        if (!egl)
            egl = dlopen("libEGL.so", RTLD_NOW | RTLD_LOCAL);
        if (!egl)
            return;
#define LOAD(type, name) \
        p_##name = reinterpret_cast<type>(dlsym(egl, #name)); \
        if (!p_##name) return;
        EGL_FUNCS(LOAD)
#undef LOAD

        dpy = p_eglGetDisplay(EGL_DEFAULT_DISPLAY);
        if (dpy == EGL_NO_DISPLAY)
            return;
        EGLint maj = 0, min = 0;
        if (!p_eglInitialize(dpy, &maj, &min))
            return;
        if (!p_eglBindAPI(EGL_OPENGL_API))
            return;
        const EGLint cfg_attr[] = {
            EGL_SURFACE_TYPE,    EGL_PBUFFER_BIT,
            EGL_RENDERABLE_TYPE, EGL_OPENGL_BIT,
            EGL_NONE
        };
        EGLConfig cfg;
        EGLint n = 0;
        if (!p_eglChooseConfig(dpy, cfg_attr, &cfg, 1, &n) || n < 1)
            return;
        const EGLint ctx_attr[] = {
            EGL_CONTEXT_MAJOR_VERSION, 4,
            EGL_CONTEXT_MINOR_VERSION, 3,
            EGL_NONE
        };
        ctx = p_eglCreateContext(dpy, cfg, EGL_NO_CONTEXT, ctx_attr);
        if (ctx == EGL_NO_CONTEXT)
            return;
        if (!p_eglMakeCurrent(dpy, EGL_NO_SURFACE, EGL_NO_SURFACE, ctx))
            return;

#define LOADGL(type, name) \
        name = reinterpret_cast<type>(p_eglGetProcAddress(#name)); \
        if (!name) return;
        GL_FUNCS(LOADGL)
#undef LOADGL

        const GLubyte* r = glGetString(GL_RENDERER);
        renderer = r ? reinterpret_cast<const char*>(r) : "(unknown)";
        ok = true;
        /*  Unbind so any thread can bind on demand  */
        p_eglMakeCurrent(dpy, EGL_NO_SURFACE, EGL_NO_SURFACE,
                         EGL_NO_CONTEXT);
    }
};

std::mutex g_lock;

Gpu& gpu()
{
    static Gpu g;
    return g;
}

////////////////////////////////////////////////////////////////////////////
//  Register rescheduling (TAPE-DESIGN Round 8).
//
//  The deck's linear-scan allocation follows graph emission order,
//  which interleaves subtrees - the merged Zeiss peaks at ~870 live
//  slots, and a per-thread array that size collapses GPU occupancy.
//  Fidget bounds its register file with LRU spilling; we hold a card
//  fidget doesn't: the export step may REORDER clauses freely (pure
//  dataflow, no side effects - values are bit-identical in any
//  topological order).  Sethi-Ullman scheduling (heavier subtree
//  first, DAG nodes emitted once at first need) drops a CSG tree's
//  true register requirement to roughly the log of its width, after
//  which a fresh linear scan packs it tight.  Constants stop being
//  registers entirely - they're thread-uniform and live in the
//  generated shader source.
//
//  Input and output are both format-v1 blobs, renumbered so that
//  slots are [0,n_const) constants, then axes, then registers.
//  Returns an empty vector if anything looks unfamiliar (caller
//  falls back to the untransformed blob).

std::vector<uint32_t> blob_reschedule(const std::vector<uint32_t>& in)
{
    struct BI
    {
        uint32_t num_slots, root;
        uint32_t n_const, n_x, n_y, n_z, n_clauses;
        uint32_t const_off, ax_off, clause_off;
    };
    BI bi;
    bi.num_slots = in[1];
    bi.root      = in[2];
    bi.n_const   = in[3];
    bi.n_x       = in[4];
    bi.n_y       = in[5];
    bi.n_z       = in[6];
    bi.n_clauses = in[7];
    bi.const_off = 8;
    bi.ax_off    = bi.const_off + 2 * bi.n_const;
    bi.clause_off = bi.ax_off + bi.n_x + bi.n_y + bi.n_z;

    const uint32_t n_axes = bi.n_x + bi.n_y + bi.n_z;
    const uint32_t n_leaves = bi.n_const + n_axes;

    /*  Value ids: [0,n_const) const leaves, [n_const,n_leaves) axis
     *  leaves, then one value per clause in original order.  */
    struct Val
    {
        uint32_t op = 0, a = 0, b = 0;   // operand VALUE ids
        uint32_t imm = 0;                // f32 bits (OP_CONST)
        uint32_t uses = 0;
        uint32_t need = 0;               // Sethi-Ullman weight
        int32_t reg = -1;
        bool emitted = false;
        uint8_t kind = 0;                // 0 leaf-const, 1 leaf-axis,
                                         // 2 clause
        uint8_t arity = 0;
    };
    std::vector<Val> vals(n_leaves + bi.n_clauses);
    for (uint32_t v = 0; v < bi.n_const; ++v)
        vals[v].kind = 0;
    for (uint32_t v = bi.n_const; v < n_leaves; ++v)
        vals[v].kind = 1;

    /*  cur[slot] = value currently living in that slot  */
    std::vector<uint32_t> cur(bi.num_slots, UINT32_MAX);
    for (uint32_t q = 0; q < bi.n_const; ++q)
        cur[in[bi.const_off + 2*q]] = q;
    for (uint32_t q = 0; q < n_axes; ++q)
        cur[in[bi.ax_off + q]] = bi.n_const + q;

    const auto op_arity = [](uint32_t op) -> uint8_t {
        switch (op)
        {
            case OP_ADD: case OP_SUB: case OP_MUL: case OP_DIV:
            case OP_MIN: case OP_MAX: case OP_POW: case OP_ATAN2:
            case OP_MOD:
                return 2;
            case OP_ABS: case OP_SQUARE: case OP_SQRT: case OP_SIN:
            case OP_COS: case OP_TAN: case OP_ASIN: case OP_ACOS:
            case OP_ATAN: case OP_NEG: case OP_EXP: case OP_FLOOR:
            case OP_LOG: case OP_COPY:
                return 1;
            case OP_CONST:
                return 0;
            default:
                return 255;   // unfamiliar - bail out
        }
    };

    for (uint32_t c = 0; c < bi.n_clauses; ++c)
    {
        const uint32_t* w = in.data() + bi.clause_off + 5*c;
        Val& v = vals[n_leaves + c];
        v.kind = 2;
        v.op = w[0];
        v.arity = op_arity(w[0]);
        if (v.arity == 255)
            return {};
        v.imm = w[4];
        if (v.arity >= 1)
        {
            if (w[2] >= bi.num_slots || cur[w[2]] == UINT32_MAX)
                return {};
            v.a = cur[w[2]];
        }
        if (v.arity >= 2)
        {
            if (w[3] >= bi.num_slots || cur[w[3]] == UINT32_MAX)
                return {};
            v.b = cur[w[3]];
        }
        cur[w[1]] = n_leaves + c;
    }
    if (bi.root >= bi.num_slots || cur[bi.root] == UINT32_MAX)
        return {};
    const uint32_t root_v = cur[bi.root];
    if (vals[root_v].kind != 2 && bi.n_clauses > 0)
        return {};

    /*  Reachability + use counts + Sethi-Ullman weights, bottom-up
     *  via an explicit post-order stack (trees can be thousands
     *  deep).  Shared DAG nodes are visited once.  */
    std::vector<uint32_t> stack = { root_v };
    std::vector<uint8_t> state(vals.size(), 0);   // 0 new 1 open 2 done
    vals[root_v].uses = 1;
    while (!stack.empty())
    {
        const uint32_t v = stack.back();
        Val& val = vals[v];
        if (state[v] == 2)
        {
            stack.pop_back();
            continue;
        }
        if (val.kind != 2)
        {
            state[v] = 2;
            stack.pop_back();
            continue;
        }
        if (state[v] == 0)
        {
            state[v] = 1;
            /*  Count uses on first reach of each edge  */
            if (val.arity >= 1)
            {
                vals[val.a].uses++;
                if (state[val.a] != 2)
                    stack.push_back(val.a);
            }
            if (val.arity >= 2)
            {
                vals[val.b].uses++;
                if (state[val.b] != 2)
                    stack.push_back(val.b);
            }
            continue;
        }
        /*  state 1: children done - compute weight  */
        const uint32_t na =
                val.arity >= 1 && vals[val.a].kind == 2
                        ? vals[val.a].need : 0;
        const uint32_t nb =
                val.arity >= 2 && vals[val.b].kind == 2
                        ? vals[val.b].need : 0;
        if (val.arity == 2)
            val.need = na == nb ? na + 1 : (na > nb ? na : nb);
        else
            val.need = na > 1 ? na : 1;
        state[v] = 2;
        stack.pop_back();
    }

    /*  Emission: DFS from the root, heavier child first, each value
     *  once.  Registers by linear scan: freed when the last use is
     *  consumed, reused LIFO.  */
    std::vector<uint32_t> order;
    order.reserve(bi.n_clauses);
    std::vector<uint32_t> remaining(vals.size());
    for (size_t v = 0; v < vals.size(); ++v)
        remaining[v] = vals[v].uses;

    std::vector<int32_t> free_regs;
    int32_t next_reg = 0;
    const auto alloc_reg = [&]() -> int32_t {
        if (!free_regs.empty())
        {
            const int32_t r = free_regs.back();
            free_regs.pop_back();
            return r;
        }
        return next_reg++;
    };
    const auto consume = [&](uint32_t v) {
        if (vals[v].kind != 2)
            return;
        if (--remaining[v] == 0 && vals[v].reg >= 0)
            free_regs.push_back(vals[v].reg);
    };

    stack = { root_v };
    std::fill(state.begin(), state.end(), 0);
    while (!stack.empty())
    {
        const uint32_t v = stack.back();
        Val& val = vals[v];
        if (val.kind != 2 || val.emitted)
        {
            stack.pop_back();
            continue;
        }
        if (state[v] == 0)
        {
            state[v] = 1;
            /*  Children, heavier first (pushed last = visited
             *  first).  Only unemitted clause children matter.  */
            uint32_t kids[2];
            int nk = 0;
            if (val.arity >= 1 && vals[val.a].kind == 2 &&
                !vals[val.a].emitted)
                kids[nk++] = val.a;
            if (val.arity >= 2 && vals[val.b].kind == 2 &&
                !vals[val.b].emitted)
                kids[nk++] = val.b;
            if (nk == 2 && vals[kids[0]].need > vals[kids[1]].need)
            {
                const uint32_t t = kids[0];
                kids[0] = kids[1];
                kids[1] = t;
            }
            for (int k = 0; k < nk; ++k)
                stack.push_back(kids[k]);
            continue;
        }
        /*  Children emitted: consume operand uses, then define  */
        if (val.arity >= 1)
            consume(val.a);
        if (val.arity >= 2)
            consume(val.b);
        val.reg = alloc_reg();
        val.emitted = true;
        order.push_back(v);
        stack.pop_back();
    }
    if (order.size() > bi.n_clauses)
        return {};

    /*  Rebuild the blob with the new numbering.  */
    const uint32_t R = uint32_t(next_reg);
    const auto slot_of = [&](uint32_t v) -> uint32_t {
        const Val& val = vals[v];
        if (val.kind == 0)
            return v;                          // const slot
        if (val.kind == 1)
            return v;                          // axis slot
        return n_leaves + uint32_t(val.reg);   // register
    };

    std::vector<uint32_t> out;
    out.reserve(in.size());
    out.push_back(1);
    out.push_back(n_leaves + R);
    out.push_back(slot_of(root_v));
    out.push_back(bi.n_const);
    out.push_back(bi.n_x);
    out.push_back(bi.n_y);
    out.push_back(bi.n_z);
    out.push_back(uint32_t(order.size()));
    for (uint32_t q = 0; q < bi.n_const; ++q)
    {
        out.push_back(q);                       // renumbered slot
        out.push_back(in[bi.const_off + 2*q + 1]);
    }
    for (uint32_t q = 0; q < n_axes; ++q)
        out.push_back(bi.n_const + q);
    for (const uint32_t v : order)
    {
        const Val& val = vals[v];
        out.push_back(val.op);
        out.push_back(n_leaves + uint32_t(val.reg));
        out.push_back(val.arity >= 1 ? slot_of(val.a) : 0);
        out.push_back(val.arity >= 2 ? slot_of(val.b) : 0);
        out.push_back(val.imm);
    }
    return out;
}

////////////////////////////////////////////////////////////////////////////
//  Shader generation (proven in tests/gpu.cpp; see that file and
//  TAPE-DESIGN Round 6 for the soundness discussion)

struct BlobInfo
{
    uint32_t num_slots, root;
    uint32_t n_const, n_x, n_y, n_z, n_clauses;
    uint32_t const_off, ax_off, clause_off;
};

BlobInfo parse_blob(const std::vector<uint32_t>& b)
{
    BlobInfo i;
    i.num_slots = b[1];
    i.root      = b[2];
    i.n_const   = b[3];
    i.n_x       = b[4];
    i.n_y       = b[5];
    i.n_z       = b[6];
    i.n_clauses = b[7];
    i.const_off = 8;
    i.ax_off    = i.const_off + 2 * i.n_const;
    i.clause_off = i.ax_off + i.n_x + i.n_y + i.n_z;
    return i;
}

std::string op_prelude()
{
    char buf[2048];
    snprintf(buf, sizeof(buf),
        "const uint OADD=%uu,OSUB=%uu,OMUL=%uu,ODIV=%uu,OMIN=%uu,"
        "OMAX=%uu,OPOW=%uu,OMOD=%uu,OABS=%uu,OSQUARE=%uu,OSQRT=%uu,"
        "OSIN=%uu,OCOS=%uu,OTAN=%uu,OASIN=%uu,OACOS=%uu,OATAN=%uu,"
        "OATAN2=%uu,ONEG=%uu,OEXP=%uu,OFLOOR=%uu,OLOG=%uu,"
        "OCONST=%uu,OCOPY=%uu;\n",
        OP_ADD, OP_SUB, OP_MUL, OP_DIV, OP_MIN, OP_MAX, OP_POW,
        OP_MOD, OP_ABS, OP_SQUARE, OP_SQRT, OP_SIN, OP_COS, OP_TAN,
        OP_ASIN, OP_ACOS, OP_ATAN, OP_ATAN2, OP_NEG, OP_EXP,
        OP_FLOOR, OP_LOG, OP_CONST, OP_COPY);
    return buf;
}

std::string blob_prelude(const BlobInfo& bi,
                         uint32_t ni, uint32_t nj, uint32_t nk)
{
    char buf[1024];
    snprintf(buf, sizeof(buf),
        "const uint NI=%uu,NJ=%uu,NK=%uu,ROOT=%uu,N_CONST=%uu,"
        "N_X=%uu,N_Y=%uu,N_Z=%uu,N_CLAUSES=%uu,CONST_OFF=%uu,"
        "AX_OFF=%uu,CLAUSE_OFF=%uu,NSLOTS=%uu;\n"
        "const uint TILE=%uu,TILES_X=%uu,TILES_Y=%uu;\n"
        "const uint SLAB=%uu,N_SLABS=%uu;\n",
        ni, nj, nk, bi.root, bi.n_const, bi.n_x, bi.n_y, bi.n_z,
        bi.n_clauses, bi.const_off, bi.ax_off, bi.clause_off,
        bi.num_slots, GPU_TILE,
        (ni + GPU_TILE - 1) / GPU_TILE, (nj + GPU_TILE - 1) / GPU_TILE,
        GPU_SLAB, (nk + GPU_SLAB - 1) / GPU_SLAB);
    return buf;
}

/*  Thread-uniform constants live in the shader source (constant
 *  memory), not in per-thread arrays: state arrays cover only
 *  [axes | registers], indexed by slot - N_CONST.  */
std::string pinned_table(const std::vector<uint32_t>& blob,
                         const BlobInfo& bi)
{
    std::string src = "const uint NSTATE = ";
    char buf[64];
    snprintf(buf, sizeof(buf), "%uu;\n",
             bi.num_slots - bi.n_const);
    src += buf;
    if (bi.n_const == 0)
    {
        src += "const float PINNED[1] = float[1](0.0);\n";
        return src;
    }
    snprintf(buf, sizeof(buf), "const float PINNED[%u] = float[%u](\n",
             bi.n_const, bi.n_const);
    src += buf;
    for (uint32_t q = 0; q < bi.n_const; ++q)
    {
        snprintf(buf, sizeof(buf), "uintBitsToFloat(%uu)%s",
                 blob[bi.const_off + 2*q + 1],
                 q + 1 < bi.n_const ? (q % 4 == 3 ? ",\n" : ",")
                                    : ");\n");
        src += buf;
    }
    return src;
}

std::string make_p1(const std::vector<uint32_t>& blob,
                    const BlobInfo& bi,
                    uint32_t ni, uint32_t nj, uint32_t nk)
{
    std::string src = "#version 430\n"
        "layout(local_size_x = 64) in;\n"
        "layout(std430, binding=0) readonly buffer Blob { uint blob[]; };\n"
        "layout(std430, binding=1) readonly buffer BX { float Xs[]; };\n"
        "layout(std430, binding=2) readonly buffer BY { float Ys[]; };\n"
        "layout(std430, binding=3) readonly buffer BZ { float Zs[]; };\n"
        "layout(std430, binding=6) writeonly buffer BC { uint choices[]; };\n"
        "layout(std430, binding=7) buffer BT { uvec2 tileinfo[]; };\n"
        "uniform uint T0;\n"
        "uniform uint NT;\n";
    src += op_prelude();
    src += blob_prelude(bi, ni, nj, nk);
    src += pinned_table(blob, bi);
    src += R"GLSL(
const float INF = uintBitsToFloat(0x7F800000u);

precise vec2 iv[NSTATE];
bool tnt[NSTATE];

vec2 ldi(uint s)
{
    return s < N_CONST ? vec2(PINNED[s]) : iv[s - N_CONST];
}
bool ldt(uint s)
{
    return s < N_CONST ? false : tnt[s - N_CONST];
}

vec2 iv_mul(vec2 A, vec2 B)
{
    precise float c1 = A.x * B.x, c2 = A.x * B.y,
                  c3 = A.y * B.x, c4 = A.y * B.y;
    float lo = min(min(c1, c2), min(c3, c4));
    float hi = max(max(c1, c2), max(c3, c4));
    if (isnan(c1) || isnan(c2) || isnan(c3) || isnan(c4))
        { lo = -INF; hi = INF; }
    return vec2(lo, hi);
}

void main()
{
    uint lt = gl_GlobalInvocationID.x;
    if (lt >= NT)
        return;
    uint t = T0 + lt;
    uint ti = t - (t / TILES_X) * TILES_X, tj = t / TILES_X;

    vec2 X = vec2(Xs[ti*TILE], Xs[min(ti*TILE + TILE, NI)]);
    vec2 Y = vec2(Ys[tj*TILE], Ys[min(tj*TILE + TILE, NJ)]);
    vec2 Z = vec2(Zs[0], Zs[NK]);
    for (uint q = 0u; q < N_X; ++q)
        { uint s = blob[AX_OFF+q] - N_CONST; iv[s] = X; tnt[s] = false; }
    for (uint q = 0u; q < N_Y; ++q)
        { uint s = blob[AX_OFF+N_X+q] - N_CONST; iv[s] = Y; tnt[s] = false; }
    for (uint q = 0u; q < N_Z; ++q)
        { uint s = blob[AX_OFF+N_X+N_Y+q] - N_CONST; iv[s] = Z; tnt[s] = false; }

    for (uint c = 0u; c < N_CLAUSES; ++c)
    {
        uint base = CLAUSE_OFF + 5u*c;
        uint op = blob[base], o = blob[base+1u];
        uint sa = blob[base+2u], sb = blob[base+3u];
        vec2 A = ldi(sa), B = ldi(sb);
        bool na = ldt(sa), nb = ldt(sb);
        precise vec2 R;
        bool u = na || nb;
        bool a_inf = isinf(A.x) || isinf(A.y);
        bool b_inf = isinf(B.x) || isinf(B.y);
        uint choice = 0u;

        if (op == OADD)      { R = vec2(A.x+B.x, A.y+B.y); u = u || a_inf || b_inf; }
        else if (op == OSUB) { R = vec2(A.x-B.y, A.y-B.x); u = u || a_inf || b_inf; }
        else if (op == OMUL) { R = iv_mul(A, B); u = u || a_inf || b_inf; }
        else if (op == ODIV) {
            if (B.x <= 0.0 && B.y >= 0.0) { R = vec2(-INF, INF); u = true; }
            else {
                precise vec2 rB = vec2(1.0/B.x, 1.0/B.y);
                R = iv_mul(A, rB);
                u = u || a_inf || b_inf;
            }
        }
        else if (op == OMIN) {
            R = vec2(A.x < B.x ? A.x : B.x, A.y < B.y ? A.y : B.y);
            if (!na && !nb) {
                if (A.y <= B.x)      choice = 1u;
                else if (B.y <= A.x) choice = 2u;
            }
        }
        else if (op == OMAX) {
            R = vec2(A.x > B.x ? A.x : B.x, A.y > B.y ? A.y : B.y);
            if (!na && !nb) {
                if (A.x >= B.y)      choice = 1u;
                else if (B.x >= A.y) choice = 2u;
            }
        }
        else if (op == OABS) {
            float lo = (A.x < 0.0) ? 0.0 : min(abs(A.x), abs(A.y));
            R = vec2(lo, max(abs(A.x), abs(A.y)));
        }
        else if (op == OSQUARE) {
            precise float l2 = A.x*A.x, u2 = A.y*A.y;
            float lo = (A.y > 0.0 && A.x < 0.0) ? 0.0 : min(u2, l2);
            R = vec2(lo, max(u2, l2));
        }
        else if (op == OSQRT) {
            R = vec2(A.x <= 0.0 ? 0.0 : sqrt(A.x),
                     A.y <= 0.0 ? 0.0 : sqrt(A.y));
            u = u || (A.x < 0.0);
        }
        else if (op == ONEG)   { R = vec2(-A.y, -A.x); }
        else if (op == OSIN || op == OCOS) { R = vec2(-1.0, 1.0); u = u || a_inf; }
        else if (op == OEXP)   { R = vec2(0.0, INF); }
        else if (op == OASIN || op == OATAN) { R = vec2(-1.5707964, 1.5707964); u = u || (op == OASIN && (A.x < -1.0 || A.y > 1.0)); }
        else if (op == OACOS)  { R = vec2(0.0, 3.1415927); u = u || A.x < -1.0 || A.y > 1.0; }
        else if (op == OATAN2) { R = vec2(-3.1415927, 3.1415927); u = true; }
        else if (op == OFLOOR) { R = vec2(floor(A.x), floor(A.y)); }
        else if (op == OCONST) { R = vec2(uintBitsToFloat(blob[base+4u])); u = false; }
        else if (op == OCOPY)  { R = A; u = na; }
        else                   { R = vec2(-INF, INF); u = true; }

        u = u || isnan(R.x) || isnan(R.y);
        iv[o - N_CONST] = R;
        tnt[o - N_CONST] = u;
        choices[lt * N_CLAUSES + c] = choice;
    }

    vec2 root = ldi(ROOT);
    bool rt = ldt(ROOT);
    uint cls = 2u;
    if (!rt && root.x > 0.0)  cls = 0u;
    else if (!rt && root.y < 0.0) cls = 1u;
    tileinfo[lt] = uvec2(cls, N_CLAUSES);
}
)GLSL";
    return src;
}

std::string make_p2(const std::vector<uint32_t>& blob,
                    const BlobInfo& bi,
                    uint32_t ni, uint32_t nj, uint32_t nk)
{
    std::string src = "#version 430\n"
        "layout(local_size_x = 64) in;\n"
        "layout(std430, binding=0) readonly buffer Blob { uint blob[]; };\n"
        "layout(std430, binding=6) readonly buffer BC { uint choices[]; };\n"
        "layout(std430, binding=7) buffer BT { uvec2 tileinfo[]; };\n"
        "layout(std430, binding=8) writeonly buffer BS { uint tapes[]; };\n"
        "uniform uint NT;\n";
    src += op_prelude();
    src += blob_prelude(bi, ni, nj, nk);
    src += pinned_table(blob, bi);
    src += R"GLSL(
bool live[NSTATE];
uint verdict[N_CLAUSES];

void main()
{
    uint t = gl_GlobalInvocationID.x;   // band-local
    if (t >= NT)
        return;
    if (tileinfo[t].x != 2u)
        return;

    for (uint s = 0u; s < NSTATE; ++s)
        live[s] = false;
    if (ROOT >= N_CONST)
        live[ROOT - N_CONST] = true;

    for (uint k = N_CLAUSES; k-- > 0u; )
    {
        uint base = CLAUSE_OFF + 5u*k;
        uint op = blob[base], o = blob[base+1u];
        uint sa = blob[base+2u], sb = blob[base+3u];
        if (!live[o - N_CONST])
        {
            verdict[k] = 0u;
            continue;
        }
        live[o - N_CONST] = false;
        uint ch = choices[t * N_CLAUSES + k];
        if ((op == OMIN || op == OMAX) && ch != 0u)
        {
            uint keep = (ch == 1u) ? sa : sb;
            verdict[k] = (keep == o) ? 0u : (ch == 1u ? 2u : 3u);
            if (keep >= N_CONST)
                live[keep - N_CONST] = true;
            continue;
        }
        verdict[k] = 1u;
        bool unary = (op == OABS || op == OSQUARE || op == OSQRT ||
                      op == OSIN || op == OCOS || op == OTAN ||
                      op == OASIN || op == OACOS || op == OATAN ||
                      op == ONEG || op == OEXP || op == OFLOOR ||
                      op == OLOG || op == OCOPY);
        bool none  = (op == OCONST);
        if (!none)
        {
            if (sa >= N_CONST)
                live[sa - N_CONST] = true;
            if (!unary && sb >= N_CONST)
                live[sb - N_CONST] = true;
        }
    }

    uint w = t * N_CLAUSES * 5u;
    uint len = 0u;
    for (uint k = 0u; k < N_CLAUSES; ++k)
    {
        if (verdict[k] == 0u)
            continue;
        uint base = CLAUSE_OFF + 5u*k;
        /*  Locals-then-store: interleaved SSBO reads and writes here
         *  miscompiled on Mesa Intel (the copy's out-slot word read
         *  back 0); same source ran correctly on NVIDIA.  */
        uint w0 = blob[base];
        uint w1 = blob[base+1u];
        uint w2 = blob[base+2u];
        uint w3 = blob[base+3u];
        uint w4 = blob[base+4u];
        if (verdict[k] != 1u)
        {
            w2 = (verdict[k] == 2u) ? w2 : w3;
            w0 = OCOPY;
            w3 = 0u;
            w4 = 0u;
        }
        tapes[w+0u] = w0;
        tapes[w+1u] = w1;
        tapes[w+2u] = w2;
        tapes[w+3u] = w3;
        tapes[w+4u] = w4;
        w += 5u;
        ++len;
    }
    tileinfo[t].y = len;
}
)GLSL";
    return src;
}

/*  Pass 1b (Round 8): the second subdivision level.  One thread per
 *  (tile, z-slab): interval-evaluate the TILE's tape over the slab
 *  box, recording min/max choices (CPU rules, taint-vetoed) in
 *  thread-local 2-bit packs; classify EMPTY/FULL/AMBIGUOUS; and for
 *  ambiguous boxes immediately replay the tape_push backward walk
 *  and emit a compacted slab tape into a shared pool via an atomic
 *  bump allocator (fidget's scheme - fixed per-box slots would need
 *  gigabytes).  Pool overflow sets a flag; the host falls back to
 *  the CPU renderer.  The interval switch mirrors pass 1 - keep
 *  them in sync.  */
std::string make_p1b(const std::vector<uint32_t>& blob,
                     const BlobInfo& bi,
                     uint32_t ni, uint32_t nj, uint32_t nk)
{
    std::string src = "#version 430\n"
        "layout(local_size_x = 64) in;\n"
        "layout(std430, binding=0) readonly buffer Blob { uint blob[]; };\n"
        "layout(std430, binding=1) readonly buffer BX { float Xs[]; };\n"
        "layout(std430, binding=2) readonly buffer BY { float Ys[]; };\n"
        "layout(std430, binding=3) readonly buffer BZ { float Zs[]; };\n"
        "layout(std430, binding=7) readonly buffer BT { uvec2 tileinfo[]; };\n"
        "layout(std430, binding=8) readonly buffer BS { uint tapes[]; };\n"
        "layout(std430, binding=9) writeonly buffer BSl { uvec4 slabinfo[]; };\n"
        "layout(std430, binding=10) writeonly buffer BP { uint pool[]; };\n"
        "layout(std430, binding=11) buffer BA { uint pool_next; uint pool_cap; uint overflow; };\n"
        "uniform uint T0;\n"
        "uniform uint NT;\n";
    src += op_prelude();
    src += blob_prelude(bi, ni, nj, nk);
    src += pinned_table(blob, bi);
    src += R"GLSL(
const float INF = uintBitsToFloat(0x7F800000u);

precise vec2 iv[NSTATE];
bool tnt[NSTATE];

/*  2-bit lanes: choices during the forward pass, verdicts during
 *  the backward pass (same clause never needs both at once).  */
uint lanes[(N_CLAUSES + 15u) / 16u];

uint get2(uint k) { return (lanes[k >> 4u] >> ((k & 15u) * 2u)) & 3u; }
void set2(uint k, uint v)
{
    uint w = k >> 4u, sh = (k & 15u) * 2u;
    lanes[w] = (lanes[w] & ~(3u << sh)) | (v << sh);
}

bool live[NSTATE];

vec2 ldi(uint s)
{
    return s < N_CONST ? vec2(PINNED[s]) : iv[s - N_CONST];
}
bool ldt(uint s)
{
    return s < N_CONST ? false : tnt[s - N_CONST];
}

vec2 iv_mul(vec2 A, vec2 B)
{
    precise float c1 = A.x * B.x, c2 = A.x * B.y,
                  c3 = A.y * B.x, c4 = A.y * B.y;
    float lo = min(min(c1, c2), min(c3, c4));
    float hi = max(max(c1, c2), max(c3, c4));
    if (isnan(c1) || isnan(c2) || isnan(c3) || isnan(c4))
        { lo = -INF; hi = INF; }
    return vec2(lo, hi);
}

void main()
{
    uint id = gl_GlobalInvocationID.x;   // band-local pair
    if (id >= NT * N_SLABS)
        return;
    uint lt = id / N_SLABS;
    uint sl = id - lt * N_SLABS;
    uint t = T0 + lt;

    uint tcls = tileinfo[lt].x;
    if (tcls != 2u)
    {
        slabinfo[id] = uvec4(tcls, 0u, 0u, 0u);
        return;
    }

    uint ti = t - (t / TILES_X) * TILES_X, tj = t / TILES_X;
    vec2 X = vec2(Xs[ti*TILE], Xs[min(ti*TILE + TILE, NI)]);
    vec2 Y = vec2(Ys[tj*TILE], Ys[min(tj*TILE + TILE, NJ)]);
    vec2 Z = vec2(Zs[sl*SLAB], Zs[min(sl*SLAB + SLAB, NK)]);
    for (uint q = 0u; q < N_X; ++q)
        { uint w = blob[AX_OFF+q] - N_CONST; iv[w] = X; tnt[w] = false; }
    for (uint q = 0u; q < N_Y; ++q)
        { uint w = blob[AX_OFF+N_X+q] - N_CONST; iv[w] = Y; tnt[w] = false; }
    for (uint q = 0u; q < N_Z; ++q)
        { uint w = blob[AX_OFF+N_X+N_Y+q] - N_CONST; iv[w] = Z; tnt[w] = false; }

    uint tbase = lt * N_CLAUSES * 5u;
    uint len = tileinfo[lt].y;
    for (uint c = 0u; c < len; ++c)
    {
        uint base = tbase + 5u*c;
        uint op = tapes[base], o = tapes[base+1u];
        uint sa = tapes[base+2u], sb = tapes[base+3u];
        vec2 A = ldi(sa), B = ldi(sb);
        bool na = ldt(sa), nb = ldt(sb);
        precise vec2 R;
        bool u = na || nb;
        bool a_inf = isinf(A.x) || isinf(A.y);
        bool b_inf = isinf(B.x) || isinf(B.y);
        uint choice = 0u;

        if (op == OADD)      { R = vec2(A.x+B.x, A.y+B.y); u = u || a_inf || b_inf; }
        else if (op == OSUB) { R = vec2(A.x-B.y, A.y-B.x); u = u || a_inf || b_inf; }
        else if (op == OMUL) { R = iv_mul(A, B); u = u || a_inf || b_inf; }
        else if (op == ODIV) {
            if (B.x <= 0.0 && B.y >= 0.0) { R = vec2(-INF, INF); u = true; }
            else {
                precise vec2 rB = vec2(1.0/B.x, 1.0/B.y);
                R = iv_mul(A, rB);
                u = u || a_inf || b_inf;
            }
        }
        else if (op == OMIN) {
            R = vec2(A.x < B.x ? A.x : B.x, A.y < B.y ? A.y : B.y);
            if (!na && !nb) {
                if (A.y <= B.x)      choice = 1u;
                else if (B.y <= A.x) choice = 2u;
            }
        }
        else if (op == OMAX) {
            R = vec2(A.x > B.x ? A.x : B.x, A.y > B.y ? A.y : B.y);
            if (!na && !nb) {
                if (A.x >= B.y)      choice = 1u;
                else if (B.x >= A.y) choice = 2u;
            }
        }
        else if (op == OABS) {
            float lo = (A.x < 0.0) ? 0.0 : min(abs(A.x), abs(A.y));
            R = vec2(lo, max(abs(A.x), abs(A.y)));
        }
        else if (op == OSQUARE) {
            precise float l2 = A.x*A.x, u2 = A.y*A.y;
            float lo = (A.y > 0.0 && A.x < 0.0) ? 0.0 : min(u2, l2);
            R = vec2(lo, max(u2, l2));
        }
        else if (op == OSQRT) {
            R = vec2(A.x <= 0.0 ? 0.0 : sqrt(A.x),
                     A.y <= 0.0 ? 0.0 : sqrt(A.y));
            u = u || (A.x < 0.0);
        }
        else if (op == ONEG)   { R = vec2(-A.y, -A.x); }
        else if (op == OSIN || op == OCOS) { R = vec2(-1.0, 1.0); u = u || a_inf; }
        else if (op == OEXP)   { R = vec2(0.0, INF); }
        else if (op == OASIN || op == OATAN) { R = vec2(-1.5707964, 1.5707964); u = u || (op == OASIN && (A.x < -1.0 || A.y > 1.0)); }
        else if (op == OACOS)  { R = vec2(0.0, 3.1415927); u = u || A.x < -1.0 || A.y > 1.0; }
        else if (op == OATAN2) { R = vec2(-3.1415927, 3.1415927); u = true; }
        else if (op == OFLOOR) { R = vec2(floor(A.x), floor(A.y)); }
        else if (op == OCONST) { R = vec2(uintBitsToFloat(tapes[base+4u])); u = false; }
        else if (op == OCOPY)  { R = A; u = na; }
        else                   { R = vec2(-INF, INF); u = true; }

        u = u || isnan(R.x) || isnan(R.y);
        iv[o - N_CONST] = R;
        tnt[o - N_CONST] = u;
        set2(c, choice);
    }

    vec2 root = ldi(ROOT);
    bool rt = ldt(ROOT);
    if (!rt && root.x > 0.0)
    {
        slabinfo[id] = uvec4(0u, 0u, 0u, 0u);
        return;
    }
    if (!rt && root.y < 0.0)
    {
        slabinfo[id] = uvec4(1u, 0u, 0u, 0u);
        return;
    }

    /*  Ambiguous: simplify the tile tape against this slab's
     *  choices (tape_push backward walk, verdicts into the same
     *  2-bit lanes: 0 dead/elide, 1 keep, 2 copy-a, 3 copy-b).  */
    for (uint s = 0u; s < NSTATE; ++s)
        live[s] = false;
    if (ROOT >= N_CONST)
        live[ROOT - N_CONST] = true;

    uint kept = 0u;
    for (uint k = len; k-- > 0u; )
    {
        uint base = tbase + 5u*k;
        uint op = tapes[base], o = tapes[base+1u];
        uint sa = tapes[base+2u], sb = tapes[base+3u];
        if (!live[o - N_CONST])
        {
            set2(k, 0u);
            continue;
        }
        live[o - N_CONST] = false;
        uint ch = get2(k);
        if ((op == OMIN || op == OMAX) && ch != 0u)
        {
            uint keep = (ch == 1u) ? sa : sb;
            if (keep == o)
                set2(k, 0u);
            else
            {
                set2(k, ch == 1u ? 2u : 3u);
                ++kept;
            }
            if (keep >= N_CONST)
                live[keep - N_CONST] = true;
            continue;
        }
        set2(k, 1u);
        ++kept;
        bool unary = (op == OABS || op == OSQUARE || op == OSQRT ||
                      op == OSIN || op == OCOS || op == OTAN ||
                      op == OASIN || op == OACOS || op == OATAN ||
                      op == ONEG || op == OEXP || op == OFLOOR ||
                      op == OLOG || op == OCOPY);
        bool none  = (op == OCONST);
        if (!none)
        {
            if (sa >= N_CONST)
                live[sa - N_CONST] = true;
            if (!unary && sb >= N_CONST)
                live[sb - N_CONST] = true;
        }
    }

    uint off = atomicAdd(pool_next, kept * 5u);
    if (off + kept * 5u > pool_cap)
    {
        overflow = 1u;
        slabinfo[id] = uvec4(2u, 0u, 0u, 0u);   // len 0: P3 must not run
        return;
    }

    uint w = off;
    for (uint k = 0u; k < len; ++k)
    {
        uint v = get2(k);
        if (v == 0u)
            continue;
        uint base = tbase + 5u*k;
        uint w0 = tapes[base];
        uint w1 = tapes[base+1u];
        uint w2 = tapes[base+2u];
        uint w3 = tapes[base+3u];
        uint w4 = tapes[base+4u];
        if (v != 1u)
        {
            w2 = (v == 2u) ? w2 : w3;
            w0 = OCOPY;
            w3 = 0u;
            w4 = 0u;
        }
        pool[w+0u] = w0;
        pool[w+1u] = w1;
        pool[w+2u] = w2;
        pool[w+3u] = w3;
        pool[w+4u] = w4;
        w += 5u;
    }
    slabinfo[id] = uvec4(2u, off, kept, 0u);
}
)GLSL";
    return src;
}

std::string make_p3(const std::vector<uint32_t>& blob,
                    const BlobInfo& bi,
                    uint32_t ni, uint32_t nj, uint32_t nk)
{
    std::string src = "#version 430\n"
        "layout(local_size_x = 8, local_size_y = 8) in;\n"
        "layout(std430, binding=0) readonly buffer Blob { uint blob[]; };\n"
        "layout(std430, binding=1) readonly buffer BX { float Xs[]; };\n"
        "layout(std430, binding=2) readonly buffer BY { float Ys[]; };\n"
        "layout(std430, binding=3) readonly buffer BZ { float Zs[]; };\n"
        "layout(std430, binding=4) readonly buffer BL { uint Ls[]; };\n"
        "layout(std430, binding=5) writeonly buffer BI { uint img[]; };\n"
        "layout(std430, binding=7) readonly buffer BT { uvec2 tileinfo[]; };\n"
        "layout(std430, binding=8) readonly buffer BS { uint tapes[]; };\n"
        "layout(std430, binding=9) readonly buffer BSl { uvec4 slabinfo[]; };\n"
        "layout(std430, binding=10) readonly buffer BP { uint pool[]; };\n"
        "uniform uint T0;\n";
    src += op_prelude();
    src += blob_prelude(bi, ni, nj, nk);
    src += pinned_table(blob, bi);
    src += R"GLSL(
uniform uint J0;   // y-offset of this dispatch slice

float regs[NSTATE];

float ld(uint s)
{
    return s < N_CONST ? PINNED[s] : regs[s - N_CONST];
}

float min_f(float A, float B)
{
    if (isnan(B)) return A;
    if (isnan(A)) return B;
    return A < B ? A : B;
}
float max_f(float A, float B)
{
    if (isnan(B)) return A;
    if (isnan(A)) return B;
    return A > B ? A : B;
}

float eval_pool(uint tbase, uint len)
{
    for (uint c = 0u; c < len; ++c)
    {
        uint base = tbase + 5u*c;
        uint op = pool[base];
        uint o  = pool[base+1u];
        float A = ld(pool[base+2u]);
        float B = ld(pool[base+3u]);
        float R = 0.0;
        if (op == OADD) R = A + B;
        else if (op == OSUB) R = A - B;
        else if (op == OMUL) R = A * B;
        else if (op == ODIV) R = A / B;
        else if (op == OMIN) R = min_f(A, B);
        else if (op == OMAX) R = max_f(A, B);
        else if (op == OPOW) R = pow(A, B);
        else if (op == OMOD) R = (B == 0.0) ? 0.0 : A - B*floor(A/B);
        else if (op == OABS) R = abs(A);
        else if (op == OSQUARE) R = A * A;
        else if (op == OSQRT) R = (A < 0.0) ? 0.0 : sqrt(A);
        else if (op == OSIN) R = sin(A);
        else if (op == OCOS) R = cos(A);
        else if (op == OTAN) R = tan(A);
        else if (op == OASIN) R = (A < -1.0) ? -1.5707963705062866
                                : (A > 1.0) ? 1.5707963705062866 : asin(A);
        else if (op == OACOS) R = (A < -1.0) ? 3.1415927410125732
                                : (A > 1.0) ? 0.0 : acos(A);
        else if (op == OATAN) R = atan(A);
        else if (op == OATAN2) R = atan(A, B);
        else if (op == ONEG) R = -A;
        else if (op == OEXP) R = exp(A);
        else if (op == OFLOOR) R = floor(A);
        else if (op == OLOG) R = log((A > 1.17549e-38) ? A : 1.17549e-38);
        else if (op == OCONST) R = uintBitsToFloat(pool[base+4u]);
        else if (op == OCOPY) R = A;
        regs[o - N_CONST] = R;
    }
    return ld(ROOT);
}

void main()
{
    uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y + J0;
    if (i >= NI || j >= NJ)
        return;
    uint t = (j / TILE) * TILES_X + (i / TILE) - T0;   // band-local
    uint cls = tileinfo[t].x;
    if (cls == 0u)          { img[j*NI+i] = 0u; return; }
    else if (cls == 1u)     { img[j*NI+i] = Ls[NK]; return; }

    for (uint q = 0u; q < N_X; ++q)
        regs[blob[AX_OFF + q] - N_CONST] = Xs[i];
    for (uint q = 0u; q < N_Y; ++q)
        regs[blob[AX_OFF + N_X + q] - N_CONST] = Ys[j];

    uint result = 0u;
    /*  March by classified z-slab on the slab's OWN pool tape:
     *  empty slabs are skipped outright, a full slab means its
     *  topmost voxel is already inside, and only ambiguous slabs
     *  pay per-voxel evaluation - on a tape pushed for exactly
     *  that box.  */
    for (int sl = int(N_SLABS) - 1; sl >= 0 && result == 0u; --sl)
    {
        uvec4 si = slabinfo[t * N_SLABS + uint(sl)];
        if (si.x == 0u)
            continue;
        int ktop = sl * int(SLAB) + int(SLAB) - 1;
        if (ktop >= int(NK))
            ktop = int(NK) - 1;
        if (si.x == 1u)
        {
            result = Ls[ktop + 1];
            break;
        }
        for (int k = ktop; k >= sl * int(SLAB); --k)
        {
            for (uint q = 0u; q < N_Z; ++q)
                regs[blob[AX_OFF + N_X + N_Y + q] - N_CONST] = Zs[k];
            if (eval_pool(si.y, si.z) < 0.0)
            {
                result = Ls[k + 1];
                break;
            }
        }
    }
    img[j*NI+i] = result;
}
)GLSL";
    return src;
}

GLuint compile_program(const std::string& src)
{
    const GLuint sh = glCreateShader(GL_COMPUTE_SHADER);
    const char* s = src.c_str();
    glShaderSource(sh, 1, &s, nullptr);
    glCompileShader(sh);
    GLint ok = 0;
    glGetShaderiv(sh, GL_COMPILE_STATUS, &ok);
    if (!ok)
    {
        char log[2048];
        glGetShaderInfoLog(sh, sizeof(log), nullptr, log);
        fprintf(stderr, "gpu_render16: shader compile failed:\n%s\n", log);
        glDeleteShader(sh);
        return 0;
    }
    const GLuint prog = glCreateProgram();
    glAttachShader(prog, sh);
    glLinkProgram(prog);
    glGetProgramiv(prog, GL_LINK_STATUS, &ok);
    glDeleteShader(sh);
    if (!ok)
    {
        fprintf(stderr, "gpu_render16: shader link failed\n");
        glDeleteProgram(prog);
        return 0;
    }
    return prog;
}

GLuint make_ssbo(GLuint binding, const void* data, size_t bytes)
{
    GLuint buf = 0;
    glGenBuffers(1, &buf);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, buf);
    glBufferData(GL_SHADER_STORAGE_BUFFER, GLsizeiptr(bytes), data,
                 GL_STATIC_DRAW);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, binding, buf);
    return buf;
}

////////////////////////////////////////////////////////////////////////////
//  Program cache: keyed by (blob content hash, dims).  The GUI
//  re-renders the same decks at a handful of refinement sizes;
//  compiles are tens of ms and must not recur per frame.

struct CacheEntry
{
    uint64_t key = 0;
    GLuint p1 = 0, p2 = 0, p1b = 0, p3 = 0;
    uint64_t stamp = 0;
};
CacheEntry g_cache[16];
uint64_t g_stamp = 0;

uint64_t blob_key(const std::vector<uint32_t>& blob,
                  uint32_t ni, uint32_t nj, uint32_t nk)
{
    uint64_t h = 0xcbf29ce484222325ull;
    for (uint32_t w : blob)
        h = (h ^ w) * 0x100000001b3ull;
    h = (h ^ ni) * 0x100000001b3ull;
    h = (h ^ nj) * 0x100000001b3ull;
    h = (h ^ nk) * 0x100000001b3ull;
    return h ? h : 1;
}

bool cache_lookup(uint64_t key, GLuint* p1, GLuint* p2,
                  GLuint* p1b, GLuint* p3)
{
    for (auto& e : g_cache)
        if (e.key == key)
        {
            e.stamp = ++g_stamp;
            *p1 = e.p1;  *p2 = e.p2;  *p1b = e.p1b;  *p3 = e.p3;
            return true;
        }
    return false;
}

void cache_insert(uint64_t key, GLuint p1, GLuint p2, GLuint p1b,
                  GLuint p3)
{
    CacheEntry* victim = &g_cache[0];
    for (auto& e : g_cache)
    {
        if (e.key == 0) { victim = &e; break; }
        if (e.stamp < victim->stamp)   victim = &e;
    }
    if (victim->key)
    {
        glDeleteProgram(victim->p1);
        glDeleteProgram(victim->p2);
        glDeleteProgram(victim->p1b);
        glDeleteProgram(victim->p3);
    }
    *victim = CacheEntry{ key, p1, p2, p1b, p3, ++g_stamp };
}

}  // namespace

////////////////////////////////////////////////////////////////////////////

extern "C" bool gpu_available(void)
{
    std::lock_guard<std::mutex> guard(g_lock);
    return gpu().ok;
}

extern "C" const char* gpu_renderer_name(void)
{
    std::lock_guard<std::mutex> guard(g_lock);
    return gpu().ok ? gpu().renderer.c_str() : "";
}

extern "C" bool gpu_render16(const Deck* deck, Region r,
                             uint16_t** img, volatile int* halt)
{
    if (r.ni == 0 || r.nj == 0 || r.nk == 0)
        return false;

    const uint32_t need = tape_export_blob(deck, deck_base(deck),
                                           nullptr, 0);
    if (need == 0)
        return false;   // OP_GRID: host-side payloads, CPU only

    std::lock_guard<std::mutex> guard(g_lock);
    Gpu& g = gpu();
    if (!g.ok)
        return false;
    if (!p_eglMakeCurrent(g.dpy, EGL_NO_SURFACE, EGL_NO_SURFACE, g.ctx))
        return false;

    std::vector<uint32_t> blob(need);
    tape_export_blob(deck, deck_base(deck), blob.data(), need);
    const uint32_t orig_slots = blob[1];

    /*  Reschedule into a compact register file (Round 8); the
     *  untransformed blob is the fallback for anything the
     *  rescheduler doesn't recognize.  */
    {
        std::vector<uint32_t> tight = blob_reschedule(blob);
        if (!tight.empty())
            blob.swap(tight);
    }
    const BlobInfo bi = parse_blob(blob);

    /*  STIBIUM_GPU_DEBUG=1: dump deck shape + opcode census once
     *  per deck (helps chase CPU-vs-GPU pixel divergence to the op
     *  responsible).  */
    static const char* dbg = getenv("STIBIUM_GPU_DEBUG");
    if (dbg && atoi(dbg) != 0)
    {
        unsigned census[LAST_OP] = {};
        for (uint32_t c = 0; c < bi.n_clauses; ++c)
            census[blob[bi.clause_off + 5*c]]++;
        fprintf(stderr,
                "[gpu] deck: slots=%u (was %u; %u const + %u axes + "
                "%u regs) clauses=%u |",
                bi.num_slots, orig_slots, bi.n_const,
                bi.n_x + bi.n_y + bi.n_z,
                bi.num_slots - bi.n_const - bi.n_x - bi.n_y - bi.n_z,
                bi.n_clauses);
        for (unsigned o = 0; o < LAST_OP; ++o)
            if (census[o])
                fprintf(stderr, " %s=%u", OPCODE_NAMES[o], census[o]);
        fprintf(stderr, "\n");
    }

    using clk = std::chrono::steady_clock;
    clk::time_point t_start = clk::now(), t_compiled, t_upload;
    (void)t_upload;

    const uint64_t key = blob_key(blob, r.ni, r.nj, r.nk);
    GLuint p1, p2, p1b, p3;
    if (!cache_lookup(key, &p1, &p2, &p1b, &p3))
    {
        p1 = compile_program(make_p1(blob, bi, r.ni, r.nj, r.nk));
        p2 = compile_program(make_p2(blob, bi, r.ni, r.nj, r.nk));
        p1b = compile_program(make_p1b(blob, bi, r.ni, r.nj, r.nk));
        p3 = compile_program(make_p3(blob, bi, r.ni, r.nj, r.nk));
        if (!p1 || !p2 || !p1b || !p3)
        {
            if (p1) glDeleteProgram(p1);
            if (p2) glDeleteProgram(p2);
            if (p1b) glDeleteProgram(p1b);
            if (p3) glDeleteProgram(p3);
            p_eglMakeCurrent(g.dpy, EGL_NO_SURFACE, EGL_NO_SURFACE,
                             EGL_NO_CONTEXT);
            return false;
        }
        cache_insert(key, p1, p2, p1b, p3);
    }
    t_compiled = clk::now();

    const uint32_t tiles_x = (r.ni + GPU_TILE - 1) / GPU_TILE;
    const uint32_t tiles_y = (r.nj + GPU_TILE - 1) / GPU_TILE;
    const uint32_t n_tiles = tiles_x * tiles_y;
    (void)n_tiles;

    /*  Band rendering: the image is processed in tile-row bands and
     *  every working buffer (choices, tile tapes, slab-tape pool) is
     *  recycled between bands - MPR never keeps a whole frame's
     *  tapes alive at once, and neither do we.  Band height is
     *  chosen so the fixed-slot tile tapes stay modest; the slab
     *  pool is a fixed budget refilled per band.  */
    uint32_t band_trows = tiles_y;
    {
        const uint64_t per_trow =
                uint64_t(tiles_x) * bi.n_clauses * 5 * 4;
        const uint64_t budget = uint64_t(128) << 20;
        if (per_trow > 0 && per_trow * band_trows > budget)
        {
            band_trows = uint32_t(budget / per_trow);
            if (band_trows == 0)
                band_trows = 1;
        }
    }
    const uint32_t band_tiles = tiles_x * band_trows;

    std::vector<uint32_t> L(r.nk + 1);
    for (uint32_t k = 0; k <= r.nk; ++k)
        L[k] = r.L[k];
    std::vector<uint32_t> out(size_t(r.ni) * r.nj);

    const uint32_t n_slabs = (r.nk + GPU_SLAB - 1) / GPU_SLAB;

    /*  Slab-tape pool (compacted, atomically bump-allocated by pass
     *  1b).  STIBIUM_GPU_POOL_MB tunes it; overflow falls back to
     *  the CPU renderer.  */
    uint32_t pool_mb = 384;
    if (const char* env = getenv("STIBIUM_GPU_POOL_MB"))
        if (atoi(env) > 0)
            pool_mb = uint32_t(atoi(env));
    const uint32_t pool_words = pool_mb * (1024u * 1024u / 4u);
    uint32_t alloc_init[3] = { 0, pool_words, 0 };

    GLuint bufs[12] = {
        make_ssbo(0, blob.data(), blob.size() * 4),
        make_ssbo(1, r.X, (r.ni + 1) * 4),
        make_ssbo(2, r.Y, (r.nj + 1) * 4),
        make_ssbo(3, r.Z, (r.nk + 1) * 4),
        make_ssbo(4, L.data(), L.size() * 4),
        make_ssbo(5, nullptr, out.size() * 4),
        make_ssbo(6, nullptr, size_t(band_tiles) * bi.n_clauses * 4),
        make_ssbo(7, nullptr, size_t(band_tiles) * 8),
        make_ssbo(8, nullptr, size_t(band_tiles) * bi.n_clauses * 5 * 4),
        make_ssbo(9, nullptr, size_t(band_tiles) * n_slabs * 16),
        make_ssbo(10, nullptr, size_t(pool_words) * 4),
        make_ssbo(11, alloc_init, sizeof(alloc_init)),
    };

    /*  STIBIUM_GPU=2: brute-force mode - skip classify/simplify and
     *  march the FULL tape per pixel.  P3 needs no changes: with the
     *  blob bound at the tapes binding, every tile pre-marked
     *  ambiguous, and len = N_CLAUSES... wait, tbase must be
     *  CLAUSE_OFF for every tile; see tileinfo fill below.  Used to
     *  bisect pipeline bugs from pointwise/orchestration ones.  */
    t_upload = clk::now();
    static const char* mode_env = getenv("STIBIUM_GPU");
    const bool brute = mode_env && atoi(mode_env) == 2;

    const GLint p1_t0 = glGetUniformLocation(p1, "T0");
    const GLint p1_nt = glGetUniformLocation(p1, "NT");
    const GLint p2_nt = glGetUniformLocation(p2, "NT");
    const GLint p1b_t0 = glGetUniformLocation(p1b, "T0");
    const GLint p1b_nt = glGetUniformLocation(p1b, "NT");
    const GLint p3_t0 = glGetUniformLocation(p3, "T0");
    const GLint p3_j0 = glGetUniformLocation(p3, "J0");

    if (brute)
    {
        /*  Brute-force diagnostic mode: every band tile pre-marked
         *  ambiguous, slab tapes = the full clause list replicated
         *  per band tile in the pool.  */
        std::vector<uint32_t> ti(size_t(band_tiles) * 2);
        for (uint32_t t = 0; t < band_tiles; ++t)
        {
            ti[2*t]   = 2u;
            ti[2*t+1] = bi.n_clauses;
        }
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufs[7]);
        glBufferData(GL_SHADER_STORAGE_BUFFER,
                     GLsizeiptr(ti.size() * 4), ti.data(),
                     GL_STATIC_DRAW);
        std::vector<uint32_t> rep(size_t(band_tiles) * bi.n_clauses * 5);
        const uint32_t* cl = blob.data() + bi.clause_off;
        for (uint32_t t = 0; t < band_tiles; ++t)
            memcpy(rep.data() + size_t(t) * bi.n_clauses * 5, cl,
                   size_t(bi.n_clauses) * 5 * 4);
        std::vector<uint32_t> si(size_t(band_tiles) * n_slabs * 4);
        for (uint32_t t = 0; t < band_tiles; ++t)
            for (uint32_t q = 0; q < n_slabs; ++q)
            {
                uint32_t* e = si.data() +
                        (size_t(t) * n_slabs + q) * 4;
                e[0] = 2u;
                e[1] = t * bi.n_clauses * 5;
                e[2] = bi.n_clauses;
                e[3] = 0u;
            }
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufs[9]);
        glBufferData(GL_SHADER_STORAGE_BUFFER,
                     GLsizeiptr(si.size() * 4), si.data(),
                     GL_STATIC_DRAW);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufs[10]);
        glBufferData(GL_SHADER_STORAGE_BUFFER,
                     GLsizeiptr(rep.size() * 4), rep.data(),
                     GL_STATIC_DRAW);
    }

    /*  Per-pixel march slice budget (watchdog protection: a killed
     *  dispatch returns normally with a half-written buffer).  */
    const uint64_t per_row =
            uint64_t(r.ni) * r.nk * (bi.n_clauses ? bi.n_clauses : 1);
    uint32_t rows_per_slice = uint32_t(
            (uint64_t(2) << 30) / (per_row ? per_row : 1));
    rows_per_slice = rows_per_slice < 8 ? 8
                   : (rows_per_slice / 8) * 8;

    double ms_p1 = 0, ms_p2 = 0, ms_p1b = 0, ms_p3 = 0;
    double pool_peak_mb = 0;
    const auto ms_since = [](clk::time_point& t) {
        const auto now = clk::now();
        const double d =
                std::chrono::duration<double, std::milli>(now - t)
                        .count();
        t = now;
        return d;
    };

    /*  The band loop: classify, simplify, and march one horizontal
     *  band of tiles at a time, recycling every working buffer.
     *  Slab-tape demand can't be known before classification, so a
     *  band that overflows the pool is simply retried at half
     *  height (down to a single tile-row before giving up).  */
    bool failed = false;
    uint32_t cur_trows = band_trows;
    for (uint32_t trow = 0; trow < tiles_y && !failed; )
    {
        if (halt && *halt)
            break;
        const uint32_t btrows =
                tiles_y - trow < cur_trows ? tiles_y - trow
                                           : cur_trows;
        const uint32_t T0 = trow * tiles_x;
        const uint32_t NT = btrows * tiles_x;
        clk::time_point tick = clk::now();

        if (!brute)
        {
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufs[11]);
            glBufferData(GL_SHADER_STORAGE_BUFFER,
                         sizeof(alloc_init), alloc_init,
                         GL_STATIC_DRAW);
            glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 11, bufs[11]);

            glUseProgram(p1);
            glUniform1ui(p1_t0, T0);
            glUniform1ui(p1_nt, NT);
            glDispatchCompute((NT + 63) / 64, 1, 1);
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
            if (dbg && atoi(dbg) != 0)
                { glFinish(); ms_p1 += ms_since(tick); }

            glUseProgram(p2);
            glUniform1ui(p2_nt, NT);
            glDispatchCompute((NT + 63) / 64, 1, 1);
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
            if (dbg && atoi(dbg) != 0)
                { glFinish(); ms_p2 += ms_since(tick); }

            glUseProgram(p1b);
            glUniform1ui(p1b_t0, T0);
            glUniform1ui(p1b_nt, NT);
            glDispatchCompute((NT * n_slabs + 63) / 64, 1, 1);
            glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
            glFinish();
            if (dbg && atoi(dbg) != 0)
                ms_p1b += ms_since(tick);

            uint32_t alloc_state[3] = {};
            glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufs[11]);
            glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0,
                               sizeof(alloc_state), alloc_state);
            if (alloc_state[2] != 0)
            {
                if (btrows > 1)
                {
                    cur_trows = btrows / 2;   // retry band narrower
                    continue;
                }
                if (dbg && atoi(dbg) != 0)
                    fprintf(stderr,
                            "[gpu] slab-tape pool overflow "
                            "(%u MB) at one tile-row - CPU "
                            "fallback\n", pool_mb);
                failed = true;
                break;
            }
            const double used =
                    double(alloc_state[0]) * 4.0 / (1024 * 1024);
            if (used > pool_peak_mb)
                pool_peak_mb = used;
        }

        /*  March this band's pixel rows in watchdog-safe slices  */
        const uint32_t band_j0 = trow * GPU_TILE;
        const uint32_t band_j1 =
                (trow + btrows) * GPU_TILE < r.nj
                        ? (trow + btrows) * GPU_TILE : r.nj;
        glUseProgram(p3);
        glUniform1ui(p3_t0, T0);
        tick = clk::now();
        for (uint32_t j0 = band_j0; j0 < band_j1;
             j0 += rows_per_slice)
        {
            if (halt && *halt)
                break;
            const uint32_t rows =
                    band_j1 - j0 < rows_per_slice ? band_j1 - j0
                                                  : rows_per_slice;
            glUniform1ui(p3_j0, j0);
            glDispatchCompute((r.ni + 7) / 8, (rows + 7) / 8, 1);
            glFinish();
        }
        if (dbg && atoi(dbg) != 0)
            ms_p3 += ms_since(tick);
        trow += btrows;
    }
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glFinish();

    if (dbg && atoi(dbg) != 0)
    {
        const auto msd = [](auto a, auto b) {
            return std::chrono::duration<double, std::milli>(b - a)
                    .count();
        };
        fprintf(stderr,
                "[gpu] compile %.1f ms | bands of %u tile-rows | "
                "p1 %.1f ms, p2 %.1f ms, p1b %.1f ms, p3 %.1f ms | "
                "pool peak %.1f MB\n",
                msd(t_start, t_compiled), band_trows,
                ms_p1, ms_p2, ms_p1b, ms_p3, pool_peak_mb);
    }
    if (failed)
    {
        glDeleteBuffers(12, bufs);
        p_eglMakeCurrent(g.dpy, EGL_NO_SURFACE, EGL_NO_SURFACE,
                         EGL_NO_CONTEXT);
        return false;
    }
    if (glGetError() != GL_NO_ERROR)
    {
        glDeleteBuffers(12, bufs);
        p_eglMakeCurrent(g.dpy, EGL_NO_SURFACE, EGL_NO_SURFACE,
                         EGL_NO_CONTEXT);
        return false;
    }

    bool wrote = false;
    if (!halt || !*halt)
    {
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufs[5]);
        glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0,
                           GLsizeiptr(out.size() * 4), out.data());

        /*  STIBIUM_GPU_DEBUG=2: render the same region on the CPU
         *  and report the depth-delta distribution (in L levels).
         *  +-1-level deltas hugging the silhouette = ulp-level
         *  pointwise divergence; anything larger = structural.  */
        if (dbg && atoi(dbg) >= 2)
        {
            std::vector<uint16_t> ref(out.size(), 0);
            std::vector<uint16_t*> rows(r.nj);
            for (uint32_t j = 0; j < r.nj; ++j)
                rows[j] = ref.data() + size_t(j) * r.ni;
            Region rr = r;
            rr.imin = 0;
            rr.jmin = 0;
            TapeCtx* ctx = tape_ctx_new(deck);
            volatile int h2 = 0;
            render16_tape(deck_base(deck), ctx, rr, rows.data(),
                          &h2, nullptr);
            tape_ctx_free(ctx);

            const uint32_t lstep = r.nk ? 65535 / r.nk : 1;
            uint64_t n_diff = 0, n_one = 0;
            int64_t max_d = 0;
            int shown = 0;
            for (size_t p = 0; p < out.size(); ++p)
            {
                const int64_t d =
                        int64_t(out[p]) - int64_t(ref[p]);
                if (d == 0)
                    continue;
                ++n_diff;
                const int64_t lv =
                        (d > 0 ? d : -d) / (lstep ? lstep : 1);
                if (lv <= 1)
                    ++n_one;
                if (lv > max_d)
                    max_d = lv;
                if (shown < 10)
                {
                    fprintf(stderr,
                            "[gpu]   px(%zu,%zu) gpu=%u cpu=%u "
                            "(%+lld levels)\n",
                            p % r.ni, p / r.ni, out[p], ref[p],
                            (long long)(d / (lstep ? lstep : 1)));
                    ++shown;
                }
            }
            fprintf(stderr,
                    "[gpu] depth diff vs CPU: %llu px of %zu "
                    "(%llu within 1 level, max %lld levels, "
                    "nk=%u lstep=%u)\n",
                    (unsigned long long)n_diff, out.size(),
                    (unsigned long long)n_one, (long long)max_d,
                    r.nk, lstep);
        }
        /*  Max-merge, exactly like render16_tape's img writes  */
        for (uint32_t j = 0; j < r.nj; ++j)
        {
            uint16_t* row = img[r.jmin + j] + r.imin;
            const uint32_t* src = out.data() + size_t(j) * r.ni;
            for (uint32_t i = 0; i < r.ni; ++i)
                if (src[i] > row[i])
                    row[i] = uint16_t(src[i]);
        }
        wrote = true;
    }

    glDeleteBuffers(12, bufs);
    p_eglMakeCurrent(g.dpy, EGL_NO_SURFACE, EGL_NO_SURFACE,
                     EGL_NO_CONTEXT);
    return wrote;
}

#else  // !STIBIUM_HAS_EGL

extern "C" bool gpu_available(void) { return false; }
extern "C" const char* gpu_renderer_name(void) { return ""; }
extern "C" bool gpu_render16(const Deck*, Region, uint16_t**,
                             volatile int*)
{
    return false;
}

#endif

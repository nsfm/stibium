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
        "const uint TILE=%uu,TILES_X=%uu,TILES_Y=%uu;\n",
        ni, nj, nk, bi.root, bi.n_const, bi.n_x, bi.n_y, bi.n_z,
        bi.n_clauses, bi.const_off, bi.ax_off, bi.clause_off,
        bi.num_slots, GPU_TILE,
        (ni + GPU_TILE - 1) / GPU_TILE, (nj + GPU_TILE - 1) / GPU_TILE);
    return buf;
}

std::string make_p1(const BlobInfo& bi,
                    uint32_t ni, uint32_t nj, uint32_t nk)
{
    std::string src = "#version 430\n"
        "layout(local_size_x = 64) in;\n"
        "layout(std430, binding=0) readonly buffer Blob { uint blob[]; };\n"
        "layout(std430, binding=1) readonly buffer BX { float Xs[]; };\n"
        "layout(std430, binding=2) readonly buffer BY { float Ys[]; };\n"
        "layout(std430, binding=3) readonly buffer BZ { float Zs[]; };\n"
        "layout(std430, binding=6) writeonly buffer BC { uint choices[]; };\n"
        "layout(std430, binding=7) buffer BT { uvec2 tileinfo[]; };\n";
    src += op_prelude();
    src += blob_prelude(bi, ni, nj, nk);
    src += R"GLSL(
const float INF = uintBitsToFloat(0x7F800000u);

precise vec2 iv[NSLOTS];
bool tnt[NSLOTS];

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
    uint t = gl_GlobalInvocationID.x;
    if (t >= TILES_X * TILES_Y)
        return;
    uint ti = t - (t / TILES_X) * TILES_X, tj = t / TILES_X;

    for (uint q = 0u; q < N_CONST; ++q)
    {
        uint s = blob[CONST_OFF + 2u*q];
        iv[s] = vec2(uintBitsToFloat(blob[CONST_OFF + 2u*q + 1u]));
        tnt[s] = false;
    }
    vec2 X = vec2(Xs[ti*TILE], Xs[min(ti*TILE + TILE, NI)]);
    vec2 Y = vec2(Ys[tj*TILE], Ys[min(tj*TILE + TILE, NJ)]);
    vec2 Z = vec2(Zs[0], Zs[NK]);
    for (uint q = 0u; q < N_X; ++q)
        { iv[blob[AX_OFF+q]] = X; tnt[blob[AX_OFF+q]] = false; }
    for (uint q = 0u; q < N_Y; ++q)
        { iv[blob[AX_OFF+N_X+q]] = Y; tnt[blob[AX_OFF+N_X+q]] = false; }
    for (uint q = 0u; q < N_Z; ++q)
        { iv[blob[AX_OFF+N_X+N_Y+q]] = Z; tnt[blob[AX_OFF+N_X+N_Y+q]] = false; }

    for (uint c = 0u; c < N_CLAUSES; ++c)
    {
        uint base = CLAUSE_OFF + 5u*c;
        uint op = blob[base], o = blob[base+1u];
        uint sa = blob[base+2u], sb = blob[base+3u];
        vec2 A = iv[sa], B = iv[sb];
        bool na = tnt[sa], nb = tnt[sb];
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
        iv[o] = R;
        tnt[o] = u;
        choices[t * N_CLAUSES + c] = choice;
    }

    vec2 root = iv[ROOT];
    bool rt = tnt[ROOT];
    uint cls = 2u;
    if (!rt && root.x > 0.0)  cls = 0u;
    else if (!rt && root.y < 0.0) cls = 1u;
    tileinfo[t] = uvec2(cls, N_CLAUSES);
}
)GLSL";
    return src;
}

std::string make_p2(const BlobInfo& bi,
                    uint32_t ni, uint32_t nj, uint32_t nk)
{
    std::string src = "#version 430\n"
        "layout(local_size_x = 64) in;\n"
        "layout(std430, binding=0) readonly buffer Blob { uint blob[]; };\n"
        "layout(std430, binding=6) readonly buffer BC { uint choices[]; };\n"
        "layout(std430, binding=7) buffer BT { uvec2 tileinfo[]; };\n"
        "layout(std430, binding=8) writeonly buffer BS { uint tapes[]; };\n";
    src += op_prelude();
    src += blob_prelude(bi, ni, nj, nk);
    src += R"GLSL(
bool live[NSLOTS];
uint verdict[N_CLAUSES];

void main()
{
    uint t = gl_GlobalInvocationID.x;
    if (t >= TILES_X * TILES_Y)
        return;
    if (tileinfo[t].x != 2u)
        return;

    for (uint s = 0u; s < NSLOTS; ++s)
        live[s] = false;
    live[ROOT] = true;

    for (uint k = N_CLAUSES; k-- > 0u; )
    {
        uint base = CLAUSE_OFF + 5u*k;
        uint op = blob[base], o = blob[base+1u];
        uint sa = blob[base+2u], sb = blob[base+3u];
        if (!live[o])
        {
            verdict[k] = 0u;
            continue;
        }
        live[o] = false;
        uint ch = choices[t * N_CLAUSES + k];
        if ((op == OMIN || op == OMAX) && ch != 0u)
        {
            uint keep = (ch == 1u) ? sa : sb;
            verdict[k] = (keep == o) ? 0u : (ch == 1u ? 2u : 3u);
            live[keep] = true;
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
            live[sa] = true;
            if (!unary)
                live[sb] = true;
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

std::string make_p3(const BlobInfo& bi,
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
        "layout(std430, binding=8) readonly buffer BS { uint tapes[]; };\n";
    src += op_prelude();
    src += blob_prelude(bi, ni, nj, nk);
    src += R"GLSL(
uniform uint J0;   // y-offset of this dispatch slice

float regs[NSLOTS];

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

float eval_tape(uint tbase, uint len)
{
    for (uint c = 0u; c < len; ++c)
    {
        uint base = tbase + 5u*c;
        uint op = tapes[base];
        uint o  = tapes[base+1u];
        float A = regs[tapes[base+2u]];
        float B = regs[tapes[base+3u]];
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
        else if (op == OCONST) R = uintBitsToFloat(tapes[base+4u]);
        else if (op == OCOPY) R = A;
        regs[o] = R;
    }
    return regs[ROOT];
}

void main()
{
    uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y + J0;
    if (i >= NI || j >= NJ)
        return;
    uint t = (j / TILE) * TILES_X + (i / TILE);
    uint cls = tileinfo[t].x;
    if (cls == 0u)          { img[j*NI+i] = 0u; return; }
    else if (cls == 1u)     { img[j*NI+i] = Ls[NK]; return; }

    for (uint q = 0u; q < N_CONST; ++q)
        regs[blob[CONST_OFF + 2u*q]] =
                uintBitsToFloat(blob[CONST_OFF + 2u*q + 1u]);
    for (uint q = 0u; q < N_X; ++q)
        regs[blob[AX_OFF + q]] = Xs[i];
    for (uint q = 0u; q < N_Y; ++q)
        regs[blob[AX_OFF + N_X + q]] = Ys[j];

    uint tbase = t * N_CLAUSES * 5u;
    uint len = tileinfo[t].y;
    uint result = 0u;
    for (int k = int(NK) - 1; k >= 0; --k)
    {
        for (uint q = 0u; q < N_Z; ++q)
            regs[blob[AX_OFF + N_X + N_Y + q]] = Zs[k];
        if (eval_tape(tbase, len) < 0.0)
        {
            result = Ls[k + 1];
            break;
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
    GLuint p1 = 0, p2 = 0, p3 = 0;
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

bool cache_lookup(uint64_t key, GLuint* p1, GLuint* p2, GLuint* p3)
{
    for (auto& e : g_cache)
        if (e.key == key)
        {
            e.stamp = ++g_stamp;
            *p1 = e.p1;  *p2 = e.p2;  *p3 = e.p3;
            return true;
        }
    return false;
}

void cache_insert(uint64_t key, GLuint p1, GLuint p2, GLuint p3)
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
        glDeleteProgram(victim->p3);
    }
    *victim = CacheEntry{ key, p1, p2, p3, ++g_stamp };
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
        fprintf(stderr, "[gpu] deck: slots=%u clauses=%u |",
                bi.num_slots, bi.n_clauses);
        for (unsigned o = 0; o < LAST_OP; ++o)
            if (census[o])
                fprintf(stderr, " %s=%u", OPCODE_NAMES[o], census[o]);
        fprintf(stderr, "\n");
    }

    const uint64_t key = blob_key(blob, r.ni, r.nj, r.nk);
    GLuint p1, p2, p3;
    if (!cache_lookup(key, &p1, &p2, &p3))
    {
        p1 = compile_program(make_p1(bi, r.ni, r.nj, r.nk));
        p2 = compile_program(make_p2(bi, r.ni, r.nj, r.nk));
        p3 = compile_program(make_p3(bi, r.ni, r.nj, r.nk));
        if (!p1 || !p2 || !p3)
        {
            if (p1) glDeleteProgram(p1);
            if (p2) glDeleteProgram(p2);
            if (p3) glDeleteProgram(p3);
            p_eglMakeCurrent(g.dpy, EGL_NO_SURFACE, EGL_NO_SURFACE,
                             EGL_NO_CONTEXT);
            return false;
        }
        cache_insert(key, p1, p2, p3);
    }

    const uint32_t tiles_x = (r.ni + GPU_TILE - 1) / GPU_TILE;
    const uint32_t tiles_y = (r.nj + GPU_TILE - 1) / GPU_TILE;
    const uint32_t n_tiles = tiles_x * tiles_y;

    /*  Per-tile tapes are fixed n_clauses-sized slots for now; big
     *  decks at big resolutions want more VRAM than is reasonable.
     *  Cap and fall back to the CPU - compacted on-device
     *  allocation (fidget's atomic-counter scheme) is the known
     *  stage-C fix (doc/TAPE-NEXT.md §4).  */
    const uint64_t tape_bytes =
            uint64_t(n_tiles) * bi.n_clauses * 5 * 4;
    if (tape_bytes > (uint64_t(512) << 20))
    {
        p_eglMakeCurrent(g.dpy, EGL_NO_SURFACE, EGL_NO_SURFACE,
                         EGL_NO_CONTEXT);
        return false;
    }

    std::vector<uint32_t> L(r.nk + 1);
    for (uint32_t k = 0; k <= r.nk; ++k)
        L[k] = r.L[k];
    std::vector<uint32_t> out(size_t(r.ni) * r.nj);

    GLuint bufs[9] = {
        make_ssbo(0, blob.data(), blob.size() * 4),
        make_ssbo(1, r.X, (r.ni + 1) * 4),
        make_ssbo(2, r.Y, (r.nj + 1) * 4),
        make_ssbo(3, r.Z, (r.nk + 1) * 4),
        make_ssbo(4, L.data(), L.size() * 4),
        make_ssbo(5, nullptr, out.size() * 4),
        make_ssbo(6, nullptr, size_t(n_tiles) * bi.n_clauses * 4),
        make_ssbo(7, nullptr, size_t(n_tiles) * 8),
        make_ssbo(8, nullptr, size_t(n_tiles) * bi.n_clauses * 5 * 4),
    };

    /*  STIBIUM_GPU=2: brute-force mode - skip classify/simplify and
     *  march the FULL tape per pixel.  P3 needs no changes: with the
     *  blob bound at the tapes binding, every tile pre-marked
     *  ambiguous, and len = N_CLAUSES... wait, tbase must be
     *  CLAUSE_OFF for every tile; see tileinfo fill below.  Used to
     *  bisect pipeline bugs from pointwise/orchestration ones.  */
    static const char* mode_env = getenv("STIBIUM_GPU");
    const bool brute = mode_env && atoi(mode_env) == 2;
    if (brute)
    {
        /*  tileinfo = AMBIG everywhere; P3 computes tbase as
         *  t*N_CLAUSES*5, so instead of binding the blob we fill
         *  the tapes buffer with the clause section replicated per
         *  tile - simpler than special-casing tbase in the shader,
         *  and this is a diagnostic path.  */
        std::vector<uint32_t> ti(size_t(n_tiles) * 2);
        for (uint32_t t = 0; t < n_tiles; ++t)
        {
            ti[2*t]   = 2u;
            ti[2*t+1] = bi.n_clauses;
        }
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufs[7]);
        glBufferData(GL_SHADER_STORAGE_BUFFER,
                     GLsizeiptr(ti.size() * 4), ti.data(),
                     GL_STATIC_DRAW);
        std::vector<uint32_t> rep(size_t(n_tiles) * bi.n_clauses * 5);
        const uint32_t* cl = blob.data() + bi.clause_off;
        for (uint32_t t = 0; t < n_tiles; ++t)
            memcpy(rep.data() + size_t(t) * bi.n_clauses * 5, cl,
                   size_t(bi.n_clauses) * 5 * 4);
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufs[8]);
        glBufferData(GL_SHADER_STORAGE_BUFFER,
                     GLsizeiptr(rep.size() * 4), rep.data(),
                     GL_STATIC_DRAW);
    }
    else
    {
        glUseProgram(p1);
        glDispatchCompute((n_tiles + 63) / 64, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
        glUseProgram(p2);
        glDispatchCompute((n_tiles + 63) / 64, 1, 1);
        glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    }
    /*  March in y-slices: one monolithic dispatch on a big model
     *  can run for minutes and the driver watchdog kills it MID-
     *  FLIGHT with no GL error - we read back a half-written buffer
     *  (learned on the merged Zeiss: the math was exact, the
     *  dispatch was dead).  Budget each slice to ~2e9 clause-evals
     *  and glFinish between slices so every submission stays far
     *  under any watchdog.  */
    glUseProgram(p3);
    const GLint j0_loc = glGetUniformLocation(p3, "J0");
    const uint64_t per_row =
            uint64_t(r.ni) * r.nk * (bi.n_clauses ? bi.n_clauses : 1);
    uint32_t rows_per_slice = uint32_t(
            (uint64_t(2) << 30) / (per_row ? per_row : 1));
    rows_per_slice = rows_per_slice < 8 ? 8
                   : (rows_per_slice / 8) * 8;
    for (uint32_t j0 = 0; j0 < r.nj; j0 += rows_per_slice)
    {
        if (halt && *halt)
            break;
        const uint32_t rows =
                r.nj - j0 < rows_per_slice ? r.nj - j0 : rows_per_slice;
        glUniform1ui(j0_loc, j0);
        glDispatchCompute((r.ni + 7) / 8, (rows + 7) / 8, 1);
        glFinish();
    }
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glFinish();
    if (glGetError() != GL_NO_ERROR)
    {
        glDeleteBuffers(9, bufs);
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

    glDeleteBuffers(9, bufs);
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

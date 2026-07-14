/*
 *  GPU tape interpreter prototype: renders a heightmap on the GPU by
 *  brute-force per-pixel ray march over the exported tape bytecode,
 *  and referees it against the CPU renderer.
 *
 *  Hidden tag "[.gpu]" - needs a working EGL + OpenGL 4.3 stack and
 *  a GPU.  Run explicitly:  SbFabTest "[gpu]"
 *  (prime-run SbFabTest "[gpu]" to use the discrete GPU.)
 *
 *  This is stage A of the GPU rung (doc/TAPE-NEXT.md §4): prove the
 *  bytecode + a GLSL interpreter reproduce the CPU renderer's exact
 *  pixels.  The MPR pipeline (interval tiles + on-device pushes,
 *  after fidget-wgpu) builds on this once the referee is green.
 *
 *  Exactness contract: the CPU renderer's BINARY pushes preserve
 *  signs exactly, and the pixel value is a function of signs only
 *  (highest k with field < 0), so a bit-exact GPU eval yields a
 *  bit-identical heightmap.  IEEE basic ops (+ - * /) match by
 *  construction; min/max mirror math_f's exact semantics below.
 *  sqrt and transcendentals are NOT guaranteed correctly-rounded in
 *  GLSL, so models are split into an exact set (algebra + min/max)
 *  and a tolerance set (sqrt/trig), and the test reports which.
 */

#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include <EGL/egl.h>
#include <EGL/eglext.h>
#include <GL/glcorearb.h>

#include "catch/catch.hpp"

#include "fab/tree/tape.h"
#include "fab/tree/gpu.h"
#include "fab/tree/tree.h"
#include "fab/tree/render.h"
#include "fab/tree/parser.h"
#include "fab/tree/node/opcodes.h"
#include "fab/util/region.h"

namespace {

////////////////////////////////////////////////////////////////////////////
//  Minimal GL loader: every entry point through eglGetProcAddress, so we
//  only link libEGL and never touch GLX / GLVND symbol questions.

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
    X(PFNGLGETPROGRAMINFOLOGPROC,     glGetProgramInfoLog) \
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
    X(PFNGLFINISHPROC,                glFinish)

#define DECLARE(type, name) type name = nullptr;
GL_FUNCS(DECLARE)
#undef DECLARE

bool load_gl()
{
#define LOAD(type, name) \
    name = reinterpret_cast<type>(eglGetProcAddress(#name)); \
    if (!name) return false;
    GL_FUNCS(LOAD)
#undef LOAD
    return true;
}

////////////////////////////////////////////////////////////////////////////
//  Headless EGL context (surfaceless).  Returns false - and the tests
//  skip with a warning - when no usable GL 4.3 stack is present.

struct GpuContext
{
    EGLDisplay dpy = EGL_NO_DISPLAY;
    EGLContext ctx = EGL_NO_CONTEXT;
    bool ok = false;
    std::string renderer;

    GpuContext()
    {
        dpy = eglGetDisplay(EGL_DEFAULT_DISPLAY);
        if (dpy == EGL_NO_DISPLAY)
            return;
        EGLint maj = 0, min = 0;
        if (!eglInitialize(dpy, &maj, &min))
            return;
        if (!eglBindAPI(EGL_OPENGL_API))
            return;

        const EGLint cfg_attr[] = {
            EGL_SURFACE_TYPE,    EGL_PBUFFER_BIT,
            EGL_RENDERABLE_TYPE, EGL_OPENGL_BIT,
            EGL_NONE
        };
        EGLConfig cfg;
        EGLint n = 0;
        if (!eglChooseConfig(dpy, cfg_attr, &cfg, 1, &n) || n < 1)
            return;

        const EGLint ctx_attr[] = {
            EGL_CONTEXT_MAJOR_VERSION, 4,
            EGL_CONTEXT_MINOR_VERSION, 3,
            EGL_NONE
        };
        ctx = eglCreateContext(dpy, cfg, EGL_NO_CONTEXT, ctx_attr);
        if (ctx == EGL_NO_CONTEXT)
            return;
        if (!eglMakeCurrent(dpy, EGL_NO_SURFACE, EGL_NO_SURFACE, ctx))
            return;
        if (!load_gl())
            return;

        const GLubyte* r = glGetString(GL_RENDERER);
        renderer = r ? reinterpret_cast<const char*>(r) : "(unknown)";
        ok = true;
    }

    ~GpuContext()
    {
        if (dpy != EGL_NO_DISPLAY)
        {
            eglMakeCurrent(dpy, EGL_NO_SURFACE, EGL_NO_SURFACE,
                           EGL_NO_CONTEXT);
            if (ctx != EGL_NO_CONTEXT)
                eglDestroyContext(dpy, ctx);
            eglTerminate(dpy);
        }
    }
};

GpuContext& gpu()
{
    static GpuContext c;
    return c;
}

////////////////////////////////////////////////////////////////////////////
//  Shader generation: the blob header is parsed on the host and every
//  offset (plus the slot count and the C opcode values) is baked into
//  the source, so the GLSL side has no protocol logic to get wrong.

struct BlobInfo
{
    uint32_t num_slots, root;
    uint32_t n_const, n_x, n_y, n_z, n_clauses;
    uint32_t const_off, ax_off, clause_off;
};

BlobInfo parse_blob(const std::vector<uint32_t>& b)
{
    BlobInfo i;
    REQUIRE(b[0] == 1);   // format version
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

std::string make_shader(const BlobInfo& bi,
                        uint32_t ni, uint32_t nj, uint32_t nk)
{
    char buf[8192];
    snprintf(buf, sizeof(buf), R"GLSL(#version 430
layout(local_size_x = 8, local_size_y = 8) in;

layout(std430, binding = 0) readonly buffer Blob { uint blob[]; };
layout(std430, binding = 1) readonly buffer BX   { float Xs[]; };
layout(std430, binding = 2) readonly buffer BY   { float Ys[]; };
layout(std430, binding = 3) readonly buffer BZ   { float Zs[]; };
layout(std430, binding = 4) readonly buffer BL   { uint  Ls[]; };
layout(std430, binding = 5) writeonly buffer BI  { uint  img[]; };

const uint NI = %uu, NJ = %uu, NK = %uu;
const uint ROOT = %uu, N_CONST = %uu, N_X = %uu, N_Y = %uu, N_Z = %uu;
const uint N_CLAUSES = %uu;
const uint CONST_OFF = %uu, AX_OFF = %uu, CLAUSE_OFF = %uu;

float regs[%u];

/*  Exact mirrors of math_f min_f/max_f (see math_f.h): ties
 *  including +-0 return B, a NaN loses to a non-NaN, NaN-vs-NaN
 *  returns A.  GLSL min()/max() leave NaN behavior undefined, so
 *  they are not used.  */
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

float eval_clauses()
{
    for (uint c = 0u; c < N_CLAUSES; ++c)
    {
        uint base = CLAUSE_OFF + 5u * c;
        uint op = blob[base];
        uint o  = blob[base + 1u];
        float A = regs[blob[base + 2u]];
        float B = regs[blob[base + 3u]];
        float R = 0.0;
        switch (int(op))
        {
            /*  Bodies mirror math_f.h exactly (clamps included)  */
            case %d: R = A + B; break;                        // ADD
            case %d: R = A - B; break;                        // SUB
            case %d: R = A * B; break;                        // MUL
            case %d: R = A / B; break;                        // DIV
            case %d: R = min_f(A, B); break;                  // MIN
            case %d: R = max_f(A, B); break;                  // MAX
            case %d: R = pow(A, B); break;                    // POW
            case %d: R = (B == 0.0) ? 0.0                     // MOD
                                    : A - B * floor(A / B); break;
            case %d: R = abs(A); break;                       // ABS
            case %d: R = A * A; break;                        // SQUARE
            case %d: R = (A < 0.0) ? 0.0 : sqrt(A); break;    // SQRT
            case %d: R = sin(A); break;                       // SIN
            case %d: R = cos(A); break;                       // COS
            case %d: R = tan(A); break;                       // TAN
            case %d: R = (A < -1.0) ? -1.5707963705062866     // ASIN
                       : (A > 1.0)  ?  1.5707963705062866
                       : asin(A); break;
            case %d: R = (A < -1.0) ? 3.1415927410125732      // ACOS
                       : (A > 1.0)  ? 0.0
                       : acos(A); break;
            case %d: R = atan(A); break;                      // ATAN
            case %d: R = atan(A, B); break;                   // ATAN2
            case %d: R = -A; break;                           // NEG
            case %d: R = exp(A); break;                       // EXP
            case %d: R = floor(A); break;                     // FLOOR
            case %d: R = log((A > 1.17549e-38) ? A            // LOG
                                               : 1.17549e-38); break;
            case %d: R = uintBitsToFloat(blob[base + 4u]); break; // CONST
            case %d: R = A; break;                            // COPY
            default: R = 0.0; break;
        }
        regs[o] = R;
    }
    return regs[ROOT];
}

void main()
{
    uint i = gl_GlobalInvocationID.x;
    uint j = gl_GlobalInvocationID.y;
    if (i >= NI || j >= NJ)
        return;

    for (uint q = 0u; q < N_CONST; ++q)
        regs[blob[CONST_OFF + 2u*q]] =
                uintBitsToFloat(blob[CONST_OFF + 2u*q + 1u]);

    for (uint q = 0u; q < N_X; ++q)
        regs[blob[AX_OFF + q]] = Xs[i];
    for (uint q = 0u; q < N_Y; ++q)
        regs[blob[AX_OFF + N_X + q]] = Ys[j];

    uint result = 0u;
    for (int k = int(NK) - 1; k >= 0; --k)
    {
        for (uint q = 0u; q < N_Z; ++q)
            regs[blob[AX_OFF + N_X + N_Y + q]] = Zs[k];
        if (eval_clauses() < 0.0)
        {
            result = Ls[k + 1];
            break;
        }
    }
    img[j * NI + i] = result;
}
)GLSL",
        ni, nj, nk,
        bi.root, bi.n_const, bi.n_x, bi.n_y, bi.n_z,
        bi.n_clauses, bi.const_off, bi.ax_off, bi.clause_off,
        bi.num_slots,
        OP_ADD, OP_SUB, OP_MUL, OP_DIV, OP_MIN, OP_MAX, OP_POW,
        OP_MOD, OP_ABS, OP_SQUARE, OP_SQRT, OP_SIN, OP_COS, OP_TAN,
        OP_ASIN, OP_ACOS, OP_ATAN, OP_ATAN2, OP_NEG, OP_EXP,
        OP_FLOOR, OP_LOG, OP_CONST, OP_COPY);
    return buf;
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
        char log[4096];
        glGetShaderInfoLog(sh, sizeof(log), nullptr, log);
        FAIL("compute shader compile failed:\n" << log);
    }
    const GLuint prog = glCreateProgram();
    glAttachShader(prog, sh);
    glLinkProgram(prog);
    glGetProgramiv(prog, GL_LINK_STATUS, &ok);
    if (!ok)
    {
        char log[4096];
        glGetProgramInfoLog(prog, sizeof(log), nullptr, log);
        FAIL("compute shader link failed:\n" << log);
    }
    glDeleteShader(sh);
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

/*  GPU render of the whole region; returns the heightmap row-major.
 *  With dispatch_ms, reports the dispatch+readback time separately
 *  from shader compilation (a per-deck cost, not per-frame).  */
std::vector<uint32_t> render_gpu(const Deck* deck, const Region& r,
                                 double* dispatch_ms = nullptr)
{
    std::vector<uint32_t> blob(
            tape_export_blob(deck, deck_base(deck), nullptr, 0));
    REQUIRE(!blob.empty());
    tape_export_blob(deck, deck_base(deck), blob.data(),
                     uint32_t(blob.size()));
    const BlobInfo bi = parse_blob(blob);

    const GLuint prog = compile_program(make_shader(bi, r.ni, r.nj, r.nk));
    glUseProgram(prog);
    const auto t_dispatch = std::chrono::steady_clock::now();

    std::vector<uint32_t> L(r.nk + 1);
    for (uint32_t k = 0; k <= r.nk; ++k)
        L[k] = r.L[k];
    std::vector<uint32_t> img(size_t(r.ni) * r.nj, 0);

    GLuint bufs[6] = {
        make_ssbo(0, blob.data(), blob.size() * 4),
        make_ssbo(1, r.X, r.ni * 4),
        make_ssbo(2, r.Y, r.nj * 4),
        make_ssbo(3, r.Z, r.nk * 4),
        make_ssbo(4, L.data(), L.size() * 4),
        make_ssbo(5, img.data(), img.size() * 4),
    };

    glDispatchCompute((r.ni + 7) / 8, (r.nj + 7) / 8, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glFinish();

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufs[5]);
    glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0,
                       GLsizeiptr(img.size() * 4), img.data());
    if (dispatch_ms)
        *dispatch_ms = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - t_dispatch).count();

    glDeleteBuffers(6, bufs);
    glDeleteProgram(prog);
    return img;
}

////////////////////////////////////////////////////////////////////////////
//  Stage B: the MPR pipeline on device, one subdivision level.
//
//  Pass 1 (per tile):  interval-evaluate the tape over the tile's
//      xy-frustum (full z), recording per-clause min/max choices with
//      the CPU's exact conditions (tape.cpp tape_eval_i): choice only
//      when neither operand is taint-flagged, <=/>= comparisons.
//      Hot interval ops are `precise` mirrors of math_i; anything
//      exotic returns [-inf,inf] + taint - CONSERVATIVE IS SOUND:
//      the GPU's choices need not match the CPU's, because any sound
//      choice set leaves values on the shortened tape exact, and a
//      pixel is a pure function of those values' signs.
//  Pass 2 (per tile):  tape simplification - the backward liveness
//      walk of tape_push (STANDARD mode: DEAD / ELIDE / COPY / KEEP),
//      emitting a compacted per-tile tape.
//  Pass 3 (per pixel): the stage-A brute-force ray march, but on the
//      tile's SHORTENED tape; EMPTY/FULL tiles fill without marching.

constexpr uint32_t GPU_TILE = 16;   // px per tile edge

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

/*  Pass 1: interval classify + record choices.  One thread per tile. */
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

/*  precise: interval endpoints must be plain rounded fp32 - a fused
 *  or reassociated bound could exclude a value pass 3 computes.  */
precise vec2 iv[NSLOTS];
bool tnt[NSLOTS];

vec2 iv_mul(vec2 A, vec2 B)
{
    precise float c1 = A.x * B.x, c2 = A.x * B.y,
                  c3 = A.y * B.x, c4 = A.y * B.y;
    float lo = min(min(c1, c2), min(c3, c4));
    float hi = max(max(c1, c2), max(c3, c4));
    if (isnan(c1) || isnan(c2) || isnan(c3) || isnan(c4))
        { lo = -INF; hi = INF; }   // conservative; taint covers it
    return vec2(lo, hi);
}

void main()
{
    uint t = gl_GlobalInvocationID.x;
    if (t >= TILES_X * TILES_Y)
        return;
    uint ti = t %% TILES_X, tj = t / TILES_X;

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
        uint choice = 0u;   // CHOICE_BOTH

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
                if (A.y <= B.x)      choice = 1u;   // CHOICE_A
                else if (B.y <= A.x) choice = 2u;   // CHOICE_B
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
        else                   { R = vec2(-INF, INF); u = true; }  // POW/MOD/LOG/TAN/...

        u = u || isnan(R.x) || isnan(R.y);
        iv[o] = R;
        tnt[o] = u;
        choices[t * N_CLAUSES + c] = choice;
    }

    vec2 root = iv[ROOT];
    bool rt = tnt[ROOT];
    uint cls = 2u;                       // AMBIG
    if (!rt && root.x > 0.0)  cls = 0u;  // EMPTY: whole frustum outside
    else if (!rt && root.y < 0.0) cls = 1u;  // FULL: topmost voxel inside
    tileinfo[t] = uvec2(cls, N_CLAUSES);
}
)GLSL";
    // GLSL has no %% escape outside printf-formatting; undo the one
    // we wrote for the C++ string (kept literal here).
    const size_t pct = src.find("%%");
    if (pct != std::string::npos)
        src.replace(pct, 2, "%");
    return src;
}

/*  Pass 2: per-tile tape simplification (tape_push STANDARD).  */
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
uint verdict[N_CLAUSES];   // 0 dead/elide, 1 keep, 2 copy-a, 3 copy-b

void main()
{
    uint t = gl_GlobalInvocationID.x;
    if (t >= TILES_X * TILES_Y)
        return;
    if (tileinfo[t].x != 2u)   // only ambiguous tiles carry a tape
        return;

    for (uint s = 0u; s < NSLOTS; ++s)
        live[s] = false;
    live[ROOT] = true;

    /*  Backward liveness, mirroring tape_push (STANDARD): a decided
     *  min/max becomes a register copy of the surviving input (or
     *  vanishes when the allocator already placed it in the output
     *  register); a dead output drops the clause entirely.  */
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
        /*  Load everything into locals BEFORE storing: interleaved
         *  SSBO reads and writes in this loop miscompiled on Mesa
         *  Intel (the copy's out-slot word read back 0), while the
         *  same source was correct on NVIDIA.  Locals-then-store
         *  sidesteps whatever aliasing analysis goes wrong.  */
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

/*  Pass 3: per-pixel ray march on the tile's shortened tape.  */
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
        "layout(std430, binding=6) readonly buffer BC { uint choices[]; };\n"
        "layout(std430, binding=7) readonly buffer BT { uvec2 tileinfo[]; };\n"
        "layout(std430, binding=8) readonly buffer BS { uint tapes[]; };\n";
    src += op_prelude();
    src += blob_prelude(bi, ni, nj, nk);
    src += R"GLSL(
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
    uint j = gl_GlobalInvocationID.y;
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

/*  Full pipeline render.  */
std::vector<uint32_t> render_gpu_tiles(const Deck* deck, const Region& r,
                                       double* dispatch_ms = nullptr)
{
    std::vector<uint32_t> blob(
            tape_export_blob(deck, deck_base(deck), nullptr, 0));
    REQUIRE(!blob.empty());
    tape_export_blob(deck, deck_base(deck), blob.data(),
                     uint32_t(blob.size()));
    const BlobInfo bi = parse_blob(blob);

    const uint32_t tiles_x = (r.ni + GPU_TILE - 1) / GPU_TILE;
    const uint32_t tiles_y = (r.nj + GPU_TILE - 1) / GPU_TILE;
    const uint32_t n_tiles = tiles_x * tiles_y;

    const GLuint p1 = compile_program(make_p1(bi, r.ni, r.nj, r.nk));
    const GLuint p2 = compile_program(make_p2(bi, r.ni, r.nj, r.nk));
    const GLuint p3 = compile_program(make_p3(bi, r.ni, r.nj, r.nk));

    std::vector<uint32_t> L(r.nk + 1);
    for (uint32_t k = 0; k <= r.nk; ++k)
        L[k] = r.L[k];
    std::vector<uint32_t> img(size_t(r.ni) * r.nj, 0);

    const auto t_dispatch = std::chrono::steady_clock::now();
    GLuint bufs[9] = {
        make_ssbo(0, blob.data(), blob.size() * 4),
        make_ssbo(1, r.X, (r.ni + 1) * 4),
        make_ssbo(2, r.Y, (r.nj + 1) * 4),
        make_ssbo(3, r.Z, (r.nk + 1) * 4),
        make_ssbo(4, L.data(), L.size() * 4),
        make_ssbo(5, img.data(), img.size() * 4),
        make_ssbo(6, nullptr, size_t(n_tiles) * bi.n_clauses * 4),
        make_ssbo(7, nullptr, size_t(n_tiles) * 8),
        make_ssbo(8, nullptr, size_t(n_tiles) * bi.n_clauses * 5 * 4),
    };

    glUseProgram(p1);
    glDispatchCompute((n_tiles + 63) / 64, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glUseProgram(p2);
    glDispatchCompute((n_tiles + 63) / 64, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glUseProgram(p3);
    glDispatchCompute((r.ni + 7) / 8, (r.nj + 7) / 8, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glFinish();

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufs[5]);
    glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0,
                       GLsizeiptr(img.size() * 4), img.data());
    if (dispatch_ms)
        *dispatch_ms = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - t_dispatch).count();

    glDeleteBuffers(9, bufs);
    glDeleteProgram(p1);
    glDeleteProgram(p2);
    glDeleteProgram(p3);
    return img;
}

/*  CPU render via the real renderer (single-threaded, full tape
 *  machinery with BINARY pushes - the production path).  */
std::vector<uint16_t> render_cpu(const Deck* deck, const Region& r)
{
    std::vector<uint16_t> store(size_t(r.ni) * r.nj, 0);
    std::vector<uint16_t*> rows(r.nj);
    for (uint32_t j = 0; j < r.nj; ++j)
        rows[j] = store.data() + size_t(j) * r.ni;

    TapeCtx* ctx = tape_ctx_new(deck);
    volatile int halt = 0;
    render16_tape(deck_base(deck), ctx, r, rows.data(), &halt, nullptr);
    tape_ctx_free(ctx);
    return store;
}

/*  Renders one model both ways and returns the mismatched-pixel count.
 *  With tiles=true, the GPU side runs the full three-pass pipeline.  */
uint64_t referee(const char* expr, uint32_t n, bool tiles = false)
{
    MathTree* tree = parse(expr);
    REQUIRE(tree != nullptr);
    Deck* deck = deck_from_tree(tree);
    free_tree(tree);

    Region r = {};
    r.ni = n;
    r.nj = n;
    r.nk = n;
    r.voxels = uint64_t(n) * n * n;
    build_arrays(&r, -1.1f, -1.1f, -1.1f, 1.1f, 1.1f, 1.1f);

    const auto g = tiles ? render_gpu_tiles(deck, r)
                         : render_gpu(deck, r);
    const auto c = render_cpu(deck, r);

    uint64_t diff = 0;
    for (size_t p = 0; p < g.size(); ++p)
        if (g[p] != c[p])
            ++diff;

    free_arrays(&r);
    deck_free(deck);
    return diff;
}

}  // namespace

////////////////////////////////////////////////////////////////////////////

TEST_CASE("GPU vs CPU render bench", "[.gpubench]")
{
    if (!gpu().ok)
    {
        WARN("no EGL / GL 4.3 stack available; skipping GPU bench");
        return;
    }
    WARN("GPU renderer: " << gpu().renderer);

    /*  Brute force GPU (full tape, every voxel, no pruning) vs the
     *  CPU carrying the whole tape arsenal, single thread.  N via
     *  STIBIUM_GPU_BENCH_N (default 256).  */
    uint32_t n = 256;
    if (const char* env = getenv("STIBIUM_GPU_BENCH_N"))
        n = uint32_t(atoi(env));

    const char* MODEL =
            "ai-r++qXqYqZf1-r++q-Xf0.5qYqZf0.7n-Zf0.2";  // nested CSG

    MathTree* tree = parse(MODEL);
    REQUIRE(tree != nullptr);
    Deck* deck = deck_from_tree(tree);
    free_tree(tree);

    Region r = {};
    r.ni = n;  r.nj = n;  r.nk = n;
    r.voxels = uint64_t(n) * n * n;
    build_arrays(&r, -1.1f, -1.1f, -1.1f, 1.1f, 1.1f, 1.1f);

    using clk = std::chrono::steady_clock;
    double dispatch = 0;
    const auto t0 = clk::now();
    const auto g = render_gpu(deck, r, &dispatch);
    const auto t1 = clk::now();
    const auto c = render_cpu(deck, r);
    const auto t2 = clk::now();

    uint64_t diff = 0;
    for (size_t p = 0; p < g.size(); ++p)
        if (g[p] != c[p])
            ++diff;

    const auto ms = [](auto a, auto b) {
        return std::chrono::duration<double, std::milli>(b - a).count();
    };
    WARN("n=" << n << ": GPU total " << ms(t0, t1)
              << " ms (dispatch+readback " << dispatch
              << " ms), CPU 1-thread " << ms(t1, t2)
              << " ms, mismatches " << diff);
    CHECK(diff == 0);

    free_arrays(&r);
    deck_free(deck);
}

TEST_CASE("GPU heightmap matches CPU renderer", "[.gpu]")
{
    if (!gpu().ok)
    {
        WARN("no EGL / GL 4.3 stack available; skipping GPU tests");
        return;
    }
    WARN("GPU renderer: " << gpu().renderer);

    /*  Algebra + min/max only: IEEE ops match bit-for-bit, so the
     *  heightmap must be IDENTICAL.  */
    SECTION("exact: algebraic CSG")
    {
        const char* MODELS[] = {
            "*+Xf2.0+Yf2.0",                    // dedup'd arithmetic sheet
            "aa-f-0.6X-Xf0.6aa-f-0.6Y-Yf0.6a-f-0.6Z-Zf0.6",  // cube
            "-+*XX*YYf0.5",                     // paraboloid (squares)
        };
        for (const char* m : MODELS)
        {
            CAPTURE(m);
            CHECK(referee(m, 128) == 0);
        }
    }

    /*  sqrt / trig are not guaranteed correctly-rounded in GLSL;
     *  mismatches sit on the surface (a voxel whose field is within
     *  an ulp of zero lands on the other side).  Report the count -
     *  bounded, not required zero.  */
    SECTION("tolerance: sqrt and trig surfaces")
    {
        const char* MODELS[] = {
            "-r++qXqYqZf1",                                  // sphere
            "i-r++qXqYqZf1-r++q-Xf0.5q-Yf0.25q-Zf0.1f0.8",   // sphere union
            "ai-r++qXqYqZf1-r++q-Xf0.5qYqZf0.7n-Zf0.2",      // nested min/max
            "+s*Xf3.0c*Yf2.0",                               // trig sheet
        };
        for (const char* m : MODELS)
        {
            const uint64_t diff = referee(m, 128);
            const double frac = double(diff) / (128.0 * 128.0);
            CAPTURE(m);
            CAPTURE(diff);
            WARN("model " << m << ": " << diff << " mismatched px ("
                          << frac * 100.0 << "%)");
            CHECK(frac < 0.005);   // < 0.5% of pixels
        }
    }
}

TEST_CASE("GPU tile pipeline matches CPU renderer", "[.gpu2]")
{
    if (!gpu().ok)
    {
        WARN("no EGL / GL 4.3 stack available; skipping GPU tests");
        return;
    }
    WARN("GPU renderer: " << gpu().renderer);

    const char* MODELS[] = {
        "*+Xf2.0+Yf2.0",
        "aa-f-0.6X-Xf0.6aa-f-0.6Y-Yf0.6a-f-0.6Z-Zf0.6",
        "-+*XX*YYf0.5",
        "-r++qXqYqZf1",
        "i-r++qXqYqZf1-r++q-Xf0.5q-Yf0.25q-Zf0.1f0.8",
        "ai-r++qXqYqZf1-r++q-Xf0.5qYqZf0.7n-Zf0.2",
        "+s*Xf3.0c*Yf2.0",
        "i-r++qXqYqZf1-r++q-Xf10qYqZf0.5",   // PRUNABLE: decided min
    };
    for (const char* m : MODELS)
    {
        CAPTURE(m);
        CHECK(referee(m, 128, true) == 0);
    }
}

TEST_CASE("GPU tile pipeline bench", "[.gpu2bench]")
{
    if (!gpu().ok)
    {
        WARN("no EGL / GL 4.3 stack available; skipping GPU bench");
        return;
    }
    WARN("GPU renderer: " << gpu().renderer);

    uint32_t n = 256;
    if (const char* env = getenv("STIBIUM_GPU_BENCH_N"))
        n = uint32_t(atoi(env));

    const char* MODEL =
            "ai-r++qXqYqZf1-r++q-Xf0.5qYqZf0.7n-Zf0.2";

    MathTree* tree = parse(MODEL);
    REQUIRE(tree != nullptr);
    Deck* deck = deck_from_tree(tree);
    free_tree(tree);

    Region r = {};
    r.ni = n;  r.nj = n;  r.nk = n;
    r.voxels = uint64_t(n) * n * n;
    build_arrays(&r, -1.1f, -1.1f, -1.1f, 1.1f, 1.1f, 1.1f);

    double d_brute = 0, d_tiles = 0;
    const auto gb = render_gpu(deck, r, &d_brute);
    const auto gt = render_gpu_tiles(deck, r, &d_tiles);

    using clk = std::chrono::steady_clock;
    const auto t0 = clk::now();
    const auto c = render_cpu(deck, r);
    const auto t1 = clk::now();
    const double d_cpu =
            std::chrono::duration<double, std::milli>(t1 - t0).count();

    uint64_t diff = 0;
    for (size_t p = 0; p < gt.size(); ++p)
        if (gt[p] != c[p])
            ++diff;

    WARN("n=" << n << ": GPU brute " << d_brute
              << " ms, GPU tile pipeline " << d_tiles
              << " ms, CPU 1-thread " << d_cpu
              << " ms, pipeline mismatches " << diff);
    CHECK(diff == 0);
    (void)gb;

    free_arrays(&r);
    deck_free(deck);
}

TEST_CASE("GPU tile pipeline debug", "[.gpudbg]")
{
    if (!gpu().ok)
        return;
    const char* MODEL = getenv("STIBIUM_GPU_DBG_MODEL");
    if (!MODEL)
        MODEL = "i-r++qXqYqZf1-r++q-Xf10qYqZf0.5";

    MathTree* tree = parse(MODEL);
    REQUIRE(tree != nullptr);
    Deck* deck = deck_from_tree(tree);
    free_tree(tree);

    Region r = {};
    r.ni = 128;  r.nj = 128;  r.nk = 128;
    r.voxels = uint64_t(128) * 128 * 128;
    build_arrays(&r, -1.1f, -1.1f, -1.1f, 1.1f, 1.1f, 1.1f);

    const auto gb = render_gpu(deck, r);
    const auto gt = render_gpu_tiles(deck, r);

    int shown = 0;
    uint64_t diff = 0;
    for (uint32_t j = 0; j < r.nj; ++j)
        for (uint32_t i = 0; i < r.ni; ++i)
        {
            const size_t p = size_t(j) * r.ni + i;
            if (gb[p] != gt[p])
            {
                ++diff;
                if (shown < 12)
                {
                    fprintf(stderr,
                            "px (%u,%u) tile (%u,%u): brute=%u tiles=%u\n",
                            i, j, i / GPU_TILE, j / GPU_TILE,
                            gb[p], gt[p]);
                    ++shown;
                }
            }
        }
    fprintf(stderr, "total brute-vs-tiles diffs: %llu\n",
            (unsigned long long)diff);

    free_arrays(&r);
    deck_free(deck);
}

TEST_CASE("GPU tile pipeline tape dump", "[.gpudump]")
{
    if (!gpu().ok)
        return;
    const char* MODEL = "i-r++qXqYqZf1-r++q-Xf10qYqZf0.5";

    MathTree* tree = parse(MODEL);
    REQUIRE(tree != nullptr);
    Deck* deck = deck_from_tree(tree);
    free_tree(tree);

    Region r = {};
    r.ni = 128;  r.nj = 128;  r.nk = 128;
    r.voxels = uint64_t(128) * 128 * 128;
    build_arrays(&r, -1.1f, -1.1f, -1.1f, 1.1f, 1.1f, 1.1f);

    std::vector<uint32_t> blob(
            tape_export_blob(deck, deck_base(deck), nullptr, 0));
    tape_export_blob(deck, deck_base(deck), blob.data(),
                     uint32_t(blob.size()));
    const BlobInfo bi = parse_blob(blob);
    fprintf(stderr, "deck: slots=%u root=%u clauses=%u\n",
            bi.num_slots, bi.root, bi.n_clauses);
    for (uint32_t c = 0; c < bi.n_clauses; ++c)
    {
        const uint32_t* w = blob.data() + bi.clause_off + 5*c;
        fprintf(stderr, "  [%2u] op=%-2u out=%u a=%u b=%u\n",
                c, w[0], w[1], w[2], w[3]);
    }

    const uint32_t tiles_x = (r.ni + GPU_TILE - 1) / GPU_TILE;
    const uint32_t tiles_y = (r.nj + GPU_TILE - 1) / GPU_TILE;
    const uint32_t n_tiles = tiles_x * tiles_y;

    const GLuint p1 = compile_program(make_p1(bi, r.ni, r.nj, r.nk));
    const GLuint p2 = compile_program(make_p2(bi, r.ni, r.nj, r.nk));

    std::vector<uint32_t> Lv(r.nk + 1);
    for (uint32_t k = 0; k <= r.nk; ++k)
        Lv[k] = r.L[k];
    std::vector<uint32_t> img(size_t(r.ni) * r.nj, 0);

    GLuint bufs[9] = {
        make_ssbo(0, blob.data(), blob.size() * 4),
        make_ssbo(1, r.X, (r.ni + 1) * 4),
        make_ssbo(2, r.Y, (r.nj + 1) * 4),
        make_ssbo(3, r.Z, (r.nk + 1) * 4),
        make_ssbo(4, Lv.data(), Lv.size() * 4),
        make_ssbo(5, img.data(), img.size() * 4),
        make_ssbo(6, nullptr, size_t(n_tiles) * bi.n_clauses * 4),
        make_ssbo(7, nullptr, size_t(n_tiles) * 8),
        make_ssbo(8, nullptr, size_t(n_tiles) * bi.n_clauses * 5 * 4),
    };

    glUseProgram(p1);
    glDispatchCompute((n_tiles + 63) / 64, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glUseProgram(p2);
    glDispatchCompute((n_tiles + 63) / 64, 1, 1);
    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    glFinish();

    std::vector<uint32_t> ti(n_tiles * 2);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufs[7]);
    glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, ti.size() * 4,
                       ti.data());
    std::vector<uint32_t> ch(size_t(n_tiles) * bi.n_clauses);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufs[6]);
    glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, ch.size() * 4,
                       ch.data());
    std::vector<uint32_t> tp(size_t(n_tiles) * bi.n_clauses * 5);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, bufs[8]);
    glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, tp.size() * 4,
                       tp.data());

    const uint32_t t = 0 * tiles_x + 3;   // tile (3,0)
    fprintf(stderr, "tile (3,0): cls=%u len=%u\n", ti[2*t], ti[2*t+1]);
    fprintf(stderr, "choices:");
    for (uint32_t c = 0; c < bi.n_clauses; ++c)
        fprintf(stderr, " %u", ch[size_t(t) * bi.n_clauses + c]);
    fprintf(stderr, "\ntile tape:\n");
    for (uint32_t c = 0; c < ti[2*t+1]; ++c)
    {
        const uint32_t* w = tp.data() + (size_t(t) * bi.n_clauses + c) * 5;
        fprintf(stderr, "  [%2u] op=%-2u out=%u a=%u b=%u\n",
                c, w[0], w[1], w[2], w[3]);
    }

    glDeleteBuffers(9, bufs);
    glDeleteProgram(p1);
    glDeleteProgram(p2);
    free_arrays(&r);
    deck_free(deck);
}

TEST_CASE("GPU division exactness", "[.gpudiv]")
{
    if (!gpu().ok)
        return;
    WARN("GPU renderer: " << gpu().renderer);
    const char* MODELS[] = {
        "-r++q/Xf1.7q/Yf1.3qZf1",       // ellipsoid via division
        "-/+qX*YYf3.1f0.2",             // rational field
        "i-r++q/Xf2q/Yf2qZf1-/Xf9f0.1", // div + decided min
    };
    for (const char* m : MODELS)
    {
        const uint64_t diff = referee(m, 128, false);
        CAPTURE(m);
        WARN("model " << m << ": " << diff << " mismatched px");
        CHECK(diff == 0);
    }
}

TEST_CASE("GPU big-deck exactness", "[.gpubig]")
{
    if (!gpu().ok)
        return;
    WARN("GPU renderer: " << gpu().renderer);

    /*  min-of-N-spheres with division transforms: linearly growing
     *  slot count, exact ops only.  Bisects "deck size breaks the
     *  GPU" from model-specific causes.  */
    for (int n_spheres : {176, 208, 256, 288})
    {
        std::string expr;
        for (int s = 0; s < n_spheres - 1; ++s)
            expr += "i";
        for (int s = 0; s < n_spheres; ++s)
        {
            char sph[160];
            const float cx = -0.9f + 1.8f * float(s) / n_spheres;
            const float cy = -0.7f + 1.4f * float((s * 7) % 13) / 13.0f;
            const float rr = 1.05f + 0.9f * float((s * 5) % 11) / 11.0f;
            snprintf(sph, sizeof(sph),
                     "-r++q-/Xf%.3ff%.3fq-/Yf%.3ff%.3fqZf%.3f",
                     rr, cx, rr, cy, 0.18f);
            expr += sph;
        }
        MathTree* tree = parse(expr.c_str());
        REQUIRE(tree != nullptr);
        Deck* deck = deck_from_tree(tree);
        free_tree(tree);
        const uint32_t slots = deck_slots(deck);
        deck_free(deck);

        const uint64_t diff = referee(expr.c_str(), 96, false);
        WARN("spheres=" << n_spheres << " slots=" << slots
                        << ": " << diff << " mismatched px");
        CHECK(diff == 0);
    }
}

/*  The PRODUCTION path (fab/tree/gpu.h - rescheduled registers,
 *  constants in shader source, z-slab tapes, band recycling), not
 *  the historical prototypes above.  This is the referee that has
 *  to stay green as the pipeline evolves.  */
TEST_CASE("Production gpu_render16 matches CPU", "[.gpulib]")
{
    if (!gpu_available())
    {
        WARN("no GPU available; skipping");
        return;
    }
    WARN("GPU renderer: " << gpu_renderer_name());

    const char* MODELS[] = {
        "*+Xf2.0+Yf2.0",
        "aa-f-0.6X-Xf0.6aa-f-0.6Y-Yf0.6a-f-0.6Z-Zf0.6",
        "-+*XX*YYf0.5",
        "-r++qXqYqZf1",
        "i-r++qXqYqZf1-r++q-Xf0.5q-Yf0.25q-Zf0.1f0.8",
        "ai-r++qXqYqZf1-r++q-Xf0.5qYqZf0.7n-Zf0.2",
        "+s*Xf3.0c*Yf2.0",
        "i-r++qXqYqZf1-r++q-Xf10qYqZf0.5",
        "-r++q/Xf1.7q/Yf1.3qZf1",
        "i-r++q/Xf2q/Yf2qZf1-/Xf9f0.1",
    };
    for (const char* m : MODELS)
    {
        MathTree* tree = parse(m);
        REQUIRE(tree != nullptr);
        Deck* deck = deck_from_tree(tree);
        free_tree(tree);

        Region r = {};
        r.ni = 160;  r.nj = 96;  r.nk = 128;   // non-square on purpose
        r.voxels = uint64_t(r.ni) * r.nj * r.nk;
        build_arrays(&r, -1.1f, -1.1f, -1.1f, 1.1f, 1.1f, 1.1f);

        std::vector<uint16_t> gstore(size_t(r.ni) * r.nj, 0);
        std::vector<uint16_t*> grows(r.nj);
        for (uint32_t j = 0; j < r.nj; ++j)
            grows[j] = gstore.data() + size_t(j) * r.ni;
        volatile int halt = 0;
        const bool ran = gpu_render16(deck, r, grows.data(), &halt);
        REQUIRE(ran);

        const auto c = render_cpu(deck, r);
        uint64_t diff = 0;
        for (size_t p = 0; p < gstore.size(); ++p)
            if (gstore[p] != c[p])
                ++diff;
        CAPTURE(m);
        CHECK(diff == 0);

        free_arrays(&r);
        deck_free(deck);
    }
}

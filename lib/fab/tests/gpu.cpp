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

/*  Renders one model both ways and returns the mismatched-pixel count. */
uint64_t referee(const char* expr, uint32_t n)
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

    const auto g = render_gpu(deck, r);
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

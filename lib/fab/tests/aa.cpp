/*
 *  Affine-arithmetic prototype (measurement phase - see
 *  doc/TAPE-NEXT.md and doc/research/2026-07-14-frep-sota.md §4).
 *
 *  Question to answer BEFORE productionizing: how much shorter do
 *  classification-driven tapes get when bounds keep correlation
 *  through transform chains?  Interval arithmetic loses x - x = 0;
 *  affine forms (Comba/Stolfi; reduced AF1 with one noise symbol
 *  per axis + one error term) keep it.
 *
 *  Design under test: every value carries BOTH arithmetics -
 *  an affine form c + kx*ex + ky*ey + kz*ez + e*[-1,1] AND the
 *  plain math_i interval; usable bounds are their intersection, so
 *  the result is never worse than today's IA, and the interval side
 *  absorbs infinities/NaN.  Ops without a hand-written affine rule
 *  fall back to an interval-only form - always sound, loses only
 *  tightness.  A decided min/max passes the WINNER's affine form
 *  through verbatim (correlation survives CSG joins).
 *
 *  Soundness invariant (spot-checked here, same standard as the
 *  pruning fuzzer): every pointwise value in the box lies inside
 *  the intersected bounds, and a push driven by AA choices leaves
 *  pointwise evaluation bit-identical.
 *
 *  Tags: "[.aa]" measurement on a dumped app blob
 *        (STIBIUM_AA_BLOB=path, default build/zeiss.blob);
 *        "[.aafuzz]" corpus soundness sweep.
 */

#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <random>
#include <string>
#include <vector>

#include "catch/catch.hpp"

#include "fab/tree/tape.h"
#include "fab/tree/tree.h"
#include "fab/tree/parser.h"
#include "fab/tree/node/opcodes.h"
#include "fab/tree/math/math_f.h"
#include "fab/tree/math/math_i.h"

namespace {

////////////////////////////////////////////////////////////////////////////
//  Blob-backed program

struct Clause5
{
    uint32_t op, out, a, b;
    float imm;
};

struct Prog
{
    uint32_t num_slots = 0, root = 0;
    uint32_t n_const = 0, n_x = 0, n_y = 0, n_z = 0;
    std::vector<std::pair<uint32_t, float>> consts;
    std::vector<uint32_t> xs, ys, zs;
    std::vector<Clause5> clauses;
};

float bits_to_f(uint32_t u)
{
    float f;
    memcpy(&f, &u, 4);
    return f;
}

Prog prog_from_blob(const std::vector<uint32_t>& b)
{
    Prog p;
    REQUIRE(b[0] == 1);
    p.num_slots = b[1];
    p.root = b[2];
    p.n_const = b[3];
    p.n_x = b[4];
    p.n_y = b[5];
    p.n_z = b[6];
    const uint32_t n_clauses = b[7];
    uint32_t w = 8;
    for (uint32_t q = 0; q < p.n_const; ++q, w += 2)
        p.consts.push_back({ b[w], bits_to_f(b[w+1]) });
    for (uint32_t q = 0; q < p.n_x; ++q)
        p.xs.push_back(b[w++]);
    for (uint32_t q = 0; q < p.n_y; ++q)
        p.ys.push_back(b[w++]);
    for (uint32_t q = 0; q < p.n_z; ++q)
        p.zs.push_back(b[w++]);
    for (uint32_t c = 0; c < n_clauses; ++c, w += 5)
        p.clauses.push_back({ b[w], b[w+1], b[w+2], b[w+3],
                              bits_to_f(b[w+4]) });
    return p;
}

Prog prog_from_tree(const char* expr)
{
    MathTree* tree = parse(expr);
    REQUIRE(tree != nullptr);
    Deck* deck = deck_from_tree(tree);
    free_tree(tree);
    std::vector<uint32_t> blob(
            tape_export_blob(deck, deck_base(deck), nullptr, 0));
    REQUIRE(!blob.empty());
    tape_export_blob(deck, deck_base(deck), blob.data(),
                     uint32_t(blob.size()));
    deck_free(deck);
    return prog_from_blob(blob);
}

////////////////////////////////////////////////////////////////////////////
//  Pointwise evaluation (mirrors math_f exactly - the ground truth)

float eval_point(const Prog& p, float x, float y, float z,
                 std::vector<float>& slots)
{
    slots.resize(p.num_slots);
    for (const auto& [s, v] : p.consts)
        slots[s] = v;
    for (uint32_t s : p.xs)   slots[s] = x;
    for (uint32_t s : p.ys)   slots[s] = y;
    for (uint32_t s : p.zs)   slots[s] = z;
    for (const Clause5& c : p.clauses)
    {
        const float A = slots[c.a], B = slots[c.b];
        float R = 0;
        switch (c.op)
        {
            case OP_ADD:    R = add_f(A, B); break;
            case OP_SUB:    R = sub_f(A, B); break;
            case OP_MUL:    R = mul_f(A, B); break;
            case OP_DIV:    R = div_f(A, B); break;
            case OP_MIN:    R = min_f(A, B); break;
            case OP_MAX:    R = max_f(A, B); break;
            case OP_POW:    R = pow_f(A, B); break;
            case OP_ATAN2:  R = atan2_f(A, B); break;
            case OP_MOD:    R = mod_f(A, B); break;
            case OP_ABS:    R = abs_f(A); break;
            case OP_SQUARE: R = square_f(A); break;
            case OP_SQRT:   R = sqrt_f(A); break;
            case OP_SIN:    R = sin_f(A); break;
            case OP_COS:    R = cos_f(A); break;
            case OP_TAN:    R = tan_f(A); break;
            case OP_ASIN:   R = asin_f(A); break;
            case OP_ACOS:   R = acos_f(A); break;
            case OP_ATAN:   R = atan_f(A); break;
            case OP_NEG:    R = neg_f(A); break;
            case OP_EXP:    R = exp_f(A); break;
            case OP_FLOOR:  R = floor_f(A); break;
            case OP_LOG:    R = log_f(A); break;
            case OP_CONST:  R = c.imm; break;
            case OP_COPY:   R = A; break;
            default:        FAIL("unsupported op in blob");
        }
        slots[c.out] = R;
    }
    return slots[p.root];
}

////////////////////////////////////////////////////////////////////////////
//  Affine form: c + k[0]*ex + k[1]*ey + k[2]*ez + e*[-1,1]

struct AAF
{
    float c = 0;
    float k[3] = { 0, 0, 0 };
    float e = 0;
};

float aaf_rad(const AAF& a)
{
    return fabsf(a.k[0]) + fabsf(a.k[1]) + fabsf(a.k[2]) + a.e;
}

Interval aaf_hull(const AAF& a)
{
    const float r = aaf_rad(a);
    return { a.c - r, a.c + r };
}

bool aaf_finite(const AAF& a)
{
    return std::isfinite(a.c) && std::isfinite(a.k[0]) &&
           std::isfinite(a.k[1]) && std::isfinite(a.k[2]) &&
           std::isfinite(a.e);
}

/*  Defensive fp smear: affine ops below use round-to-nearest, so
 *  every op widens the error term by a few ulps of the magnitude.  */
void aaf_widen(AAF& a)
{
    a.e += (fabsf(a.c) + aaf_rad(a)) * (4 * FLT_EPSILON);
}

AAF aaf_from_iv(Interval v)
{
    AAF a;
    if (!std::isfinite(v.lower) || !std::isfinite(v.upper))
    {
        a.c = 0;
        a.e = INFINITY;
        return a;
    }
    a.c = (v.lower + v.upper) * 0.5f;
    a.e = fmaxf(v.upper - a.c, a.c - v.lower);
    aaf_widen(a);
    return a;
}

AAF aaf_const(float v)
{
    AAF a;
    a.c = v;
    return a;
}

AAF aaf_add(const AAF& x, const AAF& y)
{
    AAF a;
    a.c = x.c + y.c;
    for (int i = 0; i < 3; ++i)
        a.k[i] = x.k[i] + y.k[i];
    a.e = x.e + y.e;
    aaf_widen(a);
    return a;
}

AAF aaf_neg(const AAF& x)
{
    AAF a;
    a.c = -x.c;
    for (int i = 0; i < 3; ++i)
        a.k[i] = -x.k[i];
    a.e = x.e;
    return a;
}

AAF aaf_sub(const AAF& x, const AAF& y)
{
    return aaf_add(x, aaf_neg(y));
}

AAF aaf_scale(const AAF& x, float m)
{
    AAF a;
    a.c = x.c * m;
    for (int i = 0; i < 3; ++i)
        a.k[i] = x.k[i] * m;
    a.e = x.e * fabsf(m);
    aaf_widen(a);
    return a;
}

AAF aaf_mul(const AAF& x, const AAF& y)
{
    AAF a;
    a.c = x.c * y.c;
    for (int i = 0; i < 3; ++i)
        a.k[i] = x.c * y.k[i] + y.c * x.k[i];
    const float rx = aaf_rad(x), ry = aaf_rad(y);
    a.e = fabsf(x.c) * y.e + fabsf(y.c) * x.e + rx * ry;
    aaf_widen(a);
    return a;
}

/*  Chebyshev linear approximation of sqrt over a strictly positive
 *  hull [a, b]: result = alpha*x + zeta +- delta, with the secant
 *  slope alpha = (sqrt(b)-sqrt(a))/(b-a); the deviation sqrt(x) -
 *  alpha*x peaks at x* = 1/(4 alpha^2) and equals the endpoint
 *  deviation at both ends.  Keeps the affine FORM through the sqrt,
 *  so correlations survive into min/max comparisons downstream.
 *  Boxes whose hull touches zero keep the hull fallback (sqrt_f
 *  clamps negatives - piecewise, no affine form).  */
bool aaf_sqrt(const AAF& x, Interval hull, AAF* out)
{
    const float a = hull.lower, b = hull.upper;
    if (!(a > 0) || !std::isfinite(b) || b <= a)
        return false;
    const float sa = sqrtf(a), sb = sqrtf(b);
    const float alpha = (sb - sa) / (b - a);
    if (!(alpha > 0) || !std::isfinite(alpha))
        return false;
    const float dab = sa - alpha * a;          // endpoint deviation
    const float dxs = 1.0f / (4.0f * alpha);   // peak deviation
    const float zeta = (dxs + dab) * 0.5f;
    const float delta = (dxs - dab) * 0.5f;
    if (!(delta >= 0) || !std::isfinite(delta))
        return false;
    AAF r = aaf_scale(x, alpha);
    r.c += zeta;
    r.e += delta;
    aaf_widen(r);
    *out = r;
    return true;
}

////////////////////////////////////////////////////////////////////////////
//  Classification: IA baseline or AA+IA intersection.
//
//  Choice / taint rules mirror tape.cpp's tape_eval_i exactly, run
//  on the (intersected) bounds - tighter bounds mean MORE decided
//  choices and FEWER spurious taints, both sound.

struct Cell
{
    AAF af;
    Interval iv;   // usable bounds (intersection when AA is on)
    bool tnt = false;
};

struct Classified
{
    Interval root;
    bool root_tnt = false;
    std::vector<uint8_t> choice;   // 0 both, 1 A, 2 B (per clause)
};

Interval intersect_iv(Interval a, Interval b)
{
    Interval r = { fmaxf(a.lower, b.lower), fminf(a.upper, b.upper) };
    if (r.lower > r.upper)   // fp disagreement between the two forms
        return a;            // keep the IA side (always sound)
    return r;
}

Classified classify(const Prog& p, Interval X, Interval Y, Interval Z,
                    bool use_aa)
{
    std::vector<Cell> cells(p.num_slots);
    for (const auto& [s, v] : p.consts)
        cells[s] = { aaf_const(v), { v, v }, false };

    const auto set_axis = [&](uint32_t s, Interval v, int axis) {
        AAF a;
        a.c = (v.lower + v.upper) * 0.5f;
        a.k[axis] = v.upper - a.c;
        cells[s] = { a, v, false };
    };
    for (uint32_t s : p.xs)   set_axis(s, X, 0);
    for (uint32_t s : p.ys)   set_axis(s, Y, 1);
    for (uint32_t s : p.zs)   set_axis(s, Z, 2);

    Classified out;
    out.choice.resize(p.clauses.size(), 0);

    for (size_t ci = 0; ci < p.clauses.size(); ++ci)
    {
        const Clause5& c = p.clauses[ci];
        const Cell& A = cells[c.a];
        const Cell& B = cells[c.b];
        const Interval ai = A.iv, bi = B.iv;
        const bool na = A.tnt, nb = B.tnt;

        Interval vi;      // math_i result (always computed)
        AAF af;           // affine result (or interval-only fallback)
        bool have_af = false;
        uint8_t choice = 0;

        switch (c.op)
        {
            case OP_ADD:
                vi = add_i(ai, bi);
                if (use_aa) { af = aaf_add(A.af, B.af); have_af = true; }
                break;
            case OP_SUB:
                vi = sub_i(ai, bi);
                if (use_aa) { af = aaf_sub(A.af, B.af); have_af = true; }
                break;
            case OP_NEG:
                vi = neg_i(ai);
                if (use_aa) { af = aaf_neg(A.af); have_af = true; }
                break;
            case OP_MUL:
                vi = mul_i(ai, bi);
                if (use_aa) { af = aaf_mul(A.af, B.af); have_af = true; }
                break;
            case OP_SQUARE:
                vi = square_i(ai);
                if (use_aa) { af = aaf_mul(A.af, A.af); have_af = true; }
                break;
            case OP_DIV:
                vi = div_i(ai, bi);
                /*  Division by a degenerate (constant) divisor is an
                 *  exact scalar multiply - the common case in
                 *  transform chains.  */
                if (use_aa && bi.lower == bi.upper && bi.lower != 0 &&
                    std::isfinite(1.0f / bi.lower))
                {
                    af = aaf_scale(A.af, 1.0f / bi.lower);
                    have_af = true;
                }
                break;
            case OP_MIN:
                vi = min_i(ai, bi);
                if (!na && !nb)
                {
                    if (ai.upper <= bi.lower)       choice = 1;
                    else if (bi.upper <= ai.lower)  choice = 2;
                }
                if (use_aa && choice)
                {
                    af = (choice == 1) ? A.af : B.af;   // winner verbatim
                    have_af = true;
                }
                break;
            case OP_MAX:
                vi = max_i(ai, bi);
                if (!na && !nb)
                {
                    if (ai.lower >= bi.upper)       choice = 1;
                    else if (bi.lower >= ai.upper)  choice = 2;
                }
                if (use_aa && choice)
                {
                    af = (choice == 1) ? A.af : B.af;
                    have_af = true;
                }
                break;
            case OP_ABS:
                vi = abs_i(ai);
                if (use_aa && ai.lower >= 0) { af = A.af; have_af = true; }
                else if (use_aa && ai.upper <= 0)
                    { af = aaf_neg(A.af); have_af = true; }
                break;
            case OP_SQRT:
                vi = sqrt_i(ai);
                if (use_aa)
                    have_af = aaf_sqrt(A.af, ai, &af);
                break;
            case OP_SIN:    vi = sin_i(ai); break;
            case OP_COS:    vi = cos_i(ai); break;
            case OP_TAN:    vi = tan_i(ai); break;
            case OP_ASIN:   vi = asin_i(ai); break;
            case OP_ACOS:   vi = acos_i(ai); break;
            case OP_ATAN:   vi = atan_i(ai); break;
            case OP_ATAN2:  vi = atan2_i(ai, bi); break;
            case OP_POW:    vi = pow_i(ai, bi); break;
            case OP_MOD:    vi = mod_i(ai, bi); break;
            case OP_EXP:    vi = exp_i(ai); break;
            case OP_LOG: {
                vi.lower = log_f(ai.lower);
                vi.upper = log_f(ai.upper);
                break;
            }
            case OP_FLOOR:
                vi = { floor_f(ai.lower), floor_f(ai.upper) };
                break;
            case OP_CONST:
                vi = { c.imm, c.imm };
                if (use_aa) { af = aaf_const(c.imm); have_af = true; }
                break;
            case OP_COPY:
                vi = ai;
                if (use_aa) { af = A.af; have_af = true; }
                break;
            default:
                FAIL("unsupported op in blob");
        }

        /*  Taint: mirror of tape.cpp tape_eval_i (operand taints,
         *  per-op domain conditions, infinity rules, NaN bounds).  */
        bool u = false;
        switch (c.op)
        {
            case OP_ABS: case OP_SQUARE: case OP_SQRT: case OP_SIN:
            case OP_COS: case OP_TAN: case OP_ASIN: case OP_ACOS:
            case OP_ATAN: case OP_NEG: case OP_EXP: case OP_FLOOR:
            case OP_LOG: case OP_COPY:
                u = na;
                break;
            case OP_CONST:
                u = false;
                break;
            default:
                u = na || nb;
                break;
        }
        const bool a_inf = std::isinf(ai.lower) || std::isinf(ai.upper);
        const bool b_inf = std::isinf(bi.lower) || std::isinf(bi.upper);
        switch (c.op)
        {
            case OP_SQRT:
            case OP_LOG:
                u |= ai.lower < 0;
                break;
            case OP_ASIN:
            case OP_ACOS:
                u |= ai.lower < -1 || ai.upper > 1;
                break;
            case OP_SIN:
            case OP_COS:
            case OP_TAN:
                u |= a_inf;
                break;
            case OP_ADD:
            case OP_SUB:
            case OP_MUL:
                u |= a_inf || b_inf;
                break;
            case OP_DIV:
            case OP_MOD:
                u |= (bi.lower <= 0 && bi.upper >= 0) || a_inf || b_inf;
                break;
            case OP_POW: {
                const bool int_exp = bi.lower == bi.upper &&
                                     bi.lower == (float)(int)bi.lower;
                u |= !int_exp && ai.lower < 0;
                u |= ai.lower <= 0 && ai.upper >= 0 && bi.lower <= 0;
                u |= a_inf || b_inf;
                break;
            }
            case OP_ATAN2:
                u |= ai.lower <= 0 && ai.upper >= 0 &&
                     bi.lower <= 0 && bi.upper >= 0;
                break;
            default:
                break;
        }
        u |= std::isnan(vi.lower) || std::isnan(vi.upper);

        Cell& out_cell = cells[c.out];
        if (use_aa)
        {
            if (!have_af || !aaf_finite(af))
                af = aaf_from_iv(vi);
            out_cell.af = af;
            out_cell.iv = intersect_iv(vi, aaf_hull(af));
        }
        else
        {
            out_cell.af = AAF{};
            out_cell.iv = vi;
        }
        out_cell.tnt = u;
        out.choice[ci] = choice;
    }

    out.root = cells[p.root].iv;
    out.root_tnt = cells[p.root].tnt;
    return out;
}

////////////////////////////////////////////////////////////////////////////
//  Push simulation: tape_push's STANDARD backward walk on a Prog

Prog push_prog(const Prog& p, const std::vector<uint8_t>& choice)
{
    std::vector<uint8_t> live(p.num_slots, 0);
    std::vector<uint8_t> verdict(p.clauses.size(), 0);
    live[p.root] = 1;
    for (size_t k = p.clauses.size(); k-- > 0; )
    {
        const Clause5& c = p.clauses[k];
        if (!live[c.out])
            continue;
        live[c.out] = 0;
        if ((c.op == OP_MIN || c.op == OP_MAX) && choice[k])
        {
            const uint32_t keep = choice[k] == 1 ? c.a : c.b;
            verdict[k] = (keep == c.out) ? 0 : (choice[k] == 1 ? 2 : 3);
            live[keep] = 1;
            continue;
        }
        verdict[k] = 1;
        const bool unary =
                c.op == OP_ABS || c.op == OP_SQUARE || c.op == OP_SQRT ||
                c.op == OP_SIN || c.op == OP_COS || c.op == OP_TAN ||
                c.op == OP_ASIN || c.op == OP_ACOS || c.op == OP_ATAN ||
                c.op == OP_NEG || c.op == OP_EXP || c.op == OP_FLOOR ||
                c.op == OP_LOG || c.op == OP_COPY;
        if (c.op != OP_CONST)
        {
            live[c.a] = 1;
            if (!unary)
                live[c.b] = 1;
        }
    }
    Prog out = p;
    out.clauses.clear();
    for (size_t k = 0; k < p.clauses.size(); ++k)
    {
        const Clause5& c = p.clauses[k];
        if (verdict[k] == 0)
            continue;
        if (verdict[k] == 1)
            out.clauses.push_back(c);
        else
            out.clauses.push_back({ OP_COPY, c.out,
                                    verdict[k] == 2 ? c.a : c.b, 0, 0 });
    }
    return out;
}

////////////////////////////////////////////////////////////////////////////

struct BoxStats
{
    uint64_t empty = 0, full = 0, ambig = 0;
    uint64_t len_sum = 0;
    uint64_t boxes = 0;
};

void classify_box(const Prog& p, Interval X, Interval Y, Interval Z,
                  bool use_aa, BoxStats& st, Prog* pushed_out = nullptr)
{
    const Classified c = classify(p, X, Y, Z, use_aa);
    ++st.boxes;
    if (!c.root_tnt && c.root.lower > 0)
        ++st.empty;
    else if (!c.root_tnt && c.root.upper < 0)
        ++st.full;
    else
    {
        ++st.ambig;
        Prog pushed = push_prog(p, c.choice);
        st.len_sum += pushed.clauses.size();
        if (pushed_out)
            *pushed_out = std::move(pushed);
    }
}

/*  Soundness: pointwise values must lie within the classified root
 *  bounds, and a push driven by these choices must leave pointwise
 *  evaluation bit-identical (the fuzzer's exactness standard).  */
void soundness_check(const Prog& p, Interval X, Interval Y, Interval Z,
                     std::mt19937& rng, int n_points)
{
    const Classified c = classify(p, X, Y, Z, true);
    const Prog pushed = push_prog(p, c.choice);
    std::uniform_real_distribution<float> ux(X.lower, X.upper),
            uy(Y.lower, Y.upper), uz(Z.lower, Z.upper);
    std::vector<float> slots;
    for (int i = 0; i < n_points; ++i)
    {
        const float x = ux(rng), y = uy(rng), z = uz(rng);
        const float v = eval_point(p, x, y, z, slots);
        const float vp = eval_point(pushed, x, y, z, slots);
        CAPTURE(x);
        CAPTURE(y);
        CAPTURE(z);
        CAPTURE(v);
        CAPTURE(vp);
        if (std::isnan(v))
        {
            REQUIRE(std::isnan(vp));
            continue;
        }
        REQUIRE(memcmp(&v, &vp, 4) == 0);   // push exactness, bitwise
        if (!c.root_tnt)
        {
            REQUIRE(v >= c.root.lower);
            REQUIRE(v <= c.root.upper);
        }
    }
}

const char* CORPUS[] = {
    "-r++qXqYqZf1",
    "i-r++qXqYqZf1-r++q-Xf0.5q-Yf0.25q-Zf0.1f0.8",
    "a-r++qXqYqZf1nY",
    "+s*Xf3.0c*Yf2.0",
    "ai-r++qXqYqZf1-r++q-Xf0.5qYqZf0.7n-Zf0.2",
    "*+Xf2.0+Yf2.0",
    "i-r++qXqYqZf1-r++q-Xf10qYqZf0.5",
    "-r++q/Xf1.7q/Yf1.3qZf1",
    "aa-f-0.6X-Xf0.6aa-f-0.6Y-Yf0.6a-f-0.6Z-Zf0.6",
};

}  // namespace

////////////////////////////////////////////////////////////////////////////

TEST_CASE("AA soundness sweep on corpus", "[.aafuzz]")
{
    std::mt19937 rng(0xAFF17E5);
    std::uniform_real_distribution<float> uc(-1.4f, 1.4f),
            uw(0.01f, 1.2f);
    int boxes = 400;
    if (const char* env = getenv("STIBIUM_AA_BOXES"))
        boxes = atoi(env);

    for (const char* expr : CORPUS)
    {
        const Prog p = prog_from_tree(expr);
        CAPTURE(expr);
        for (int b = 0; b < boxes; ++b)
        {
            const float cx = uc(rng), cy = uc(rng), cz = uc(rng);
            const float wx = uw(rng), wy = uw(rng), wz = uw(rng);
            soundness_check(p,
                            { cx - wx, cx + wx },
                            { cy - wy, cy + wy },
                            { cz - wz, cz + wz }, rng, 8);
        }
    }
}

TEST_CASE("AA vs IA on a dumped app deck", "[.aa]")
{
    const char* path = getenv("STIBIUM_AA_BLOB");
    if (!path)
        path = "zeiss.blob";
    FILE* f = fopen(path, "rb");
    if (!f)
    {
        WARN("no blob at " << path
             << " (dump one with STIBIUM_GPU=1 STIBIUM_GPU_BLOB_DUMP=...)");
        return;
    }
    uint32_t hdr[9];
    REQUIRE(fread(hdr, 4, 9, f) == 9);
    std::vector<uint32_t> blob;
    uint32_t w;
    while (fread(&w, 4, 1, f) == 1)
        blob.push_back(w);
    fclose(f);

    const Prog p = prog_from_blob(blob);
    const uint32_t ni = hdr[0], nj = hdr[1];
    const float xmin = bits_to_f(hdr[3]), xmax = bits_to_f(hdr[4]);
    const float ymin = bits_to_f(hdr[5]), ymax = bits_to_f(hdr[6]);
    const float zmin = bits_to_f(hdr[7]), zmax = bits_to_f(hdr[8]);
    WARN("deck: " << p.clauses.size() << " clauses, "
                  << p.num_slots << " slots; region " << ni << "x" << nj);

    constexpr uint32_t TILE = 16;
    const uint32_t tx = (ni + TILE - 1) / TILE;
    const uint32_t ty = (nj + TILE - 1) / TILE;
    const float sx = (xmax - xmin) / ni, sy = (ymax - ymin) / nj;

    /*  Level 1: (16px tile x full z), both arithmetics.  */
    BoxStats ia, aa;
    std::vector<std::pair<Prog, Prog>> ambig_tiles;   // (ia, aa) pushed
    std::vector<std::pair<Interval, Interval>> ambig_xy;
    for (uint32_t tj = 0; tj < ty; ++tj)
        for (uint32_t ti = 0; ti < tx; ++ti)
        {
            const Interval X = { xmin + sx * (ti * TILE),
                                 xmin + sx * fminf((ti + 1) * TILE, ni) };
            const Interval Y = { ymin + sy * (tj * TILE),
                                 ymin + sy * fminf((tj + 1) * TILE, nj) };
            const Interval Z = { zmin, zmax };
            Prog pushed_ia, pushed_aa;
            classify_box(p, X, Y, Z, false, ia, &pushed_ia);
            classify_box(p, X, Y, Z, true, aa, &pushed_aa);
            if (!pushed_ia.clauses.empty() && !pushed_aa.clauses.empty()
                && ambig_tiles.size() < 64 && (ti + tj) % 7 == 0)
            {
                ambig_tiles.push_back({ std::move(pushed_ia),
                                        std::move(pushed_aa) });
                ambig_xy.push_back({ X, Y });
            }
        }

    WARN("L1 tiles (" << ia.boxes << "): IA "
         << ia.empty << "E/" << ia.full << "F/" << ia.ambig
         << "A avg-tape " << (ia.ambig ? ia.len_sum / ia.ambig : 0)
         << "  |  AA "
         << aa.empty << "E/" << aa.full << "F/" << aa.ambig
         << "A avg-tape " << (aa.ambig ? aa.len_sum / aa.ambig : 0));

    /*  Level 2: (tile x 16-voxel z-slab) on the level-1 pushed tapes
     *  (sampled subset of ambiguous tiles).  */
    const uint32_t nk = hdr[2];
    const float sz = (zmax - zmin) / nk;
    BoxStats ia2, aa2;
    for (size_t t = 0; t < ambig_tiles.size(); ++t)
        for (uint32_t sl = 0; sl * TILE < nk; ++sl)
        {
            const Interval Z = { zmin + sz * (sl * TILE),
                                 zmin + sz * fminf((sl + 1) * TILE, nk) };
            classify_box(ambig_tiles[t].first, ambig_xy[t].first,
                         ambig_xy[t].second, Z, false, ia2);
            classify_box(ambig_tiles[t].second, ambig_xy[t].first,
                         ambig_xy[t].second, Z, true, aa2);
        }

    WARN("L2 slabs (" << ia2.boxes << ", from " << ambig_tiles.size()
         << " sampled tiles): IA "
         << ia2.empty << "E/" << ia2.full << "F/" << ia2.ambig
         << "A avg-tape " << (ia2.ambig ? ia2.len_sum / ia2.ambig : 0)
         << "  |  AA "
         << aa2.empty << "E/" << aa2.full << "F/" << aa2.ambig
         << "A avg-tape " << (aa2.ambig ? aa2.len_sum / aa2.ambig : 0));

    /*  Soundness spot-check on real-deck boxes.  */
    std::mt19937 rng(0x5EED);
    for (uint32_t t = 0; t < 24; ++t)
    {
        const uint32_t ti = (t * 7) % tx, tj = (t * 5) % ty;
        const Interval X = { xmin + sx * (ti * TILE),
                             xmin + sx * fminf((ti + 1) * TILE, ni) };
        const Interval Y = { ymin + sy * (tj * TILE),
                             ymin + sy * fminf((tj + 1) * TILE, nj) };
        soundness_check(p, X, Y, { zmin, zmax }, rng, 6);
    }
}

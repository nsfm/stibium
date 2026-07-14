#include <algorithm>
#include <atomic>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <utility>
#include <vector>

#include "fab/tree/tape.h"

#include "fab/tree/tree.h"
#include "fab/tree/node/node.h"
#include "fab/tree/grid.h"

#include "fab/tree/math/math_f.h"
#include "fab/tree/math/math_i.h"
#include "fab/tree/math/math_r.h"
#include "fab/tree/math/math_g.h"

#include "fab/util/switches.h"

namespace {

/*  One operation over slot-indexed values.  out/a/b/m are slot ids;
 *  m and payload are used only by OP_GRID.  OP_CONST clauses (imm
 *  writes) exist only in binary-pushed tapes - base-tape constants
 *  are pre-filled in the ctx and never appear as clauses. */
struct Clause
{
    Opcode op;
    uint32_t out;
    uint32_t a, b, m;
    float imm;
    void* payload;
};

enum Arity { ARITY_NONE, ARITY_UNARY, ARITY_BINARY, ARITY_GRID };

/*  "no remap" sentinel: slot 0 is a real slot (the first constant or
 *  axis), so 0 cannot mean "none" here  */
constexpr uint32_t NO_REMAP = UINT32_MAX;

Arity arity(Opcode op)
{
    switch (op)
    {
        case OP_ADD: case OP_SUB: case OP_MUL: case OP_DIV:
        case OP_MIN: case OP_MAX: case OP_POW: case OP_ATAN2:
        case OP_MOD:
            return ARITY_BINARY;
        case OP_ABS: case OP_SQUARE: case OP_SQRT: case OP_SIN:
        case OP_COS: case OP_TAN: case OP_ASIN: case OP_ACOS:
        case OP_ATAN: case OP_NEG: case OP_EXP: case OP_FLOOR:
        case OP_LOG:
            return ARITY_UNARY;
        case OP_GRID:
            return ARITY_GRID;
        default:
            return ARITY_NONE;
    }
}

}  // namespace

struct Tape_
{
    std::vector<Clause> clauses;
    uint32_t root;

    /*  Valid when this tape was pushed for an interval region  */
    Interval X, Y, Z;
    bool has_bounds = false;

    /*  No min/max clauses left; pushing is a no-op  */
    bool terminal = false;

    Tape_* parent = nullptr;   // strong reference
    std::atomic<int> refs { 1 };
};

struct Deck_
{
    Tape_* base = nullptr;
    unsigned num_slots = 0;

    /*  (slot, value) pairs, filled once per ctx  */
    std::vector<std::pair<uint32_t, float>> constants;

    /*  Slots the eval drivers fill with coordinate data.  Normally
     *  at most one slot per axis, but kept as vectors so a tree with
     *  duplicate axis nodes would still evaluate correctly. */
    std::vector<uint32_t> xs, ys, zs;
    uint8_t axes = 0;

    /*  Grid payloads referenced by clauses, retained for our lifetime  */
    std::vector<MeshGrid*> grids;
};

struct TapeCtx_
{
    unsigned num_slots = 0;
    const Deck_* deck = nullptr;

    float* f = nullptr;        // [num_slots]
    Interval* i = nullptr;     // [num_slots]
    float* r = nullptr;        // [num_slots * MIN_VOLUME]

    /*  Whether constant r-rows currently hold floats or derivatives
     *  (the g rows alias the r rows, exactly like Results.r)  */
    bool g_mode = false;

    /*  tape_push scratch, all [num_slots]  */
    std::vector<uint8_t> live;
    std::vector<uint32_t> remap;
    std::vector<uint8_t> boolean_;
    std::vector<uint8_t> has_imm;
    std::vector<float> immv;
};

////////////////////////////////////////////////////////////////////////////////

extern "C" Deck* deck_from_tree(MathTree* tree)
{
    Deck_* deck = new Deck_;

    unsigned total = tree->num_constants;
    for (unsigned level = 0; level < tree->num_levels; ++level)
        total += tree->active[level];

    std::unordered_map<const Node*, uint32_t> slot;
    slot.reserve(total);
    uint32_t next = 0;

    for (unsigned c = 0; c < tree->num_constants; ++c)
    {
        const Node* n = tree->constants[c];
        slot[n] = next;
        deck->constants.emplace_back(next, n->results.f);
        next++;
    }

    std::vector<Clause> clauses;
    for (unsigned level = 0; level < tree->num_levels; ++level)
    {
        for (unsigned n = 0; n < tree->active[level]; ++n)
        {
            const Node* node = tree->nodes[level][n];
            const uint32_t out = next++;
            slot[node] = out;

            switch (node->opcode)
            {
                case OP_X:  deck->xs.push_back(out);
                            deck->axes |= (1 << 2);
                            continue;
                case OP_Y:  deck->ys.push_back(out);
                            deck->axes |= (1 << 1);
                            continue;
                case OP_Z:  deck->zs.push_back(out);
                            deck->axes |= (1 << 0);
                            continue;
                default:    break;
            }

            Clause c = {};
            c.op = node->opcode;
            c.out = out;
            c.a = node->lhs ? slot.at(node->lhs) : 0;
            c.b = node->rhs ? slot.at(node->rhs) : 0;
            c.m = node->mhs ? slot.at(node->mhs) : 0;
            c.payload = node->payload;
            if (node->opcode == OP_GRID)
            {
                MeshGrid* g = static_cast<MeshGrid*>(node->payload);
                grid_retain(g);
                deck->grids.push_back(g);
            }
            clauses.push_back(c);
        }
    }

    deck->num_slots = next;

    Tape_* base = new Tape_;
    base->clauses = std::move(clauses);
    base->root = slot.at(tree->head);
    base->terminal = true;
    for (const auto& c : base->clauses)
    {
        if (c.op == OP_MIN || c.op == OP_MAX)
        {
            base->terminal = false;
            break;
        }
    }
    deck->base = base;

    return deck;
}

extern "C" void deck_free(Deck* deck)
{
    if (deck == nullptr)    return;
    tape_release(deck->base);
    for (MeshGrid* g : deck->grids)
        grid_release(g);
    delete deck;
}

extern "C" Tape* deck_base(const Deck* deck)        { return deck->base; }
extern "C" unsigned deck_slots(const Deck* deck)    { return deck->num_slots; }
extern "C" uint8_t deck_active_axes(const Deck* deck) { return deck->axes; }

////////////////////////////////////////////////////////////////////////////////

static void fill_const_rows(TapeCtx_* ctx, bool g_mode)
{
    for (const auto& cv : ctx->deck->constants)
    {
        float* row = ctx->r + size_t(cv.first) * MIN_VOLUME;
        if (g_mode)
        {
            derivative* d = reinterpret_cast<derivative*>(row);
            for (int q = 0; q < MIN_VOLUME / 4; ++q)
                d[q] = derivative{ cv.second, 0, 0, 0 };
        }
        else
        {
            for (int q = 0; q < MIN_VOLUME; ++q)
                row[q] = cv.second;
        }
    }
    ctx->g_mode = g_mode;
}

extern "C" TapeCtx* tape_ctx_new(const Deck* deck)
{
    TapeCtx_* ctx = new TapeCtx_;
    ctx->deck = deck;
    ctx->num_slots = deck->num_slots;

    ctx->f = static_cast<float*>(calloc(deck->num_slots, sizeof(float)));
    ctx->i = static_cast<Interval*>(calloc(deck->num_slots, sizeof(Interval)));
    ctx->r = static_cast<float*>(
            calloc(size_t(deck->num_slots) * MIN_VOLUME, sizeof(float)));

    ctx->live.resize(deck->num_slots);
    ctx->remap.resize(deck->num_slots);
    ctx->boolean_.resize(deck->num_slots);
    ctx->has_imm.resize(deck->num_slots);
    ctx->immv.resize(deck->num_slots);

    for (const auto& cv : deck->constants)
    {
        ctx->f[cv.first] = cv.second;
        ctx->i[cv.first] = Interval{ cv.second, cv.second };
    }
    fill_const_rows(ctx, false);

    return ctx;
}

extern "C" void tape_ctx_free(TapeCtx* ctx)
{
    if (ctx == nullptr) return;
    free(ctx->f);
    free(ctx->i);
    free(ctx->r);
    delete ctx;
}

////////////////////////////////////////////////////////////////////////////////

extern "C" float tape_eval_f(const Tape* tape, TapeCtx* ctx,
                             const float x, const float y, const float z)
{
    float* v = ctx->f;
    for (uint32_t s : ctx->deck->xs)  v[s] = X_f(x);
    for (uint32_t s : ctx->deck->ys)  v[s] = Y_f(y);
    for (uint32_t s : ctx->deck->zs)  v[s] = Z_f(z);

    for (const auto& c : tape->clauses)
    {
        const float A = v[c.a], B = v[c.b];
        switch (c.op)
        {
            case OP_ADD:    v[c.out] = add_f(A, B); break;
            case OP_SUB:    v[c.out] = sub_f(A, B); break;
            case OP_MUL:    v[c.out] = mul_f(A, B); break;
            case OP_DIV:    v[c.out] = div_f(A, B); break;
            case OP_MIN:    v[c.out] = min_f(A, B); break;
            case OP_MAX:    v[c.out] = max_f(A, B); break;
            case OP_POW:    v[c.out] = pow_f(A, B); break;
            case OP_ATAN2:  v[c.out] = atan2_f(A, B); break;
            case OP_MOD:    v[c.out] = mod_f(A, B); break;

            case OP_ABS:    v[c.out] = abs_f(A); break;
            case OP_SQUARE: v[c.out] = square_f(A); break;
            case OP_SQRT:   v[c.out] = sqrt_f(A); break;
            case OP_SIN:    v[c.out] = sin_f(A); break;
            case OP_COS:    v[c.out] = cos_f(A); break;
            case OP_TAN:    v[c.out] = tan_f(A); break;
            case OP_ASIN:   v[c.out] = asin_f(A); break;
            case OP_ACOS:   v[c.out] = acos_f(A); break;
            case OP_ATAN:   v[c.out] = atan_f(A); break;
            case OP_NEG:    v[c.out] = neg_f(A); break;
            case OP_EXP:    v[c.out] = exp_f(A); break;
            case OP_FLOOR:  v[c.out] = floor_f(A); break;
            case OP_LOG:    v[c.out] = log_f(A); break;

            case OP_GRID:
                v[c.out] = grid_eval_f(
                        static_cast<const MeshGrid*>(c.payload),
                        A, B, v[c.m]);
                break;

            case OP_CONST:  v[c.out] = c.imm; break;
            default: ;
        }
    }
    return v[tape->root];
}

extern "C" Interval tape_eval_i(const Tape* tape, TapeCtx* ctx,
                                const Interval X, const Interval Y,
                                const Interval Z)
{
    Interval* v = ctx->i;
    for (uint32_t s : ctx->deck->xs)  v[s] = X_i(X);
    for (uint32_t s : ctx->deck->ys)  v[s] = Y_i(Y);
    for (uint32_t s : ctx->deck->zs)  v[s] = Z_i(Z);

    for (const auto& c : tape->clauses)
    {
        const Interval A = v[c.a], B = v[c.b];
        switch (c.op)
        {
            case OP_ADD:    v[c.out] = add_i(A, B); break;
            case OP_SUB:    v[c.out] = sub_i(A, B); break;
            case OP_MUL:    v[c.out] = mul_i(A, B); break;
            case OP_DIV:    v[c.out] = div_i(A, B); break;
            case OP_MIN:    v[c.out] = min_i(A, B); break;
            case OP_MAX:    v[c.out] = max_i(A, B); break;
            case OP_POW:    v[c.out] = pow_i(A, B); break;
            case OP_ATAN2:  v[c.out] = atan2_i(A, B); break;
            case OP_MOD:    v[c.out] = mod_i(A, B); break;

            case OP_ABS:    v[c.out] = abs_i(A); break;
            case OP_SQUARE: v[c.out] = square_i(A); break;
            case OP_SQRT:   v[c.out] = sqrt_i(A); break;
            case OP_SIN:    v[c.out] = sin_i(A); break;
            case OP_COS:    v[c.out] = cos_i(A); break;
            case OP_TAN:    v[c.out] = tan_i(A); break;
            case OP_ASIN:   v[c.out] = asin_i(A); break;
            case OP_ACOS:   v[c.out] = acos_i(A); break;
            case OP_ATAN:   v[c.out] = atan_i(A); break;
            case OP_NEG:    v[c.out] = neg_i(A); break;
            case OP_EXP:    v[c.out] = exp_i(A); break;
            case OP_FLOOR:  v[c.out] = floor_i(A); break;
            case OP_LOG:    v[c.out] = log_i(A); break;

            case OP_GRID:
                v[c.out] = grid_eval_i(
                        static_cast<const MeshGrid*>(c.payload),
                        A, B, v[c.m]);
                break;

            case OP_CONST:  v[c.out] = Interval{ c.imm, c.imm }; break;
            default: ;
        }
    }
    return v[tape->root];
}

extern "C" const float* tape_eval_r(const Tape* tape, TapeCtx* ctx, Region r)
{
    if (ctx->g_mode)
        fill_const_rows(ctx, false);

    const int c_ = r.voxels;
    float* base = ctx->r;
#define ROW(s)  (base + size_t(s) * MIN_VOLUME)

    for (uint32_t s : ctx->deck->xs)  X_r(r.X, ROW(s), c_);
    for (uint32_t s : ctx->deck->ys)  Y_r(r.Y, ROW(s), c_);
    for (uint32_t s : ctx->deck->zs)  Z_r(r.Z, ROW(s), c_);

    for (const auto& c : tape->clauses)
    {
        float *A = ROW(c.a), *B = ROW(c.b), *R = ROW(c.out);
        switch (c.op)
        {
            case OP_ADD:    add_r(A, B, R, c_); break;
            case OP_SUB:    sub_r(A, B, R, c_); break;
            case OP_MUL:    mul_r(A, B, R, c_); break;
            case OP_DIV:    div_r(A, B, R, c_); break;
            case OP_MIN:    min_r(A, B, R, c_); break;
            case OP_MAX:    max_r(A, B, R, c_); break;
            case OP_POW:    pow_r(A, B, R, c_); break;
            case OP_ATAN2:  atan2_r(A, B, R, c_); break;
            case OP_MOD:    mod_r(A, B, R, c_); break;

            case OP_ABS:    abs_r(A, R, c_); break;
            case OP_SQUARE: square_r(A, R, c_); break;
            case OP_SQRT:   sqrt_r(A, R, c_); break;
            case OP_SIN:    sin_r(A, R, c_); break;
            case OP_COS:    cos_r(A, R, c_); break;
            case OP_TAN:    tan_r(A, R, c_); break;
            case OP_ASIN:   asin_r(A, R, c_); break;
            case OP_ACOS:   acos_r(A, R, c_); break;
            case OP_ATAN:   atan_r(A, R, c_); break;
            case OP_NEG:    neg_r(A, R, c_); break;
            case OP_EXP:    exp_r(A, R, c_); break;
            case OP_FLOOR:  floor_r(A, R, c_); break;
            case OP_LOG:    log_r(A, R, c_); break;

            case OP_GRID: {
                const MeshGrid* grid =
                        static_cast<const MeshGrid*>(c.payload);
                const float* M = ROW(c.m);
                for (int q = 0; q < c_; ++q)
                    R[q] = grid_eval_f(grid, A[q], B[q], M[q]);
                break;
            }

            case OP_CONST:
                for (int q = 0; q < c_; ++q)
                    R[q] = c.imm;
                break;
            default: ;
        }
    }
    return ROW(tape->root);
#undef ROW
}

extern "C" const derivative* tape_eval_g(const Tape* tape, TapeCtx* ctx,
                                         Region r)
{
    if (!ctx->g_mode)
        fill_const_rows(ctx, true);

    const int c_ = r.voxels;
    float* base = ctx->r;
#define ROW(s)  (reinterpret_cast<derivative*>(base + size_t(s) * MIN_VOLUME))

    for (uint32_t s : ctx->deck->xs)  X_g(r.X, ROW(s), c_);
    for (uint32_t s : ctx->deck->ys)  Y_g(r.Y, ROW(s), c_);
    for (uint32_t s : ctx->deck->zs)  Z_g(r.Z, ROW(s), c_);

    for (const auto& c : tape->clauses)
    {
        derivative *A = ROW(c.a), *B = ROW(c.b), *R = ROW(c.out);
        switch (c.op)
        {
            case OP_ADD:    add_g(A, B, R, c_); break;
            case OP_SUB:    sub_g(A, B, R, c_); break;
            case OP_MUL:    mul_g(A, B, R, c_); break;
            case OP_DIV:    div_g(A, B, R, c_); break;
            case OP_MIN:    min_g(A, B, R, c_); break;
            case OP_MAX:    max_g(A, B, R, c_); break;
            case OP_POW:    pow_g(A, B, R, c_); break;
            case OP_ATAN2:  atan2_g(A, B, R, c_); break;
            case OP_MOD:    mod_g(A, B, R, c_); break;

            case OP_ABS:    abs_g(A, R, c_); break;
            case OP_SQUARE: square_g(A, R, c_); break;
            case OP_SQRT:   sqrt_g(A, R, c_); break;
            case OP_SIN:    sin_g(A, R, c_); break;
            case OP_COS:    cos_g(A, R, c_); break;
            case OP_TAN:    tan_g(A, R, c_); break;
            case OP_ASIN:   asin_g(A, R, c_); break;
            case OP_ACOS:   acos_g(A, R, c_); break;
            case OP_ATAN:   atan_g(A, R, c_); break;
            case OP_NEG:    neg_g(A, R, c_); break;
            case OP_EXP:    exp_g(A, R, c_); break;
            case OP_FLOOR:  floor_g(A, R, c_); break;
            case OP_LOG:    log_g(A, R, c_); break;

            case OP_GRID: {
                // Chain rule through the (possibly remapped)
                // sample coordinates in A, B, and M
                const MeshGrid* grid =
                        static_cast<const MeshGrid*>(c.payload);
                const derivative* M = ROW(c.m);
                for (int q = 0; q < c_; ++q)
                {
                    float v, gx, gy, gz;
                    grid_eval_g(grid, A[q].v, B[q].v, M[q].v,
                                &v, &gx, &gy, &gz);
                    R[q].v = v;
                    R[q].dx = gx*A[q].dx + gy*B[q].dx + gz*M[q].dx;
                    R[q].dy = gx*A[q].dy + gy*B[q].dy + gz*M[q].dy;
                    R[q].dz = gx*A[q].dz + gy*B[q].dz + gz*M[q].dz;
                }
                break;
            }

            case OP_CONST:
                for (int q = 0; q < c_; ++q)
                    R[q] = derivative{ c.imm, 0, 0, 0 };
                break;
            default: ;
        }
    }
    return ROW(tape->root);
#undef ROW
}

////////////////////////////////////////////////////////////////////////////////

extern "C" Tape* tape_push(Tape* tape, TapeCtx* ctx,
                           const Interval X, const Interval Y,
                           const Interval Z, const TapePushMode mode)
{
    if (tape->terminal)
        return tape_retain(tape);

    auto& live = ctx->live;
    auto& remap = ctx->remap;
    auto& boolean_ = ctx->boolean_;
    auto& has_imm = ctx->has_imm;
    const Interval* iv = ctx->i;

    std::fill(live.begin(), live.end(), 0);
    std::fill(remap.begin(), remap.end(), NO_REMAP);
    std::fill(boolean_.begin(), boolean_.end(), 1);
    std::fill(has_imm.begin(), has_imm.end(), 0);

    live[tape->root] = 1;
    bool changed = false;

    /*  Backward = reverse-topological, so every parent is processed
     *  before its children and boolean flags are final on arrival
     *  (the old top-down level walk, over clauses).  */
    for (auto it = tape->clauses.rbegin(); it != tape->clauses.rend(); ++it)
    {
        const Clause& c = *it;
        if (!live[c.out])
            continue;

        /*  disable_nodes_binary: a sign-determined value in boolean
         *  (min/max/neg) context collapses to its upper bound, and
         *  its subtree dies unless referenced elsewhere.  */
        if (mode == TAPE_PUSH_BINARY && boolean_[c.out] &&
            (iv[c.out].lower >= 0 || iv[c.out].upper < 0))
        {
            has_imm[c.out] = 1;
            ctx->immv[c.out] = iv[c.out].upper;
            changed = true;
            continue;
        }

        /*  disable_nodes: decided min/max branches drop the clause
         *  and remap readers to the survivor (exact same value the
         *  min/max would return, since the comparison guarantees
         *  which side wins pointwise across the region).  */
        if (c.op == OP_MAX &&
                (iv[c.a].lower >= iv[c.b].upper ||
                 iv[c.b].lower >= iv[c.a].upper))
        {
            const uint32_t keep =
                    (iv[c.a].lower >= iv[c.b].upper) ? c.a : c.b;
            remap[c.out] = keep;
            live[keep] = 1;
            boolean_[keep] &= boolean_[c.out];
            changed = true;
            continue;
        }
        if (c.op == OP_MIN &&
                (iv[c.a].upper <= iv[c.b].lower ||
                 iv[c.b].upper <= iv[c.a].lower))
        {
            const uint32_t keep =
                    (iv[c.a].upper <= iv[c.b].lower) ? c.a : c.b;
            remap[c.out] = keep;
            live[keep] = 1;
            boolean_[keep] &= boolean_[c.out];
            changed = true;
            continue;
        }

        /*  Kept clause: mark children live and propagate boolean
         *  context (min/max/neg are sign-transparent; every other
         *  opcode breaks the chain).  */
        const bool transparent =
                (c.op == OP_MIN || c.op == OP_MAX || c.op == OP_NEG);
        switch (arity(c.op))
        {
            case ARITY_GRID:
                live[c.m] = 1;
                boolean_[c.m] &= transparent ? boolean_[c.out] : 0;
                // fallthrough
            case ARITY_BINARY:
                live[c.b] = 1;
                boolean_[c.b] &= transparent ? boolean_[c.out] : 0;
                // fallthrough
            case ARITY_UNARY:
                live[c.a] = 1;
                boolean_[c.a] &= transparent ? boolean_[c.out] : 0;
                break;
            case ARITY_NONE:
                break;
        }
    }

    if (!changed)
        return tape_retain(tape);

    Tape_* out = new Tape_;
    out->clauses.reserve(tape->clauses.size());

    auto resolve = [&remap](uint32_t s) {
        while (remap[s] != NO_REMAP)
            s = remap[s];
        return s;
    };

    bool terminal = true;
    for (const Clause& c : tape->clauses)
    {
        if (!live[c.out] || remap[c.out] != NO_REMAP)
            continue;
        if (has_imm[c.out])
        {
            Clause k = {};
            k.op = OP_CONST;
            k.out = c.out;
            k.imm = ctx->immv[c.out];
            out->clauses.push_back(k);
            continue;
        }
        Clause k = c;
        switch (arity(c.op))
        {
            case ARITY_GRID:    k.m = resolve(c.m);  // fallthrough
            case ARITY_BINARY:  k.b = resolve(c.b);  // fallthrough
            case ARITY_UNARY:   k.a = resolve(c.a);  break;
            case ARITY_NONE:    break;
        }
        if (k.op == OP_MIN || k.op == OP_MAX)
            terminal = false;
        out->clauses.push_back(k);
    }

    out->root = resolve(tape->root);
    out->terminal = terminal;
    out->X = X;
    out->Y = Y;
    out->Z = Z;
    out->has_bounds = true;
    out->parent = tape_retain(tape);

    return out;
}

////////////////////////////////////////////////////////////////////////////////

extern "C" Tape* tape_retain(Tape* tape)
{
    tape->refs.fetch_add(1, std::memory_order_relaxed);
    return tape;
}

extern "C" void tape_release(Tape* tape)
{
    while (tape &&
           tape->refs.fetch_sub(1, std::memory_order_acq_rel) == 1)
    {
        Tape_* parent = tape->parent;
        delete tape;
        tape = parent;
    }
}

extern "C" Tape* tape_base_for_point(Tape* tape,
                                     const float x, const float y,
                                     const float z)
{
    while (tape->parent)
    {
        if (tape->has_bounds &&
            x >= tape->X.lower && x <= tape->X.upper &&
            y >= tape->Y.lower && y <= tape->Y.upper &&
            z >= tape->Z.lower && z <= tape->Z.upper)
        {
            break;
        }
        tape = tape->parent;
    }
    return tape;
}

extern "C" unsigned tape_length(const Tape* tape)
{
    return tape->clauses.size();
}

extern "C" bool tape_is_terminal(const Tape* tape)
{
    return tape->terminal;
}

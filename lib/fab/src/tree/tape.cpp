#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
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
        case OP_LOG: case OP_COPY:
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

    /*  Push depth from the base tape (diagnostics)  */
    unsigned depth = 0;

    Tape_* parent = nullptr;   // strong reference
    std::atomic<int> refs { 1 };
};

namespace {

/*  STIBIUM_TAPE_STATS=1: accumulate tape length by push depth and
 *  dump a shrinkage curve to stderr at exit (the MPR paper's
 *  clauses-per-level table, for our models).  */
struct TapeStats
{
    static constexpr unsigned DEPTHS = 48;
    std::atomic<uint64_t> count[DEPTHS] {};
    std::atomic<uint64_t> len_sum[DEPTHS] {};
    std::atomic<uint64_t> len_min[DEPTHS] {};
    std::atomic<uint64_t> len_max[DEPTHS] {};
    const bool enabled;

    TapeStats() : enabled(getenv("STIBIUM_TAPE_STATS") != nullptr)
    {
        for (auto& m : len_min)
            m.store(UINT64_MAX, std::memory_order_relaxed);
    }

    void record(unsigned depth, uint64_t len)
    {
        if (!enabled)   return;
        const unsigned d = depth < DEPTHS ? depth : DEPTHS - 1;
        count[d].fetch_add(1, std::memory_order_relaxed);
        len_sum[d].fetch_add(len, std::memory_order_relaxed);
        uint64_t cur = len_min[d].load(std::memory_order_relaxed);
        while (len < cur &&
               !len_min[d].compare_exchange_weak(cur, len)) {}
        cur = len_max[d].load(std::memory_order_relaxed);
        while (len > cur &&
               !len_max[d].compare_exchange_weak(cur, len)) {}
    }

    ~TapeStats()
    {
        if (!enabled)   return;
        fprintf(stderr,
                "[tape-stats] depth      tapes    avg-len    min    max\n");
        for (unsigned d = 0; d < DEPTHS; ++d)
        {
            const uint64_t n = count[d].load();
            if (!n) continue;
            fprintf(stderr, "[tape-stats] %5u %10llu %10.1f %6llu %6llu\n",
                    d, (unsigned long long)n,
                    double(len_sum[d].load()) / n,
                    (unsigned long long)len_min[d].load(),
                    (unsigned long long)len_max[d].load());
        }
    }
};

TapeStats g_tape_stats;

/*  Pushed tapes are created and released on the same thread
 *  (push/release pairs bracketing recursion), so a thread-local
 *  freelist recycles Tape objects - and their clause vectors'
 *  capacity - without locks.  A big model pushes millions of times
 *  per mesh (see STIBIUM_TAPE_STATS), so this matters.  */
struct SpareTapes
{
    std::vector<Tape_*> v;
    ~SpareTapes();
};
thread_local SpareTapes tape_spares;
constexpr size_t SPARE_CAP = 32;

Tape_* tape_alloc()
{
    if (!tape_spares.v.empty())
    {
        Tape_* t = tape_spares.v.back();
        tape_spares.v.pop_back();
        t->clauses.clear();            // keeps capacity
        t->has_bounds = false;
        t->terminal = false;
        t->depth = 0;
        t->parent = nullptr;
        t->refs.store(1, std::memory_order_relaxed);
        return t;
    }
    return new Tape_;
}

void tape_dealloc(Tape_* t)
{
    if (t->parent != nullptr && tape_spares.v.size() < SPARE_CAP)
        tape_spares.v.push_back(t);    // pushed tape: recycle
    else
        delete t;                      // base tape or full pool
}

SpareTapes::~SpareTapes()
{
    for (Tape_* t : v)
        delete t;
}

}  // namespace

/*  Per-clause records from the most recent tape_eval_i, indexed by
 *  clause position in the evaluated tape.  Registers are reused, so
 *  the register file only holds each value until its last read -
 *  push decisions need the values AS OF each clause, which is
 *  exactly what these are (the MPR paper records its min/max choices
 *  during the interval pass for the same reason).
 *
 *  nan is the maybe-NaN taint: set when a domain error (sqrt/log of
 *  negatives, asin/acos out of range, 0/0, inf-inf, ...) could make
 *  pointwise values NaN inside the region even though the bounds
 *  look clean.  Pushes refuse to act on tainted records - domain-
 *  error bounds are exactly the garbage that makes pruning unsound
 *  (libfive tracks the same flag inside its Interval type).
 *  Interval VALUES are never altered, so culling is unchanged.  */
enum Choice : uint8_t { CHOICE_BOTH = 0, CHOICE_A, CHOICE_B };
struct ClauseIv
{
    Interval iv;
    uint8_t nan;
    uint8_t choice;
};

struct Deck_
{
    Tape_* base = nullptr;

    /*  Register-file layout: [constants | axes | registers].  The
     *  pinned prefix is written once (constants) or per-eval (axes);
     *  clause outputs live in num_regs linear-scan-allocated
     *  registers, so workspaces scale with peak liveness instead of
     *  tree size.  */
    unsigned num_pinned = 0;
    unsigned num_regs = 0;
    unsigned num_slots = 0;   // num_pinned + num_regs

    /*  (pinned slot, value) pairs, filled once per ctx  */
    std::vector<std::pair<uint32_t, float>> constants;

    /*  Pinned slots the eval drivers fill with coordinate data.
     *  Normally at most one per axis, but kept as vectors so a tree
     *  with duplicate axis nodes would still evaluate correctly. */
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

    /*  Register-taint mirror of i[] during interval evaluation  */
    std::vector<uint8_t> nan;

    /*  Per-clause records from the most recent tape_eval_i (sized to
     *  the base tape; the evaluated tape is never longer)  */
    std::vector<ClauseIv> civ;

    /*  Whether constant r-rows currently hold floats or derivatives
     *  (the g rows alias the r rows, exactly like Results.r)  */
    bool g_mode = false;

    /*  tape_push scratch: per-register liveness/boolean-context for
     *  the backward pass, per-clause verdicts for emission  */
    std::vector<uint8_t> live;
    std::vector<uint8_t> boolean_;
    std::vector<uint8_t> verdict;
};

////////////////////////////////////////////////////////////////////////////////

extern "C" Deck* deck_from_tree(MathTree* tree)
{
    Deck_* deck = new Deck_;

    unsigned total = tree->num_constants;
    for (unsigned level = 0; level < tree->num_levels; ++level)
        total += tree->active[level];

    /*  Pass 1: virtual slots.  Constants and axes take the pinned
     *  prefix; every clause output gets a unique virtual id above
     *  it (SSA-style), to be register-allocated in pass 2.  */
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
            switch (node->opcode)
            {
                case OP_X:  slot[node] = next;
                            deck->xs.push_back(next++);
                            deck->axes |= (1 << 2);
                            continue;
                case OP_Y:  slot[node] = next;
                            deck->ys.push_back(next++);
                            deck->axes |= (1 << 1);
                            continue;
                case OP_Z:  slot[node] = next;
                            deck->zs.push_back(next++);
                            deck->axes |= (1 << 0);
                            continue;
                default:    break;
            }
        }
    }
    uint32_t pinned = next;

    for (unsigned level = 0; level < tree->num_levels; ++level)
    {
        for (unsigned n = 0; n < tree->active[level]; ++n)
        {
            const Node* node = tree->nodes[level][n];
            if (node->opcode == OP_X || node->opcode == OP_Y ||
                node->opcode == OP_Z)
                continue;

            Clause c = {};
            c.op = node->opcode;
            c.out = next++;
            slot[node] = c.out;
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

    uint32_t root_virtual = slot.at(tree->head);

    /*  Affine-chain measurement (stats only): how many add/sub/mul-
     *  by-constant clauses feed another one - i.e. how much an
     *  affine-collapse pass (libfive tree.cpp:612) could fold.  */
    if (g_tape_stats.enabled)
    {
        const uint32_t n_const = deck->constants.size();
        std::vector<uint8_t> aff(next, 0);
        uint32_t links = 0;
        for (const Clause& c : clauses)
        {
            if (c.op != OP_ADD && c.op != OP_SUB && c.op != OP_MUL)
                continue;
            const bool a_const = c.a < n_const, b_const = c.b < n_const;
            if (!a_const && !b_const)
                continue;
            aff[c.out] = 1;
            if ((c.a < next && aff[c.a]) || (c.b < next && aff[c.b]))
                links++;
        }
        fprintf(stderr,
                "[tape-stats] affine: %u foldable chain links of %zu clauses\n",
                links, clauses.size());
    }

    /*  Affine collapse (STIBIUM_AFFINE=1, off by default): fold
     *  chains of add/sub-by-constant into one add/sub, and chains of
     *  mul-by-constant into one mul - the arithmetic every stacked
     *  move/scale transform emits (measured at ~35% of the merged
     *  Zeiss).  Opt-in because reassociating constants changes float
     *  rounding: results are equal in exact arithmetic but not bit-
     *  identical, and perceptual/print fidelity defaults belong to
     *  the person with the printer.  */
    static const bool affine_enabled = getenv("STIBIUM_AFFINE") != nullptr;
    if (affine_enabled && !clauses.empty())
    {
        const uint32_t n_const0 = deck->constants.size();

        /*  Constant values by virtual id (originals + ones we mint) */
        std::vector<float> cval(n_const0);
        for (uint32_t ci = 0; ci < n_const0; ++ci)
            cval[ci] = deck->constants[ci].second;
        std::unordered_map<uint32_t, uint32_t> cbits;  // bits -> id
        for (uint32_t ci = 0; ci < n_const0; ++ci)
        {
            uint32_t bits;
            memcpy(&bits, &cval[ci], 4);
            cbits.emplace(bits, ci);
        }
        std::vector<std::pair<uint32_t, float>> new_consts;
        auto is_const = [&](uint32_t v) { return v < n_const0 ||
                (v >= next && v - next < new_consts.size()); };
        auto const_of = [&](uint32_t v) -> float {
            return v < n_const0 ? cval[v] : new_consts[v - next].second; };
        auto mint = [&](float value) -> uint32_t {
            uint32_t bits;
            memcpy(&bits, &value, 4);
            auto found = cbits.find(bits);
            if (found != cbits.end())
                return found->second;
            const uint32_t id = next + new_consts.size();
            new_consts.emplace_back(id, value);
            cbits.emplace(bits, id);
            return id;
        };

        /*  value(out) = A*w + B with A in {1,-1} for add/sub chains,
         *  or value(out) = A*w (B unused) for mul chains  */
        enum AffClass : uint8_t { AFF_NONE, AFF_ADDSUB, AFF_MUL };
        struct Aff { uint8_t cls; float A, B; uint32_t w; };
        std::vector<Aff> aff(next, Aff{ AFF_NONE, 0, 0, 0 });

        std::vector<uint32_t> readers(next, 0);
        for (const Clause& c : clauses)
        {
            switch (arity(c.op))
            {
                case ARITY_GRID:    readers[c.m]++;  // fallthrough
                case ARITY_BINARY:  readers[c.b]++;  // fallthrough
                case ARITY_UNARY:   readers[c.a]++;  break;
                case ARITY_NONE:    break;
            }
        }
        readers[root_virtual] += 2;   // the root always has a reader

        uint32_t folded = 0;
        for (Clause& c : clauses)
        {
            if (c.op != OP_ADD && c.op != OP_SUB && c.op != OP_MUL)
                continue;
            const bool a_const = is_const(c.a), b_const = is_const(c.b);
            if (a_const == b_const)   // need exactly one constant side
                continue;
            const uint32_t x = a_const ? c.b : c.a;
            const float k = const_of(a_const ? c.a : c.b);

            /*  This clause as an affine form over x  */
            Aff f;
            if (c.op == OP_MUL)
                f = Aff{ AFF_MUL, k, 0, x };
            else if (c.op == OP_ADD)
                f = Aff{ AFF_ADDSUB, 1, k, x };
            else if (a_const)   // SUB(k, x)
                f = Aff{ AFF_ADDSUB, -1, k, x };
            else                // SUB(x, k)
                f = Aff{ AFF_ADDSUB, 1, -k, x };

            /*  Fold through x when x is a same-class link.  Shared
             *  links stay alive for their other readers (the sweep
             *  reaps them if every reader folds past); this clause
             *  stays a single clause either way, so the count never
             *  grows and dependency chains only get shorter.  */
            const Aff& fx = aff[x];
            if (fx.cls == f.cls && fx.cls != AFF_NONE)
            {
                if (f.cls == AFF_MUL)
                    f = Aff{ AFF_MUL, f.A * fx.A, 0, fx.w };
                else
                    f = Aff{ AFF_ADDSUB, f.A * fx.A,
                             f.A * fx.B + f.B, fx.w };
                readers[x]--;
                readers[fx.w]++;
                folded++;

                /*  Rewrite this clause to compute the folded form  */
                if (f.cls == AFF_MUL)
                {
                    c.op = OP_MUL;
                    c.a = f.w;
                    c.b = mint(f.A);
                }
                else if (f.A > 0 && f.B >= 0)
                {
                    c.op = OP_ADD;
                    c.a = f.w;
                    c.b = mint(f.B);
                }
                else if (f.A > 0)
                {
                    c.op = OP_SUB;
                    c.a = f.w;
                    c.b = mint(-f.B);
                }
                else
                {
                    c.op = OP_SUB;
                    c.a = mint(f.B);
                    c.b = f.w;
                }
            }
            aff[c.out] = f;
        }

        if (folded)
        {
            /*  Dead-code sweep: folded-through links lost their only
             *  reader.  Backward so chains die whole.  (Resized so
             *  minted-constant operand ids index safely.)  */
            readers.resize(next + new_consts.size(), 1);
            std::vector<Clause> live_clauses;
            std::vector<uint8_t> dead(clauses.size(), 0);
            for (uint32_t ki = clauses.size(); ki-- > 0;)
            {
                Clause& c = clauses[ki];
                if (readers[c.out] == 0)
                {
                    dead[ki] = 1;
                    switch (arity(c.op))
                    {
                        case ARITY_GRID:    readers[c.m]--;  // fallthrough
                        case ARITY_BINARY:  readers[c.b]--;  // fallthrough
                        case ARITY_UNARY:   readers[c.a]--;  break;
                        case ARITY_NONE:    break;
                    }
                }
            }
            for (uint32_t ki = 0; ki < clauses.size(); ++ki)
                if (!dead[ki])
                    live_clauses.push_back(clauses[ki]);
            clauses = std::move(live_clauses);

            /*  Renumber so minted constants join the pinned prefix:
             *  [old consts | new consts | axes | clause outs]  */
            std::unordered_map<uint32_t, uint32_t> renum;
            uint32_t nn = 0;
            for (uint32_t ci = 0; ci < n_const0; ++ci)
                renum[ci] = nn++;
            for (auto& nc : new_consts)
            {
                renum[nc.first] = nn;
                deck->constants.emplace_back(nn, nc.second);
                nn++;
            }
            auto shift = [&](std::vector<uint32_t>& ids) {
                for (auto& s : ids) { renum[s] = nn++; s = renum[s]; }
            };
            shift(deck->xs);
            shift(deck->ys);
            shift(deck->zs);
            for (uint32_t ci = 0; ci < n_const0; ++ci)
                deck->constants[ci].first = renum[ci];
            for (Clause& c : clauses)
            {
                renum.emplace(c.out, nn++);   // first definition wins
                c.out = renum.at(c.out);
                c.a = renum.count(c.a) ? renum.at(c.a) : c.a;
                c.b = renum.count(c.b) ? renum.at(c.b) : c.b;
                c.m = renum.count(c.m) ? renum.at(c.m) : c.m;
            }
            /* inputs always defined before use, so all were renumbered */
            root_virtual = renum.at(root_virtual);
            pinned = n_const0 + new_consts.size() + uint32_t(
                            deck->xs.size() + deck->ys.size() +
                            deck->zs.size());
            next = nn;
            if (g_tape_stats.enabled)
                fprintf(stderr, "[tape-stats] affine: folded %u links, "
                        "%zu clauses remain\n", folded, clauses.size());
        }
    }

    /*  Pass 2: linear-scan register allocation over the clause
     *  outputs.  last_use per virtual; a register frees the moment
     *  its value is read for the last time, so an input dying at a
     *  clause can donate its register to that clause's output
     *  (all evaluators are elementwise, so in-place is safe).
     *  Workspaces then scale with peak liveness, not tree size.  */
    const uint32_t N_INPUTS = 3;
    std::vector<uint32_t> last_use(next, 0);
    last_use[root_virtual] = UINT32_MAX;
    for (uint32_t k = 0; k < clauses.size(); ++k)
    {
        const uint32_t in[N_INPUTS] = { clauses[k].a, clauses[k].b,
                                        clauses[k].m };
        for (uint32_t j = 0; j < N_INPUTS; ++j)
            if (in[j] >= pinned && last_use[in[j]] != UINT32_MAX)
                last_use[in[j]] = k;
    }

    std::vector<uint32_t> reg(next, UINT32_MAX);
    std::vector<uint32_t> free_regs;
    uint32_t num_regs = 0;
    for (uint32_t k = 0; k < clauses.size(); ++k)
    {
        Clause& c = clauses[k];

        // Rewrite inputs first (they were allocated earlier), then
        // release the ones dying here so the output can reuse them
        uint32_t* in[N_INPUTS] = { &c.a, &c.b, &c.m };
        uint32_t dying[N_INPUTS];
        uint32_t n_dying = 0;
        for (uint32_t j = 0; j < N_INPUTS; ++j)
        {
            const uint32_t v = *in[j];
            if (v < pinned)
                continue;
            if (last_use[v] == k)
            {
                bool dup = false;
                for (uint32_t d = 0; d < n_dying; ++d)
                    dup |= (dying[d] == v);
                if (!dup)
                    dying[n_dying++] = v;
            }
            *in[j] = pinned + reg[v];
        }
        for (uint32_t d = 0; d < n_dying; ++d)
            free_regs.push_back(reg[dying[d]]);

        const uint32_t out_virtual = c.out;
        if (!free_regs.empty())
        {
            reg[out_virtual] = free_regs.back();
            free_regs.pop_back();
        }
        else
        {
            reg[out_virtual] = num_regs++;
        }
        c.out = pinned + reg[out_virtual];
    }

    deck->num_pinned = pinned;
    deck->num_regs = num_regs;
    deck->num_slots = pinned + num_regs;

    Tape_* base = new Tape_;
    base->clauses = std::move(clauses);
    base->root = root_virtual < pinned ? root_virtual
                                       : pinned + reg[root_virtual];
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
    g_tape_stats.record(0, base->clauses.size());
    if (g_tape_stats.enabled)
        fprintf(stderr, "[tape-stats] deck: %zu clauses, %u pinned, "
                "%u regs (workspace %u slots vs %zu unallocated)\n",
                base->clauses.size(), deck->num_pinned, deck->num_regs,
                deck->num_slots, base->clauses.size() + deck->num_pinned);

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
    ctx->boolean_.resize(deck->num_slots);
    ctx->nan.resize(deck->num_slots);
    ctx->civ.resize(deck->base->clauses.size());
    ctx->verdict.resize(deck->base->clauses.size());

    for (const auto& cv : deck->constants)
    {
        ctx->f[cv.first] = cv.second;
        ctx->i[cv.first] = Interval{ cv.second, cv.second };
        ctx->nan[cv.first] = std::isnan(cv.second);
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
            case OP_COPY:   v[c.out] = A; break;
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
    uint8_t* nan = ctx->nan.data();
    for (uint32_t s : ctx->deck->xs)
        { v[s] = X_i(X); nan[s] = std::isnan(X.lower) || std::isnan(X.upper); }
    for (uint32_t s : ctx->deck->ys)
        { v[s] = Y_i(Y); nan[s] = std::isnan(Y.lower) || std::isnan(Y.upper); }
    for (uint32_t s : ctx->deck->zs)
        { v[s] = Z_i(Z); nan[s] = std::isnan(Z.lower) || std::isnan(Z.upper); }

    ClauseIv* civ = ctx->civ.data();
    uint32_t k = 0;
    for (const auto& c : tape->clauses)
    {
        /*  Operand values and taints are captured before the write:
         *  registers are reused, and c.out may alias an input.  */
        const Interval A = v[c.a], B = v[c.b];
        const uint8_t na = nan[c.a], nb = nan[c.b];
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
            case OP_COPY:   v[c.out] = A; break;
            default: ;
        }

        /*  Maybe-NaN taint (see ClauseIv).  Conservative: any tainted
         *  input taints the output (even min/max, whose too-tight
         *  bounds can't be trusted when a side may be NaN), plus the
         *  per-op domain-error conditions from libfive's
         *  eval/interval.hpp, plus NaN leaking into the bounds
         *  themselves (inf-inf and friends).  */
        uint8_t u = na;
        switch (arity(c.op))
        {
            case ARITY_GRID:    u |= nan[c.m];  // fallthrough
            case ARITY_BINARY:  u |= nb;        break;
            case ARITY_NONE:    u = 0;          break;
            default:            break;
        }
        /*  Infinite bounds are NaN factories one op later (inf-inf,
         *  0*inf, inf/inf, sin/cos/tan of inf), and the bounds
         *  arithmetic hides it ([a,inf]-[b,inf] = [-inf,inf], no NaN
         *  in sight).  Taint conservatively.  */
        const bool a_inf = std::isinf(A.lower) || std::isinf(A.upper);
        const bool b_inf = std::isinf(B.lower) || std::isinf(B.upper);
        switch (c.op)
        {
            case OP_SQRT:
            case OP_LOG:
                u |= A.lower < 0;
                break;
            case OP_ASIN:
            case OP_ACOS:
                u |= A.lower < -1 || A.upper > 1;
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
                u |= (B.lower <= 0 && B.upper >= 0) || a_inf || b_inf;
                break;
            case OP_POW: {
                // Exact integer exponents are safe with negative
                // bases; real exponents make them NaN.  A base
                // touching zero with a non-positive exponent is
                // 0^0 / 0^negative either way.
                const bool int_exp = B.lower == B.upper &&
                                     B.lower == (float)(int)B.lower;
                u |= !int_exp && A.lower < 0;
                u |= A.lower <= 0 && A.upper >= 0 && B.lower <= 0;
                u |= a_inf || b_inf;
                break;
            }
            case OP_ATAN2:
                u |= A.lower <= 0 && A.upper >= 0 &&
                     B.lower <= 0 && B.upper >= 0;
                break;
            default:
                break;
        }
        u |= std::isnan(v[c.out].lower) || std::isnan(v[c.out].upper);
        nan[c.out] = u;

        /*  Per-clause record for tape_push: bounds, taint, and (for
         *  min/max) which side wins - decided here, while the operand
         *  values are still live in their registers.  Tainted
         *  operands veto the decision: a domain error's bounds are
         *  garbage, and fmin/fmax pass a lone NaN through to the
         *  OTHER side, so no branch can be dropped.  */
        uint8_t choice = CHOICE_BOTH;
        if (c.op == OP_MIN && !na && !nb)
        {
            if (A.upper <= B.lower)         choice = CHOICE_A;
            else if (B.upper <= A.lower)    choice = CHOICE_B;
        }
        else if (c.op == OP_MAX && !na && !nb)
        {
            if (A.lower >= B.upper)         choice = CHOICE_A;
            else if (B.lower >= A.upper)    choice = CHOICE_B;
        }
        civ[k].iv = v[c.out];
        civ[k].nan = u;
        civ[k].choice = choice;
        k++;
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
            case OP_COPY:
                if (R != A)
                    memcpy(R, A, c_ * sizeof(float));
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
            case OP_COPY:
                if (R != A)
                    memcpy(R, A, c_ * sizeof(derivative));
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
    auto& boolean_ = ctx->boolean_;
    const ClauseIv* civ = ctx->civ.data();
    uint8_t* verdict = ctx->verdict.data();

    enum : uint8_t { V_DEAD, V_KEEP, V_COPY_A, V_COPY_B, V_ELIDE, V_IMM };

    std::fill(live.begin(), live.end(), 0);
    std::fill(boolean_.begin(), boolean_.end(), 1);

    live[tape->root] = 1;
    bool changed = false;

    /*  Backward pass = reverse-topological liveness over the
     *  register file: live[r] means "the value the nearest preceding
     *  writer of r produces is still needed"; boolean_[r] carries
     *  the same pending-readers semantics for sign-only context.
     *  When a clause (a definition) is reached, both flags are
     *  consumed and reset for the next-earlier writer of that
     *  register.  Per-clause bounds/taint/choices come from the
     *  ClauseIv records of the most recent tape_eval_i (registers
     *  are reused, so the register file itself only holds each
     *  value until its last read).  */
    const uint32_t len = tape->clauses.size();
    for (uint32_t k = len; k-- > 0;)
    {
        const Clause& c = tape->clauses[k];
        if (!live[c.out])
        {
            verdict[k] = V_DEAD;
            continue;
        }
        const uint8_t boolctx = boolean_[c.out];
        live[c.out] = 0;
        boolean_[c.out] = 1;

        /*  disable_nodes_binary: a sign-determined value in boolean
         *  (min/max/neg) context collapses to its upper bound, and
         *  its subtree dies unless referenced elsewhere.  Tainted
         *  bounds (possible NaN inside the region) prove nothing.
         *  Note the STRICT lower test: a bound touching zero isn't a
         *  single sign class under negation (neg(+0) is still >= 0,
         *  but neg(upper) would go negative - the old disable_nodes_
         *  binary used >= and could flip exactly-zero fields).  */
        if (mode == TAPE_PUSH_BINARY && boolctx && !civ[k].nan &&
            (civ[k].iv.lower > 0 || civ[k].iv.upper < 0))
        {
            verdict[k] = V_IMM;
            changed = true;
            continue;
        }

        /*  disable_nodes: a decided min/max becomes a register copy
         *  of the surviving input (the exact value the min/max would
         *  return; readers can't be remapped to the survivor because
         *  its register may be redefined in between).  The choice
         *  was recorded during interval evaluation, taint-vetoed.
         *  When the allocator already placed the survivor in the
         *  output register, even the copy disappears.  */
        if ((c.op == OP_MIN || c.op == OP_MAX) &&
            civ[k].choice != CHOICE_BOTH)
        {
            const uint32_t keep = civ[k].choice == CHOICE_A ? c.a : c.b;
            verdict[k] = keep == c.out
                    ? V_ELIDE
                    : (civ[k].choice == CHOICE_A ? V_COPY_A : V_COPY_B);
            live[keep] = 1;
            boolean_[keep] &= boolctx;
            changed = true;
            continue;
        }

        /*  Kept clause: mark inputs live and propagate boolean
         *  context (min/max/neg/copy are sign-transparent; every
         *  other opcode breaks the chain).  */
        verdict[k] = V_KEEP;
        const bool transparent =
                (c.op == OP_MIN || c.op == OP_MAX || c.op == OP_NEG ||
                 c.op == OP_COPY);
        switch (arity(c.op))
        {
            case ARITY_GRID:
                live[c.m] = 1;
                boolean_[c.m] &= transparent ? boolctx : 0;
                // fallthrough
            case ARITY_BINARY:
                live[c.b] = 1;
                boolean_[c.b] &= transparent ? boolctx : 0;
                // fallthrough
            case ARITY_UNARY:
                live[c.a] = 1;
                boolean_[c.a] &= transparent ? boolctx : 0;
                break;
            case ARITY_NONE:
                break;
        }
    }

    if (!changed)
        return tape_retain(tape);

    Tape_* out = tape_alloc();
    out->clauses.reserve(tape->clauses.size());

    bool terminal = true;
    for (uint32_t k = 0; k < len; ++k)
    {
        const Clause& c = tape->clauses[k];
        switch (verdict[k])
        {
            case V_DEAD:
            case V_ELIDE:
                break;
            case V_IMM: {
                Clause n = {};
                n.op = OP_CONST;
                n.out = c.out;
                n.imm = civ[k].iv.upper;
                out->clauses.push_back(n);
                break;
            }
            case V_COPY_A:
            case V_COPY_B: {
                Clause n = {};
                n.op = OP_COPY;
                n.out = c.out;
                n.a = verdict[k] == V_COPY_A ? c.a : c.b;
                out->clauses.push_back(n);
                break;
            }
            case V_KEEP:
                if (c.op == OP_MIN || c.op == OP_MAX)
                    terminal = false;
                out->clauses.push_back(c);
                break;
        }
    }

    out->root = tape->root;
    out->terminal = terminal;
    out->X = X;
    out->Y = Y;
    out->Z = Z;
    out->has_bounds = true;
    out->depth = tape->depth + 1;
    out->parent = tape_retain(tape);

    g_tape_stats.record(out->depth, out->clauses.size());

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
        tape_dealloc(tape);
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

extern "C" Tape* tape_base_for_region(Tape* tape,
                                      const float xmin, const float xmax,
                                      const float ymin, const float ymax,
                                      const float zmin, const float zmax)
{
    while (tape->parent)
    {
        if (tape->has_bounds &&
            xmin >= tape->X.lower && xmax <= tape->X.upper &&
            ymin >= tape->Y.lower && ymax <= tape->Y.upper &&
            zmin >= tape->Z.lower && zmax <= tape->Z.upper)
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

extern "C" void tape_dump(const Tape* tape, const TapeCtx* ctx)
{
    fprintf(stderr, "tape %p: %zu clauses, root s%u%s\n",
            (const void*)tape, tape->clauses.size(), tape->root,
            tape->terminal ? " (terminal)" : "");
    for (const Clause& c : tape->clauses)
    {
        fprintf(stderr, "  s%-4u = %-8s", c.out, OPCODE_NAMES[c.op]);
        switch (arity(c.op))
        {
            case ARITY_GRID:
                fprintf(stderr, " s%u s%u s%u", c.a, c.b, c.m); break;
            case ARITY_BINARY:
                fprintf(stderr, " s%u s%u", c.a, c.b); break;
            case ARITY_UNARY:
                fprintf(stderr, " s%u", c.a); break;
            case ARITY_NONE:
                fprintf(stderr, " imm=%g", c.imm); break;
        }
        if (ctx)
            fprintf(stderr, "   i=[%g,%g]%s", ctx->i[c.out].lower,
                    ctx->i[c.out].upper,
                    ctx->nan[c.out] ? " NAN?" : "");
        fprintf(stderr, "\n");
    }
}

extern "C" bool tape_is_terminal(const Tape* tape)
{
    return tape->terminal;
}

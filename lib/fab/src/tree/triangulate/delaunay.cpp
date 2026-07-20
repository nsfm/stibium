/*
 *  Adaptive-Delaunay meshing, stages A-C (doc/MESH-NEXT.md).
 *  Blueprint: Keeter, "Please Steal My Meshing Algorithm Idea"
 *  (2026-07-03), annotated with this kernel's machinery: the octree
 *  descent runs on STANDARD-pushed tapes (pointwise-exact by the
 *  fuzzer's standard), interval culls hand the far field over as
 *  proven-sign box corners, and bisection batches through eval_r.
 */

#include <array>
#include <atomic>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstring>
#include <functional>
#include <iostream>
#include <mutex>
#include <sstream>
#include <thread>
#include <unordered_map>
#include <unordered_set>

#include <Eigen/Dense>
#include <Eigen/SVD>

#include "fab/tree/triangulate/delaunay.h"
#include "fab/util/switches.h"

namespace {

////////////////////////////////////////////////////////////////////////////
//  Stage A: point soup

/*  Bit-packed coordinate key: lattice coordinates come from the
 *  region's shared X/Y/Z arrays, so shared corners are bitwise
 *  identical across blocks (same trick as the mesher's vertex
 *  interning; -0 normalized).  */
uint64_t coord_hash(float x, float y, float z)
{
    if (x == 0)   x = 0;
    if (y == 0)   y = 0;
    if (z == 0)   z = 0;
    uint32_t u[3];
    memcpy(&u[0], &x, 4);
    memcpy(&u[1], &y, 4);
    memcpy(&u[2], &z, 4);
    uint64_t h = 0xcbf29ce484222325ull;
    for (uint32_t w : u)
        h = (h ^ w) * 0x100000001b3ull;
    return h ? h : 1;
}

/*  Canonical vertex ids under EXACT coordinate identity: the
 *  manifold pass deliberately mints index-distinct coincident
 *  vertices (one per sheet), and edges between two split sites
 *  read "open" in index space while the surface is sealed.  A
 *  hole is a geometric fact - count it on these welded ids.
 *  (Zeiss autod22 taught this the expensive way: "4 open" = 0
 *  geometric, and the retreat loop vandalized a healthy strip
 *  band chasing them - 26 real opens minted from 0.)  */
static std::vector<uint32_t> weld_ids(const std::vector<float>& V)
{
    const size_t nv = V.size() / 3;
    std::vector<uint32_t> remap(nv);
    std::unordered_map<uint64_t, uint32_t> canon;
    canon.reserve(nv);
    for (size_t v = 0; v < nv; ++v)
        remap[v] = canon.emplace(
                coord_hash(V[3*v], V[3*v + 1], V[3*v + 2]),
                uint32_t(v)).first->second;
    return remap;
}

struct PendingEdge
{
    /*  Bracket endpoints: a is strictly inside (f < 0), b is not  */
    float ax, ay, az;
    float bx, by, bz;
};

/*  A lattice cell that might straddle a sharp feature: its bounds
 *  and the indices (into the edge/surface-point arrays) of its
 *  crossing edges.  Collected during block sampling, judged after
 *  bisection when positions and normals exist.  */
struct FeatCell
{
    float lo[3], hi[3];
    uint64_t leaf_key = 0;   // coord_hash of the owning leaf's lo
                             // corner (stage-D drill-down address)
    uint32_t pts[12];
    uint8_t n = 0;
};

struct Collector
{
    const Deck* deck;
    TapeCtx* ctx;
    volatile int* halt;

    DSoup soup;
    std::unordered_set<uint64_t> seen_samples;
    std::unordered_map<uint64_t, uint32_t> edge_index;
    std::vector<PendingEdge> edges;
    std::vector<FeatCell> cells;
    float spacing = 0;   // finest lattice pitch (uniform grid)
    struct Census* census = nullptr;    // stage-D leaf census (opt-in)

    /*  Stage-D drill-down: pass 2 reads dense_map (leaf lo-corner
     *  hash -> extra levels); pass 1 writes want_dense from the
     *  QEF-residual and hidden-candidate triggers.  */
    const std::unordered_map<uint64_t, int>* dense_map = nullptr;
    std::unordered_map<uint64_t, int> want_dense;
    uint64_t phantom_rejected = 0;   // QEF corners the surface vetoed

    /*  Weld rollback: sites where a previous attempt's weld seam
     *  minted an open edge - no welding within 1 cell of these.  */
    const std::vector<std::array<float, 3>>* noweld = nullptr;

    /*  Every crease-suspect leaf's box, for flag dilation: a blend
     *  band crossing a leaf boundary must not change pitch ON the
     *  crease (the seam T-junctions read as chips), so flags spread
     *  to touching crease-suspect leaves before pass 2.  live is
     *  the leaf's surviving-pair count - the tangle gate reads it.  */
    /*  mindot: min pairwise dot of the 8 corner full-field
     *  gradients.  Anti-parallel corners (~ -1) mean opposing
     *  sheets in the leaf - but opposition alone is not danger:
     *  sep estimates the DISTANCE between the opposing sheets
     *  (Newton each corner of an anti-parallel pair onto its own
     *  sheet, project the difference onto the shared normal, keep
     *  the minimum).  A 2 mm bore-to-wall gap takes quarter-cell
     *  cores safely; a 0.3 mm web pinches.  HUGE_VALF when no
     *  anti-parallel pair exists.  */
    struct LeafBox
    {
        float lo[3], hi[3];
        unsigned live;
        float mindot;
        float sep;
    };
    std::unordered_map<uint64_t, LeafBox> crease_leaves;
    /*  Smooth (live = 0) surface leaves: recorded so the density
     *  dilation can SEE them - they are invisible to the flood
     *  (which walks crease_leaves only), leaving level-0 pockets
     *  inside promoted neighborhoods (the wireframe splotch
     *  class).  Box only; no survey.  */
    std::unordered_map<uint64_t, LeafBox> smooth_leaves;
    /*  Leaves the RESIDUAL formula flagged (flag_leaf) - the
     *  crowding grant's second signal; want_dense alone cannot
     *  distinguish residual-driven 2s from crowding-driven 2s.  */
    std::unordered_set<uint64_t> resid_hot;
    uint64_t tangle_suppressed = 0;
    /*  Hidden-feature oracle verdicts + (env-gated) verdict
     *  boxes for the STL dump.  */
    uint64_t hidden_feature = 0, hidden_graze = 0;
    uint64_t hidden_contact = 0;     // reach-collapse markers
    /*  Hidden-CELL witnesses (campaign 3/3): certified thin
     *  cells inside mixed leaves + proof-carrying bisection
     *  segments minted from them.  */
    uint64_t hidden_cell_feat = 0, hidden_cell_wit = 0;
    std::vector<std::array<float, 6>> hidden_feat_boxes;
    std::vector<std::array<float, 6>> hidden_graze_boxes;
    std::vector<std::array<float, 6>> hidden_contact_boxes;
    uint64_t curve_flagged = 0;      // chainless-curvature trigger
    uint64_t curve_seen = 0, curve_cross = 0;
    float curve_theta_max = 0.f;     // worst facet angle seen (deg)

    void add_sample(float x, float y, float z, bool inside,
                    bool on_surface = false)
    {
        if (seen_samples.insert(coord_hash(x, y, z)).second)
            soup.samples.push_back({ x, y, z, inside, on_surface });
    }

    static uint64_t edge_hash(float ax, float ay, float az,
                              float bx, float by, float bz)
    {
        return coord_hash(ax, ay, az) ^ (coord_hash(bx, by, bz) * 31);
    }

    /*  Returns the crossing's index, or UINT32_MAX if this edge has
     *  no sign change recorded.  */
    uint32_t find_edge(float ax, float ay, float az,
                       float bx, float by, float bz) const
    {
        const auto it = edge_index.find(edge_hash(ax, ay, az,
                                                  bx, by, bz));
        return it == edge_index.end() ? UINT32_MAX : it->second;
    }

    void add_edge(float ax, float ay, float az, bool a_in,
                  float bx, float by, float bz)
    {
        const uint64_t h = edge_hash(ax, ay, az, bx, by, bz);
        if (edge_index.count(h))
            return;
        edge_index.emplace(h, uint32_t(edges.size()));
        if (a_in)
            edges.push_back({ ax, ay, az, bx, by, bz });
        else
            edges.push_back({ bx, by, bz, ax, ay, az });
    }
};

/*  Corners of a culled region, sign proven by the interval bound -
 *  the far field costs nothing pointwise.  */
void add_box_corners(Collector& c, const Region& r, bool inside)
{
    for (int k = 0; k < 2; ++k)
        for (int j = 0; j < 2; ++j)
            for (int i = 0; i < 2; ++i)
                c.add_sample(r.X[i ? r.ni : 0],
                             r.Y[j ? r.nj : 0],
                             r.Z[k ? r.nk : 0], inside);
}

/*  Field-eval accounting (TIME=2): eval_points is the inner
 *  currency of every stage - sampling, oracle, referees, snap,
 *  repair.  One global tally answers "how much of the wall is
 *  the FIELD" in a single line.  */
struct EvalStats
{
    double secs = 0;
    uint64_t pts = 0, calls = 0;
};
static EvalStats g_eval;
static std::mutex g_eval_mx;

/*  Eval-side threading (P5, perf round 2, 2026-07-18): the
 *  TIME=2 map read 73 s of 221 s - one third of the r2 bino wall
 *  - as single-threaded field evaluation.  Octree subtrees are
 *  independent; tape_ctx-per-thread is designed-for.  Default
 *  hardware-2 (leave the OS and the renderer a lane);
 *  STIBIUM_DMESH_THREADS=1 reverts to serial.  */
int mesh_threads()
{
    static const char* env = getenv("STIBIUM_DMESH_THREADS");
    if (env)
        return std::max(1, atoi(env));
    /*  hw-1: leave the UI thread a lane (Nate's call, app
     *  integration 2026-07-18).  Env-only - not a UI knob.  */
    const unsigned hw = std::thread::hardware_concurrency();
    return hw > 1 ? int(hw - 1) : 1;
}

/*  Evaluate an arbitrary point list in MIN_VOLUME-sized batches on
 *  the given tape (same dummy-Region pattern as get_normals).  */
static void eval_range(const Tape* tape, TapeCtx* ctx,
                       const std::vector<float>& xs,
                       const std::vector<float>& ys,
                       const std::vector<float>& zs,
                       std::vector<float>& out,
                       size_t lo, size_t hi)
{
    for (size_t at = lo; at < hi; at += MIN_VOLUME)
    {
        const unsigned count =
                unsigned(hi - at < MIN_VOLUME ? hi - at
                                              : MIN_VOLUME);
        Region dummy = {};
        dummy.voxels = count;
        dummy.X = const_cast<float*>(xs.data() + at);
        dummy.Y = const_cast<float*>(ys.data() + at);
        dummy.Z = const_cast<float*>(zs.data() + at);
        const float* v = tape_eval_r(tape, ctx, dummy);
        memcpy(out.data() + at, v, count * sizeof(float));
    }
}

void eval_points(const Tape* tape, TapeCtx* ctx,
                 const std::vector<float>& xs,
                 const std::vector<float>& ys,
                 const std::vector<float>& zs,
                 std::vector<float>& out)
{
    const auto e0 = std::chrono::steady_clock::now();
    const size_t n = xs.size();
    out.resize(n);
    eval_range(tape, ctx, xs, ys, zs, out, 0, n);
    {
        std::lock_guard<std::mutex> lk(g_eval_mx);
        g_eval.secs += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - e0).count();
        g_eval.pts += n;
        ++g_eval.calls;
    }
}

int mesh_threads();

}  // namespace

/*  Progress sink (see delaunay.h).  Stage weights are the
 *  measured r1/r2 bino phase ratios (2026-07-18 TIME=2 map),
 *  normalized - an estimate by design; the app may also derive
 *  its own from the stable stage ids.  */
namespace {
struct ProgStage
{
    const char* name;
    float weight;
};
/*  Names are USER-FACING (the export dialog's status line reads
 *  them); the TIME instrument keeps its own terse pt.mark labels
 *  separately.  */
constexpr ProgStage PROG_STAGES[] = {
    { "starting",           0.00f },
    { "sampling",           0.08f },
    { "sampling detail",    0.10f },
    { "tracing features",   0.09f },
    { "verifying features", 0.04f },
    { "building lattice",   0.09f },
    { "placing surface",    0.17f },
    { "constraining",       0.01f },
    { "refining",           0.07f },
    { "repairing",          0.24f },
    { "splitting pinches",  0.01f },
    { "snapping creases",   0.03f },
    { "optimizing",         0.03f },
    { "finishing",          0.04f },
};
constexpr int PROG_N =
        int(sizeof(PROG_STAGES) / sizeof(PROG_STAGES[0]));
}  // namespace

DMeshProgress* dmesh_progress()
{
    static DMeshProgress p;
    return &p;
}

const char* dmesh_stage_name(int stage)
{
    return (stage >= 0 && stage < PROG_N)
            ? PROG_STAGES[stage].name : "?";
}

int dmesh_stage_count()
{
    return PROG_N;
}

namespace {

/*  Cumulative weight below a stage.  */
static float prog_base(int s)
{
    float b = 0;
    for (int i = 0; i < s && i < PROG_N; ++i)
        b += PROG_STAGES[i].weight;
    return b;
}

static void prog_emit()
{
    static const char* env = getenv("STIBIUM_DMESH_PROGRESS");
    if (!env || atoi(env) == 0)
        return;
    /*  Throttle to ~1 Hz; racy last-time is fine (worst case an
     *  extra line).  */
    static std::atomic<int64_t> last{ 0 };
    const int64_t now = std::chrono::duration_cast<
            std::chrono::milliseconds>(
            std::chrono::steady_clock::now().time_since_epoch())
            .count();
    const int64_t prev = last.load(std::memory_order_relaxed);
    if (now - prev < 1000)
        return;
    last.store(now, std::memory_order_relaxed);
    DMeshProgress* p = dmesh_progress();
    fprintf(stderr, "PROGRESS %s %.1f%% %.0f%%\n",
            dmesh_stage_name(p->stage.load()),
            100.f * p->overall.load(),
            100.f * p->in_stage.load());
}

static void prog_stage(int s)
{
    DMeshProgress* p = dmesh_progress();
    p->stage.store(s, std::memory_order_relaxed);
    p->in_stage.store(0, std::memory_order_relaxed);
    p->overall.store(s >= PROG_N - 1 ? 1.f : prog_base(s),
                     std::memory_order_relaxed);
    prog_emit();
}

static void prog_frac(float f)
{
    DMeshProgress* p = dmesh_progress();
    const int s = p->stage.load(std::memory_order_relaxed);
    f = f < 0 ? 0 : f > 1 ? 1 : f;
    p->in_stage.store(f, std::memory_order_relaxed);
    if (s < PROG_N)
        p->overall.store(prog_base(s) +
                         PROG_STAGES[s].weight * f,
                         std::memory_order_relaxed);
    prog_emit();
}

/*  Threaded batch eval on the BASE tape (perf round 2 increment
 *  2): the repair oracle, chip detection, segment referee, snap
 *  stash and Newton projection all push 100K+ point batches
 *  through eval_points serially.  Base-tape evals are pure per
 *  point, so ranges split freely; each worker gets its own ctx
 *  (pushed tapes belong to their owning ctx's arena and stay on
 *  the serial path).  Ranges align to MIN_VOLUME so batch
 *  boundaries match the serial run bit-for-bit.  */
void eval_points_mt(const Deck* deck, TapeCtx* ctx,
                    const std::vector<float>& xs,
                    const std::vector<float>& ys,
                    const std::vector<float>& zs,
                    std::vector<float>& out)
{
    const size_t n = xs.size();
    const size_t nt = std::min<size_t>(
            size_t(mesh_threads()),
            std::max<size_t>(1, n / (4 * MIN_VOLUME)));
    if (nt <= 1)
    {
        eval_points(deck_base(deck), ctx, xs, ys, zs, out);
        return;
    }
    const auto e0 = std::chrono::steady_clock::now();
    out.resize(n);
    const size_t chunks = (n + MIN_VOLUME - 1) / MIN_VOLUME;
    std::vector<std::thread> pool;
    for (size_t t = 0; t < nt; ++t)
        pool.emplace_back([&, t]() {
            TapeCtx* wctx = tape_ctx_new(deck);
            const size_t lo = std::min(n,
                    chunks * t / nt * MIN_VOLUME);
            const size_t hi = std::min(n,
                    chunks * (t + 1) / nt * MIN_VOLUME);
            if (lo < hi)
                eval_range(deck_base(deck), wctx, xs, ys, zs,
                           out, lo, hi);
            tape_ctx_free(wctx);
        });
    for (auto& th : pool)
        th.join();
    {
        std::lock_guard<std::mutex> lk(g_eval_mx);
        g_eval.secs += std::chrono::duration<double>(
                std::chrono::steady_clock::now() - e0).count();
        g_eval.pts += n;
        ++g_eval.calls;
    }
}

/*  Leaf block: evaluate the corner lattice, tag signs, and queue
 *  every sign-change lattice edge for bisection.  leaf_key is the
 *  leaf's drill-down address; crease_leaf says the pushed tape kept
 *  a live crease pair (gates the hidden-feature trigger - smooth
 *  tangent grazes read hidden too, and they never need density).  */
float box_dense_factor(const DSoup& soup, float x, float y, float z);

/*  Interval-subdivision hidden-feature oracle (cull -> certify
 *  move 1, doc/reviews/2026-07-18-interval-certify-reach-locus.md
 *  sec 5): at a leaf whose root interval straddles zero but whose
 *  lattice samples are sign-unanimous, subdivide with the SOUND
 *  bound the mesher already trusts for culling.  Interval
 *  arithmetic is conservative: a decided sign is a PROOF, so
 *  "every sub-box decided" proves no surface here (a tangent
 *  graze / loose bound - the class that killed the old blind
 *  HIDDEN flag), while a box still ambiguous at half the lattice
 *  pitch is a real sub-lattice feature and earns density.  The
 *  VERDICT RULE (refined by the first referee round: the
 *  stacked-cylinders tangent rings - the exact 18-leaf
 *  population that killed the blind flag - stay ambiguous at
 *  EVERY depth, because an internal tangency touches zero from
 *  one side and the interval can never decide it): a hidden
 *  feature means, by definition, that the lattice missed
 *  OPPOSITE-SIGN volume.  So FEATURE_PRESENT iff subdivision
 *  PROVES a sub-box of the sign opposite the lattice's
 *  unanimous read - a proof, so a flag is never wasted on
 *  contact loci.  Persistent single-sign ambiguity (tangency
 *  shells, loose bounds, walls far below any achievable
 *  density) resolves EMPTY: nothing meshable at the densities
 *  the trigger can buy.  Miss floor: features thinner than
 *  pitch/8 read EMPTY - 8x below lattice is unrecoverable at
 *  level 2 anyway.  Batched through tape_eval_i_batch (64
 *  boxes/pass; batch mode records no push state, so mid-
 *  pipeline use is safe); decided octants terminate
 *  immediately, so cost tracks the ambiguous shell.  */
bool autodense();
/*  contact out-param: true when ambiguity PERSISTED to the
 *  subdivision floor without deciding - a loose bound resolves
 *  once subdivided, a (near-)tangency never does.  Persistent
 *  ambiguity kissing the surface is a REACH-COLLAPSE marker
 *  (Nate's meshlab find, 2026-07-19: the one close-touching
 *  graze box on the off-axis bench sat on a 5.7x sliver
 *  hotspot).  */
static bool certify_hidden(const Tape* sub, TapeCtx* ctx,
                           const float lo[3], const float hi[3],
                           float pitch, bool lat_inside,
                           bool* contact,
                           std::vector<std::array<float, 3>>*
                                   proven = nullptr)
{
    struct HB { float lo[3], hi[3]; };
    std::vector<HB> work;
    work.push_back({ { lo[0], lo[1], lo[2] },
                     { hi[0], hi[1], hi[2] } });
    /*  Box budget: past it, give up WITHOUT a feature verdict -
     *  an opposite-sign proof never arrived, and unbounded
     *  digging in a tangency shell buys nothing.  */
    size_t budget = 4096;
    bool found = false;
    Interval Xs[TAPE_BATCH], Ys[TAPE_BATCH], Zs[TAPE_BATCH],
             out[TAPE_BATCH];
    while (!work.empty())
    {
        const int n = int(std::min(work.size(),
                                   size_t(TAPE_BATCH)));
        for (int i = 0; i < n; ++i)
        {
            const HB& b = work[work.size() - n + i];
            Xs[i] = { b.lo[0], b.hi[0] };
            Ys[i] = { b.lo[1], b.hi[1] };
            Zs[i] = { b.lo[2], b.hi[2] };
        }
        std::vector<HB> batch(work.end() - n, work.end());
        work.resize(work.size() - n);
        tape_eval_i_batch(sub, ctx, Xs, Ys, Zs, out, n);
        for (int i = 0; i < n; ++i)
        {
            /*  With a proven-box collector the search continues
             *  past the first proof (cap 4): the collector's
             *  centers seed proof-carrying WITNESS bisections -
             *  see the hidden-cell pass.  */
            if (out[i].lower > 0)       /* proven AIR */
            {
                if (lat_inside)
                {
                    if (!proven)
                        return true;    /* opposite sign: FEATURE */
                    found = true;
                    const HB& b = batch[i];
                    proven->push_back({
                            0.5f * (b.lo[0] + b.hi[0]),
                            0.5f * (b.lo[1] + b.hi[1]),
                            0.5f * (b.lo[2] + b.hi[2]) });
                    if (proven->size() >= 4)
                        return true;
                }
                continue;
            }
            if (out[i].upper < 0)       /* proven MATERIAL */
            {
                if (!lat_inside)
                {
                    if (!proven)
                        return true;    /* opposite sign: FEATURE */
                    found = true;
                    const HB& b = batch[i];
                    proven->push_back({
                            0.5f * (b.lo[0] + b.hi[0]),
                            0.5f * (b.lo[1] + b.hi[1]),
                            0.5f * (b.lo[2] + b.hi[2]) });
                    if (proven->size() >= 4)
                        return true;
                }
                continue;
            }
            const HB& b = batch[i];
            const float ex = b.hi[0] - b.lo[0];
            if (ex <= 0.125f * pitch || budget < 8)
            {
                *contact = true;        /* persisted to the floor */
                continue;
            }
            budget -= 8;
            for (int o = 0; o < 8; ++o)
            {
                HB h2;
                for (int a = 0; a < 3; ++a)
                {
                    const float mid = 0.5f * (b.lo[a] +
                                              b.hi[a]);
                    h2.lo[a] = (o >> a) & 1 ? mid : b.lo[a];
                    h2.hi[a] = (o >> a) & 1 ? b.hi[a] : mid;
                }
                work.push_back(h2);
            }
        }
    }
    return found;      /* (collector mode: any proof counts)
                          no opposite-sign volume: EMPTY  */
}

void sample_block(Collector& c, const Region& r, Tape* tape,
                  uint64_t leaf_key, bool crease_leaf)
{
    const uint32_t ni = r.ni, nj = r.nj, nk = r.nk;
    const uint32_t cx = ni + 1, cy = nj + 1, cz = nk + 1;
    const size_t n = size_t(cx) * cy * cz;

    std::vector<float> xs(n), ys(n), zs(n), vals;
    size_t q = 0;
    for (uint32_t k = 0; k < cz; ++k)
        for (uint32_t j = 0; j < cy; ++j)
            for (uint32_t i = 0; i < cx; ++i, ++q)
            {
                xs[q] = r.X[i];
                ys[q] = r.Y[j];
                zs[q] = r.Z[k];
            }
    eval_points(tape, c.ctx, xs, ys, zs, vals);

    /*  Clearance weld (sheet-separation campaign, 2026-07-17):
     *  the pinch anatomy read 160 of 183 non-manifold sites as AIR
     *  gaps at or below the first probe step - near-tangent
     *  contacts between assembled parts, far below the 0.1 mm
     *  product bar ("features >= 0.1 mm retained").  A sample
     *  whose air gap is thinner than the bar is flipped to solid,
     *  so the parts fuse cleanly at extraction - exactly what a
     *  printer would do with the clearance.  STIBIUM_DMESH_WELD
     *  sets the bar in model units; 0 disables.  */
    /*  Default OFF until the weld-induced fan tear is cured: with
     *  welding on, the bino referee trades nm 178 -> 13 (93% of
     *  the pinch class) but mints 4 open edges at ONE constrained
     *  vertex (-9.07, 47.64, 51.94) via repair dynamics - not
     *  spatially weld-adjacent, so the no-weld rollback can't
     *  reach it.  Next session's opening autopsy.  */
    static const char* weld_env = getenv("STIBIUM_DMESH_WELD");
    const float weld = weld_env ? float(atof(weld_env)) : 0.f;
    if (weld > 0)
    {
        std::vector<size_t> candv;
        for (size_t p = 0; p < n; ++p)
            if (vals[p] > 0 && vals[p] < weld)
                candv.push_back(p);
        if (!candv.empty())
        {
            const float h = 0.01f * c.spacing > 0
                    ? 0.01f * c.spacing : 0.01f;
            std::vector<float> px, py, pz, pv;
            for (const size_t p : candv)
                for (int q = 0; q < 7; ++q)
                {
                    px.push_back(xs[p] + (q==1?h:q==2?-h:0));
                    py.push_back(ys[p] + (q==3?h:q==4?-h:0));
                    pz.push_back(zs[p] + (q==5?h:q==6?-h:0));
                }
            eval_points(tape, c.ctx, px, py, pz, pv);
            std::vector<std::array<float, 3>> nd(candv.size());
            for (size_t i = 0; i < candv.size(); ++i)
            {
                float gx = pv[i*7+1] - pv[i*7+2];
                float gy = pv[i*7+3] - pv[i*7+4];
                float gz = pv[i*7+5] - pv[i*7+6];
                const float l = sqrtf(gx*gx + gy*gy + gz*gz);
                nd[i] = l > 0
                    ? std::array<float,3>{ gx/l, gy/l, gz/l }
                    : std::array<float,3>{ 0, 0, 0 };
            }
            /*  March both ways along the gap normal; solid on BOTH
             *  sides within the bar means the sample sits in a
             *  sub-bar clearance.  4 steps of weld/2 each way.  */
            constexpr int WK = 4;
            px.clear(); py.clear(); pz.clear();
            for (size_t i = 0; i < candv.size(); ++i)
                for (int s2 = 0; s2 < 2; ++s2)
                    for (int j = 1; j <= WK; ++j)
                    {
                        const float t = (s2 ? 1.f : -1.f) *
                                0.5f * weld * j;
                        const size_t p = candv[i];
                        px.push_back(xs[p] + t*nd[i][0]);
                        py.push_back(ys[p] + t*nd[i][1]);
                        pz.push_back(zs[p] + t*nd[i][2]);
                    }
            eval_points(tape, c.ctx, px, py, pz, pv);
            /*  Weld rollback (measured: an enclosure pre-rule
             *  strangled the weld to 1 sample - most gaps are one
             *  sample thin).  Flip aggressively; if a weld seam
             *  mints an open edge, the retreat loop feeds the
             *  site back as a no-weld zone and re-runs.  */
            for (size_t i = 0; i < candv.size(); ++i)
            {
                bool neg_lo = false, neg_hi = false;
                for (int j = 0; j < WK; ++j)
                {
                    if (pv[i*2*WK + j] < 0)
                        neg_lo = true;
                    if (pv[i*2*WK + WK + j] < 0)
                        neg_hi = true;
                }
                if (!(neg_lo && neg_hi))
                    continue;
                const size_t p = candv[i];
                bool banned = false;
                if (c.noweld)
                    for (const auto& q : *c.noweld)
                    {
                        const float dx = xs[p] - q[0];
                        const float dy = ys[p] - q[1];
                        const float dz = zs[p] - q[2];
                        if (dx*dx + dy*dy + dz*dz <
                            9.f * c.spacing * c.spacing)
                        {
                            banned = true;
                            break;
                        }
                    }
                if (banned)
                    continue;
                vals[p] = -vals[p];
                ++c.soup.welded;
            }
        }
    }

    const auto idx = [&](uint32_t i, uint32_t j, uint32_t k) {
        return (size_t(k) * cy + j) * cx + i;
    };
    const auto inside = [&](size_t p) { return vals[p] < 0; };

    bool any_in = false, any_out = false;
    for (size_t p = 0; p < n; ++p)
        (inside(p) ? any_in : any_out) = true;

    /*  Sample thinning (perf front P2, 2026-07-17): the DT only
     *  needs sign witnesses where signs can DISAGREE - refinement
     *  bisects inside<->outside tet edges, and the interval-culled
     *  far field already witnesses at 8 corners per box.  Deep
     *  lattice interiors contribute nothing but insert time and
     *  RAM (bino: insert 56% of wall).  Keep the sign-change band
     *  (both endpoints of every crossing edge) dilated by
     *  STIBIUM_DMESH_THIN rings of 26-neighbourhood shell, plus
     *  the 8 block corners everywhere (the culled-region witness
     *  granularity - an unwitnessed air gap between parts would
     *  let an all-inside tet span the air).  Exact-zero samples
     *  (on_surface) always survive.  -1 disables (dense legacy
     *  lattice).  DEFAULT 1 (Nate's eyes, 2026-07-18: "I couldn't
     *  spot any visible difference.  Worth it" - bino 170 -> 59 s,
     *  r2 unlocked at 198 s / 1.77 GB; the snap damage/churn
     *  referees are the enabling guard).  Dropping is per-block;
     *  keeping is the UNION
     *  across blocks (shared faces dedup in add_sample), so a
     *  crossing just over a block boundary keeps its shell from
     *  the block that sees it.  */
    static const char* thin_env = getenv("STIBIUM_DMESH_THIN");
    const int shell = thin_env ? atoi(thin_env) : 1;
    std::vector<uint8_t> keep;
    if (shell >= 0)
    {
        keep.assign(n, 0);
        if (any_in && any_out)
        {
            for (uint32_t k = 0; k < cz; ++k)
                for (uint32_t j = 0; j < cy; ++j)
                    for (uint32_t i = 0; i < cx; ++i)
                    {
                        const size_t p = idx(i, j, k);
                        const size_t nb[3] = {
                            i + 1 < cx ? idx(i + 1, j, k) : p,
                            j + 1 < cy ? idx(i, j + 1, k) : p,
                            k + 1 < cz ? idx(i, j, k + 1) : p,
                        };
                        for (const size_t pn : nb)
                            if (pn != p && inside(p) != inside(pn))
                                keep[p] = keep[pn] = 1;
                    }
            const auto near_kept = [&](uint32_t i, uint32_t j,
                                       uint32_t k) {
                const uint32_t i0 = i ? i - 1 : i,
                        i1 = i + 1 < cx ? i + 1 : i,
                        j0 = j ? j - 1 : j,
                        j1 = j + 1 < cy ? j + 1 : j,
                        k0 = k ? k - 1 : k,
                        k1 = k + 1 < cz ? k + 1 : k;
                for (uint32_t kk = k0; kk <= k1; ++kk)
                    for (uint32_t jj = j0; jj <= j1; ++jj)
                        for (uint32_t ii = i0; ii <= i1; ++ii)
                            if (keep[idx(ii, jj, kk)])
                                return true;
                return false;
            };
            std::vector<uint8_t> next;
            for (int ring = 0; ring < shell; ++ring)
            {
                next = keep;
                for (uint32_t k = 0; k < cz; ++k)
                    for (uint32_t j = 0; j < cy; ++j)
                        for (uint32_t i = 0; i < cx; ++i)
                            if (!keep[idx(i, j, k)] &&
                                near_kept(i, j, k))
                                next[idx(i, j, k)] = 1;
                keep.swap(next);
            }
        }
        for (int k = 0; k < 2; ++k)
            for (int j = 0; j < 2; ++j)
                for (int i = 0; i < 2; ++i)
                    keep[idx(i ? ni : 0, j ? nj : 0,
                             k ? nk : 0)] = 1;
    }

    /*  Thinning autopsy (STIBIUM_DMESH_THIN_DEBUG="x,y,z,r"):
     *  every lattice decision within r of the point, with class -
     *  the witness-level FPROBE.  */
    static const char* tdb_env = getenv("STIBIUM_DMESH_THIN_DEBUG");
    static float tdb[4] = { 0, 0, 0, -1 };
    if (tdb_env && tdb[3] < 0)
        sscanf(tdb_env, "%f,%f,%f,%f",
               tdb, tdb + 1, tdb + 2, tdb + 3);

    for (uint32_t k = 0; k < cz; ++k)
        for (uint32_t j = 0; j < cy; ++j)
            for (uint32_t i = 0; i < cx; ++i)
            {
                const size_t p = idx(i, j, k);
                if (tdb_env && tdb[3] > 0)
                {
                    const float ddx = xs[p] - tdb[0],
                            ddy = ys[p] - tdb[1],
                            ddz = zs[p] - tdb[2];
                    if (ddx*ddx + ddy*ddy + ddz*ddz <
                        tdb[3] * tdb[3])
                        fprintf(stderr, "THINDBG (%.4f, %.4f, "
                                "%.4f) f=%.5g %s block=%s%s\n",
                                xs[p], ys[p], zs[p], vals[p],
                                shell < 0 ? "legacy"
                                : (keep[p] || vals[p] == 0)
                                        ? "KEEP" : "DROP",
                                any_in && any_out ? "mixed"
                                : any_in ? "uni-in" : "uni-out",
                                box_dense_factor(c.soup, xs[p],
                                        ys[p], zs[p]) > 1
                                        ? "-dense" : "");
                }
                if (shell < 0 || keep[p] || vals[p] == 0)
                    c.add_sample(xs[p], ys[p], zs[p], inside(p),
                                 vals[p] == 0);
                else
                    ++c.soup.thinned;
                const size_t nb[3] = {
                    i + 1 < cx ? idx(i + 1, j, k) : p,
                    j + 1 < cy ? idx(i, j + 1, k) : p,
                    k + 1 < cz ? idx(i, j, k + 1) : p,
                };
                for (const size_t pn : nb)
                    if (pn != p && inside(p) != inside(pn))
                        c.add_edge(xs[p], ys[p], zs[p], inside(p),
                                   xs[pn], ys[pn], zs[pn]);
            }

    ++c.soup.leaf_blocks;
    if (!(any_in && any_out))
    {
        ++c.soup.hidden_candidates;   // interval said maybe; samples say no
        /*  Hidden trigger, now WITH the graze-vs-feature oracle
         *  it always needed (certify_hidden above; the 2026-07-16
         *  measurement killed the blind version because tangent
         *  grazes fired it - 18 leaves, 4x tris, zero depth).
         *  STIBIUM_DMESH_HIDDEN: 2 = interval oracle (DEFAULT),
         *  1 = legacy blind promote (the old opt-in), 0 = off.  */
        const char* hid_env = getenv("STIBIUM_DMESH_HIDDEN");
        const int hid = hid_env ? atoi(hid_env) : 2;
        if (hid == 1 && crease_leaf)
        {
            int& lv = c.want_dense[leaf_key];
            lv = std::max(lv, 2);
        }
        else if (hid >= 2 && !c.dense_map && autodense())
        {
            const float blo[3] = { r.X[0], r.Y[0], r.Z[0] };
            const float bhi[3] = { r.X[ni], r.Y[nj], r.Z[nk] };
            bool contact = false;
            const bool feat = certify_hidden(tape, c.ctx, blo,
                                             bhi, r.X[1] - r.X[0],
                                             any_in, &contact);
            if (feat)
            {
                ++c.hidden_feature;
                int& lv = c.want_dense[leaf_key];
                lv = std::max(lv, 2);
            }
            else if (contact)
            {
                ++c.hidden_contact;
                /*  Reach-collapse density (STIBIUM_DMESH_HIDDEN_
                 *  CONTACT=<level>, DEFAULT 2 - Nate's call
                 *  2026-07-19: sliver metrics 4.2 -> 1.2% at the
                 *  off-axis tangency, visually marginal but
                 *  cheap, "leave it enabled for the sake of our
                 *  metrics").  Density damps tangency damage;
                 *  the cure is strategy-doc #4.  0 disables.  */
                static const char* hc_env =
                        getenv("STIBIUM_DMESH_HIDDEN_CONTACT");
                const int hc = hc_env ? atoi(hc_env) : 2;
                if (hc > 0)
                {
                    int& lv = c.want_dense[leaf_key];
                    lv = std::max(lv, hc);
                }
            }
            else
                ++c.hidden_graze;
            /*  Contact boxes are ALWAYS recorded - they seed the
             *  #4a contact tracer; feature/graze boxes only under
             *  the dump env.  */
            if (contact && !feat)
                c.hidden_contact_boxes.push_back(
                        { blo[0], blo[1], blo[2],
                          bhi[0], bhi[1], bhi[2] });
            else if (getenv("STIBIUM_DMESH_DUMP_HIDDEN"))
                (feat ? c.hidden_feat_boxes
                      : c.hidden_graze_boxes).push_back(
                        { blo[0], blo[1], blo[2],
                          bhi[0], bhi[1], bhi[2] });
        }
    }

    /*  HIDDEN-CELL WITNESSES (density campaign 3 of 3, the
     *  dense-details autopsy 2026-07-20): sub-lattice-thin
     *  material hides INSIDE mixed leaves too - the screws head
     *  rim reads 136/136 lattice points OUT across its whole
     *  band even at level 3 (thickness < the local pitch), so
     *  the leaf-level hidden branch (which needs an all-quiet
     *  leaf) never sees it, and the surface there hangs from law
     *  chains alone: the chips-and-air-chords-where-details-
     *  crowd class.  The lattice cannot see the material, but
     *  the INTERVAL can prove it: for sign-unanimous cells near
     *  crossing cells (2-cell shell - the thin band always
     *  borders visible surface), certify with the subdivision
     *  oracle, COLLECTING proven opposite-sign sub-boxes; each
     *  proven center paired with a cell corner is a guaranteed
     *  sign-change segment, and the ordinary bisection queue
     *  turns those into exact surface witnesses.  Proof-carrying
     *  witnesses: never minted on grazes (the certify verdict
     *  rule), never guessed.  STIBIUM_DMESH_HIDDEN_CELLS=0
     *  disables.  Pass-2 dense leaves only (level >= 2): the
     *  class lives where density was already granted and the
     *  fine pitch makes the certify floor (pitch/8) deep.  */
    {
        static const char* hcw_env =
                getenv("STIBIUM_DMESH_HIDDEN_CELLS");
        const int hcw = hcw_env ? atoi(hcw_env) : 1;
        int dlevel = 0;
        if (hcw && c.dense_map)
        {
            const auto dit = c.dense_map->find(leaf_key);
            if (dit != c.dense_map->end())
                dlevel = dit->second;
        }
        if (hcw && dlevel >= 2 && (any_in || any_out))
        {
            const uint32_t mi = ni, mj = nj, mk = nk;
            const auto cidx = [&](uint32_t i, uint32_t j,
                                  uint32_t k) {
                return (size_t(k) * mj + j) * mi + i;
            };
            /*  Cell classes: mixed (some corner signs differ)
             *  vs unanimous; then a 2-cell dilation of mixed.  */
            std::vector<uint8_t> mixed(size_t(mi) * mj * mk, 0);
            for (uint32_t k = 0; k < mk; ++k)
                for (uint32_t j = 0; j < mj; ++j)
                    for (uint32_t i = 0; i < mi; ++i)
                    {
                        const bool s0 = inside(idx(i, j, k));
                        bool mix = false;
                        for (int o = 1; o < 8 && !mix; ++o)
                            mix = inside(idx(i + (o & 1),
                                             j + ((o >> 1) & 1),
                                             k + (o >> 2)))
                                  != s0;
                        mixed[cidx(i, j, k)] = mix;
                    }
            std::vector<std::array<uint32_t, 3>> shell;
            for (uint32_t k = 0; k < mk; ++k)
                for (uint32_t j = 0; j < mj; ++j)
                    for (uint32_t i = 0; i < mi; ++i)
                    {
                        if (mixed[cidx(i, j, k)])
                            continue;
                        bool near = false;
                        for (int dk2 = -2; dk2 <= 2 && !near;
                             ++dk2)
                            for (int dj2 = -2; dj2 <= 2 && !near;
                                 ++dj2)
                                for (int di2 = -2;
                                     di2 <= 2 && !near; ++di2)
                                {
                                    const int qi = int(i) + di2,
                                            qj = int(j) + dj2,
                                            qk = int(k) + dk2;
                                    if (qi < 0 || qj < 0 ||
                                        qk < 0 ||
                                        qi >= int(mi) ||
                                        qj >= int(mj) ||
                                        qk >= int(mk))
                                        continue;
                                    near = mixed[cidx(qi, qj,
                                                      qk)];
                                }
                        if (near)
                            shell.push_back({ i, j, k });
                    }
            /*  Batch-interval the shell; certify straddlers.  */
            Interval bXs[TAPE_BATCH], bYs[TAPE_BATCH],
                     bZs[TAPE_BATCH], bout[TAPE_BATCH];
            for (size_t s0 = 0; s0 < shell.size();
                 s0 += TAPE_BATCH)
            {
                const int n2 = int(std::min(
                        shell.size() - s0, size_t(TAPE_BATCH)));
                for (int q = 0; q < n2; ++q)
                {
                    const auto& s = shell[s0 + q];
                    bXs[q] = { xs[idx(s[0], s[1], s[2])],
                               xs[idx(s[0] + 1, s[1], s[2])] };
                    bYs[q] = { ys[idx(s[0], s[1], s[2])],
                               ys[idx(s[0], s[1] + 1, s[2])] };
                    bZs[q] = { zs[idx(s[0], s[1], s[2])],
                               zs[idx(s[0], s[1], s[2] + 1)] };
                }
                tape_eval_i_batch(tape, c.ctx, bXs, bYs, bZs,
                                  bout, n2);
                for (int q = 0; q < n2; ++q)
                {
                    if (!(bout[q].lower < 0 &&
                          bout[q].upper > 0))
                        continue;
                    const auto& s = shell[s0 + q];
                    const float clo[3] = { bXs[q].lower,
                                           bYs[q].lower,
                                           bZs[q].lower };
                    const float chi[3] = { bXs[q].upper,
                                           bYs[q].upper,
                                           bZs[q].upper };
                    const bool cin =
                            inside(idx(s[0], s[1], s[2]));
                    /*  MATERIAL ONLY (bino referee, 2026-07-20:
                     *  the symmetric version minted 127K
                     *  bisections and dragged worst 0.091 ->
                     *  0.372): a proven-air cell inside solid is
                     *  a near-tangent assembly GAP - the class
                     *  the clearance weld fuses on purpose
                     *  (sheet-separation verdict: pinches are
                     *  thin air; resolving them mints the
                     *  pinches back).  Hidden MATERIAL in air is
                     *  the rim/thin-wall class witnesses cure.  */
                    if (cin)
                        continue;
                    bool contact2 = false;
                    std::vector<std::array<float, 3>> proven;
                    if (!certify_hidden(tape, c.ctx, clo, chi,
                                        chi[0] - clo[0], cin,
                                        &contact2, &proven) ||
                        proven.empty())
                        continue;
                    ++c.hidden_cell_feat;
                    /*  Four tetrahedral corners: spread the
                     *  bisection directions around each proof.  */
                    static const int TC[4] = { 0, 3, 5, 6 };
                    for (const auto& P : proven)
                        for (const int o : TC)
                        {
                            const uint32_t pi = s[0] + (o & 1),
                                    pj = s[1] + ((o >> 1) & 1),
                                    pk = s[2] + (o >> 2);
                            const size_t cp = idx(pi, pj, pk);
                            c.add_edge(P[0], P[1], P[2], !cin,
                                       xs[cp], ys[cp], zs[cp]);
                            ++c.hidden_cell_wit;
                        }
                }
            }
        }
    }

    /*  Candidate feature cells: any voxel with >= 3 crossing edges.
     *  Judged later, once bisected positions and normals exist.  */
    for (uint32_t k = 0; k < nk; ++k)
        for (uint32_t j = 0; j < nj; ++j)
            for (uint32_t i = 0; i < ni; ++i)
            {
                FeatCell fc;
                fc.leaf_key = leaf_key;
                fc.lo[0] = r.X[i];    fc.hi[0] = r.X[i + 1];
                fc.lo[1] = r.Y[j];    fc.hi[1] = r.Y[j + 1];
                fc.lo[2] = r.Z[k];    fc.hi[2] = r.Z[k + 1];
                /*  12 edges: 4 per axis  */
                for (int axis = 0; axis < 3; ++axis)
                    for (int u = 0; u < 2; ++u)
                        for (int v = 0; v < 2; ++v)
                        {
                            float a[3] = { fc.lo[0], fc.lo[1], fc.lo[2] };
                            float b[3];
                            const int u_ax = (axis + 1) % 3;
                            const int v_ax = (axis + 2) % 3;
                            a[u_ax] = u ? fc.hi[u_ax] : fc.lo[u_ax];
                            a[v_ax] = v ? fc.hi[v_ax] : fc.lo[v_ax];
                            b[0] = a[0];  b[1] = a[1];  b[2] = a[2];
                            b[axis] = fc.hi[axis];
                            const uint32_t e = c.find_edge(
                                    a[0], a[1], a[2], b[0], b[1], b[2]);
                            if (e != UINT32_MAX && fc.n < 12)
                                fc.pts[fc.n++] = e;
                        }
                if (fc.n >= 3)
                    c.cells.push_back(fc);
            }
}

/*  Phase wall-clock under STIBIUM_DMESH_TIME=1: the profile that
 *  aims the performance pass (measure before optimizing - the
 *  house rule applies to speed too).  TIME=2 adds sub-stage
 *  lines (indented TIME2) - repair-round anatomy, snap
 *  internals, the fix-stage trio - without disturbing the
 *  level-1 phase deltas.  */
struct PhaseTimer
{
    const int level = [] {
        const char* e = getenv("STIBIUM_DMESH_TIME");
        return e ? std::max(1, atoi(e)) : 0;
    }();
    std::chrono::steady_clock::time_point t0 =
            std::chrono::steady_clock::now();
    std::chrono::steady_clock::time_point s0 = t0;
    void mark(const char* name)
    {
        if (!level)
            return;
        const auto t1 = std::chrono::steady_clock::now();
        fprintf(stderr, "TIME %-20s %8.2f s\n", name,
                std::chrono::duration<double>(t1 - t0).count());
        t0 = s0 = t1;
    }
    void sub(const char* name)
    {
        if (level < 2)
            return;
        const auto t1 = std::chrono::steady_clock::now();
        fprintf(stderr, "TIME2   %-24s %8.2f s\n", name,
                std::chrono::duration<double>(t1 - s0).count());
        s0 = t1;
    }
};


/*  Crease-band density (MESH-NEXT density round, measured
 *  2026-07-15): extra midpoint-refinement levels applied to
 *  crease-suspect leaves, and the factor by which the crease-local
 *  length scales (redundant-vertex corridor, chip repair keep-out)
 *  shrink with them.  A quality KNOB, default off: raw chip depth
 *  tracks the band pitch (0.500/0.375/0.192/0.098 sp at 1-8x), but
 *  the crease-snap pass zeroes the metric at EVERY density, and
 *  the band costs 2.2x triangles on crease-heavy models (csg 47K
 *  vs 22K at 96^3) for sub-threshold smoothness only.  Insert-
 *  repair can NOT substitute for the snap pass at any density: the
 *  blocked-repair residue is self-similar (keep-out and chip
 *  target both scale with local pitch - plateau 0.196/0.198/0.192
 *  sp at 1/2/4x, both repair strategies).  */
int dense_levels()
{
    /*  Read per call (cheap: once per leaf block) so tests can pin
     *  a density with setenv.  */
    const char* env = getenv("STIBIUM_DMESH_DENSE");
    const int lv = env ? atoi(env) : 0;
    /*  Cap 8 (256x): correctness review #6 - the level feeds
     *  1 << lv, and an extreme knob value is UB, not density.  */
    return lv > 0 ? std::min(lv, 8) : 0;
}

float dense_factor()
{
    return float(1 << dense_levels());
}

/*  Stage-D auto-density (STIBIUM_DMESH_AUTODENSE, default ON;
 *  STIBIUM_DMESH_DENSE > 0 overrides with the uniform manual knob).
 *  Trigger bar: QEF residual in cell units.  0.03 is the torus-lip
 *  referee bar itself - for near-parallel normals the clamped
 *  solve degenerates to a plane fit, so the residual IS the
 *  surface sagitta across the cell.  Level from magnitude:
 *  curvature sagitta quarters per level, so one level per factor
 *  of 4 over the bar (measured: sub-lattice lips read nr 0.2-0.3
 *  and need level 2; gentle blends read 0.03-0.12 and need 1).  */
constexpr float AUTOD_BAR = 0.03f;

/*  Level cap (STIBIUM_DMESH_AUTODENSE_MAX, default 2): the
 *  discriminator knob for damage bisection - cap 1 removes every
 *  level-2 tile (and with it both 4x-pitch seams and 4x-resolved
 *  thin walls) in one move.  */
int autod_max_level()
{
    /*  Default 3 was tried and REFUTED (2026-07-17): local
     *  eighth-cell cores on zeiss = 1,147 nm (7x), 2x triangles,
     *  3x wall-clock - the global-8x disaster in local form.
     *  Level 3 remains an opt-in knob for models whose sub-lattice
     *  detail justifies it (e.g. ~0.2 mm emboss at 1 vox/mm).  */
    const char* env = getenv("STIBIUM_DMESH_AUTODENSE_MAX");
    const int lv = env ? atoi(env) : 2;
    /*  Same 1 << lv foot-gun cap as dense_levels (review #6).  */
    return lv < 1 ? 1 : std::min(lv, 8);
}

/*  Crowding-grant ceiling (STIBIUM_DMESH_CROWD_MAX, default 3):
 *  the MAX=3 referee round (2026-07-19) split the verdict by
 *  TRIGGER - blanket level 3 wrecks organic models through the
 *  residual formula (bino: 1,555 of 1,665 leaves @3, worst
 *  0.093 -> 0.372, nm 6.4x) while the extreme-crowding grant is
 *  exactly the screws-class cure (26 leaves @3: tilt 4.03 ->
 *  3.04% of area, weighted tilt 0.947 -> 0.692 deg, worst
 *  held).  So the ceilings are SEPARATE: the residual formula
 *  stays behind AUTODENSE_MAX (default 2), and only live-pair
 *  crowding (>= 4x the bar) may take level 3.  */
int crowd_max_level()
{
    const char* env = getenv("STIBIUM_DMESH_CROWD_MAX");
    const int lv = env ? atoi(env)
                       : std::max(3, autod_max_level());
    return lv < 1 ? 1 : std::min(lv, 8);
}

/*  Tangle gate (zeiss autopsy, 2026-07-16): every open/non-manifold
 *  site in the auto-density zeiss sat in a leaf with 9-201 live
 *  pairs - thin-wall CSG tangles where a denser lattice resolves
 *  pinch geometry faster than it fixes divots (level 2 there:
 *  nm 241; level 1 patchy: 8 holes).  Every referee model that
 *  NEEDS density maxes out at 4-5 live pairs.  Leaves above the
 *  gate are not densified at all; their divots are the untraced-
 *  blend class (tracer-coverage work, not density work).  */
/*  STIBIUM_DMESH_TANGLE_DOT: leaves whose corner gradients contain
 *  a pair more anti-parallel than this are thin-wall tangles - no
 *  level-2 core (density resolves opposing sheets into pinches
 *  faster than it fixes divots).  Live-pair count was the first
 *  gate tried and REFUTED as the discriminator (zeiss collars and
 *  damage sites both read live 55-201; the real variable is
 *  opposing sheets, which fillet bands don't have).  -1 disables
 *  the gate.  */
float autod_tangle_dot()
{
    const char* env = getenv("STIBIUM_DMESH_TANGLE_DOT");
    return env ? float(atof(env)) : -0.5f;
}

/*  Wall-gap bar for the level-2 core gate, in cells (see the
 *  tangle comment at the flag site).  */
float autod_sep_bar()
{
    const char* env = getenv("STIBIUM_DMESH_SEP");
    return env ? float(atof(env)) : 0.7f;
}

/*  Crowding trigger bar (STIBIUM_DMESH_LIVE, live min/max pairs
 *  per leaf; 0 disables): leaves at or above it take level 2 from
 *  the survey directly - the signal for crowded sub-lattice
 *  geometry that the QEF residual cannot see.  Default 16: bench
 *  referees max 4-5 live, screw/eyepiece fuzz reads 51-75, a 10x
 *  gap.  */
int autod_live_bar()
{
    const char* env = getenv("STIBIUM_DMESH_LIVE");
    const int v = env ? atoi(env) : 16;
    return v <= 0 ? INT_MAX : v;
}

/*  Shallow-seed channel (STIBIUM_DMESH_SHALLOW, a normal dot;
 *  default 0.97 ~ 14 degrees, 0 disables): cells whose normal
 *  spread falls between feature grade (~25 degrees) and this bar
 *  mint TRACER-ONLY seeds.  The tracer verifies every seed
 *  against the oracle before tracing, so shallow candidates are
 *  free to offer - but they must never enter the feature tail:
 *  widening SPREAD_DOT itself to 0.95/0.97 poisons crossing
 *  suppression and the fallback radius graph ([.dchain] csg loop
 *  breaks; measured twice, 2026-07-17).  Below ~14 degrees the
 *  SSI corrector is chasing near-tangent intersections the det
 *  gate rejects anyway.  */
float shallow_seed_dot()
{
    const char* env = getenv("STIBIUM_DMESH_SHALLOW");
    return env ? float(atof(env)) : 0.97f;
}

bool autodense()
{
    const char* env = getenv("STIBIUM_DMESH_AUTODENSE");
    return (!env || atoi(env) != 0) && dense_levels() == 0;
}

/*  Chainless-curvature trigger bar (STIBIUM_DMESH_CURVEBAR,
 *  facet angle per CELL in degrees; 0 disables): the visible-
 *  quilting threshold.  Leaves whose projected level-set
 *  curvature turns more than the bar per cell promote level 1,
 *  more than twice the bar level 2 (ceiling 25 deg/cell =
 *  crease grade, law territory).  Density campaign trigger 2
 *  of 3.  */
float curve_bar_deg()
{
    const char* env = getenv("STIBIUM_DMESH_CURVEBAR");
    return env ? float(atof(env)) : 3.f;
}

/*  Local band factor at a point: the densest recorded drill-down
 *  box containing it (1 outside all bands).  Crease-local radii
 *  divide by max(global knob, this) - the density round measured
 *  that corridor and keep-out must track the LOCAL pitch or the
 *  band starves repair (281/281 blocked at fixed radii).  */
float box_dense_factor(const DSoup& soup, float x, float y, float z)
{
    float f = 1.f;
    for (const auto& b : soup.dense_boxes)
        if (x >= b.lo[0] && x <= b.hi[0] &&
            y >= b.lo[1] && y <= b.hi[1] &&
            z >= b.lo[2] && z <= b.hi[2])
            f = std::max(f, float(1 << b.level));
    return f;
}

/*  Midpoint-refine one lattice axis: n cells -> 2n, original corner
 *  coordinates preserved bit-exact (shared faces of mixed-density
 *  neighbours must dedup in add_sample).  */
void refine_axis(std::vector<float>& a)
{
    const size_t n = a.size() - 1;
    std::vector<float> b(2 * n + 1);
    for (size_t i = 0; i < n; ++i)
    {
        b[2 * i] = a[i];
        b[2 * i + 1] = 0.5f * (a[i] + a[i + 1]);
    }
    b[2 * n] = a[n];
    a = std::move(b);
}

/*  Stage-D leaf census (STIBIUM_DMESH_CENSUS=1 for a summary,
 *  =path for a per-leaf dump): the signals available at the
 *  density decision point, measured so the auto-density triggers
 *  are chosen from populations, not guesses.  Per crease-suspect
 *  leaf (live pair in the pushed tape):
 *    - SHEET-active: the pair's kink sheet {f_A = f_B} (or {g = 0}
 *      for abs) changes sign across the 8 leaf corners;
 *    - CREASE-active: sheet crossing AND the crease locus itself is
 *      near (min over corners of max(|f_A|,|f_B|) under the leaf
 *      diagonal - both fields small means the curve is close);
 *    - corner-gradient spread: min pairwise dot of the 8 normalized
 *      full-field gradients (a blend band curls the normal inside
 *      one leaf; a flat or gently curved patch does not).  */
struct Census
{
    bool on = false;
    FILE* dump = nullptr;
    uint64_t leaves = 0;        // crease-suspect (live pairs > 0)
    uint64_t by_sheet[4] = {};  // 0 / 1 / 2 / 3+ sheet-active pairs
    uint64_t by_crease[4] = {};
    uint64_t dot_hist[12] = {}; // min-dot: [-1,1] in 12 buckets,
                                // leaves with >= 1 crease-active pair

    /*  Per-cell QEF residual (feature_points): RMS distance of the
     *  cell's crossing-planes from the minimizer, in cell units.
     *  For near-parallel normals the clamped solve is a plane fit,
     *  so the residual IS the surface sagitta across the cell - the
     *  same units as the torus-lip referee bar (0.03 sp).  */
    static constexpr float NR_EDGES[17] = {
        0.005f, 0.01f, 0.02f, 0.03f, 0.05f, 0.075f, 0.1f, 0.15f,
        0.2f, 0.3f, 0.5f, 0.75f, 1.f, 1.5f, 2.f, 3.f, 5.f };
    uint64_t nr_hist[18] = {};
};

void leaf_census(Collector& c, const Region& r, Tape* sub)
{
    Census& cs = *c.census;
    const unsigned nmm = tape_pairs(sub, nullptr, 0);
    const unsigned nabs = tape_abs_pairs(sub, nullptr, 0);
    const unsigned np = nmm + nabs;
    if (!np)
        return;
    ++cs.leaves;

    std::vector<TapePair> pairs(np);
    tape_pairs(sub, pairs.data(), nmm);
    tape_abs_pairs(sub, pairs.data() + nmm, nabs);
    std::sort(pairs.begin(), pairs.end(),
              [](const TapePair& a, const TapePair& b) {
                  return a.clause < b.clause;
              });

    const float x0 = r.X[0], x1 = r.X[r.ni];
    const float y0 = r.Y[0], y1 = r.Y[r.nj];
    const float z0 = r.Z[0], z1 = r.Z[r.nk];
    const float diag = sqrtf((x1 - x0) * (x1 - x0) +
                             (y1 - y0) * (y1 - y0) +
                             (z1 - z0) * (z1 - z0));

    /*  Pair fields at the 8 corners: one recording walk each.  */
    std::vector<float> fa(8 * size_t(np)), fb(8 * size_t(np));
    float cx[8], cy[8], cz[8];
    for (int q = 0; q < 8; ++q)
    {
        cx[q] = (q & 1) ? x1 : x0;
        cy[q] = (q & 2) ? y1 : y0;
        cz[q] = (q & 4) ? z1 : z0;
        tape_eval_f_pairs(sub, c.ctx, cx[q], cy[q], cz[q],
                          pairs.data(), np,
                          fa.data() + size_t(q) * np,
                          fb.data() + size_t(q) * np);
    }

    unsigned nsheet = 0, ncrease = 0;
    float leaf_close = HUGE_VALF;
    for (unsigned i = 0; i < np; ++i)
    {
        bool pos = false, neg = false;
        float close = HUGE_VALF;
        for (int q = 0; q < 8; ++q)
        {
            const float A = fa[size_t(q) * np + i];
            const float B = fb[size_t(q) * np + i];
            if (!std::isfinite(A) || !std::isfinite(B))
                continue;
            const float d = pairs[i].is_max == 2 ? A : A - B;
            ((d < 0) ? neg : pos) = true;
            close = std::min(close, std::max(fabsf(A), fabsf(B)));
        }
        const bool sheet = pos && neg;
        const bool crease = sheet && close < diag;
        nsheet += sheet;
        ncrease += crease;
        if (crease)
            leaf_close = std::min(leaf_close, close);
    }
    ++cs.by_sheet[std::min(nsheet, 3u)];
    ++cs.by_crease[std::min(ncrease, 3u)];

    /*  Corner gradients of the full field (analytic g-mode walk):
     *  the spread says how far the normal turns inside one leaf.  */
    Region dummy = {};
    dummy.voxels = 8;
    dummy.X = cx;
    dummy.Y = cy;
    dummy.Z = cz;
    const derivative* g = tape_eval_g(sub, c.ctx, dummy);
    float mindot = 1.f;
    for (int a = 0; a < 8; ++a)
        for (int b = a + 1; b < 8; ++b)
        {
            const float la = sqrtf(g[a].dx * g[a].dx +
                                   g[a].dy * g[a].dy +
                                   g[a].dz * g[a].dz);
            const float lb = sqrtf(g[b].dx * g[b].dx +
                                   g[b].dy * g[b].dy +
                                   g[b].dz * g[b].dz);
            if (la < 1e-12f || lb < 1e-12f ||
                !std::isfinite(la) || !std::isfinite(lb))
                continue;
            const float dot = (g[a].dx * g[b].dx +
                               g[a].dy * g[b].dy +
                               g[a].dz * g[b].dz) / (la * lb);
            mindot = std::min(mindot, dot);
        }
    if (ncrease > 0)
    {
        const int bucket = std::min(11, std::max(0,
                int((mindot + 1.f) * 6.f)));
        ++cs.dot_hist[bucket];
    }

    if (cs.dump)
    {
        const float cell = r.X[1] - r.X[0];
        fprintf(cs.dump,
                "LEAF %.6g %.6g %.6g live=%u sheet=%u crease=%u "
                "mindot=%.4f close=%.4f\n",
                0.5f * (x0 + x1), 0.5f * (y0 + y1),
                0.5f * (z0 + z1), np, nsheet, ncrease, mindot,
                leaf_close == HUGE_VALF ? -1.f : leaf_close / cell);
    }
}

/*  Octree descent, mirroring the production mesher's structure:
 *  scalar interval eval + STANDARD push per level, cull on decided
 *  sign, leaf blocks small enough that their corner lattice fits a
 *  few eval_r batches.  */
constexpr uint64_t LEAF_VOXELS = 64;

/*  Principal curvatures of the level set at a point, from a
 *  19-point second-difference stencil on the field oracle (the
 *  shared dependency under valley snap, the churn-gate curvature
 *  carve-out, and the curvature-density trigger - built once per
 *  quality-research P1).  k = kappa_min with its tangent
 *  eigendirection (tx,ty,tz); kmax = kappa_max.  Sign
 *  convention: gradient is outward, so convex reads positive and
 *  a concave blend throat reads kappa_min < 0.  h is the stencil
 *  pitch - second differences divide by h^2, so it must sit well
 *  above eval noise (0.1 sp; 0.01 sp measured as pure static).
 *  False on degenerate gradient / non-finite readings.  */
struct KOut
{
    float k, kmax;
    float tx, ty, tz;   // kappa_min tangent eigendirection
    float sx, sy, sz;   // kappa_max tangent eigendirection
};
static const int K_SO[19][3] = {
    { 0, 0, 0 },
    { 1, 0, 0 }, { -1, 0, 0 }, { 0, 1, 0 }, { 0, -1, 0 },
    { 0, 0, 1 }, { 0, 0, -1 },
    { 1, 1, 0 }, { 1, -1, 0 }, { -1, 1, 0 }, { -1, -1, 0 },
    { 1, 0, 1 }, { 1, 0, -1 }, { -1, 0, 1 }, { -1, 0, -1 },
    { 0, 1, 1 }, { 0, 1, -1 }, { 0, -1, 1 }, { 0, -1, -1 },
};
static bool kout_from_vals(const float* v, float h, KOut* o);

/*  Facet angle per CELL at a leaf, from a PROJECTED Hessian
 *  probe (density campaign 2 of 3): leaf-level corner-normal
 *  spread cannot read one-sided curvature on thin geometry -
 *  a shade wall thinner than the leaf puts both sheets in every
 *  leaf and mindot pegs at -1 (measured on lamp).  Project the
 *  leaf center onto the surface (base tape - the point may exit
 *  the leaf box) and read the level set's own curvature there;
 *  kappa * cell_pitch is the facet angle the eye sees.
 *  Single-threaded evals only: descend runs inside worker
 *  threads.  */
static bool leaf_curve_theta(Collector& c, const Region& r,
                             float* deg)
{
    const Tape* base = deck_base(c.deck);
    const float cell = r.X[1] - r.X[0];
    float px = 0.5f * (r.X[0] + r.X[r.ni]),
          py = 0.5f * (r.Y[0] + r.Y[r.nj]),
          pz = 0.5f * (r.Z[0] + r.Z[r.nk]);
    const float clamp = 1.5f * cell * float(r.ni);
    const float h = 0.01f * c.spacing;
    for (int it = 0; it < 2; ++it)
    {
        std::vector<float> xs(7, px), ys(7, py), zs(7, pz), v;
        xs[1] += h;  xs[2] -= h;
        ys[3] += h;  ys[4] -= h;
        zs[5] += h;  zs[6] -= h;
        eval_points(base, c.ctx, xs, ys, zs, v);
        const float f = v[0];
        const float gx = (v[1]-v[2]) / (2*h),
                    gy = (v[3]-v[4]) / (2*h),
                    gz = (v[5]-v[6]) / (2*h);
        const float g2 = gx*gx + gy*gy + gz*gz;
        if (!(g2 > 0) || !std::isfinite(g2) || !std::isfinite(f))
            return false;
        float sx2 = -f * gx / g2, sy2 = -f * gy / g2,
              sz2 = -f * gz / g2;
        const float sl = sqrtf(sx2*sx2 + sy2*sy2 + sz2*sz2);
        if (sl > clamp)
            return false;   /* projection left the neighborhood */
        px += sx2;  py += sy2;  pz += sz2;
    }
    const float hh = 0.1f * c.spacing;
    std::vector<float> xs(19), ys(19), zs(19), v;
    for (int q = 0; q < 19; ++q)
    {
        xs[q] = px + K_SO[q][0]*hh;
        ys[q] = py + K_SO[q][1]*hh;
        zs[q] = pz + K_SO[q][2]*hh;
    }
    eval_points(base, c.ctx, xs, ys, zs, v);
    KOut kd;
    if (!kout_from_vals(v.data(), hh, &kd))
        return false;
    const float kmag = std::max(fabsf(kd.k), fabsf(kd.kmax));
    *deg = kmag * cell * (180.f / 3.14159265f);
    return true;
}

void descend(Collector& c, const Region& r, Tape* tape)
{
    if (*c.halt)
        return;

    const Interval X = { r.X[0], r.X[r.ni] },
                   Y = { r.Y[0], r.Y[r.nj] },
                   Z = { r.Z[0], r.Z[r.nk] };
    const Interval v = tape_eval_i(tape, c.ctx, X, Y, Z);
    if (v.lower > 0)
    {
        ++c.soup.culled_empty;
        add_box_corners(c, r, false);
        return;
    }
    if (v.upper < 0)
    {
        ++c.soup.culled_full;
        add_box_corners(c, r, true);
        return;
    }

    Tape* sub = tape_push(tape, c.ctx, X, Y, Z, TAPE_PUSH_STANDARD);
    if (r.voxels <= LEAF_VOXELS)
    {
        if (c.spacing == 0)
            c.spacing = r.X[1] - r.X[0];

        if (c.census)
            leaf_census(c, r, sub);

        const uint64_t key = coord_hash(r.X[0], r.Y[0], r.Z[0]);
        const unsigned live = tape_pairs(sub, nullptr, 0);
        const bool crease_leaf = live > 0;
        if (crease_leaf && !c.dense_map)
        {
            float cgx[8], cgy[8], cgz[8];
            for (int q = 0; q < 8; ++q)
            {
                cgx[q] = (q & 1) ? r.X[r.ni] : r.X[0];
                cgy[q] = (q & 2) ? r.Y[r.nj] : r.Y[0];
                cgz[q] = (q & 4) ? r.Z[r.nk] : r.Z[0];
            }
            Region dummy = {};
            dummy.voxels = 8;
            dummy.X = cgx;
            dummy.Y = cgy;
            dummy.Z = cgz;
            const derivative* g = tape_eval_g(sub, c.ctx, dummy);
            float mindot = 1.f;
            float sep = HUGE_VALF;
            for (int a2 = 0; a2 < 8; ++a2)
                for (int b2 = a2 + 1; b2 < 8; ++b2)
                {
                    const float la = sqrtf(g[a2].dx * g[a2].dx +
                                           g[a2].dy * g[a2].dy +
                                           g[a2].dz * g[a2].dz);
                    const float lb = sqrtf(g[b2].dx * g[b2].dx +
                                           g[b2].dy * g[b2].dy +
                                           g[b2].dz * g[b2].dz);
                    if (la < 1e-12f || lb < 1e-12f ||
                        !std::isfinite(la) || !std::isfinite(lb))
                        continue;
                    const float dot = (g[a2].dx * g[b2].dx +
                                       g[a2].dy * g[b2].dy +
                                       g[a2].dz * g[b2].dz) /
                                      (la * lb);
                    mindot = std::min(mindot, dot);
                    if (dot > autod_tangle_dot() ||
                        !std::isfinite(g[a2].v) ||
                        !std::isfinite(g[b2].v))
                        continue;
                    /*  Newton each corner onto its own sheet,
                     *  project the difference onto A's OUTWARD
                     *  normal.  The SIGN separates the two
                     *  anti-parallel configurations (they read
                     *  identically otherwise - measured): a thin
                     *  WALL has sheets facing away (gradients
                     *  diverge, signed gap negative) and pinches
                     *  under density; a concave GROOVE has flanks
                     *  facing each other (positive) and is exactly
                     *  what density cures.  Only wall gaps count.  */
                    const float ka = g[a2].v / (la * la);
                    const float kb = g[b2].v / (lb * lb);
                    const float ax = cgx[a2] - ka * g[a2].dx;
                    const float ay = cgy[a2] - ka * g[a2].dy;
                    const float az = cgz[a2] - ka * g[a2].dz;
                    const float bx = cgx[b2] - kb * g[b2].dx;
                    const float by = cgy[b2] - kb * g[b2].dy;
                    const float bz = cgz[b2] - kb * g[b2].dz;
                    const float s =
                            ((bx - ax) * g[a2].dx +
                             (by - ay) * g[a2].dy +
                             (bz - az) * g[a2].dz) / la;
                    if (s < 0)
                        sep = std::min(sep, -s);
                }
            c.crease_leaves[key] = {
                    { r.X[0], r.Y[0], r.Z[0] },
                    { r.X[r.ni], r.Y[r.nj], r.Z[r.nk] },
                    live, mindot, sep };
            if (c.census && c.census->dump)
                fprintf(c.census->dump,
                        "GATE %.6g %.6g %.6g live=%u mindot=%.4f "
                        "sep=%.4f\n",
                        0.5f * (r.X[0] + r.X[r.ni]),
                        0.5f * (r.Y[0] + r.Y[r.nj]),
                        0.5f * (r.Z[0] + r.Z[r.nk]),
                        live, mindot,
                        sep == HUGE_VALF ? -1.f : sep);
            /*  Crowding trigger (density campaign, 2026-07-18):
             *  the QEF residual is BLIND to crowded sub-lattice
             *  geometry - each cell's crossings self-explain
             *  (screw census: 284/292 cells silent below .005
             *  while the leaves carry live 51-75; bench referees
             *  max 4-5).  Live-pair count is the loudest
             *  separable signal the survey already computes.
             *  Level 2, same tangle demotion as every flag (the
             *  old live-count fear predates law + graduated
             *  retreat + geometric holes; the screws densify to
             *  0 open / 0 nm at r2 - measured).  */
            if (autodense() && live >= autod_live_bar())
            {
                const bool tangle =
                        sep < autod_sep_bar() * c.spacing;
                /*  Crowding asks for level 2 here; the level-3
                 *  grant is decided LATER (post-residual-survey,
                 *  see the crowding-grant block in
                 *  delaunay_sample) because the discriminator
                 *  needs both signals at once.  */
                int& lv = c.want_dense[key];
                lv = std::max(lv, tangle ? 1 : 2);
            }
            /*  Chainless-curvature trigger (density campaign 2
             *  of 3, 2026-07-18): cones quilt and curved walls
             *  carry random level-0 holes because gentle
             *  curvature is invisible to every other trigger -
             *  the QEF residual is near-zero on developable
             *  walls and live-pair crowding keys on clause
             *  count, not bending.  Leaf-level normal spread
             *  cannot see it either (thin shade walls put both
             *  sheets in one leaf, mindot pegs at -1 -
             *  measured); the honest signal is the level set's
             *  own curvature at a PROJECTED probe point.
             *  BAND: below the bar is visually flat; above 25
             *  degrees/cell is crease-grade (sub-cell fillets,
             *  law territory - the probe cannot tell a crease
             *  from a tight blend and both are handled
             *  elsewhere).  Tangles keep their demotion; the
             *  surface must CROSS the leaf.  */
            if (autodense() && curve_bar_deg() > 0)
            {
                bool pos = false, neg = false;
                for (int a2 = 0; a2 < 8; ++a2)
                    if (std::isfinite(g[a2].v))
                        ((g[a2].v < 0) ? neg : pos) = true;
                ++c.curve_seen;
                const bool tangle2 =
                        sep < autod_sep_bar() * c.spacing;
                if (pos && neg && !tangle2)
                {
                    ++c.curve_cross;
                    float th = 0;
                    if (leaf_curve_theta(c, r, &th))
                    {
                        c.curve_theta_max =
                                std::max(c.curve_theta_max, th);
                        if (th > curve_bar_deg() && th < 25.f)
                        {
                            const int want2 =
                                    th > 2 * curve_bar_deg()
                                            ? 2 : 1;
                            int& lv = c.want_dense[key];
                            lv = std::max(lv, want2);
                            ++c.curve_flagged;
                        }
                    }
                }
            }
        }

        else if (!crease_leaf && !c.dense_map && autodense())
        {
            /*  Box only - so the smooth-pocket fill can see
             *  leaves the flood cannot walk.  */
            c.smooth_leaves[key] = {
                    { r.X[0], r.Y[0], r.Z[0] },
                    { r.X[r.ni], r.Y[r.nj], r.Z[r.nk] },
                    0, 1.f, HUGE_VALF };
            /*  Chainless-curvature trigger, SMOOTH branch (the
             *  V-quilt autopsy, 2026-07-20): the trigger built
             *  for "curved walls carry random level-0 holes"
             *  only probed crease-suspect leaves - but a tilted
             *  cylinder wall is a live = 0 SMOOTH expanse, and
             *  its interior leaves are beyond the smooth-fill's
             *  reach too (the fill needs >= 3 touching promoted
             *  CREASE leaves; a big smooth region's interior
             *  touches none - measured at (-3.775, 32.14,
             *  40.98): 9 leaves, all smooth, 7 at level 0, the
             *  quilted V exactly).  Same probe, same band, same
             *  magnitude ladder as the crease branch;
             *  leaf_curve_theta already refuses leaves whose
             *  surface is out of reach.  */
            if (curve_bar_deg() > 0)
            {
                ++c.curve_seen;
                float th = 0;
                if (leaf_curve_theta(c, r, &th))
                {
                    ++c.curve_cross;
                    c.curve_theta_max =
                            std::max(c.curve_theta_max, th);
                    if (th > curve_bar_deg() && th < 25.f)
                    {
                        const int want2 =
                                th > 2 * curve_bar_deg() ? 2 : 1;
                        int& lv = c.want_dense[key];
                        lv = std::max(lv, want2);
                        ++c.curve_flagged;
                    }
                }
            }
        }

        /*  Crease-band density: the push rewrites every DECIDED
         *  min/max clause to a copy, so a min/max clause surviving
         *  in the leaf's pushed tape means this box straddles an
         *  unresolved crease choice.  STIBIUM_DMESH_DENSE > 0
         *  refines every such leaf uniformly (the manual knob);
         *  otherwise stage-D's pass-2 drill-down map picks the
         *  level per leaf from pass-1 evidence.  */
        int levels = 0;
        if (dense_levels() > 0)
            levels = crease_leaf ? dense_levels() : 0;
        else if (c.dense_map)
        {
            const auto it = c.dense_map->find(key);
            if (it != c.dense_map->end())
                levels = it->second;
        }

        if (levels > 0)
        {
            ++c.soup.dense_blocks;
            c.soup.dense_boxes.push_back({
                    { r.X[0], r.Y[0], r.Z[0] },
                    { r.X[r.ni], r.Y[r.nj], r.Z[r.nk] },
                    levels, key });
            std::vector<float> RX(r.X, r.X + r.ni + 1),
                               RY(r.Y, r.Y + r.nj + 1),
                               RZ(r.Z, r.Z + r.nk + 1);
            for (int lv = 0; lv < levels; ++lv)
            {
                refine_axis(RX);
                refine_axis(RY);
                refine_axis(RZ);
            }
            Region rd = r;
            rd.ni = uint32_t(RX.size() - 1);
            rd.nj = uint32_t(RY.size() - 1);
            rd.nk = uint32_t(RZ.size() - 1);
            rd.voxels = uint64_t(rd.ni) * rd.nj * rd.nk;
            rd.X = RX.data();
            rd.Y = RY.data();
            rd.Z = RZ.data();
            sample_block(c, rd, sub, key, crease_leaf);
        }
        else
            sample_block(c, r, sub, key, crease_leaf);
    }
    else
    {
        Region octants[8];
        const uint8_t split = octsect(r, octants);
        for (int i = 0; i < 8; ++i)
            if (split & (1 << i))
                descend(c, octants[i], sub);
    }
    tape_release(sub);
}

/*  Bisect all queued sign-change edges in lockstep batches: 16
 *  rounds, every midpoint evaluated on the base tape (pointwise
 *  values are pushed-tape-exact, so any covering tape gives the
 *  same bits; base keeps the batching trivial).  Invariant per
 *  round: f(a) < 0 <= f(b).  */
constexpr int BISECT_ROUNDS = 16;

static void bisect_range(const Deck* deck, TapeCtx* ctx,
                         volatile int* halt,
                         std::vector<PendingEdge>& edges,
                         size_t lo, size_t hi)
{
    const size_t n = hi - lo;
    Tape* base = deck_base(deck);
    std::vector<float> xs(n), ys(n), zs(n), vals;
    for (int round = 0; round < BISECT_ROUNDS; ++round)
    {
        if (*halt)
            return;
        for (size_t e = 0; e < n; ++e)
        {
            const PendingEdge& ed = edges[lo + e];
            xs[e] = (ed.ax + ed.bx) * 0.5f;
            ys[e] = (ed.ay + ed.by) * 0.5f;
            zs[e] = (ed.az + ed.bz) * 0.5f;
        }
        eval_points(base, ctx, xs, ys, zs, vals);
        for (size_t e = 0; e < n; ++e)
        {
            PendingEdge& ed = edges[lo + e];
            if (vals[e] < 0)
            {
                ed.ax = xs[e];
                ed.ay = ys[e];
                ed.az = zs[e];
            }
            else
            {
                ed.bx = xs[e];
                ed.by = ys[e];
                ed.bz = zs[e];
            }
        }
    }
}

void bisect_edges(Collector& c)
{
    const size_t n = c.edges.size();
    if (n == 0)
        return;
    /*  Disjoint ranges, own ctx each, byte-identical to serial
     *  (each edge's bisection depends on nothing but itself).  */
    const size_t nt = std::min<size_t>(
            size_t(mesh_threads()), std::max<size_t>(1, n / 4096));
    if (nt > 1)
    {
        std::vector<std::thread> pool;
        for (size_t t = 0; t < nt; ++t)
            pool.emplace_back([&, t]() {
                TapeCtx* ctx = tape_ctx_new(c.deck);
                bisect_range(c.deck, ctx, c.halt, c.edges,
                             n * t / nt, n * (t + 1) / nt);
                tape_ctx_free(ctx);
            });
        for (auto& th : pool)
            th.join();
    }
    else
        bisect_range(c.deck, c.ctx, c.halt, c.edges, 0, n);
    if (*c.halt)
        return;
    c.soup.surface.reserve(n);
    for (const PendingEdge& ed : c.edges)
        c.soup.surface.push_back({ (ed.ax + ed.bx) * 0.5f,
                                   (ed.ay + ed.by) * 0.5f,
                                   (ed.az + ed.bz) * 0.5f });
}

}  // namespace

namespace {

/*  Batched Newton projection onto the surface: p -= f * grad /
 *  |grad|^2, two rounds, gradients by 6-tap central differences.
 *  Steps are clamped so a bad gradient can't launch a point.  */
void project_points_impl(const Deck* deck, TapeCtx* ctx,
                    std::vector<float>& px, std::vector<float>& py,
                    std::vector<float>& pz, const std::vector<float>& hs,
                    const std::vector<float>& clamp)
{
    const size_t n = px.size();
    Tape* base = deck_base(deck);
    std::vector<float> xs, ys, zs, vals;
    for (int round = 0; round < 2; ++round)
    {
        xs.resize(n * 7);
        ys.resize(n * 7);
        zs.resize(n * 7);
        for (size_t i = 0; i < n; ++i)
        {
            const float h = hs[i];
            for (int q = 0; q < 7; ++q)
            {
                xs[i*7 + q] = px[i];
                ys[i*7 + q] = py[i];
                zs[i*7 + q] = pz[i];
            }
            xs[i*7 + 1] += h;   xs[i*7 + 2] -= h;
            ys[i*7 + 3] += h;   ys[i*7 + 4] -= h;
            zs[i*7 + 5] += h;   zs[i*7 + 6] -= h;
        }
        eval_points_mt(deck, ctx, xs, ys, zs, vals);
        for (size_t i = 0; i < n; ++i)
        {
            const float f = vals[i*7];
            const float gx = vals[i*7 + 1] - vals[i*7 + 2];
            const float gy = vals[i*7 + 3] - vals[i*7 + 4];
            const float gz = vals[i*7 + 5] - vals[i*7 + 6];
            const float h2 = 2 * hs[i];
            const float g2 = (gx*gx + gy*gy + gz*gz) / (h2 * h2);
            if (!(g2 > 0) || !std::isfinite(g2) || !std::isfinite(f))
                continue;
            float sx = -f * (gx / h2) / g2;
            float sy = -f * (gy / h2) / g2;
            float sz = -f * (gz / h2) / g2;
            const float sl = sqrtf(sx*sx + sy*sy + sz*sz);
            if (sl > clamp[i])
            {
                const float k = clamp[i] / sl;
                sx *= k;  sy *= k;  sz *= k;
            }
            px[i] += sx;
            py[i] += sy;
            pz[i] += sz;
        }
    }
}

/*  Hessian curvature-valley projection (quality P1, doc/reviews/
 *  2026-07-18-quality-research.md): the FAR snap skips are divots
 *  on UNTRACED tangent-blend seams - the clause system
 *  {f_A = f_B, f = 0} has a double root at G1 contact and cannot
 *  be marched, but the blend throat is still a sharp negative
 *  extremum of the minimum principal curvature, and curvature
 *  stays finite straight through the contact [Ohtake-Belyaev-
 *  Seidel 2004].  Project the chip midpoint onto that valley
 *  line: surface-project, then walk ACROSS the valley along the
 *  minimum-curvature eigendirection to the curvature extremum
 *  (parabola fit on kappa_min - the OBS04 extremality condition
 *  by differencing), re-project, iterate.  The result is a snap
 *  TARGET only - no constraint, no CCDT contact ("pure output
 *  surgery", same contract as SNAP_SURF).
 *  Returns false where no credible valley lives: kappa_min must
 *  read concave beyond a curvature floor (third-order quantities
 *  are noise on flats - the floor keeps the pass silent there).  */
/*  Raw Hessian + gradient from a 19-point stencil's values.
 *  H packed [xx, yy, zz, xy, xz, yz].  */
static bool hg_from_vals(const float* v, float h,
                         float H[6], float g[3])
{
    g[0] = (v[1] - v[2]) / (2*h);
    g[1] = (v[3] - v[4]) / (2*h);
    g[2] = (v[5] - v[6]) / (2*h);
    const float gl = sqrtf(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
    if (!(gl > 0) || !std::isfinite(gl))
        return false;
    H[0] = (v[1] - 2*v[0] + v[2]) / (h*h);
    H[1] = (v[3] - 2*v[0] + v[4]) / (h*h);
    H[2] = (v[5] - 2*v[0] + v[6]) / (h*h);
    H[3] = (v[7] - v[8] - v[9] + v[10]) / (4*h*h);
    H[4] = (v[11] - v[12] - v[13] + v[14]) / (4*h*h);
    H[5] = (v[15] - v[16] - v[17] + v[18]) / (4*h*h);
    return std::isfinite(H[0] + H[1] + H[2] + H[3] + H[4] +
                         H[5]);
}

/*  Normal curvature along a unit tangent direction c:
 *  kappa_c = (c^T H c) / |grad|.  */
static float kappa_along(const float H[6], const float g[3],
                         const float c[3])
{
    const float gl = sqrtf(g[0]*g[0] + g[1]*g[1] + g[2]*g[2]);
    const float q = c[0]*c[0]*H[0] + c[1]*c[1]*H[1] +
                    c[2]*c[2]*H[2] +
                    2*(c[0]*c[1]*H[3] + c[0]*c[2]*H[4] +
                       c[1]*c[2]*H[5]);
    return q / (gl > 0 ? gl : 1.f);
}

static bool kout_from_vals(const float* v, float h, KOut* o)
{
    float H[6], g[3];
    if (!hg_from_vals(v, h, H, g))
        return false;
    const float gx = g[0], gy = g[1], gz = g[2];
    const float gl = sqrtf(gx*gx + gy*gy + gz*gz);
    const float hxx = H[0], hyy = H[1], hzz = H[2],
                hxy = H[3], hxz = H[4], hyz = H[5];
    const float nx = gx/gl, ny = gy/gl, nz = gz/gl;
    /*  Tangent basis (u, w) perpendicular to the normal.  */
    float ux, uy, uz;
    if (fabsf(nx) <= fabsf(ny) && fabsf(nx) <= fabsf(nz))
        { ux = 0; uy = -nz; uz = ny; }
    else if (fabsf(ny) <= fabsf(nz))
        { ux = nz; uy = 0; uz = -nx; }
    else
        { ux = -ny; uy = nx; uz = 0; }
    const float ul = sqrtf(ux*ux + uy*uy + uz*uz);
    ux /= ul; uy /= ul; uz /= ul;
    const float wx = ny*uz - nz*uy,
                wy = nz*ux - nx*uz,
                wz = nx*uy - ny*ux;
    /*  Shape operator restricted to the tangent plane:
     *  S2 = [u w]^T H [u w] / |grad|.  */
    const float Hux = hxx*ux + hxy*uy + hxz*uz,
                Huy = hxy*ux + hyy*uy + hyz*uz,
                Huz = hxz*ux + hyz*uy + hzz*uz;
    const float Hwx = hxx*wx + hxy*wy + hxz*wz,
                Hwy = hxy*wx + hyy*wy + hyz*wz,
                Hwz = hxz*wx + hyz*wy + hzz*wz;
    const float a = (ux*Hux + uy*Huy + uz*Huz) / gl,
                b = (wx*Hux + wy*Huy + wz*Huz) / gl,
                c = (wx*Hwx + wy*Hwy + wz*Hwz) / gl;
    const float mean = 0.5f * (a + c);
    const float r = sqrtf(0.25f*(a - c)*(a - c) + b*b);
    o->k = mean - r;
    o->kmax = mean + r;
    /*  Eigenvector of the smaller eigenvalue, in (u, w)
     *  coordinates - pick the better-conditioned form.  */
    float e1 = b, e2 = o->k - a;
    if (fabsf(o->k - c) > fabsf(e2))
        { e1 = o->k - c; e2 = b; }
    const float el = sqrtf(e1*e1 + e2*e2);
    if (el > 0)
        { e1 /= el; e2 /= el; }
    else
        { e1 = 1; e2 = 0; }  /* umbilic: any direction */
    o->tx = e1*ux + e2*wx;
    o->ty = e1*uy + e2*wy;
    o->tz = e1*uz + e2*wz;
    /*  kappa_max's eigendirection: the tangent-plane
     *  perpendicular of kappa_min's.  */
    o->sx = -e2*ux + e1*wx;
    o->sy = -e2*uy + e1*wy;
    o->sz = -e2*uz + e1*wz;
    return true;
}

static bool curvature_probe(const Deck* deck, TapeCtx* ctx,
                            float h, float X, float Y, float Z,
                            KOut* o)
{
    std::vector<float> xs(19), ys(19), zs(19), v;
    for (int q = 0; q < 19; ++q)
    {
        xs[q] = X + K_SO[q][0]*h;
        ys[q] = Y + K_SO[q][1]*h;
        zs[q] = Z + K_SO[q][2]*h;
    }
    eval_points_mt(deck, ctx, xs, ys, zs, v);
    return kout_from_vals(v.data(), h, o);
}

static uint64_t g_vly_try = 0, g_vly_proj = 0, g_vly_floor = 0;
static float g_vly_kmin = 1e30f;
/*  ridge = false: follow kappa_min CONCAVE extremum lines (blend
 *  valleys).  ridge = true: follow kappa_max CONVEX extremum
 *  lines - the fading-crease class (a rim meeting a face at
 *  grazing incidence is a RIDGE whose sharpness fades; it has no
 *  valley to find).  */
static bool valley_project(const Deck* deck, TapeCtx* ctx,
                           float sp, float lsp,
                           float& X, float& Y, float& Z,
                           float kcoef = -1.f,
                           bool ridge = false,
                           float rprom = 0.f)
{
    ++g_vly_try;
    /*  Curvature floor in 1/sp units: fire only where the blend
     *  radius is under ~1/bar cells - beyond that the surface is
     *  visually flat and a tent is churn, not cure.  kcoef >= 0
     *  overrides the env coefficient (the contact tracer follows
     *  FADING creases - at a grazing-incidence rim the effective
     *  blend radius grows exactly where the defects live, and
     *  the default floor is deaf there).  */
    static const char* kb_env = getenv("STIBIUM_DMESH_VALLEY_K");
    const float kbar = (kcoef >= 0 ? kcoef
                        : kb_env ? float(atof(kb_env)) : 0.02f)
                       / std::max(sp, 1e-20f);
    /*  Second-difference stencils divide by h^2: at 0.01 sp the
     *  denominator sits at float-eval noise and kappa reads pure
     *  static (measured: 0 of 993 valley attempts survived).
     *  0.1 sp trades locality (well under any visible blend
     *  radius) for ~100x signal-to-noise.  */
    const float h = 0.1f * sp;
    float px2 = X, py2 = Y, pz2 = Z;
    for (int it = 0; it < 3; ++it)
    {
        /*  Pull onto the surface first (the chip midpoint starts
         *  inside the divot; later iterations re-project after
         *  the cross-valley step).  */
        std::vector<float> sx{ px2 }, sy{ py2 }, sz{ pz2 },
                hs{ 0.01f * sp }, cl{ 0.75f * sp };
        project_points_impl(deck, ctx, sx, sy, sz, hs, cl);
        if (!std::isfinite(sx[0] + sy[0] + sz[0]))
            { ++g_vly_proj; return false; }
        px2 = sx[0]; py2 = sy[0]; pz2 = sz[0];
        KOut k0;
        if (!curvature_probe(deck, ctx, h, px2, py2, pz2, &k0))
            { ++g_vly_proj; return false; }
        /*  Floor gate on the CENTER reading: concave for
         *  valleys, convex for ridges.  */
        if (it == 0 && !ridge)
            g_vly_kmin = std::min(g_vly_kmin, k0.k * sp);
        const float v0 = ridge ? k0.kmax : k0.k;
        if (ridge ? !(v0 > kbar) : !(v0 < -kbar))
            { ++g_vly_floor; return false; }
        const float cdx = ridge ? k0.sx : k0.tx,
                    cdy = ridge ? k0.sy : k0.ty,
                    cdz = ridge ? k0.sz : k0.tz;
        /*  Two cross-line probes; the span rides the local
         *  pitch - the extremum line is at most a chip-edge
         *  away.  */
        const float delta = 0.35f * lsp;
        KOut km, kp;
        if (!curvature_probe(deck, ctx, h,
                px2 - delta*cdx, py2 - delta*cdy,
                pz2 - delta*cdz, &km) ||
            !curvature_probe(deck, ctx, h,
                px2 + delta*cdx, py2 + delta*cdy,
                pz2 + delta*cdz, &kp))
            return false;
        /*  Parabola through (-d, vm) (0, v0) (+d, vp2):
         *  stationary point at d/2 (vm - vp2) / denom.  The
         *  denominator sign says whether the stationary point is
         *  the extremum this mode seeks (minimum of kappa_min
         *  for valleys, maximum of kappa_max for ridges);
         *  otherwise walk toward it.  */
        const float vm = ridge ? km.kmax : km.k;
        const float vp2 = ridge ? kp.kmax : kp.k;
        const float denom = vm - 2*v0 + vp2;
        float step;
        if (ridge)
        {
            /*  STRICT: a genuine ridge is a local MAXIMUM of
             *  kappa_max across the line.  A uniform cylinder
             *  wall passes the magnitude floor (kappa_max = 1/r
             *  everywhere) but has no extremum - the old
             *  walk-uphill branch chased probe noise across it
             *  (Nate's "curvy chaotic lines", measured).  No
             *  max, no ridge, march ends.  rprom > 0 (the tseed
             *  rung) additionally demands PROMINENCE over both
             *  shoulders: near a smeared curvature STEP the
             *  stencil overshoots into a shoulder bump barely
             *  above its plateau (measured on the off-axis lip -
             *  a wobbling false rail 0.5 sp inside the fillet);
             *  a true shallow-crossing spike clears its plateaus
             *  by the kink angle over the stencil width.  */
            if (!(v0 > vm + rprom && v0 > vp2 + rprom) ||
                !(denom < 0))
                return false;
            step = 0.5f * delta * (vm - vp2) / denom;
            step = step < -delta ? -delta
                 : step >  delta ?  delta : step;
        }
        else if (denom > 0)
        {
            step = 0.5f * delta * (vm - vp2) / denom;
            step = step < -delta ? -delta
                 : step >  delta ?  delta : step;
        }
        else
            step = (vp2 < vm ? 0.5f : -0.5f) * delta;
        px2 += step * cdx;
        py2 += step * cdy;
        pz2 += step * cdz;
    }
    /*  Final pull onto the surface.  */
    std::vector<float> sx{ px2 }, sy{ py2 }, sz{ pz2 },
            hs{ 0.01f * sp }, cl{ 0.75f * sp };
    project_points_impl(deck, ctx, sx, sy, sz, hs, cl);
    if (!std::isfinite(sx[0] + sy[0] + sz[0]))
        { ++g_vly_proj; return false; }
    X = sx[0]; Y = sy[0]; Z = sz[0];
    return true;
}

/*  Curvature-STEP corrector (tracer round 3 - the primitive
 *  verdict): a fillet boundary is two constant-curvature
 *  plateaus joined by a JUMP the stencil smears into a ramp -
 *  no curvature peak exists anywhere on the honest geometry
 *  (on-axis lip control: ridge mode 0/28).  The junction line
 *  is the ramp's STEEPEST point: sample normal curvature at
 *  five stations across the march frame (c = n x dir - never
 *  eigen ordering, which tilt rotates), take the derivative at
 *  three, and center on the strict peak of |d kappa / ds| via
 *  the parabola vertex.  Prominence bar: the plateau-to-plateau
 *  jump |kappa(+2d) - kappa(-2d)| must beat STIBIUM_DMESH_
 *  STEP_K / sp (default 0.05: fillets to ~20 cells).  */
static uint64_t g_stp[6] = { 0, 0, 0, 0, 0, 0 };
/*  g_stp: 0 proj, 1 hess, 2 cdir, 3 prominence, 4 peak, 5 ok  */
static bool step_project(const Deck* deck, TapeCtx* ctx,
                         float sp, const float dirv[3],
                         float& X, float& Y, float& Z)
{
    static const char* sk_env = getenv("STIBIUM_DMESH_STEP_K");
    const float kstep = (sk_env ? float(atof(sk_env)) : 0.05f)
                        / std::max(sp, 1e-20f);
    const float hh = 0.1f * sp;
    const float delta = 0.35f * sp;
    float px = X, py = Y, pz = Z;
    /*  Newton pull onto the surface.  */
    const auto project = [&](float& x, float& y, float& z) {
        const float h = 0.01f * sp;
        for (int it = 0; it < 2; ++it)
        {
            std::vector<float> xs(7, x), ys(7, y), zs(7, z), v;
            xs[1] += h;  xs[2] -= h;
            ys[3] += h;  ys[4] -= h;
            zs[5] += h;  zs[6] -= h;
            eval_points(deck_base(deck), ctx, xs, ys, zs, v);
            const float f = v[0];
            const float gx = (v[1]-v[2]) / (2*h),
                        gy = (v[3]-v[4]) / (2*h),
                        gz = (v[5]-v[6]) / (2*h);
            const float g2 = gx*gx + gy*gy + gz*gz;
            if (!(g2 > 0) || !std::isfinite(g2) ||
                !std::isfinite(f))
                return false;
            float sx2 = -f*gx/g2, sy2 = -f*gy/g2, sz2 = -f*gz/g2;
            const float sl = sqrtf(sx2*sx2 + sy2*sy2 + sz2*sz2);
            if (sl > 0.75f * sp)
                return false;
            x += sx2;  y += sy2;  z += sz2;
        }
        return true;
    };
    if (!project(px, py, pz))
        { ++g_stp[0]; return false; }
    /*  Center stencil: normal + the march-frame cross
     *  direction.  */
    std::vector<float> xs(5 * 19), ys(5 * 19), zs(5 * 19), v;
    for (int q = 0; q < 19; ++q)
    {
        xs[2*19 + q] = px + K_SO[q][0]*hh;
        ys[2*19 + q] = py + K_SO[q][1]*hh;
        zs[2*19 + q] = pz + K_SO[q][2]*hh;
    }
    {
        std::vector<float> cx(xs.begin() + 2*19,
                              xs.begin() + 3*19),
                cy(ys.begin() + 2*19, ys.begin() + 3*19),
                cz(zs.begin() + 2*19, zs.begin() + 3*19), cv;
        eval_points(deck_base(deck), ctx, cx, cy, cz, cv);
        float H[6], g[3];
        if (!hg_from_vals(cv.data(), hh, H, g))
            { ++g_stp[1]; return false; }
        const float gl = sqrtf(g[0]*g[0] + g[1]*g[1] +
                               g[2]*g[2]);
        const float n[3] = { g[0]/gl, g[1]/gl, g[2]/gl };
        float c[3] = { n[1]*dirv[2] - n[2]*dirv[1],
                       n[2]*dirv[0] - n[0]*dirv[2],
                       n[0]*dirv[1] - n[1]*dirv[0] };
        const float cl = sqrtf(c[0]*c[0] + c[1]*c[1] +
                               c[2]*c[2]);
        if (!(cl > 0.5f))
            { ++g_stp[2]; return false; }
        c[0] /= cl;  c[1] /= cl;  c[2] /= cl;
        /*  Five stations across; all stencils in one batch.  */
        for (int s2 = 0; s2 < 5; ++s2)
        {
            if (s2 == 2)
                continue;
            const float off = (s2 - 2) * delta;
            for (int q = 0; q < 19; ++q)
            {
                xs[s2*19 + q] = px + off*c[0] + K_SO[q][0]*hh;
                ys[s2*19 + q] = py + off*c[1] + K_SO[q][1]*hh;
                zs[s2*19 + q] = pz + off*c[2] + K_SO[q][2]*hh;
            }
        }
        eval_points(deck_base(deck), ctx, xs, ys, zs, v);
        float kap[5];
        for (int s2 = 0; s2 < 5; ++s2)
        {
            float H2[6], g2[3];
            if (!hg_from_vals(v.data() + s2*19, hh, H2, g2))
                { ++g_stp[1]; return false; }
            kap[s2] = kappa_along(H2, g2, c);
        }
        /*  Prominence: a real step jumps plateau to plateau.  */
        if (fabsf(kap[4] - kap[0]) < kstep)
            { ++g_stp[3]; return false; }
        /*  Derivative stations; sign-normalize; strict peak.  */
        const float gm = (kap[2] - kap[0]) / (2*delta),
                    g0 = (kap[3] - kap[1]) / (2*delta),
                    gp = (kap[4] - kap[2]) / (2*delta);
        const float sgn = g0 >= 0 ? 1.f : -1.f;
        const float Gm = gm*sgn, G0 = g0*sgn, Gp = gp*sgn;
        const float denom = Gm - 2*G0 + Gp;
        if (!(G0 > Gm) || !(G0 > Gp) || !(denom < 0))
            { ++g_stp[4]; return false; }
        float st = 0.5f * delta * (Gm - Gp) / denom;
        st = st < -delta ? -delta : st > delta ? delta : st;
        px += st * c[0];
        py += st * c[1];
        pz += st * c[2];
    }
    if (!project(px, py, pz))
        { ++g_stp[0]; return false; }
    ++g_stp[5];
    X = px;  Y = py;  Z = pz;
    return true;
}

/*  Shared tracer probes (contact tracer + tseed step tracer).
 *  Claim cells are 0.5 sp; a claimed 27-neighborhood means some
 *  law (crease, contact, or step) already owns the spot.  */
static inline uint64_t cc_key(float cell, float x, float y,
                              float z)
{
    return coord_hash(floorf(x / cell), floorf(y / cell),
                      floorf(z / cell));
}

static bool cc_claimed(const std::unordered_set<uint64_t>& s,
                       float cell, float x, float y, float z)
{
    for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz)
                if (s.count(cc_key(cell, x + dx*cell,
                                   y + dy*cell, z + dz*cell)))
                    return true;
    return false;
}

/*  Outward normal via a 7-point gradient probe.  */
static bool probe_normal7(const Deck* deck, TapeCtx* ctx,
                          float sp, float x, float y, float z,
                          float n[3])
{
    const float h = 0.01f * sp;
    std::vector<float> xs(7, x), ys(7, y), zs(7, z), v;
    xs[1] += h;  xs[2] -= h;
    ys[3] += h;  ys[4] -= h;
    zs[5] += h;  zs[6] -= h;
    eval_points(deck_base(deck), ctx, xs, ys, zs, v);
    n[0] = v[1] - v[2];
    n[1] = v[3] - v[4];
    n[2] = v[5] - v[6];
    const float l = sqrtf(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (!(l > 0) || !std::isfinite(l))
        return false;
    n[0] /= l;  n[1] /= l;  n[2] /= l;
    return true;
}

/*  Sign-cross filter (the buried-sheet lesson, 2026-07-19:
 *  coincident CSG faces leave INTERNAL f = 0 sheets - the
 *  off-axis bench carries a whole buried disk where its two
 *  cylinders meet, and the tracer faithfully marched it
 *  instead of the defect rim).  A REAL surface point has
 *  opposite field signs a step along +/- its normal; an
 *  interior zero-sheet touches zero from one side only.  */
static bool probe_real_surface(const Deck* deck, TapeCtx* ctx,
                               float sp, float x, float y,
                               float z)
{
    static int dbg_sheet = 0;
    float n[3];
    if (!probe_normal7(deck, ctx, sp, x, y, z, n))
        return false;
    const float h = 0.25f * sp;
    std::vector<float> xs{ x + h*n[0], x - h*n[0] },
            ys{ y + h*n[1], y - h*n[1] },
            zs{ z + h*n[2], z - h*n[2] }, v;
    eval_points(deck_base(deck), ctx, xs, ys, zs, v);
    const bool ok = std::isfinite(v[0]) &&
            std::isfinite(v[1]) && v[0] > 0 && v[1] < 0;
    if (!ok && dbg_sheet < 6 &&
        getenv("STIBIUM_DMESH_CHIP_DEBUG"))
    {
        ++dbg_sheet;
        fprintf(stderr, "  SHEETDBG land=(%.3f,%.3f,%.3f) "
                "n=(%.2f,%.2f,%.2f) f(+h)=%.4f f(-h)=%.4f\n",
                x, y, z, n[0], n[1], n[2], v[0], v[1]);
    }
    return ok;
}

/*  Two Newton pulls onto f = 0 (seed tightening).  */
static bool newton7(const Deck* deck, TapeCtx* ctx, float sp,
                    float& x, float& y, float& z)
{
    const float h = 0.01f * sp;
    for (int it = 0; it < 2; ++it)
    {
        std::vector<float> xs(7, x), ys(7, y), zs(7, z), v;
        xs[1] += h;  xs[2] -= h;
        ys[3] += h;  ys[4] -= h;
        zs[5] += h;  zs[6] -= h;
        eval_points(deck_base(deck), ctx, xs, ys, zs, v);
        const float f = v[0];
        const float gx = (v[1]-v[2]) / (2*h),
                    gy = (v[3]-v[4]) / (2*h),
                    gz = (v[5]-v[6]) / (2*h);
        const float g2 = gx*gx + gy*gy + gz*gz;
        if (!(g2 > 0) || !std::isfinite(g2) ||
            !std::isfinite(f))
            return false;
        x -= f*gx/g2;
        y -= f*gy/g2;
        z -= f*gz/g2;
    }
    return true;
}

/*  Tangential gradient of MEAN curvature at an on-surface
 *  point: mean curvature is basis-free (trace of the shape
 *  operator - immune to the eigen-ordering swap that killed
 *  ridge mode under tilt), and its surface gradient points
 *  ACROSS a curvature step, waking exactly where two plateaus
 *  meet.  n and c0 (unnormalized, magnitude in c0l) out.  */
static bool mc_tgrad(const Deck* deck, TapeCtx* ctx, float sp,
                     float x, float y, float z,
                     float n[3], float c0[3], float& c0l)
{
    KOut kc, ku, kw;
    if (!probe_normal7(deck, ctx, sp, x, y, z, n) ||
        !curvature_probe(deck, ctx, 0.1f * sp, x, y, z, &kc))
        return false;
    float u[3];
    if (fabsf(n[0]) <= fabsf(n[1]) &&
        fabsf(n[0]) <= fabsf(n[2]))
        { u[0] = 0; u[1] = -n[2]; u[2] = n[1]; }
    else if (fabsf(n[1]) <= fabsf(n[2]))
        { u[0] = n[2]; u[1] = 0; u[2] = -n[0]; }
    else
        { u[0] = -n[1]; u[1] = n[0]; u[2] = 0; }
    const float ul = sqrtf(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
    u[0] /= ul;  u[1] /= ul;  u[2] /= ul;
    const float w[3] = { n[1]*u[2] - n[2]*u[1],
                         n[2]*u[0] - n[0]*u[2],
                         n[0]*u[1] - n[1]*u[0] };
    const float dp = 0.35f * sp;
    if (!curvature_probe(deck, ctx, 0.1f * sp,
            x + dp*u[0], y + dp*u[1], z + dp*u[2], &ku) ||
        !curvature_probe(deck, ctx, 0.1f * sp,
            x + dp*w[0], y + dp*w[1], z + dp*w[2], &kw))
        return false;
    const float m0 = 0.5f * (kc.k + kc.kmax);
    const float du = 0.5f * (ku.k + ku.kmax) - m0;
    const float dw = 0.5f * (kw.k + kw.kmax) - m0;
    c0[0] = du*u[0] + dw*w[0];
    c0[1] = du*u[1] + dw*w[1];
    c0[2] = du*u[2] + dw*w[2];
    c0l = sqrtf(c0[0]*c0[0] + c0[1]*c0[1] + c0[2]*c0[2]);
    return std::isfinite(c0l);
}

/*  March an extremum line (mode 0 valley, 1 curvature step,
 *  2 ridge) from (px,py,pz) along (dx,dy,dz), claiming cells as
 *  it goes.  Returns true when the march CLOSED a loop (came
 *  back within a step of the start after real progress) - the
 *  closure TODO from round 2; without it a ring march laps
 *  itself.  */
static bool march_xline(const Deck* deck, TapeCtx* ctx,
        float sp, float px, float py, float pz,
        float dx, float dy, float dz, int mode,
        std::unordered_set<uint64_t>* claimed,
        std::vector<std::array<float, 3>>* pts,
        float rprom = 0.f)
{
    const float step = 0.5f * sp;
    const float cell = 0.5f * sp;
    const float sx0 = px, sy0 = py, sz0 = pz;
    for (int i = 0; i < 300; ++i)
    {
        float nx2 = px + step*dx, ny2 = py + step*dy,
              nz2 = pz + step*dz;
        const float d0[3] = { dx, dy, dz };
        /*  0.002: contact traces follow FADING creases -
         *  the default valley floor is deaf exactly at
         *  grazing incidence.  */
        const bool ok2 = mode == 1
                ? step_project(deck, ctx, sp, d0,
                               nx2, ny2, nz2)
                : valley_project(deck, ctx, sp, sp,
                                 nx2, ny2, nz2, 0.002f,
                                 mode == 2, rprom);
        if (!ok2)
            break;
        if (!probe_real_surface(deck, ctx, sp, nx2, ny2, nz2))
            break;      /* wandered onto a buried sheet */
        const float mx = nx2 - px, my2 = ny2 - py,
                    mz2 = nz2 - pz;
        const float ml = sqrtf(mx*mx + my2*my2 + mz2*mz2);
        /*  Stalled, jumped, or reversed: the extremum line
         *  ended (or the corrector left the curve).  */
        if (ml < 0.05f * sp || ml > 2.f * sp ||
            mx*dx + my2*dy + mz2*dz < 0)
            break;
        dx = mx / ml;  dy = my2 / ml;  dz = mz2 / ml;
        px = nx2;  py = ny2;  pz = nz2;
        pts->push_back({ px, py, pz });
        claimed->insert(cc_key(cell, px, py, pz));
        const float lx = px - sx0, ly = py - sy0, lz = pz - sz0;
        if (i > 8 &&
            lx*lx + ly*ly + lz*lz < 0.75f*0.75f * step*step)
            return true;
    }
    return false;
}

/*  Strategy-doc #4a: the contact-curve tracer.  The clause
 *  system {f_A = f_B, f = 0} has a DOUBLE ROOT at tangent
 *  contact (graveyard: Newton cannot march it) - but the contact
 *  locus is a curvature-extremum line, and curvature stays
 *  finite straight through the tangency.  So: seed from the
 *  hidden oracle's contact verdicts (persistent interval
 *  ambiguity = reach collapse), pull each seed onto the extremum
 *  line with valley_project, then march with a predictor-
 *  corrector - predictor steps along n x t_min (the along-curve
 *  tangent), corrector is valley_project itself.  Output:
 *  world-space polylines, SNAP TARGETS ONLY (never CCDT
 *  constraints - the crossing-suppression tombstones stay
 *  honored; if the marcher proves itself the full #4 constraint
 *  version is the next rung).  STIBIUM_DMESH_CONTACT_TRACE=0
 *  disables.  */
static void trace_contact_chains(const Deck* deck, float sp,
        const std::vector<std::array<float, 6>>& boxes,
        std::vector<std::vector<std::array<float, 3>>>* out)
{
    const char* env = getenv("STIBIUM_DMESH_CONTACT_TRACE");
    if ((env && atoi(env) == 0) || boxes.empty() || sp <= 0)
        return;
    /*  Own context: the collector's per-worker ctxs are freed
     *  after the descent merge (measured the hard way).  */
    TapeCtx* ctx = tape_ctx_new(deck);
    const float cell = 0.5f * sp;
    std::unordered_set<uint64_t> claimed;
    size_t nchains = 0, npts = 0;
    size_t r_vp = 0, r_sheet = 0, r_claim = 0, r_probe = 0,
           r_tang = 0, r_short = 0;
    for (const auto& b : boxes)
    {
        if (nchains >= 256)
            break;
        float sx = 0.5f * (b[0] + b[3]),
              sy = 0.5f * (b[1] + b[4]),
              sz = 0.5f * (b[2] + b[5]);
        /*  Valley first (blend seams), then the curvature-STEP
         *  fallback (fillet boundaries: plateaus joined by a
         *  jump - no peak exists; seed direction from the
         *  tangential gradient of MEAN curvature, basis-free).  */
        int mode = 0;
        float ax = 0, ay = 0, az = 0;
        if (!valley_project(deck, ctx, sp, sp, sx, sy, sz,
                            0.002f))
        {
            sx = 0.5f * (b[0] + b[3]);
            sy = 0.5f * (b[1] + b[4]);
            sz = 0.5f * (b[2] + b[5]);
            mode = 1;
            /*  Project the seed, then measure the tangential
             *  gradient of mean curvature.  Seeds routinely land
             *  on a featureless flat 1-2 cells from the junction
             *  (Newton takes the shortest path, not the
             *  interesting one - measured: every on-axis seed
             *  lands on the z=15 top flat, gradient EXACTLY
             *  zero); when that happens, WALK toward the seed's
             *  own contact box (tangent-projected) until the
             *  gradient wakes, then center on the step.  */
            float tx2 = sx, ty2 = sy, tz2 = sz;
            float n2[3], c0[3];
            float c0l = 0;
            if (!newton7(deck, ctx, sp, tx2, ty2, tz2))
                { ++r_vp; continue; }
            if (!mc_tgrad(deck, ctx, sp, tx2, ty2, tz2,
                          n2, c0, c0l))
                { ++r_vp; continue; }
            /*  Wake bar: a tenth of the step-prominence bar.  */
            const float cbar = 0.005f / sp;
            bool alive = c0l > cbar;
            for (int wk = 0; wk < 8 && !alive; ++wk)
            {
                float bx2 = 0.5f * (b[0] + b[3]) - tx2,
                      by2 = 0.5f * (b[1] + b[4]) - ty2,
                      bz2 = 0.5f * (b[2] + b[5]) - tz2;
                const float bn = bx2*n2[0] + by2*n2[1] +
                                 bz2*n2[2];
                bx2 -= bn*n2[0];
                by2 -= bn*n2[1];
                bz2 -= bn*n2[2];
                const float bl = sqrtf(bx2*bx2 + by2*by2 +
                                       bz2*bz2);
                if (!(bl > 1e-3f * sp))
                    break;
                tx2 += 0.5f * sp * bx2 / bl;
                ty2 += 0.5f * sp * by2 / bl;
                tz2 += 0.5f * sp * bz2 / bl;
                if (!newton7(deck, ctx, sp, tx2, ty2, tz2) ||
                    !mc_tgrad(deck, ctx, sp, tx2, ty2, tz2,
                              n2, c0, c0l))
                    break;
                alive = c0l > cbar;
            }
            if (!alive)
                { ++r_vp; continue; }
            c0[0] /= c0l;  c0[1] /= c0l;  c0[2] /= c0l;
            const float al0[3] = {
                    n2[1]*c0[2] - n2[2]*c0[1],
                    n2[2]*c0[0] - n2[0]*c0[2],
                    n2[0]*c0[1] - n2[1]*c0[0] };
            sx = tx2;  sy = ty2;  sz = tz2;
            if (!step_project(deck, ctx, sp, al0, sx, sy, sz))
                { ++r_vp; continue; }
        }
        if (!probe_real_surface(deck, ctx, sp, sx, sy, sz))
            { ++r_sheet; continue; }
        if (cc_claimed(claimed, cell, sx, sy, sz))
            { ++r_claim; continue; }
        KOut k0;
        if (!curvature_probe(deck, ctx, 0.1f * sp,
                             sx, sy, sz, &k0))
            { ++r_probe; continue; }
        float n[3];
        if (!probe_normal7(deck, ctx, sp, sx, sy, sz, n))
            { ++r_probe; continue; }
        float cdx, cdy, cdz;
        if (mode == 0)
            { cdx = k0.tx;  cdy = k0.ty;  cdz = k0.tz; }
        else
        {
            /*  Step mode: cross = steepest mean-curvature
             *  direction at the CENTERED point.  */
            float cd[3], cl2 = 0;
            if (!mc_tgrad(deck, ctx, sp, sx, sy, sz,
                          n, cd, cl2))
                { ++r_probe; continue; }
            if (!(cl2 > 0) || !std::isfinite(cl2))
                { ++r_tang; continue; }
            cdx = cd[0] / cl2;
            cdy = cd[1] / cl2;
            cdz = cd[2] / cl2;
        }
        ax = n[1]*cdz - n[2]*cdy;
        ay = n[2]*cdx - n[0]*cdz;
        az = n[0]*cdy - n[1]*cdx;
        const float al = sqrtf(ax*ax + ay*ay + az*az);
        if (!(al > 0.5f))
            { ++r_tang; continue; }
        ax /= al;  ay /= al;  az /= al;
        std::vector<std::array<float, 3>> fwd, bwd;
        fwd.push_back({ sx, sy, sz });
        claimed.insert(cc_key(cell, sx, sy, sz));
        const bool closed = march_xline(deck, ctx, sp,
                sx, sy, sz, ax, ay, az, mode, &claimed, &fwd);
        if (!closed)
            march_xline(deck, ctx, sp, sx, sy, sz,
                        -ax, -ay, -az, mode, &claimed, &bwd);
        if (fwd.size() + bwd.size() < 4)
            { ++r_short; continue; }
        std::vector<std::array<float, 3>> chain(bwd.rbegin(),
                                                bwd.rend());
        chain.insert(chain.end(), fwd.begin(), fwd.end());
        npts += chain.size();
        out->push_back(std::move(chain));
        ++nchains;
    }
    if (getenv("STIBIUM_DMESH_TIME") ||
        getenv("STIBIUM_DMESH_CHIP_DEBUG"))
    {
        fprintf(stderr, "CONTACT: %zu chains traced (%zu points) "
                "from %zu contact boxes (rejects: %zu vproj, "
                "%zu sheet, %zu claimed, %zu probe, %zu tangent, "
                "%zu short)\n",
                nchains, npts, boxes.size(), r_vp, r_sheet,
                r_claim, r_probe, r_tang, r_short);
        fprintf(stderr, "STEP: proj %llu, hess %llu, cdir %llu, "
                "prominence %llu, peak %llu, ok %llu\n",
                (unsigned long long)g_stp[0],
                (unsigned long long)g_stp[1],
                (unsigned long long)g_stp[2],
                (unsigned long long)g_stp[3],
                (unsigned long long)g_stp[4],
                (unsigned long long)g_stp[5]);
    }
    tape_ctx_free(ctx);
}

/*  #4a round 4 - TSEED STEP-TRACING (2026-07-19).  Fillet
 *  boundaries are curvature STEPS on ordinary sampled surface:
 *  they never become hidden candidates, so the contact boxes are
 *  the wrong seed door for them (measured on the on-axis lip -
 *  every buried-ring seed Newton-projects onto a featureless
 *  flat, gradient exactly zero).  The right door already exists:
 *  soup.tseeds, the shallow channel - on-surface QEF points from
 *  crease-leaf cells whose normal spread sits in exactly the
 *  sub-feature band (~14-25 deg) a fillet boundary occupies.
 *  Runs strictly AFTER delaunay_trace: crease chains claim their
 *  cells first, so crease law and step law never double-cover;
 *  contact chains from the sample stage claim theirs too.
 *  Output joins soup.contact_chains - SNAP TARGETS ONLY, same
 *  contract as the contact tracer (never CCDT constraints, never
 *  surface vertices).  STIBIUM_DMESH_STEP_TRACE=0 disables.  */
static void trace_step_seams(const Deck* deck, DSoup* soup,
                             volatile int* halt)
{
    const char* env = getenv("STIBIUM_DMESH_STEP_TRACE");
    if ((env && atoi(env) == 0) || soup->tseeds.empty() ||
        soup->spacing <= 0 || (halt && *halt))
        return;
    const float sp = soup->spacing;
    const float cell = 0.5f * sp;
    TapeCtx* ctx = tape_ctx_new(deck);
    std::unordered_set<uint64_t> claimed;
    for (const auto& ch : soup->tchains)
        for (const uint32_t i : ch)
            claimed.insert(cc_key(cell, soup->surface[i].x,
                                  soup->surface[i].y,
                                  soup->surface[i].z));
    for (const auto& ch : soup->contact_chains)
        for (const auto& p : ch)
            claimed.insert(cc_key(cell, p[0], p[1], p[2]));
    size_t nchains = 0, npts = 0, nclosed = 0;
    size_t r_claim = 0, r_proj = 0, r_flat = 0, r_probe = 0,
           r_step = 0, r_sheet = 0, r_tang = 0, r_short = 0;
    uint64_t stp0[6];
    memcpy(stp0, g_stp, sizeof(stp0));
    int dbg_rej = 0;
    /*  Wake bar: a tenth of the step-prominence bar (the same
     *  calibration the contact tracer measured).  */
    const float cbar = 0.005f / sp;
    /*  Ridge-rung prominence bar (1/sp units, default 0 = the
     *  plain strict rule): measured 2026-07-19, prominence
     *  cannot separate crossing spikes from step shoulders at
     *  any bar - the FOLD-EXCESS referee below does that; the
     *  knob stays for experiments.  */
    static const char* rp_env =
            getenv("STIBIUM_DMESH_RIDGE_PROM");
    const float rprom = (rp_env ? float(atof(rp_env)) : 0.f)
                        / std::max(sp, 1e-20f);
    for (const auto& s : soup->tseeds)
    {
        if (nchains >= 256 || (halt && *halt))
            break;
        float x = s.x, y = s.y, z = s.z;
        if (cc_claimed(claimed, cell, x, y, z))
            { ++r_claim; continue; }
        if (!newton7(deck, ctx, sp, x, y, z))
            { ++r_proj; continue; }
        float n[3], c0[3];
        float c0l = 0;
        if (!mc_tgrad(deck, ctx, sp, x, y, z, n, c0, c0l))
            { ++r_probe; continue; }
        /*  Quiet surface: most tseeds sit on blend interiors
         *  where mean curvature is locally constant - only the
         *  ramp near a junction wakes the gradient.  */
        if (!(c0l > cbar))
            { ++r_flat; continue; }
        c0[0] /= c0l;  c0[1] /= c0l;  c0[2] /= c0l;
        const float al0[3] = { n[1]*c0[2] - n[2]*c0[1],
                               n[2]*c0[0] - n[0]*c0[2],
                               n[0]*c0[1] - n[1]*c0[0] };
        /*  Step first (fillet boundaries - curvature RAMPS), then
         *  the valley corrector (shallow-crossing seams - the
         *  cone-torus lip meets at ~13 deg, below the crease
         *  tracer's reach, and its smeared kappa SPIKE flips the
         *  derivative sign mid-window: the strict peak gate
         *  refuses ramp-wise by design).  Same pairing as the
         *  contact seed path, priority reversed.  */
        int mode = 1;
        if (!step_project(deck, ctx, sp, al0, x, y, z))
        {
            /*  OPT-IN spike rungs (STIBIUM_DMESH_TSEED_RIDGE=1):
             *  a shallow CROSSING seam (the cone-torus lip meets
             *  at ~13 deg, below the crease tracer's reach) is a
             *  kappa spike, not a ramp - valley catches concave
             *  crossings, ridge convex ones.  NOT DEFAULT yet:
             *  measured 2026-07-19, the ridge rule cannot
             *  separate a true crossing spike from the stencil's
             *  overshoot SHOULDER beside a smeared G1 step (the
             *  off-axis lip mints a wobbling false rail 0.5 sp
             *  inside the fillet at every prominence bar that
             *  keeps the cone seam - 0.02/0.05/0.1 all fail one
             *  side).  Next primitive: normal-fold excess over
             *  the curvature-explained rotation.  */
            static const char* tr_env =
                    getenv("STIBIUM_DMESH_TSEED_RIDGE");
            bool got = false;
            if (tr_env && atoi(tr_env) != 0)
            {
                mode = 0;
                x = s.x;  y = s.y;  z = s.z;
                got = valley_project(deck, ctx, sp, sp,
                                     x, y, z, 0.002f);
                if (!got)
                {
                    x = s.x;  y = s.y;  z = s.z;
                    got = valley_project(deck, ctx, sp, sp,
                                         x, y, z, 0.002f, true,
                                         rprom);
                    /*  SCALE-DIVERGENCE referee: neither
                     *  prominence nor fold-excess separates a
                     *  crossing KINK from a smeared G1 STEP
                     *  (both measured failing, see MESH-WAR -
                     *  total normal turn is conserved under
                     *  smearing, so every step "folds").  What
                     *  separates them is how kappa RESPONDS TO
                     *  THE STENCIL: a kink is a delta - its
                     *  reading grows ~theta/2h as h shrinks; a
                     *  C1 step saturates at its taller plateau
                     *  (curvature between two smooth sheets
                     *  never exceeds either side).  Probe kmax
                     *  at h and h/2: divergent = kink, keep;
                     *  saturated = shoulder, drop.  */
                    if (got)
                    {
                        static const char* sd_env = getenv(
                                "STIBIUM_DMESH_SCALE_DIV");
                        const float sbar = sd_env
                                ? float(atof(sd_env)) : 1.4f;
                        /*  The projected candidate sits at the
                         *  SMEARED kappa peak, which the step
                         *  asymmetry biases off the kink line
                         *  (measured: ratio 1.000 at the
                         *  candidate itself) - sweep stations
                         *  across; some station lands within
                         *  h/2 of the kink and diverges.  */
                        /*  Station pitch must be finer than the
                         *  SMALL stencil's capture radius (h/2
                         *  = 0.025 sp) or every station misses
                         *  the kink and the small probe reads
                         *  LOWER than the big one (measured:
                         *  candidates sit ~0.035 sp off the
                         *  kink, 0.1-pitch stations all read
                         *  1.001).  */
                        KOut kc0;
                        float best = 0;
                        got = curvature_probe(deck, ctx,
                                      0.1f * sp, x, y, z, &kc0);
                        for (int st = -5; got && st <= 5; ++st)
                        {
                            const float qx =
                                    x + 0.02f*sp*st*kc0.sx,
                                    qy = y + 0.02f*sp*st*kc0.sy,
                                    qz = z + 0.02f*sp*st*kc0.sz;
                            KOut kh, kh2;
                            if (curvature_probe(deck, ctx,
                                        0.1f * sp, qx, qy, qz,
                                        &kh) &&
                                curvature_probe(deck, ctx,
                                        0.05f * sp, qx, qy, qz,
                                        &kh2) &&
                                kh.kmax > 0)
                                best = std::max(best,
                                        kh2.kmax / kh.kmax);
                        }
                        if (got)
                            got = best > sbar;
                        if (dbg_rej < 8 && getenv(
                                "STIBIUM_DMESH_CHIP_DEBUG"))
                        {
                            ++dbg_rej;
                            fprintf(stderr, "  TSEEDDBG ridge "
                                    "cand (%.3f,%.3f,%.3f) "
                                    "scale-ratio=%.3f %s\n",
                                    x, y, z, best,
                                    got ? "KEEP" : "drop");
                        }
                    }
                    if (got)
                        mode = 2;
                }
            }
            if (!got)
            {
                ++r_step;
                if (dbg_rej < 5 &&
                    getenv("STIBIUM_DMESH_CHIP_DEBUG"))
                {
                    ++dbg_rej;
                    fprintf(stderr, "  TSEEDDBG step-reject at "
                            "(%.3f,%.3f,%.3f) |mc-grad|=%.4f\n",
                            s.x, s.y, s.z, c0l);
                }
                continue;
            }
        }
        /*  Re-check at the CENTERED point: the first seed on a
         *  ring traces and claims the whole line; every later
         *  tseed centers onto the claimed cells and stops here.  */
        if (cc_claimed(claimed, cell, x, y, z))
            { ++r_claim; continue; }
        if (!probe_real_surface(deck, ctx, sp, x, y, z))
            { ++r_sheet; continue; }
        float cd[3];
        if (mode == 1)
        {
            float cl = 0;
            if (!mc_tgrad(deck, ctx, sp, x, y, z, n, cd, cl) ||
                !(cl > 0) || !std::isfinite(cl))
                { ++r_probe; continue; }
            cd[0] /= cl;  cd[1] /= cl;  cd[2] /= cl;
        }
        else
        {
            KOut k0;
            if (!curvature_probe(deck, ctx, 0.1f * sp,
                                 x, y, z, &k0) ||
                !probe_normal7(deck, ctx, sp, x, y, z, n))
                { ++r_probe; continue; }
            /*  Cross = the eigendirection the mode extremizes:
             *  t (kappa_min) for valleys, s (kappa_max) for
             *  ridges.  */
            if (mode == 0)
                { cd[0] = k0.tx;  cd[1] = k0.ty;  cd[2] = k0.tz; }
            else
                { cd[0] = k0.sx;  cd[1] = k0.sy;  cd[2] = k0.sz; }
        }
        float ax = n[1]*cd[2] - n[2]*cd[1],
              ay = n[2]*cd[0] - n[0]*cd[2],
              az = n[0]*cd[1] - n[1]*cd[0];
        const float al = sqrtf(ax*ax + ay*ay + az*az);
        if (!(al > 0.5f))
            { ++r_tang; continue; }
        ax /= al;  ay /= al;  az /= al;
        std::vector<std::array<float, 3>> fwd, bwd;
        fwd.push_back({ x, y, z });
        claimed.insert(cc_key(cell, x, y, z));
        const bool closed = march_xline(deck, ctx, sp,
                x, y, z, ax, ay, az, mode, &claimed, &fwd,
                rprom);
        if (!closed)
            march_xline(deck, ctx, sp, x, y, z,
                        -ax, -ay, -az, mode, &claimed, &bwd,
                        rprom);
        if (fwd.size() + bwd.size() < 4)
            { ++r_short; continue; }
        std::vector<std::array<float, 3>> chain(bwd.rbegin(),
                                                bwd.rend());
        chain.insert(chain.end(), fwd.begin(), fwd.end());
        npts += chain.size();
        nclosed += closed ? 1 : 0;
        soup->contact_chains.push_back(std::move(chain));
        ++nchains;
    }
    if (getenv("STIBIUM_DMESH_TIME") ||
        getenv("STIBIUM_DMESH_CHIP_DEBUG"))
    {
        fprintf(stderr, "TSEED: %zu chains (%zu closed, %zu "
                "points) from %zu shallow seeds (rejects: %zu "
                "claimed, %zu proj, %zu flat, %zu probe, %zu "
                "step, %zu sheet, %zu tangent, %zu short)\n",
                nchains, nclosed, npts, soup->tseeds.size(),
                r_claim, r_proj, r_flat, r_probe, r_step,
                r_sheet, r_tang, r_short);
        fprintf(stderr, "TSEED step gates: proj %llu, hess %llu, "
                "cdir %llu, prominence %llu, peak %llu, ok %llu\n",
                (unsigned long long)(g_stp[0] - stp0[0]),
                (unsigned long long)(g_stp[1] - stp0[1]),
                (unsigned long long)(g_stp[2] - stp0[2]),
                (unsigned long long)(g_stp[3] - stp0[3]),
                (unsigned long long)(g_stp[4] - stp0[4]),
                (unsigned long long)(g_stp[5] - stp0[5]));
    }
    tape_ctx_free(ctx);
}

/*  Keeter step 3: sharp-feature points.  For every candidate cell,
 *  probe normals at its crossing points (batched central
 *  differences); if the normals disagree enough, the cell straddles
 *  a crease - solve the QEF (SVD with singular-value clamping,
 *  regularized toward the centroid, the classic robust DC solve)
 *  and add the minimizer as a surface point, rejecting solutions
 *  that escape the cell.  */
void feature_points(Collector& c)
{
    if (c.cells.empty() || c.soup.surface.empty())
        return;

    /*  Normals at every surface point (many cells share points) -
     *  six probes per point, batched.  */
    const size_t np = c.soup.surface.size();
    const float h = c.spacing * 0.01f;
    std::vector<float> xs(np * 6), ys(np * 6), zs(np * 6), vals;
    for (size_t i = 0; i < np; ++i)
    {
        const DSurfPoint& p = c.soup.surface[i];
        for (int q = 0; q < 6; ++q)
        {
            xs[i * 6 + q] = p.x;
            ys[i * 6 + q] = p.y;
            zs[i * 6 + q] = p.z;
        }
        xs[i * 6 + 0] += h;   xs[i * 6 + 1] -= h;
        ys[i * 6 + 2] += h;   ys[i * 6 + 3] -= h;
        zs[i * 6 + 4] += h;   zs[i * 6 + 5] -= h;
    }
    eval_points(deck_base(c.deck), c.ctx, xs, ys, zs, vals);

    std::vector<Eigen::Vector3f> normals(np);
    for (size_t i = 0; i < np; ++i)
    {
        Eigen::Vector3f n(vals[i * 6 + 0] - vals[i * 6 + 1],
                          vals[i * 6 + 2] - vals[i * 6 + 3],
                          vals[i * 6 + 4] - vals[i * 6 + 5]);
        const float len = n.norm();
        normals[i] = len > 0 ? Eigen::Vector3f(n / len)
                             : Eigen::Vector3f::Zero();
    }

    constexpr float SPREAD_DOT = 0.9f;   // ~25 degrees

    /*  Stage-D flag with the standing gates (bar, level-from-
     *  magnitude, tangle demotion - zeiss: level 2 in thin-wall
     *  tangles is nm 241 vs 48).  */
    const auto flag_leaf = [&c](uint64_t key, float nrv) {
        if (!autodense() || nrv < AUTOD_BAR)
            return;
        /*  Tangle = a thin WALL within pinch reach of the level-2
         *  pitch (STIBIUM_DMESH_SEP, in cells, default 0.7 =
         *  2.8x the quarter-cell pitch).  Anti-parallelism alone
         *  was the previous gate and over-blocked: concave grooves
         *  (the thing density CURES) read anti-parallel too; only
         *  the signed wall gap separates them (measured - the
         *  populations are inseparable without the sign).  */
        const auto lit = c.crease_leaves.find(key);
        const bool tangle = lit != c.crease_leaves.end() &&
                lit->second.sep < autod_sep_bar() * c.spacing;
        if (tangle)
            ++c.tangle_suppressed;
        /*  One level per factor of 4 over the bar (sagitta
         *  quarters per level): 0.12 cells wants 2, 0.48 wants 3.
         *  Level 3 is LOCAL eighth-cells in grooves - safe only
         *  because the retreat loop measures and rolls back any
         *  hole it opens (global 8x remains a proven disaster).  */
        int lvl = 1;
        if (!tangle)
            for (float t = 4 * AUTOD_BAR;
                 nrv >= t && lvl < autod_max_level(); t *= 4)
                ++lvl;
        int& lv = c.want_dense[key];
        lv = std::max(lv, lvl);
        c.resid_hot.insert(key);
    };

    /*  QEF candidates awaiting the phantom oracle (below).  */
    struct FeatCand
    {
        const FeatCell* fc;
        float x[3];
        bool shallow;   // tracer-only seed: no tail, no suppression
    };
    std::vector<FeatCand> cands;
    uint64_t added = 0;
    /*  A QEF vertex REPLACES the crossings it consumed (dual-
     *  contouring semantics): leaving both in the soup makes the
     *  triangulation alternate between on-crease and off-crease
     *  vertices - a beautifully regular sawtooth (caught by eyeball
     *  review, as usual).  */
    std::vector<uint8_t> suppress(np, 0);
    for (const FeatCell& fc : c.cells)
    {
        /*  Spread test: does any normal pair disagree?  */
        float min_dot = 1;
        for (uint8_t a = 0; a < fc.n; ++a)
            for (uint8_t b = uint8_t(a + 1); b < fc.n; ++b)
            {
                const float d = normals[fc.pts[a]]
                                        .dot(normals[fc.pts[b]]);
                if (d < min_dot)
                    min_dot = d;
            }
        const bool feat = min_dot <= SPREAD_DOT;
        /*  Shallow channel: too flat for feature status, but a
         *  live crease pair runs through this leaf - the crease
         *  may simply be shallower than the 25-degree gate (the
         *  under-polylined crowded bands in Nate's chain overlay,
         *  2026-07-17).  Offer the tracer a seed; it verifies
         *  everything.  The crease-leaf test keeps smooth curved
         *  shells (whose spread is pure curvature) out of the
         *  pool.  */
        const bool shallow = !feat &&
                min_dot <= shallow_seed_dot() &&
                c.crease_leaves.count(fc.leaf_key) != 0;
        if (!feat && !shallow && !c.census && !autodense())
            continue;

        Eigen::MatrixXf A(fc.n, 3);
        Eigen::VectorXf b(fc.n);
        Eigen::Vector3f centroid = Eigen::Vector3f::Zero();
        for (uint8_t q = 0; q < fc.n; ++q)
        {
            const DSurfPoint& p = c.soup.surface[fc.pts[q]];
            centroid += Eigen::Vector3f(p.x, p.y, p.z);
        }
        centroid /= float(fc.n);
        for (uint8_t q = 0; q < fc.n; ++q)
        {
            const DSurfPoint& p = c.soup.surface[fc.pts[q]];
            A.row(q) = normals[fc.pts[q]];
            b(q) = normals[fc.pts[q]].dot(
                    Eigen::Vector3f(p.x, p.y, p.z) - centroid);
        }
        Eigen::JacobiSVD<Eigen::MatrixXf> svd(
                A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        Eigen::Vector3f sv = svd.singularValues();
        const float cutoff = sv(0) * 0.1f;
        Eigen::Vector3f inv = Eigen::Vector3f::Zero();
        for (int q = 0; q < 3; ++q)
            if (sv(q) > cutoff)
                inv(q) = 1.0f / sv(q);
        const Eigen::Vector3f x = centroid +
                svd.matrixV() *
                        inv.asDiagonal() *
                        (svd.matrixU().transpose() * b);

        /*  QEF residual in units of the cell's own edge (pass-2
         *  cells are smaller; the bar is scale-free).  This is the
         *  stage-D trigger: sharp creases solve to ~0 at any
         *  resolution, while blend aliasing, crease crowding, and
         *  sub-lattice churn all leave a residual the minimizer
         *  cannot explain (census 2026-07-16: solved models are
         *  silent above 0.01; every problem model lights up).  */
        const float nr = (A * (x - centroid) - b).norm() /
                (sqrtf(float(fc.n)) * (fc.hi[0] - fc.lo[0]));
        flag_leaf(fc.leaf_key, nr);
        if (c.census)
        {
            int bucket = 17;
            for (int q = 0; q < 17; ++q)
                if (nr < Census::NR_EDGES[q])
                {
                    bucket = q;
                    break;
                }
            ++c.census->nr_hist[bucket];
            if (c.census->dump)
                fprintf(c.census->dump,
                        "CELL %.6g %.6g %.6g n=%u mindot=%.4f "
                        "nr=%.5f feat=%d\n",
                        0.5f * (fc.lo[0] + fc.hi[0]),
                        0.5f * (fc.lo[1] + fc.hi[1]),
                        0.5f * (fc.lo[2] + fc.hi[2]),
                        unsigned(fc.n), min_dot, nr, int(feat));
        }
        if (!feat && !shallow)
            continue;

        /*  Reject escapes (with a small margin): a QEF that leaves
         *  its cell is extrapolating noise.  */
        const float margin = c.spacing * 0.25f;
        bool ok = true;
        for (int q = 0; q < 3; ++q)
            if (x(q) < fc.lo[q] - margin || x(q) > fc.hi[q] + margin)
                ok = false;
        if (!ok)
            continue;

        cands.push_back({ &fc, { x(0), x(1), x(2) }, shallow });
    }

    /*  Phantom oracle (2026-07-16, the zeiss "knurled collar"
     *  autopsy): a sub-lattice fillet turns entirely inside one
     *  cell, so the cell's crossings read the two FLANK normals
     *  and the QEF reconstructs the sharp corner the model rounds
     *  off - a vertex up to half a cell OFF the surface, invisible
     *  to the residual (a corner explains the crossings; the
     *  surface just disagrees).  One question the pipeline never
     *  asked: is the minimizer ON the surface?  Real corners read
     *  |f|/|grad| ~ 1e-3 cells; phantoms read the fillet depth.
     *  Phantoms are rejected (their true crossings stay), and the
     *  miss distance feeds stage-D - it is exactly the sagitta the
     *  cell hid from the residual.  */
    if (!cands.empty())
    {
        const size_t nc = cands.size();
        const float h = c.spacing * 0.01f;
        std::vector<float> qx(nc * 7), qy(nc * 7), qz(nc * 7), qv;
        for (size_t i = 0; i < nc; ++i)
        {
            for (int q = 0; q < 7; ++q)
            {
                qx[i*7 + q] = cands[i].x[0];
                qy[i*7 + q] = cands[i].x[1];
                qz[i*7 + q] = cands[i].x[2];
            }
            qx[i*7 + 1] += h;   qx[i*7 + 2] -= h;
            qy[i*7 + 3] += h;   qy[i*7 + 4] -= h;
            qz[i*7 + 5] += h;   qz[i*7 + 6] -= h;
        }
        eval_points(deck_base(c.deck), c.ctx, qx, qy, qz, qv);

        const auto fdist_of = [&](size_t i) {
            const float f = qv[i*7];
            const float gx = (qv[i*7 + 1] - qv[i*7 + 2]);
            const float gy = (qv[i*7 + 3] - qv[i*7 + 4]);
            const float gz = (qv[i*7 + 5] - qv[i*7 + 6]);
            const float gl = sqrtf(gx*gx + gy*gy + gz*gz) / (2 * h);
            const float f2 = fabsf(f) / (gl > 1e-12f ? gl : 1e30f);
            return (std::isfinite(f) && std::isfinite(gl)) ? f2
                                                           : 0.f;
        };

        /*  Off-surface minimizers are PROJECTED back onto f = 0,
         *  not discarded: a tightly curved TRUE crease leaves the
         *  plane-fit corner a small sagitta off the surface, and
         *  rejecting those punches gaps in the fallback chains
         *  (csg loop -> 3 open chains, [.dchain]).  A fillet
         *  phantom projects onto the mid-band - an on-surface
         *  vertex where the knurl used to bulge.  Every
         *  off-surface minimizer flags stage-D with its
         *  PRE-projection miss - that distance is the
         *  corner-rounding depth, exactly the groove sagitta the
         *  cell hides from the residual (the first projection
         *  round only flagged UNprojectable points, silently
         *  starving projected groove cells of their cores).  */
        std::vector<size_t> off;
        for (size_t i = 0; i < nc; ++i)
        {
            const FeatCell& fc = *cands[i].fc;
            const float cell = fc.hi[0] - fc.lo[0];
            const float fd = fdist_of(i);
            if (fd > 0.05f * cell)
            {
                flag_leaf(fc.leaf_key, fd / cell);
                off.push_back(i);
            }
        }
        if (!off.empty())
        {
            std::vector<float> px2(off.size()), py2(off.size()),
                    pz2(off.size()),
                    hs2(off.size(), h),
                    cl2(off.size(), 0.75f * c.spacing);
            for (size_t k = 0; k < off.size(); ++k)
            {
                px2[k] = cands[off[k]].x[0];
                py2[k] = cands[off[k]].x[1];
                pz2[k] = cands[off[k]].x[2];
            }
            project_points_impl(c.deck, c.ctx, px2, py2, pz2,
                                hs2, cl2);
            for (size_t k = 0; k < off.size(); ++k)
            {
                cands[off[k]].x[0] = px2[k];
                cands[off[k]].x[1] = py2[k];
                cands[off[k]].x[2] = pz2[k];
            }
            /*  Re-measure the projected points (batched).  */
            for (size_t k = 0; k < off.size(); ++k)
            {
                const size_t i = off[k];
                for (int q = 0; q < 7; ++q)
                {
                    qx[i*7 + q] = cands[i].x[0];
                    qy[i*7 + q] = cands[i].x[1];
                    qz[i*7 + q] = cands[i].x[2];
                }
                qx[i*7 + 1] += h;   qx[i*7 + 2] -= h;
                qy[i*7 + 3] += h;   qy[i*7 + 4] -= h;
                qz[i*7 + 5] += h;   qz[i*7 + 6] -= h;
            }
            eval_points(deck_base(c.deck), c.ctx, qx, qy, qz, qv);
        }

        for (size_t i = 0; i < nc; ++i)
        {
            const FeatCell& fc = *cands[i].fc;
            const float cell = fc.hi[0] - fc.lo[0];
            const float fdist = fdist_of(i);
            if (fdist > 0.05f * cell)
            {
                ++c.phantom_rejected;
                flag_leaf(fc.leaf_key, fdist / cell);
                if (c.census && c.census->dump)
                    fprintf(c.census->dump,
                            "PHANTOM %.6g %.6g %.6g fdist=%.5f\n",
                            cands[i].x[0], cands[i].x[1],
                            cands[i].x[2], fdist / cell);
                continue;
            }
            if (cands[i].shallow)
            {
                /*  Tracer-only: the seed pool sees it; the
                 *  feature tail, crossing suppression, and the
                 *  fallback chain graph never do.  */
                c.soup.tseeds.push_back({ cands[i].x[0],
                                          cands[i].x[1],
                                          cands[i].x[2] });
                continue;
            }
            c.soup.surface.push_back({ cands[i].x[0],
                                       cands[i].x[1],
                                       cands[i].x[2] });
            ++added;
            for (uint8_t q = 0; q < fc.n; ++q)
                suppress[fc.pts[q]] = 1;
        }
    }
    c.soup.feature_points = added;

    /*  Chain mediation (the 'outliers' Nate spotted): a feature
     *  point whose two nearest siblings run roughly THROUGH it
     *  (chain-interior; a cube corner's neighbors turn 90 degrees,
     *  so corners are immune) but which deviates from their segment
     *  is a bad QEF solve - replace it with the siblings' midpoint
     *  projected onto the surface.  Straight chains measure zero
     *  deviation and are untouched.  */
    if (added > 2)
    {
        const size_t f0 = np;   // features start here (pre-compact)
        const size_t nf = added;
        std::vector<uint32_t> med_idx;
        std::vector<float> mx, my, mz, mh, mc;
        for (size_t i = 0; i < nf; ++i)
        {
            const DSurfPoint& F = c.soup.surface[f0 + i];
            /*  two nearest siblings within 2.5 cells  */
            float d1 = 1e30f, d2 = 1e30f;
            size_t i1 = SIZE_MAX, i2 = SIZE_MAX;
            const float lim = c.spacing * 2.5f;
            for (size_t j = 0; j < nf; ++j)
            {
                if (j == i)
                    continue;
                const DSurfPoint& G = c.soup.surface[f0 + j];
                const float dx = G.x - F.x, dy = G.y - F.y,
                            dz = G.z - F.z;
                const float d = sqrtf(dx*dx + dy*dy + dz*dz);
                if (d > lim)
                    continue;
                if (d < d1)
                    { d2 = d1; i2 = i1; d1 = d; i1 = j; }
                else if (d < d2)
                    { d2 = d; i2 = j; }
            }
            if (i2 == SIZE_MAX)
                continue;
            const DSurfPoint& A = c.soup.surface[f0 + i1];
            const DSurfPoint& B = c.soup.surface[f0 + i2];
            /*  chain-interior test: A-F and F-B roughly opposed  */
            const float ax = F.x - A.x, ay = F.y - A.y, az = F.z - A.z;
            const float bx = B.x - F.x, by = B.y - F.y, bz = B.z - F.z;
            const float la = sqrtf(ax*ax + ay*ay + az*az);
            const float lb = sqrtf(bx*bx + by*by + bz*bz);
            if (!(la > 0) || !(lb > 0))
                continue;
            const float cosang =
                    (ax*bx + ay*by + az*bz) / (la * lb);
            if (cosang < 0.7f)
                continue;   // corner or junction: leave it alone
            /*  deviation of F from segment AB  */
            const float abx = B.x - A.x, aby = B.y - A.y,
                        abz = B.z - A.z;
            const float ab2 = abx*abx + aby*aby + abz*abz;
            if (!(ab2 > 0))
                continue;
            float t = ((F.x - A.x) * abx + (F.y - A.y) * aby +
                       (F.z - A.z) * abz) / ab2;
            t = t < 0 ? 0 : (t > 1 ? 1 : t);
            const float qx = A.x + t * abx - F.x;
            const float qy = A.y + t * aby - F.y;
            const float qz = A.z + t * abz - F.z;
            const float dev = sqrtf(qx*qx + qy*qy + qz*qz);
            if (dev <= c.spacing * 0.35f)
                continue;
            med_idx.push_back(uint32_t(f0 + i));
            mx.push_back((A.x + B.x) * 0.5f);
            my.push_back((A.y + B.y) * 0.5f);
            mz.push_back((A.z + B.z) * 0.5f);
            mh.push_back(c.spacing * 0.01f);
            mc.push_back(c.spacing);
        }
        if (!med_idx.empty())
        {
            project_points_impl(c.deck, c.ctx, mx, my, mz, mh, mc);
            for (size_t q = 0; q < med_idx.size(); ++q)
                c.soup.surface[med_idx[q]] =
                        { mx[q], my[q], mz[q] };
            c.soup.mediated = med_idx.size();
        }
    }

    if (added)
    {
        /*  Compact: drop suppressed crossings, keep the appended
         *  feature points (they live past index np).  */
        std::vector<DSurfPoint> kept;
        kept.reserve(c.soup.surface.size());
        for (size_t i = 0; i < c.soup.surface.size(); ++i)
            if (i >= np || !suppress[i])
                kept.push_back(c.soup.surface[i]);
        c.soup.suppressed = c.soup.surface.size() - kept.size();
        c.soup.surface.swap(kept);
    }
}

}  // namespace

DChains delaunay_chains(const DSoup& soup)
{
    DChains out;
    const size_t nf_raw = size_t(soup.feature_points);
    if (nf_raw < 2)
    {
        out.stray = nf_raw;
        return out;
    }
    const size_t f0 = soup.surface.size() - nf_raw;

    /*  1. Merge near-duplicates: adjacent cells\' QEFs converge to
     *  nearly the same crease point on curved creases, clumping the
     *  graph into junction storms.  Representatives only, within
     *  half a cell.  */
    const float mr = soup.spacing * 0.5f;
    const float mr2 = mr * mr;
    std::vector<uint32_t> rep;   // indices into soup.surface
    for (size_t i = 0; i < nf_raw; ++i)
    {
        const DSurfPoint& P = soup.surface[f0 + i];
        bool dup = false;
        for (const uint32_t r : rep)
        {
            const DSurfPoint& Q = soup.surface[r];
            const float dx = Q.x - P.x, dy = Q.y - P.y,
                        dz = Q.z - P.z;
            if (dx*dx + dy*dy + dz*dz < mr2)
            {
                dup = true;
                break;
            }
        }
        if (!dup)
            rep.push_back(uint32_t(f0 + i));
    }
    const size_t nf = rep.size();
    out.reps = nf;
    if (nf < 2)
    {
        out.stray = nf;
        return out;
    }

    /*  1b. Junction classification (the junction-split walk): a
     *  crease-interior point's neighborhood is a LINE (one dominant
     *  covariance eigenvalue); where crease curves meet or cross
     *  (cube corners, csg's seam-meets-circle points) it is not.
     *  Junction reps are cut out before linking - the walk cannot
     *  hop between curves through a crossing it never sees - and
     *  reattached afterward as SHARED chain endpoints.  (Gumhold
     *  2001 / Pauly 2003 covariance classification; see
     *  MESH-NEXT.)  */
    const float lim = soup.spacing * 2.0f;
    const float lim2 = lim * lim;
    std::vector<uint8_t> junction(nf, 0);
    for (size_t i = 0; i < nf; ++i)
    {
        const DSurfPoint& P = soup.surface[rep[i]];
        float c[6] = { 0, 0, 0, 0, 0, 0 };   // xx xy xz yy yz zz
        int nn = 0;
        for (size_t j = 0; j < nf; ++j)
        {
            if (j == i)
                continue;
            const DSurfPoint& Q = soup.surface[rep[j]];
            const float dx = Q.x - P.x, dy = Q.y - P.y,
                        dz = Q.z - P.z;
            if (dx*dx + dy*dy + dz*dz >= lim2)
                continue;
            c[0] += dx*dx;  c[1] += dx*dy;  c[2] += dx*dz;
            c[3] += dy*dy;  c[4] += dy*dz;  c[5] += dz*dz;
            ++nn;
        }
        if (nn < 3)
            continue;
        Eigen::Matrix3f m;
        m << c[0], c[1], c[2],
             c[1], c[3], c[4],
             c[2], c[4], c[5];
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> es(m);
        const auto& ev = es.eigenvalues();   // ascending
        if (ev[2] > 0 && ev[1] / ev[2] > 0.2f)
        {
            junction[i] = 1;
            ++out.junctions;
        }
    }

    /*  2. Two-sided linking: each point takes its nearest neighbor,
     *  then its nearest neighbor in the OPPOSING direction (dot <
     *  -0.3).  Chains get at most two links per node; corners end
     *  chains naturally (the second edge turns 90 degrees and fails
     *  the opposition test).  Junction reps do not link.  */
    std::vector<std::array<uint32_t, 2>> link(
            nf, { UINT32_MAX, UINT32_MAX });
    for (size_t i = 0; i < nf; ++i)
    {
        if (junction[i])
            continue;
        const DSurfPoint& P = soup.surface[rep[i]];
        float d1 = 1e30f;
        size_t n1 = SIZE_MAX;
        for (size_t j = 0; j < nf; ++j)
        {
            if (j == i || junction[j])
                continue;
            const DSurfPoint& Q = soup.surface[rep[j]];
            const float dx = Q.x - P.x, dy = Q.y - P.y,
                        dz = Q.z - P.z;
            const float d = dx*dx + dy*dy + dz*dz;
            if (d < d1 && d < lim2)
            {
                d1 = d;
                n1 = j;
            }
        }
        if (n1 == SIZE_MAX)
            continue;
        link[i][0] = uint32_t(n1);
        const DSurfPoint& N1 = soup.surface[rep[n1]];
        const float ax = N1.x - P.x, ay = N1.y - P.y,
                    az = N1.z - P.z;
        const float la = sqrtf(ax*ax + ay*ay + az*az);
        float d2 = 1e30f;
        size_t n2 = SIZE_MAX;
        for (size_t j = 0; j < nf; ++j)
        {
            if (j == i || j == n1 || junction[j])
                continue;
            const DSurfPoint& Q = soup.surface[rep[j]];
            const float bx = Q.x - P.x, by = Q.y - P.y,
                        bz = Q.z - P.z;
            const float d = bx*bx + by*by + bz*bz;
            if (d >= d2 || d >= lim2)
                continue;
            const float lb = sqrtf(d);
            if (!(la > 0) || !(lb > 0))
                continue;
            const float cosang =
                    (ax*bx + ay*by + az*bz) / (la * lb);
            if (cosang < -0.3f)
            {
                d2 = d;
                n2 = j;
            }
        }
        if (n2 != SIZE_MAX)
            link[i][1] = uint32_t(n2);
    }

    /*  3. Symmetrize: an edge exists only if EITHER endpoint links
     *  the other (one-sided is enough - curved spacing is uneven),
     *  capped at degree 2 by construction of the walk.  */
    const auto linked = [&](size_t i, uint32_t j) {
        return link[i][0] == j || link[i][1] == j;
    };
    const auto next_of = [&](size_t cur, uint32_t prev,
                             std::vector<uint8_t>& used) -> uint32_t {
        for (const uint32_t cand : link[cur])
        {
            if (cand == UINT32_MAX || cand == prev || used[cand])
                continue;
            /*  One-sided links are accepted (correctness review
             *  #4: the symmetric-link test here was dead code via
             *  "|| true" through every refereed rev - kept
             *  one-sided EXPLICITLY; curved spacing is uneven and
             *  the walk caps degree at 2 regardless).  */
            return cand;
        }
        return UINT32_MAX;
    };

    std::vector<uint8_t> used(nf, 0);
    const auto walk = [&](size_t start) {
        std::vector<uint32_t> chain;
        /*  Walk both directions from start  */
        std::vector<uint32_t> fwd, bwd;
        uint32_t prev = UINT32_MAX, cur = uint32_t(start);
        used[start] = 1;
        fwd.push_back(uint32_t(start));
        while (true)
        {
            const uint32_t nx = next_of(cur, prev, used);
            if (nx == UINT32_MAX)
                break;
            prev = cur;
            cur = nx;
            used[cur] = 1;
            fwd.push_back(cur);
        }
        prev = fwd.size() > 1 ? fwd[1] : UINT32_MAX;
        cur = uint32_t(start);
        while (true)
        {
            const uint32_t nx = next_of(cur, prev, used);
            if (nx == UINT32_MAX)
                break;
            prev = cur;
            cur = nx;
            used[cur] = 1;
            bwd.push_back(cur);
        }
        std::vector<uint32_t> chain_idx;
        for (auto it = bwd.rbegin(); it != bwd.rend(); ++it)
            chain_idx.push_back(*it);
        for (const uint32_t v : fwd)
            chain_idx.push_back(v);
        return chain_idx;
    };

    /*  Reattach: chains end where junctions were cut out; each
     *  open end takes the nearest junction rep within reach as a
     *  SHARED endpoint - three cube edges meet at one corner
     *  vertex, crossing curves share the crossing.  A curve
     *  passing through a single junction comes back as an open
     *  chain with the junction at both ends (the segments still
     *  close the loop through it).  */
    const auto nearest_junction = [&](uint32_t v,
                                      uint32_t exclude) -> uint32_t {
        const DSurfPoint& P = soup.surface[rep[v]];
        float best = lim2;
        uint32_t bj = UINT32_MAX;
        for (size_t j = 0; j < nf; ++j)
        {
            if (!junction[j] || uint32_t(j) == exclude)
                continue;
            const DSurfPoint& Q = soup.surface[rep[j]];
            const float dx = Q.x - P.x, dy = Q.y - P.y,
                        dz = Q.z - P.z;
            const float d = dx*dx + dy*dy + dz*dz;
            if (d < best)
            {
                best = d;
                bj = uint32_t(j);
            }
        }
        return bj;
    };

    for (size_t i = 0; i < nf; ++i)
    {
        if (used[i] || junction[i])
            continue;
        auto local = walk(i);
        if (local.size() < 2)
        {
            /*  A lone rep between two junctions is still an arm:
             *  bridge it.  */
            if (local.size() == 1)
            {
                const uint32_t j1 =
                        nearest_junction(local[0], UINT32_MAX);
                const uint32_t j2 = j1 == UINT32_MAX ? UINT32_MAX
                        : nearest_junction(local[0], j1);
                if (j1 != UINT32_MAX && j2 != UINT32_MAX)
                {
                    out.chains.push_back(
                            { rep[j1], rep[local[0]], rep[j2] });
                    out.closed.push_back(0);
                    continue;
                }
            }
            ++out.stray;
            continue;
        }
        /*  Closed if the two ends link each other  */
        const bool closed =
                local.size() >= 3 &&
                (linked(local.front(), local.back()) ||
                 linked(local.back(), local.front()));
        std::vector<uint32_t> chain;
        chain.reserve(local.size() + 2);
        for (const uint32_t v : local)
            chain.push_back(rep[v]);
        if (!closed)
        {
            const uint32_t jf =
                    nearest_junction(local.front(), UINT32_MAX);
            const uint32_t jb =
                    nearest_junction(local.back(), UINT32_MAX);
            if (jf != UINT32_MAX)
                chain.insert(chain.begin(), rep[jf]);
            if (jb != UINT32_MAX)
                chain.push_back(rep[jb]);
        }
        out.chains.push_back(std::move(chain));
        out.closed.push_back(closed ? 1 : 0);
    }
    return out;
}

static void census_print(const Census& census, const DSoup& soup,
                         const char* tag)
{
    fprintf(stderr, "%s %llu crease leaves of %llu\n", tag,
            (unsigned long long)census.leaves,
            (unsigned long long)soup.leaf_blocks);
    fprintf(stderr, "%s sheet-active  0/1/2/3+: "
            "%llu %llu %llu %llu\n", tag,
            (unsigned long long)census.by_sheet[0],
            (unsigned long long)census.by_sheet[1],
            (unsigned long long)census.by_sheet[2],
            (unsigned long long)census.by_sheet[3]);
    fprintf(stderr, "%s crease-active 0/1/2/3+: "
            "%llu %llu %llu %llu\n", tag,
            (unsigned long long)census.by_crease[0],
            (unsigned long long)census.by_crease[1],
            (unsigned long long)census.by_crease[2],
            (unsigned long long)census.by_crease[3]);
    fprintf(stderr, "%s min-dot hist [-1..1]:", tag);
    for (int i = 0; i < 12; ++i)
        fprintf(stderr, " %llu",
                (unsigned long long)census.dot_hist[i]);
    fprintf(stderr, "\n");
    fprintf(stderr, "%s cell nr hist (edges .005 .01 .02 "
            ".03 .05 .075 .1 .15 .2 .3 .5 .75 1 1.5 2 3 5):", tag);
    for (int i = 0; i < 18; ++i)
        fprintf(stderr, " %llu",
                (unsigned long long)census.nr_hist[i]);
    fprintf(stderr, "\n");
}

/*  One stage-A pass: descend + bisect + QEF.  dense_map is the
 *  drill-down input (null on pass 1); the pass's own drill-down
 *  requests land in c.want_dense.  */
/*  Depth-limited geometric split for the parallel descend: up to
 *  8^depth task regions whose X/Y/Z pointers alias the ROOT
 *  coordinate arrays (octsect preserves storage), so tasks are
 *  copies-by-value with shared read-only backing.  No interval
 *  culling here - each worker culls from the root tape itself
 *  (two extra interval evals per task; the leaf tapes are
 *  box-idempotent so leaf state matches the serial chain).  */
static void collect_tasks(const Region& r, int depth,
                          std::vector<Region>& out)
{
    if (depth <= 0 || r.voxels <= LEAF_VOXELS)
    {
        out.push_back(r);
        return;
    }
    Region oct[8];
    const uint8_t split = octsect(r, oct);
    for (int i = 0; i < 8; ++i)
        if (split & (1 << i))
            collect_tasks(oct[i], depth - 1, out);
}

/*  Merge a worker's Collector into the shared one, in FIXED task
 *  order (determinism: results must not depend on thread
 *  scheduling).  Samples and edges dedup through the shared
 *  hash indices exactly as add_sample/add_edge would; feature
 *  cells remap their crossing indices through the edge merge.  */
static void merge_collector(Collector& c, Collector& w)
{
    for (const DSample& s : w.soup.samples)
        c.add_sample(s.x, s.y, s.z, s.inside, s.on_surface);
    std::vector<uint64_t> h_of(w.edges.size());
    for (const auto& [h, i] : w.edge_index)
        h_of[i] = h;
    std::vector<uint32_t> remap(w.edges.size());
    for (size_t i = 0; i < w.edges.size(); ++i)
    {
        const auto it = c.edge_index.find(h_of[i]);
        if (it != c.edge_index.end())
            remap[i] = it->second;
        else
        {
            remap[i] = uint32_t(c.edges.size());
            c.edge_index.emplace(h_of[i], remap[i]);
            c.edges.push_back(w.edges[i]);
        }
    }
    for (FeatCell fc : w.cells)
    {
        for (uint8_t q = 0; q < fc.n; ++q)
            fc.pts[q] = remap[fc.pts[q]];
        c.cells.push_back(fc);
    }
    for (auto& [k, v] : w.crease_leaves)
        c.crease_leaves.emplace(k, v);
    for (const auto& [k, v] : w.want_dense)
    {
        int& lv = c.want_dense[k];
        lv = std::max(lv, v);
    }
    for (const auto& b : w.soup.dense_boxes)
        c.soup.dense_boxes.push_back(b);
    c.soup.leaf_blocks += w.soup.leaf_blocks;
    c.soup.culled_empty += w.soup.culled_empty;
    c.soup.culled_full += w.soup.culled_full;
    c.soup.hidden_candidates += w.soup.hidden_candidates;
    c.soup.dense_blocks += w.soup.dense_blocks;
    c.soup.welded += w.soup.welded;
    c.soup.thinned += w.soup.thinned;
    c.tangle_suppressed += w.tangle_suppressed;
    c.hidden_feature += w.hidden_feature;
    c.hidden_graze += w.hidden_graze;
    c.hidden_contact += w.hidden_contact;
    c.hidden_cell_feat += w.hidden_cell_feat;
    c.hidden_cell_wit += w.hidden_cell_wit;
    for (const auto& b : w.hidden_feat_boxes)
        c.hidden_feat_boxes.push_back(b);
    for (const auto& b : w.hidden_graze_boxes)
        c.hidden_graze_boxes.push_back(b);
    for (const auto& b : w.hidden_contact_boxes)
        c.hidden_contact_boxes.push_back(b);
    c.curve_flagged += w.curve_flagged;
    c.curve_seen += w.curve_seen;
    c.curve_cross += w.curve_cross;
    c.curve_theta_max = std::max(c.curve_theta_max,
                                 w.curve_theta_max);
    for (const auto& [k, b] : w.smooth_leaves)
        c.smooth_leaves.emplace(k, b);
    if (c.spacing == 0)
        c.spacing = w.spacing;
}

static void sample_pass(Collector& c, const Deck* deck, Region r,
                        volatile int* halt,
                        const std::unordered_map<uint64_t, int>* map,
                        Census* census,
                        const std::vector<std::array<float, 3>>*
                                noweld = nullptr)
{
    c.deck = deck;
    c.halt = halt;
    c.dense_map = map;
    c.census = census;
    c.noweld = noweld;
    /*  Parallel descend: census (and the THIN_DEBUG autopsy) are
     *  single-threaded instruments - they force the serial path.  */
    const int nt = mesh_threads();
    if (nt > 1 && !census && !getenv("STIBIUM_DMESH_THIN_DEBUG"))
    {
        std::vector<Region> tasks;
        collect_tasks(r, 2, tasks);
        std::vector<Collector> ws(tasks.size());
        std::atomic<size_t> next{ 0 };
        std::vector<std::thread> pool;
        const size_t use = std::min<size_t>(size_t(nt),
                                            tasks.size());
        for (size_t t = 0; t < use; ++t)
            pool.emplace_back([&]() {
                TapeCtx* ctx = tape_ctx_new(deck);
                for (;;)
                {
                    const size_t i = next.fetch_add(1);
                    prog_frac(float(next.load()) /
                              float(tasks.size()));
                    if (i >= tasks.size() || *halt)
                        break;
                    Collector& wc = ws[i];
                    wc.deck = deck;
                    wc.ctx = ctx;
                    wc.halt = halt;
                    wc.dense_map = map;
                    wc.census = nullptr;
                    wc.noweld = noweld;
                    descend(wc, tasks[i], deck_base(deck));
                }
                tape_ctx_free(ctx);
            });
        for (auto& th : pool)
            th.join();
        for (size_t i = 0; i < ws.size(); ++i)
        {
            merge_collector(c, ws[i]);
            ws[i] = Collector();   // release per-worker soup early
        }
    }
    else
    {
        c.ctx = tape_ctx_new(deck);
        descend(c, r, deck_base(deck));
        tape_ctx_free(c.ctx);
        c.ctx = nullptr;
    }
    c.ctx = tape_ctx_new(deck);
    bisect_edges(c);
    feature_points(c);
    c.soup.spacing = c.spacing;
    tape_ctx_free(c.ctx);
}

/*  Field probe (STIBIUM_DMESH_FPROBE="x,y,z;x,y,z;..."): print f
 *  and |f|/|grad| at listed points and continue.  The instrument
 *  that settles is-it-really-air arguments with the oracle
 *  instead of a fitted surface model (plinth teeth, 2026-07-18).  */
static void field_probe(const Deck* deck)
{
    const char* pe = getenv("STIBIUM_DMESH_FPROBE");
    if (!pe || !*pe)
        return;
    TapeCtx* ctx = tape_ctx_new(deck);
    std::vector<float> px, py, pz;
    float x, y, z;
    const char* s = pe;
    while (sscanf(s, "%f,%f,%f", &x, &y, &z) == 3)
    {
        px.push_back(x);
        py.push_back(y);
        pz.push_back(z);
        const char* n = strchr(s, ';');
        if (!n)
            break;
        s = n + 1;
    }
    const float h = 0.001f;
    std::vector<float> qx, qy, qz, qv;
    for (size_t i = 0; i < px.size(); ++i)
        for (int q = 0; q < 7; ++q)
        {
            qx.push_back(px[i] + (q==1 ? h : q==2 ? -h : 0));
            qy.push_back(py[i] + (q==3 ? h : q==4 ? -h : 0));
            qz.push_back(pz[i] + (q==5 ? h : q==6 ? -h : 0));
        }
    eval_points_mt(deck, ctx, qx, qy, qz, qv);
    for (size_t i = 0; i < px.size(); ++i)
    {
        const float f = qv[i*7];
        const float gx = qv[i*7+1] - qv[i*7+2];
        const float gy = qv[i*7+3] - qv[i*7+4];
        const float gz = qv[i*7+5] - qv[i*7+6];
        const float gl = sqrtf(gx*gx + gy*gy + gz*gz) / (2 * h);
        fprintf(stderr, "FPROBE (%.4f, %.4f, %.4f): f=%.6g "
                "|f|/|g|=%.6g\n", px[i], py[i], pz[i], f,
                gl > 1e-12f ? fabsf(f) / gl : -1.f);
    }
    tape_ctx_free(ctx);
}

DSoup delaunay_sample(const Deck* deck, Region r, volatile int* halt,
                      const std::unordered_map<uint64_t, int>* demote,
                      const std::vector<std::array<float, 3>>*
                              noweld,
                      const std::unordered_map<uint64_t, int>*
                              promote)
{
    Census census;
    const char* cenv = getenv("STIBIUM_DMESH_CENSUS");
    if (cenv && *cenv)
    {
        census.on = true;
        if (strcmp(cenv, "1") != 0)
            census.dump = fopen(cenv, "w");
    }

    Collector c;
    field_probe(deck);
    prog_stage(1);
    sample_pass(c, deck, r, halt, nullptr,
                census.on ? &census : nullptr, noweld);
    if (census.on)
        census_print(census, c.soup, "CENSUS");
    if (c.soup.welded && (census.on ||
                          getenv("STIBIUM_DMESH_TIME") ||
                          getenv("STIBIUM_DMESH_CHIP_DEBUG")))
        fprintf(stderr, "WELD: %llu air samples fused (sub-bar "
                "clearances)\n",
                (unsigned long long)c.soup.welded);
    if (c.soup.thinned && (census.on ||
                           getenv("STIBIUM_DMESH_TIME")))
        fprintf(stderr, "THIN: %llu witnesses dropped, %zu kept\n",
                (unsigned long long)c.soup.thinned,
                c.soup.samples.size());

    /*  Stage-D drill-down: pass 1 is also the survey.  When any
     *  leaf triggered (QEF residual over the bar, or a hidden
     *  crease-suspect leaf), re-run stage A with the map; flagged
     *  leaves sample a midpoint-refined lattice at their chosen
     *  level.  One drill-down round (the level formula already
     *  predicts the need from the residual magnitude); leftover
     *  pass-2 requests are reported, not chased.  */
    /*  Hidden-oracle verdicts: counter line + (env) the verdict
     *  boxes as STL - certified sub-lattice features in one
     *  file, proven grazes in the other, viewable next to the
     *  mesh exactly like the chain tubes
     *  (STIBIUM_DMESH_DUMP_HIDDEN=prefix ->
     *  prefix_feature.stl / prefix_graze.stl).  */
    if ((c.hidden_feature || c.hidden_graze ||
         c.hidden_contact || c.hidden_cell_feat) &&
        (census.on || getenv("STIBIUM_DMESH_TIME") ||
         getenv("STIBIUM_DMESH_CHIP_DEBUG")))
    {
        fprintf(stderr, "HIDDEN: %llu candidates -> %llu "
                "certified features, %llu contact (reach-"
                "collapse), %llu proven grazes\n",
                (unsigned long long)c.soup.hidden_candidates,
                (unsigned long long)c.hidden_feature,
                (unsigned long long)c.hidden_contact,
                (unsigned long long)c.hidden_graze);
        if (c.hidden_cell_feat)
            fprintf(stderr, "HIDDEN cells: %llu thin cells "
                    "certified in mixed leaves, %llu witness "
                    "bisections minted\n",
                    (unsigned long long)c.hidden_cell_feat,
                    (unsigned long long)c.hidden_cell_wit);
    }
    if (const char* hp = getenv("STIBIUM_DMESH_DUMP_HIDDEN"))
    {
        const auto wboxes = [](const std::string& path,
                const std::vector<std::array<float, 6>>& bs) {
            FILE* f = fopen(path.c_str(), "wb");
            if (!f)
                return;
            char hdr[80] = { 0 };
            fwrite(hdr, 1, 80, f);
            const uint32_t nt = uint32_t(bs.size() * 12);
            fwrite(&nt, 4, 1, f);
            static const int F[12][3][3] = {
                /* corner bit order: x=1, y=2, z=4 */
                { {0,0,0},{0,1,1},{0,0,1} }, { {0,0,0},{0,1,0},{0,1,1} },
                { {1,0,0},{1,0,1},{1,1,1} }, { {1,0,0},{1,1,1},{1,1,0} },
                { {0,0,0},{1,0,1},{1,0,0} }, { {0,0,0},{0,0,1},{1,0,1} },
                { {0,1,0},{1,1,0},{1,1,1} }, { {0,1,0},{1,1,1},{0,1,1} },
                { {0,0,0},{1,0,0},{1,1,0} }, { {0,0,0},{1,1,0},{0,1,0} },
                { {0,0,1},{0,1,1},{1,1,1} }, { {0,0,1},{1,1,1},{1,0,1} },
            };
            for (const auto& b : bs)
                for (int t = 0; t < 12; ++t)
                {
                    float rec[12] = { 0, 0, 0 };
                    for (int v = 0; v < 3; ++v)
                    {
                        rec[3 + v*3]     = F[t][v][0] ? b[3] : b[0];
                        rec[3 + v*3 + 1] = F[t][v][1] ? b[4] : b[1];
                        rec[3 + v*3 + 2] = F[t][v][2] ? b[5] : b[2];
                    }
                    fwrite(rec, 4, 12, f);
                    const uint16_t attr = 0;
                    fwrite(&attr, 2, 1, f);
                }
            fclose(f);
        };
        wboxes(std::string(hp) + "_feature.stl",
               c.hidden_feat_boxes);
        wboxes(std::string(hp) + "_contact.stl",
               c.hidden_contact_boxes);
        wboxes(std::string(hp) + "_graze.stl",
               c.hidden_graze_boxes);
        fprintf(stderr, "HIDDEN dump: %zu feature + %zu contact "
                "+ %zu graze boxes -> "
                "%s_{feature,contact,graze}.stl\n",
                c.hidden_feat_boxes.size(),
                c.hidden_contact_boxes.size(),
                c.hidden_graze_boxes.size(), hp);
    }
    /*  #4a: march the contact curves (snap targets + the vertex
     *  rail injected later, after delaunay_trace - injecting
     *  here would land inside the crease tracer's seed-tail
     *  window and die with pass-1's soup anyway).  */
    if (!(*halt))
        trace_contact_chains(c.deck, c.soup.spacing,
                             c.hidden_contact_boxes,
                             &c.soup.contact_chains);
    /*  Level-3 crowding grant, TWO-SIGNAL rule (MAX=3 referee
     *  round, 2026-07-19): blanket level 3 wrecks organic models
     *  (bino: worst 0.093 -> 0.372, nm 6.4x) and live-count
     *  magnitude alone cannot target it - bino's blend bands
     *  honestly read live >= 64 in 1,288 leaves.  The measured
     *  discriminator is the COINCIDENCE of signals: screws-class
     *  sub-lattice slots are extreme-crowded AND residual-hot
     *  (92% of screws grants), while assembly collars are
     *  crowded but residual-silent (86% of bino's false grants).
     *  So level 3 requires both: live >= 4x the bar AND the
     *  residual formula already flagged the leaf >= 2.  Ceiling
     *  crowd_max_level() (default 3), refereed on screws (tilt
     *  4.03 -> 3.04% of area, worst held, 0 geometric opens).  */
    if (autodense() && crowd_max_level() >= 3)
        for (auto& [key, lv] : c.want_dense)
        {
            if (lv < 2 || lv >= 3 || !c.resid_hot.count(key))
                continue;
            const auto it = c.crease_leaves.find(key);
            if (it != c.crease_leaves.end() &&
                it->second.live >= 4u * autod_live_bar() &&
                !(it->second.sep <
                  autod_sep_bar() * c.spacing))
                lv = 3;
        }
    if (autodense() && !c.want_dense.empty() && !(*halt))
    {
        const bool dbg = census.on ||
                getenv("STIBIUM_DMESH_TIME") != nullptr;

        /*  Graduated coverage (every configuration measured
         *  2026-07-16 taught the same lesson: coverage must be
         *  CONTIGUOUS - pitch-band borders may only sit on quiet
         *  geometry).  Level 1 floods the whole connected
         *  crease-suspect component containing any flag - exactly
         *  the global-d1 treatment, but only for components that
         *  asked; sharp models raise no flag and pay nothing.
         *  Level-2 cores then dilate one ring into touching
         *  non-tangle crease leaves, so every seam is at most a
         *  2x jump and never lands mid-tangle.  */
        const float eps = 0.5f * c.soup.spacing;
        const auto touches = [&](const Collector::LeafBox& a,
                                 const Collector::LeafBox& b) {
            return !(b.lo[0] > a.hi[0] + eps ||
                     b.hi[0] < a.lo[0] - eps ||
                     b.lo[1] > a.hi[1] + eps ||
                     b.hi[1] < a.lo[1] - eps ||
                     b.lo[2] > a.hi[2] + eps ||
                     b.hi[2] < a.lo[2] - eps);
        };

        std::unordered_map<uint64_t, int> dilated = c.want_dense;
        /*  BFS flood at level 1 through box adjacency.  */
        std::vector<uint64_t> wave;
        wave.reserve(c.want_dense.size());
        for (const auto& [key, lvl] : c.want_dense)
            wave.push_back(key);
        std::unordered_set<uint64_t> visited(wave.begin(),
                                             wave.end());
        while (!wave.empty())
        {
            std::vector<uint64_t> next;
            for (const auto& [k2, b2] : c.crease_leaves)
            {
                if (visited.count(k2))
                    continue;
                for (const uint64_t k : wave)
                {
                    const auto it = c.crease_leaves.find(k);
                    if (it != c.crease_leaves.end() &&
                        touches(it->second, b2))
                    {
                        visited.insert(k2);
                        next.push_back(k2);
                        int& lv = dilated[k2];
                        lv = std::max(lv, 1);
                        break;
                    }
                }
            }
            wave = std::move(next);
        }
        /*  Level-2 core dilation (one ring, non-tangle only).  */
        for (const auto& [key, lvl] : c.want_dense)
        {
            if (lvl < 2)
                continue;
            const auto it = c.crease_leaves.find(key);
            if (it == c.crease_leaves.end())
                continue;
            for (const auto& [k2, b2] : c.crease_leaves)
            {
                if (b2.sep < autod_sep_bar() * c.soup.spacing ||
                    !touches(it->second, b2))
                    continue;
                int& l2 = dilated[k2];
                l2 = std::max(l2, lvl);
            }
        }
        /*  Smooth-pocket fill (the wireframe splotch class,
         *  2026-07-18): a live = 0 leaf is invisible to the
         *  flood - it walks crease_leaves only - so smooth
         *  pockets inside promoted neighborhoods stay level 0
         *  and render as coarse polygonal splotches on curved
         *  walls.  Promote any smooth leaf touching >= 3
         *  promoted leaves to the MIN touching level (contrast
         *  removal, not extra detail); boundary strips of big
         *  smooth regions touch fewer and stay untouched.  */
        {
            std::vector<std::pair<uint64_t, int>> fills;
            for (const auto& [k2, b2] : c.smooth_leaves)
            {
                int cnt = 0, lvmin = INT_MAX;
                for (const auto& [k3, l3v] : dilated)
                {
                    const auto it3 = c.crease_leaves.find(k3);
                    if (it3 == c.crease_leaves.end() ||
                        !touches(it3->second, b2))
                        continue;
                    ++cnt;
                    lvmin = std::min(lvmin, std::max(1, l3v));
                }
                if (cnt >= 3)
                    fills.push_back({ k2, lvmin });
            }
            for (const auto& [k2, lv2] : fills)
                dilated[k2] = lv2;
            if (dbg && !fills.empty())
                fprintf(stderr, "AUTOD smooth-fill: %zu of %zu "
                        "level-0 pockets promoted\n",
                        fills.size(), c.smooth_leaves.size());
        }
        /*  Tangle-demotion DESPECKLE (the wireframe splotch
         *  class, cracked 2026-07-18): bino is blanket level-2
         *  yet walls carry random polygonal 2x-coarse patches -
         *  they are tangle-demoted islands (404 on bino), minted
         *  wherever the per-leaf sep reading (a min over corner
         *  pairs - noisy) flickers across the bar on a thin
         *  wall.  A REAL tangle band is contiguous: its
         *  neighbors demote too, and a majority filter leaves it
         *  alone.  Only isolated speckles - level-1 leaves with
         *  >= 5 touching level-2 neighbors - re-promote.  Runs
         *  before the retreat demote so convictions still land
         *  on top.  */
        {
            std::vector<uint64_t> speckle;
            for (const auto& [k2, l2v] : dilated)
            {
                if (l2v >= 2)
                    continue;
                const auto it2 = c.crease_leaves.find(k2);
                if (it2 == c.crease_leaves.end())
                    continue;
                int hi2 = 0;
                for (const auto& [k3, l3v] : dilated)
                {
                    if (l3v < 2 || k3 == k2)
                        continue;
                    const auto it3 = c.crease_leaves.find(k3);
                    if (it3 != c.crease_leaves.end() &&
                        touches(it3->second, it2->second))
                        ++hi2;
                }
                if (hi2 >= 5)
                    speckle.push_back(k2);
            }
            for (const uint64_t k2 : speckle)
                dilated[k2] = 2;
            if (dbg && !speckle.empty())
                fprintf(stderr, "AUTOD despeckle: %zu demoted "
                        "islands re-promoted\n", speckle.size());
        }
        /*  Close-ring strip promotion (level 3 between paired
         *  rims) lands before the demote cap so retreat can still
         *  roll a torn strip back.  */
        if (promote)
            for (const auto& [k, v] : *promote)
            {
                int& lv = dilated[k];
                lv = std::max(lv, v);
            }
        /*  Retreat rollback, GRADUATED (2026-07-18): each
         *  conviction steps the leaf down ONE level, not to the
         *  floor - a level-3 strip leaf dropped straight to flood
         *  put an 8x pitch cliff mid-band and turned 4 opens into
         *  12 (zeiss autod21 attempt 3).  Contiguity is law for
         *  rollback too; repeat convictions keep stepping.  */
        if (demote)
            for (auto& [k, l] : dilated)
            {
                const auto it = demote->find(k);
                if (it != demote->end() && l >= 2)
                    l = std::max(1, l - it->second);
            }
        if (dbg)
        {
            size_t l1 = 0, l2 = 0, l3 = 0;
            for (const auto& [k, l] : dilated)
                ++(l >= 3 ? l3 : l >= 2 ? l2 : l1);
            fprintf(stderr, "CURVE: %llu seen, %llu probed, "
                    "%llu flagged, max theta %.2f deg/cell "
                    "(bar %.2f)\n",
                    (unsigned long long)c.curve_seen,
                    (unsigned long long)c.curve_cross,
                    (unsigned long long)c.curve_flagged,
                    c.curve_theta_max, curve_bar_deg());
            fprintf(stderr, "AUTOD drill-down: %zu leaves flagged "
                    "(%llu curve), "
                    "%zu after flood+cores (%zu @1, %zu @2, "
                    "%zu @3), %llu tangle-demoted, "
                    "%llu phantom QEF\n",
                    c.want_dense.size(),
                    (unsigned long long)c.curve_flagged,
                    dilated.size(), l1, l2,
                    l3, (unsigned long long)c.tangle_suppressed,
                    (unsigned long long)c.phantom_rejected);
        }
        /*  Final density map (census dump): one DENS row per
         *  granted leaf - the offline referee for density-shadow
         *  hunts (quilted patches are LEVEL HOLES; this is the
         *  map that says which trigger never fired).  Crease
         *  leaves not in `dilated` are base-level by omission.  */
        if (census.dump)
        {
            for (const auto& [k, l] : dilated)
            {
                const auto it = c.crease_leaves.find(k);
                if (it == c.crease_leaves.end())
                    continue;
                fprintf(census.dump, "DENS %.6g %.6g %.6g %d C\n",
                        0.5f * (it->second.lo[0] +
                                it->second.hi[0]),
                        0.5f * (it->second.lo[1] +
                                it->second.hi[1]),
                        0.5f * (it->second.lo[2] +
                                it->second.hi[2]), l);
            }
            for (const auto& [k, b] : c.smooth_leaves)
            {
                const auto it = dilated.find(k);
                fprintf(census.dump, "DENS %.6g %.6g %.6g %d S\n",
                        0.5f * (b.lo[0] + b.hi[0]),
                        0.5f * (b.lo[1] + b.hi[1]),
                        0.5f * (b.lo[2] + b.hi[2]),
                        it == dilated.end() ? 0 : it->second);
            }
        }

        Census census2;
        census2.on = census.on;
        /*  Perf review, the cheap first win: pass 1's soup and
         *  index structures are DEAD once the drill-down map is
         *  built (pass 2 re-samples everything), but they stayed
         *  resident through pass 2 - a full duplicate of the
         *  biggest transients at peak.  Release before c2
         *  samples.  */
        c.soup.samples.clear();
        c.soup.samples.shrink_to_fit();
        c.soup.surface.clear();
        c.soup.surface.shrink_to_fit();
        c.seen_samples.clear();
        c.seen_samples.rehash(0);
        c.edge_index.clear();
        c.edge_index.rehash(0);
        c.edges.clear();
        c.edges.shrink_to_fit();
        c.cells.clear();
        c.cells.shrink_to_fit();
        Collector c2;
        /*  Pass 2 skips the leaf gradient survey (dense_map set),
         *  but the cell judges - the shallow-seed gate and the
         *  tangle flag - read the crease-leaf map, and the soup
         *  that SHIPS is pass 2's.  Same leaf boxes, same keys:
         *  hand pass 1's survey over (without this the shallow
         *  channel judged against an empty map and minted zero
         *  seeds - measured, bino 2026-07-17).  */
        c2.crease_leaves = std::move(c.crease_leaves);
        prog_stage(2);
        sample_pass(c2, deck, r, halt, &dilated,
                    census2.on ? &census2 : nullptr, noweld);
        if (census2.on)
            census_print(census2, c2.soup, "CENSUS2");
        if (c2.soup.thinned && dbg)
            fprintf(stderr, "THIN2: %llu witnesses dropped, "
                    "%zu kept\n",
                    (unsigned long long)c2.soup.thinned,
                    c2.soup.samples.size());
        if (c2.hidden_cell_feat &&
            (dbg || getenv("STIBIUM_DMESH_CHIP_DEBUG")))
            fprintf(stderr, "HIDDEN cells: %llu thin cells "
                    "certified in mixed leaves, %llu witness "
                    "bisections minted\n",
                    (unsigned long long)c2.hidden_cell_feat,
                    (unsigned long long)c2.hidden_cell_wit);
        if (dbg && !c2.want_dense.empty())
            fprintf(stderr, "AUTOD residue: %zu leaves still hot "
                    "after drill-down\n", c2.want_dense.size());
        if (census.dump)
            fclose(census.dump);
        /*  Contact chains were traced on the pass-1 collector;
         *  the drill-down soup is the one that ships.  */
        c2.soup.contact_chains = std::move(c.soup.contact_chains);
        return std::move(c2.soup);
    }

    if (census.dump)
        fclose(census.dump);
    return std::move(c.soup);
}

////////////////////////////////////////////////////////////////////////////
//  The crease tracer: CSG-driven surface-surface intersection marching
//  (doc/MESH-NEXT.md crease-tracing round; research in
//  doc/research/2026-07-15-junction-extraction.md).
//
//  Every min/max clause in the tape is a crease generator: the
//  surface kinks exactly where its two operands are equal (and both
//  zero, and the branch is active up to the root).  The tape's
//  prefix evaluation hands us f_A and f_B - with every upstream
//  coordinate transform already applied - so the crease is the
//  solution set {f_A = 0, f_B = 0}, marched with the classical SSI
//  predictor-corrector: tangent = grad A x grad B, corrector = 2-eq
//  Newton onto both zero sets.  The full oracle trims the curve to
//  the actual boundary (past a junction the crease continues but a
//  third branch owns the surface, so |f| lifts off zero - the trim
//  boundary IS the junction, and we bisect onto it so incident
//  edges share their corner).

namespace {

struct CreaseTracer
{
    const Deck* deck;
    TapeCtx* ctx;
    Tape* base;
    TapePair tp;
    float sp;                        // lattice spacing
    float step;                      // march step (0.5 sp / dense)
    float bmin[3], bmax[3];          // region bounds + margin

    /*  Law delivery (density campaign, 2026-07-18): the march
     *  step was GLOBAL, so lattice density outran constraint
     *  density inside stage-D bands - a level-2 screw bore
     *  walks in six half-cell strides while its lattice is
     *  quarter-cell, and the chain dies where r2's survives
     *  (screws: 339 constraints at r1+level2 vs 681 at r2;
     *  churn triples to fill the law gap).  The 2026-07-15
     *  "step scaling refuted" verdict was GLOBAL halving on
     *  at-scale features; local band scaling is a new claim,
     *  re-refereed.  Dense leaves quantized on the uniform leaf
     *  grid; lookup is O(1) per march step.
     *  STIBIUM_DMESH_TRACE_LOCAL=0 disables.  */
    std::unordered_map<uint64_t, int> band;
    float leaf_edge = 0;
    float borigin[3] = { 0, 0, 0 };

    static uint64_t band_key(int qx, int qy, int qz)
    {
        return (uint64_t(uint32_t(qx + (1 << 20))) << 42) |
               (uint64_t(uint32_t(qy + (1 << 20))) << 21) |
                uint64_t(uint32_t(qz + (1 << 20)));
    }

    float local_step(const float p[3]) const
    {
        if (band.empty() || !(leaf_edge > 0))
            return step;
        const auto it = band.find(band_key(
                int(std::floor((p[0] - borigin[0]) / leaf_edge)),
                int(std::floor((p[1] - borigin[1]) / leaf_edge)),
                int(std::floor((p[2] - borigin[2]) / leaf_edge))));
        return it == band.end()
                ? step : step / float(1 << it->second);
    }

    /*  7-tap probe at p: pair values/gradients (prefix eval) and
     *  full-surface distance |f|/|grad f| (full eval).  */
    struct Probe
    {
        float fa, fb;
        float ga[3], gb[3];
        float fdist;                 // estimated distance to surface
        bool finite;
    };

    Probe probe(const float p[3]) const
    {
        const float h = 0.01f * sp;
        float xs[7], ys[7], zs[7];
        for (int q = 0; q < 7; ++q)
        {
            xs[q] = p[0];
            ys[q] = p[1];
            zs[q] = p[2];
        }
        xs[1] += h;  xs[2] -= h;
        ys[3] += h;  ys[4] -= h;
        zs[5] += h;  zs[6] -= h;

        Region dummy = {};
        dummy.voxels = 7;
        dummy.X = xs;
        dummy.Y = ys;
        dummy.Z = zs;

        Probe out;
        tape_eval_r_prefix(base, ctx, dummy, tp.clause);
        const float* ra = tape_ctx_r_row(ctx, tp.slot_a);
        out.fa = ra[0];
        const float h2 = 2 * h;
        out.ga[0] = (ra[1] - ra[2]) / h2;
        out.ga[1] = (ra[3] - ra[4]) / h2;
        out.ga[2] = (ra[5] - ra[6]) / h2;
        float fa7[7], fb7[7];
        memcpy(fa7, ra, sizeof(fa7));
        if (tp.is_max != 2)
        {
            /*  min/max: field B is the sibling operand  */
            const float* rb = tape_ctx_r_row(ctx, tp.slot_b);
            out.fb = rb[0];
            out.gb[0] = (rb[1] - rb[2]) / h2;
            out.gb[1] = (rb[3] - rb[4]) / h2;
            out.gb[2] = (rb[5] - rb[6]) / h2;
            memcpy(fb7, rb, sizeof(fb7));
        }

        const float* rf = tape_eval_r(base, ctx, dummy);
        const float f = rf[0];
        const float gfx = (rf[1] - rf[2]) / h2;
        const float gfy = (rf[3] - rf[4]) / h2;
        const float gfz = (rf[5] - rf[6]) / h2;
        const float gfl = sqrtf(gfx*gfx + gfy*gfy + gfz*gfz);
        out.fdist = gfl > 0 ? fabsf(f) / gfl : 1e30f;
        if (tp.is_max == 2)
        {
            /*  abs(g): the kink locus is {g = 0} INTERSECT the
             *  surface, so field B is the FULL oracle - already
             *  evaluated for the trim distance.  */
            out.fb = f;
            out.gb[0] = gfx;
            out.gb[1] = gfy;
            out.gb[2] = gfz;
            memcpy(fb7, rf, sizeof(fb7));
        }
        out.finite = std::isfinite(out.fa) && std::isfinite(out.fb) &&
                std::isfinite(f);
        for (int q = 0; q < 7 && out.finite; ++q)
            if (!std::isfinite(fa7[q]) || !std::isfinite(fb7[q]))
                out.finite = false;
        return out;
    }

    /*  Kink-activity test for abs creases: the |f| trim can never
     *  fire on {g = 0, f = 0} (the curve is ON the surface by
     *  construction), so the trim analog asks whether the surface
     *  is actually creased here - one-sided full-oracle gradients
     *  across the g = 0 sheet must disagree.  A smooth crossing
     *  (the abs branch pruned away or dominated elsewhere) reads
     *  parallel gradients and the march stops.  */
    bool kinked(const float p[3], const Probe& pr) const
    {
        const float gal = sqrtf(pr.ga[0]*pr.ga[0] +
                pr.ga[1]*pr.ga[1] + pr.ga[2]*pr.ga[2]);
        if (!(gal > 0))
            return false;
        const float e = 0.05f * sp;
        float gl3[3], gr3[3];
        for (int side = 0; side < 2; ++side)
        {
            const float s = side ? e : -e;
            const float q[3] = { p[0] + s * pr.ga[0] / gal,
                                 p[1] + s * pr.ga[1] / gal,
                                 p[2] + s * pr.ga[2] / gal };
            float xs[7], ys[7], zs[7];
            const float h = 0.01f * sp;
            for (int u = 0; u < 7; ++u)
            {
                xs[u] = q[0];
                ys[u] = q[1];
                zs[u] = q[2];
            }
            xs[1] += h;  xs[2] -= h;
            ys[3] += h;  ys[4] -= h;
            zs[5] += h;  zs[6] -= h;
            Region dummy = {};
            dummy.voxels = 7;
            dummy.X = xs;
            dummy.Y = ys;
            dummy.Z = zs;
            const float* rf = tape_eval_r(base, ctx, dummy);
            float* g3 = side ? gr3 : gl3;
            const float h2 = 2 * h;
            g3[0] = (rf[1] - rf[2]) / h2;
            g3[1] = (rf[3] - rf[4]) / h2;
            g3[2] = (rf[5] - rf[6]) / h2;
            const float l = sqrtf(g3[0]*g3[0] + g3[1]*g3[1] +
                                  g3[2]*g3[2]);
            if (!(l > 0) || !std::isfinite(l))
                return false;
            g3[0] /= l;
            g3[1] /= l;
            g3[2] /= l;
        }
        return gl3[0]*gr3[0] + gl3[1]*gr3[1] + gl3[2]*gr3[2]
                < 0.9f;
    }

    /*  Newton onto {f_A = 0, f_B = 0}: minimal-norm step
     *  dx = J^T (J J^T)^-1 (-F).  False on a singular J (the
     *  gradients are parallel: tangency) or non-convergence.  */
    bool correct(float p[3], Probe* last = nullptr) const
    {
        for (int it = 0; it < 6; ++it)
        {
            const Probe pr = probe(p);
            if (!pr.finite)
                return false;
            const float a = pr.ga[0]*pr.ga[0] + pr.ga[1]*pr.ga[1] +
                            pr.ga[2]*pr.ga[2];
            const float b = pr.gb[0]*pr.gb[0] + pr.gb[1]*pr.gb[1] +
                            pr.gb[2]*pr.gb[2];
            const float c = pr.ga[0]*pr.gb[0] + pr.ga[1]*pr.gb[1] +
                            pr.ga[2]*pr.gb[2];
            /*  Tangency gate with teeth (2026-07-17, the additive-
             *  joint autopsy): coincident-surface seams - the same
             *  wall computed through two arithmetic paths, e.g. an
             *  eyepiece fitted exactly into a bore - degenerate the
             *  kink locus into a 2D sheet.  Rounding noise keeps
             *  the operand gradients JUST non-parallel (the old
             *  1e-6 bar = accepts pairs meeting at ~0.06 deg), so
             *  Newton converges anywhere on the sheet and the
             *  march wanders along noise, hallucinating constraint
             *  polylines that become the visible interference
             *  rings.  det > 0.01*a*b demands the surfaces meet at
             *  >~5.7 deg - real creases are 30-90 deg, rounding
             *  ghosts are thousandths.  */
            const float det = a*b - c*c;
            if (!(det > 1e-2f * a * b) || !(a > 0) || !(b > 0))
                return false;                    // tangency/ghost
            const float u = (-pr.fa * b + pr.fb * c) / det;
            const float v = (-pr.fb * a + pr.fa * c) / det;
            float dx = pr.ga[0]*u + pr.gb[0]*v;
            float dy = pr.ga[1]*u + pr.gb[1]*v;
            float dz = pr.ga[2]*u + pr.gb[2]*v;
            const float dl = sqrtf(dx*dx + dy*dy + dz*dz);
            const float cap = 0.5f * sp;
            if (dl > cap)
            {
                dx *= cap / dl;
                dy *= cap / dl;
                dz *= cap / dl;
            }
            p[0] += dx;
            p[1] += dy;
            p[2] += dz;
            if (dl < 1e-4f * sp)
            {
                if (last)
                    *last = pr;
                return true;
            }
        }
        const Probe pr = probe(p);
        if (last)
            *last = pr;
        const float gal = sqrtf(pr.ga[0]*pr.ga[0] +
                pr.ga[1]*pr.ga[1] + pr.ga[2]*pr.ga[2]);
        const float gbl = sqrtf(pr.gb[0]*pr.gb[0] +
                pr.gb[1]*pr.gb[1] + pr.gb[2]*pr.gb[2]);
        return pr.finite && gal > 0 && gbl > 0 &&
               fabsf(pr.fa) / gal < 0.02f * sp &&
               fabsf(pr.fb) / gbl < 0.02f * sp;
    }

    bool in_bounds(const float p[3]) const
    {
        for (int q = 0; q < 3; ++q)
            if (p[q] < bmin[q] || p[q] > bmax[q])
                return false;
        return true;
    }

    static void tangent(const Probe& pr, float t[3])
    {
        t[0] = pr.ga[1]*pr.gb[2] - pr.ga[2]*pr.gb[1];
        t[1] = pr.ga[2]*pr.gb[0] - pr.ga[0]*pr.gb[2];
        t[2] = pr.ga[0]*pr.gb[1] - pr.ga[1]*pr.gb[0];
        const float l = sqrtf(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
        if (l > 0)
        {
            t[0] /= l;
            t[1] /= l;
            t[2] /= l;
        }
    }

    /*  March one direction from p0; appends corrected points.
     *  Returns true when the walk closed back onto p0.  */
    bool march(const float p0[3], const float t0[3], int dir,
               std::vector<DSurfPoint>& out_pts,
               volatile int* halt) const
    {
        float p[3] = { p0[0], p0[1], p0[2] };
        float tprev[3] = { dir * t0[0], dir * t0[1], dir * t0[2] };
        for (int steps = 0; steps < 16384 && !*halt; ++steps)
        {
            const float step = local_step(p);
            /*  The trim is per-STRIDE sensitivity: a coarse
             *  stride leaps past a junction's lift-off zone and
             *  bisects cleanly, but a fine stride CRAWLS through
             *  it - with a fixed trim it emits points with fdist
             *  just under the bar and never fires the bisect
             *  ([.dtrace] csg: worst |f| 1.9e-3, the same number
             *  the 2026-07-15 global experiment measured).  Scale
             *  the trim with the stride and fine bands tighten
             *  their own tolerance.  */
            const float trim = 0.1f * step;
            float pt[3] = { p[0] + step * tprev[0],
                            p[1] + step * tprev[1],
                            p[2] + step * tprev[2] };
            Probe pr;
            if (!correct(pt, &pr) || !in_bounds(pt))
                break;
            /*  Minimum progress: near a tangency the corrector
             *  pulls every predictor step back to the same point -
             *  emitting float-noise micro-segments (measured: three
             *  "points" 1e-15 apart on csg's seam).  No progress =
             *  the curve ends here.  */
            {
                const float px2 = pt[0] - p[0], py2 = pt[1] - p[1],
                            pz2 = pt[2] - p[2];
                if (px2*px2 + py2*py2 + pz2*pz2 <
                    0.01f * step * step)
                    break;
            }
            /*  abs trim analog: fdist can't lift (field B IS the
             *  oracle) - the crease ends where the surface stops
             *  being creased.  No corner bisection yet: the
             *  overshoot is at most one step of on-surface
             *  polyline, which constrains harmlessly.  */
            if (tp.is_max == 2 && !kinked(pt, pr))
                break;
            if (pr.fdist > trim)
            {
                /*  Left the boundary: a third branch owns the
                 *  surface past here - the junction.  Bisect [p,
                 *  pt] for the corner so every incident crease
                 *  shares it.  On the boundary side fdist is the
                 *  Newton residual (~0); past the corner it grows
                 *  smoothly - so the bisection target must be far
                 *  TIGHTER than the marching trim, or endpoints
                 *  stop trim-short of the corner (measured: worst
                 *  |f| exactly 0.05*sp).  */
                const float ctol = 5e-4f * sp;
                float lo = 0, hi = 1;
                for (int r2 = 0; r2 < 16; ++r2)
                {
                    const float mid = (lo + hi) * 0.5f;
                    float pm[3] = { p[0] + mid * (pt[0] - p[0]),
                                    p[1] + mid * (pt[1] - p[1]),
                                    p[2] + mid * (pt[2] - p[2]) };
                    Probe pm2;
                    if (correct(pm, &pm2) && pm2.fdist <= ctol)
                        lo = mid;
                    else
                        hi = mid;
                }
                if (lo > 0.05f)
                {
                    float pe[3] = { p[0] + lo * (pt[0] - p[0]),
                                    p[1] + lo * (pt[1] - p[1]),
                                    p[2] + lo * (pt[2] - p[2]) };
                    /*  Femto-segment guard (cascade autopsy,
                     *  2026-07-17): when the march has already
                     *  reached the corner, the bisected endpoint
                     *  lands within float noise of the previous
                     *  point - the chain ends in a ~1e-14
                     *  segment, and the CCDT throws 'constraint
                     *  passes through a vertex' on it.  ONE such
                     *  segment cost a full 660 s rebuild on the
                     *  zeiss (twice - its mirror twin too).  The
                     *  corner must still be reached, but never as
                     *  a degenerate hop.  */
                    const float ddx = pe[0] - p[0],
                                ddy = pe[1] - p[1],
                                ddz = pe[2] - p[2];
                    if (correct(pe) &&
                        ddx*ddx + ddy*ddy + ddz*ddz >
                                1e-10f * sp * sp)
                        out_pts.push_back({ pe[0], pe[1], pe[2] });
                }
                break;
            }
            float t[3];
            tangent(pr, t);
            if (t[0]*tprev[0] + t[1]*tprev[1] + t[2]*tprev[2] < 0)
            {
                t[0] = -t[0];
                t[1] = -t[1];
                t[2] = -t[2];
            }
            if (t[0]*tprev[0] + t[1]*tprev[1] + t[2]*tprev[2] < 0.5f)
                break;                            // hard turn: corner
            const float cx = pt[0] - p0[0], cy = pt[1] - p0[1],
                        cz = pt[2] - p0[2];
            if (steps > 3 &&
                cx*cx + cy*cy + cz*cz < 0.5f * step * step)
                return true;                      // closed the loop
            out_pts.push_back({ pt[0], pt[1], pt[2] });
            p[0] = pt[0];
            p[1] = pt[1];
            p[2] = pt[2];
            tprev[0] = t[0];
            tprev[1] = t[1];
            tprev[2] = t[2];
        }
        return false;
    }
};

}  // namespace

static void write_stl_raw(const char* path,
                          const std::vector<float>& V,
                          const std::vector<uint32_t>& T);
static void dump_chains_stl(const char* path, const DSoup& soup);
static void emit_tube(std::vector<float>& V,
                      std::vector<uint32_t>& T,
                      const float A[3], const float B[3],
                      float rr);

bool delaunay_trace(const Deck* deck, Region r, DSoup* soup,
                    volatile int* halt)
{
    Tape* base = deck_base(deck);
    const unsigned nmm = tape_pairs(base, nullptr, 0);
    const unsigned nabs = tape_abs_pairs(base, nullptr, 0);
    const unsigned npairs = nmm + nabs;
    if (!npairs ||
        (soup->feature_points == 0 && soup->tseeds.empty()))
        return false;
    std::vector<TapePair> pairs(npairs);
    tape_pairs(base, pairs.data(), nmm);
    tape_abs_pairs(base, pairs.data() + nmm, nabs);
    if (getenv("STIBIUM_DMESH_TRACE_DEBUG"))
        fprintf(stderr, "TRACE generators: %u min/max, %u abs\n",
                nmm, nabs);

    TapeCtx* ctx = tape_ctx_new(deck);
    const float sp = soup->spacing;
    const size_t f0 = soup->surface.size() -
            size_t(soup->feature_points);

    /*  Seed pool: the feature tail plus the shallow channel's
     *  tracer-only seeds (sub-25-degree creases in crease leaves;
     *  they exist nowhere else in the soup).  */
    std::vector<DSurfPoint> seeds(soup->surface.begin() + f0,
                                  soup->surface.end());
    seeds.insert(seeds.end(), soup->tseeds.begin(),
                 soup->tseeds.end());
    if (getenv("STIBIUM_DMESH_TRACE_DEBUG"))
        fprintf(stderr, "TRACE seeds: %llu feature + %zu shallow\n",
                (unsigned long long)soup->feature_points,
                soup->tseeds.size());

    CreaseTracer tr;
    tr.deck = deck;
    tr.ctx = ctx;
    tr.base = base;
    tr.sp = sp;
    /*  NOT scaled by the dense round (measured 2026-07-15): halving
     *  the march step left the chip chord distribution byte-
     *  identical (chords track the LATTICE pitch, not the polyline
     *  pitch), doubled the constraint count, and pushed the kink-
     *  straddle noise past the [.dtrace] bound (1.9e-3 vs 1.5e-3).
     *  The snap pass needs only the polyline GEOMETRY, whose chord
     *  sagitta at 0.5 sp is already ~1e-3 sp.  */
    tr.step = 0.5f * sp;
    tr.bmin[0] = r.X[0] - sp;
    tr.bmax[0] = r.X[r.ni] + sp;
    tr.bmin[1] = r.Y[0] - sp;
    tr.bmax[1] = r.Y[r.nj] + sp;
    tr.bmin[2] = r.Z[0] - sp;
    tr.bmax[2] = r.Z[r.nk] + sp;
    /*  Law delivery: march at the LOCAL pitch inside stage-D
     *  bands (see the CreaseTracer comment).  */
    static const char* tl_env = getenv("STIBIUM_DMESH_TRACE_LOCAL");
    if ((!tl_env || atoi(tl_env) != 0) && !soup->dense_boxes.empty())
    {
        tr.leaf_edge = soup->dense_boxes[0].hi[0] -
                       soup->dense_boxes[0].lo[0];
        tr.borigin[0] = r.X[0];
        tr.borigin[1] = r.Y[0];
        tr.borigin[2] = r.Z[0];
        if (tr.leaf_edge > 0)
            for (const auto& b : soup->dense_boxes)
            {
                if (b.level < 1)
                    continue;
                const uint64_t k = CreaseTracer::band_key(
                        int(std::floor((b.lo[0] - r.X[0]) /
                                       tr.leaf_edge + 0.5f)),
                        int(std::floor((b.lo[1] - r.Y[0]) /
                                       tr.leaf_edge + 0.5f)),
                        int(std::floor((b.lo[2] - r.Z[0]) /
                                       tr.leaf_edge + 0.5f)));
                int& lv = tr.band[k];
                lv = std::max(lv, b.level);
            }
        if (getenv("STIBIUM_DMESH_TRACE_DEBUG"))
            fprintf(stderr, "TRACE local-step band: %zu dense "
                    "leaves, edge %.3f\n", tr.band.size(),
                    tr.leaf_edge);
    }

    std::vector<std::vector<DSurfPoint>> polys;
    std::vector<uint8_t> closed;
    std::vector<uint32_t> poly_pair;

    /*  Seed gate (perf round, 2026-07-15: the tracer was 69% of a
     *  zeiss export - 957 s - because every pair ran full Newton
     *  from every seed, and almost every (pair, seed) match is
     *  impossible).  One recording walk per seed captures every
     *  pair's operand values; Newton then runs only where BOTH
     *  fields are already small.  The gate is generous (Newton
     *  moves seeds up to 1.5 cells and verifies everything) and
     *  STIBIUM_DMESH_SEEDGATE=0 restores the exhaustive loop.  */
    static const char* sg_env = getenv("STIBIUM_DMESH_SEEDGATE");
    const bool seedgate = !sg_env || atoi(sg_env) != 0;
    std::vector<std::vector<uint32_t>> pair_seeds(npairs);
    if (seedgate)
    {
        std::vector<TapePair> sorted(pairs);
        std::vector<uint32_t> order(npairs);
        for (unsigned i = 0; i < npairs; ++i)
            order[i] = i;
        std::sort(order.begin(), order.end(),
                  [&](uint32_t l, uint32_t r2) {
                      return pairs[l].clause < pairs[r2].clause;
                  });
        for (unsigned i = 0; i < npairs; ++i)
            sorted[i] = pairs[order[i]];
        std::vector<float> fa(npairs), fb(npairs);
        const float gate = 8.f * sp;
        for (size_t s = 0; s < seeds.size(); ++s)
        {
            const DSurfPoint& seed = seeds[s];
            tape_eval_f_pairs(base, ctx, seed.x, seed.y, seed.z,
                              sorted.data(), npairs,
                              fa.data(), fb.data());
            for (unsigned i = 0; i < npairs; ++i)
                if (std::isfinite(fa[i]) && std::isfinite(fb[i]) &&
                    fabsf(fa[i]) < gate && fabsf(fb[i]) < gate)
                    pair_seeds[order[i]].push_back(uint32_t(s));
        }
        if (getenv("STIBIUM_DMESH_TRACE_DEBUG"))
        {
            size_t total = 0;
            for (const auto& v : pair_seeds)
                total += v.size();
            fprintf(stderr, "TRACE seed gate: %zu candidate "
                    "(pair, seed) matches of %zu exhaustive\n",
                    total, size_t(npairs) * seeds.size());
        }
    }

    std::vector<uint8_t> consumed(seeds.size(), 0);
    for (unsigned pi = 0; pi < npairs; ++pi)
    {
        const TapePair& tp = pairs[pi];
        if (*halt)
            break;
        tr.tp = tp;
        const size_t nseeds = seedgate ? pair_seeds[pi].size()
                                       : consumed.size();
        for (size_t si = 0; si < nseeds; ++si)
        {
            const size_t s = seedgate ? size_t(pair_seeds[pi][si])
                                      : si;
            if (consumed[s] || *halt)
                continue;
            const DSurfPoint& seed = seeds[s];
            float p0[3] = { seed.x, seed.y, seed.z };
            CreaseTracer::Probe pr;
            if (!tr.correct(p0, &pr))
                continue;
            const float mx = p0[0] - seed.x, my = p0[1] - seed.y,
                        mz = p0[2] - seed.z;
            if (mx*mx + my*my + mz*mz > 2.25f * sp * sp)
                continue;                // converged elsewhere
            if (pr.fdist > 0.05f * sp)
                continue;                // crease not on the surface
            /*  abs creases can't lift off the surface (field B IS
             *  the oracle), so the on-surface gate always passes -
             *  the activity question is whether the surface is
             *  actually creased here.  */
            if (tp.is_max == 2 && !tr.kinked(p0, pr))
                continue;
            /*  Duplicate-trace guard: a seed can converge onto a
             *  curve this pair already traced.  Identity by
             *  DISTANCE-TO-SEGMENT: p0 is post-correct, ON the
             *  crease, so a true duplicate sits within chord sag
             *  of the existing polyline (<= ~0.06 sp at the
             *  tightest traceable curvature; Newton noise 1e-3),
             *  while the closest real neighbour rings measured
             *  0.22 sp (the additive-joint autopsy).  0.1 sp
             *  splits the populations with 2x margin both ways.
             *  The old 0.5 sp point-to-VERTEX radius was a second
             *  ring-eater: same-clause loops in sub-half-cell
             *  gaps lost every seed (Nate's overlay read,
             *  2026-07-17: crowded bands under-polylined).  */
            {
                bool dup = false;
                const float dr2 = 0.1f * sp * 0.1f * sp;
                for (size_t c = 0; c < polys.size() && !dup; ++c)
                {
                    if (poly_pair[c] != tp.clause)
                        continue;
                    const std::vector<DSurfPoint>& pl = polys[c];
                    for (size_t qi = 0; qi + 1 < pl.size(); ++qi)
                    {
                        const DSurfPoint& A = pl[qi];
                        const DSurfPoint& B = pl[qi + 1];
                        const float bx = B.x - A.x,
                                    by = B.y - A.y,
                                    bz = B.z - A.z;
                        const float px = p0[0] - A.x,
                                    py = p0[1] - A.y,
                                    pz = p0[2] - A.z;
                        const float bb = bx*bx + by*by + bz*bz;
                        float t = bb > 0
                                ? (px*bx + py*by + pz*bz) / bb
                                : 0;
                        t = t < 0 ? 0 : t > 1 ? 1 : t;
                        const float dx = px - t*bx,
                                    dy = py - t*by,
                                    dz = pz - t*bz;
                        if (dx*dx + dy*dy + dz*dz < dr2)
                        {
                            dup = true;
                            break;
                        }
                    }
                }
                if (dup)
                {
                    consumed[s] = 1;
                    continue;
                }
            }
            float t0[3];
            CreaseTracer::tangent(pr, t0);
            if (!(t0[0] || t0[1] || t0[2]))
                continue;                // tangency at the seed

            std::vector<DSurfPoint> fwd, bwd;
            const bool loop = tr.march(p0, t0, +1, fwd, halt);
            if (!loop)
                tr.march(p0, t0, -1, bwd, halt);
            std::vector<DSurfPoint> poly;
            poly.reserve(bwd.size() + 1 + fwd.size());
            for (auto it = bwd.rbegin(); it != bwd.rend(); ++it)
                poly.push_back(*it);
            poly.push_back({ p0[0], p0[1], p0[2] });
            for (const DSurfPoint& q : fwd)
                poly.push_back(q);

            /*  Consume every feature near the traced curve so the
             *  next unconsumed seed is a NEW component.  Identity
             *  by DISTANCE-TO-SEGMENT at a tight radius (0.05 sp):
             *  own-curve seeds sit at Newton noise (~1e-3 sp) from
             *  the polyline SEGMENTS, while the closest real
             *  neighbour rings measured 0.22 sp away (the additive-
             *  joint autopsy).  The old point-to-vertex test needed
             *  a FULL-CELL radius to cover the chord spacing - and
             *  ate entire neighbouring rings in crowded bands
             *  (Nate's eyeball read: "averaged, merged, or
             *  skipped").  */
            const float cr2 = 0.05f * sp * 0.05f * sp;
            for (size_t s2 = 0; s2 < consumed.size(); ++s2)
            {
                if (consumed[s2])
                    continue;
                const DSurfPoint& F = seeds[s2];
                for (size_t qi = 0; qi + 1 < poly.size(); ++qi)
                {
                    const DSurfPoint& A = poly[qi];
                    const DSurfPoint& B = poly[qi + 1];
                    const float bx = B.x - A.x, by = B.y - A.y,
                                bz = B.z - A.z;
                    const float px = F.x - A.x, py = F.y - A.y,
                                pz = F.z - A.z;
                    const float bb = bx*bx + by*by + bz*bz;
                    float t = bb > 0
                            ? (px*bx + py*by + pz*bz) / bb : 0;
                    t = t < 0 ? 0 : t > 1 ? 1 : t;
                    const float dx = px - t*bx, dy = py - t*by,
                                dz = pz - t*bz;
                    if (dx*dx + dy*dy + dz*dz < cr2)
                    {
                        consumed[s2] = 1;
                        break;
                    }
                }
            }
            consumed[s] = 1;
            if (poly.size() >= 3)
            {
                polys.push_back(std::move(poly));
                closed.push_back(loop ? 1 : 0);
                poly_pair.push_back(tp.clause);
            }
        }
    }
    if (polys.empty())
    {
        tape_ctx_free(ctx);
        return false;
    }

    if (getenv("STIBIUM_DMESH_TRACE_DEBUG"))
        for (size_t c = 0; c < polys.size(); ++c)
        {
            float mind = 1e30f;
            for (size_t q = 1; q < polys[c].size(); ++q)
            {
                const float dx = polys[c][q].x - polys[c][q-1].x;
                const float dy = polys[c][q].y - polys[c][q-1].y;
                const float dz = polys[c][q].z - polys[c][q-1].z;
                mind = std::min(mind,
                        sqrtf(dx*dx + dy*dy + dz*dz));
            }
            float minx = 1e30f;
            for (size_t c2 = 0; c2 < polys.size(); ++c2)
            {
                if (c2 == c)
                    continue;
                for (const DSurfPoint& A : polys[c])
                    for (const DSurfPoint& B : polys[c2])
                    {
                        const float dx = A.x - B.x,
                                    dy = A.y - B.y,
                                    dz = A.z - B.z;
                        minx = std::min(minx, sqrtf(
                                dx*dx + dy*dy + dz*dz));
                    }
            }
            fprintf(stderr, "TRACE poly %zu (pair@%u): %zu pts, "
                    "closed %d, min-consec %.2e, min-cross %.2e, "
                    "start (%.4f,%.4f,%.4f)\n", c, poly_pair[c],
                    polys[c].size(), int(closed[c]), mind, minx,
                    polys[c][0].x, polys[c][0].y, polys[c][0].z);
        }

    /*  Junction sharing: open-chain endpoints within reach cluster
     *  to their mean, so incident creases carry the IDENTICAL
     *  corner coordinates and the triangulation's coincidence
     *  guard makes them one shared constrained vertex.  */
    {
        struct End { uint32_t poly; int end; };
        std::vector<End> ends;
        for (uint32_t c = 0; c < polys.size(); ++c)
            if (!closed[c])
            {
                ends.push_back({ c, 0 });
                ends.push_back({ c, 1 });
            }
        const auto pt_of = [&](const End& e) -> DSurfPoint& {
            return e.end ? polys[e.poly].back()
                         : polys[e.poly].front();
        };
        std::vector<uint8_t> used(ends.size(), 0);
        const float jr2 = 0.75f * 0.75f * sp * sp;
        for (size_t i = 0; i < ends.size(); ++i)
        {
            if (used[i])
                continue;
            std::vector<size_t> cl{ i };
            const DSurfPoint& A = pt_of(ends[i]);
            for (size_t j = i + 1; j < ends.size(); ++j)
            {
                if (used[j])
                    continue;
                const DSurfPoint& B = pt_of(ends[j]);
                const float dx = B.x - A.x, dy = B.y - A.y,
                            dz = B.z - A.z;
                if (dx*dx + dy*dy + dz*dz < jr2)
                    cl.push_back(j);
            }
            if (cl.size() < 2)
                continue;
            float mx = 0, my = 0, mz = 0;
            for (const size_t q : cl)
            {
                mx += pt_of(ends[q]).x;
                my += pt_of(ends[q]).y;
                mz += pt_of(ends[q]).z;
            }
            mx /= cl.size();
            my /= cl.size();
            mz /= cl.size();
            /*  The mean of near-corner endpoints sits off the
             *  surface by their spread (~2.5e-3 measured on cube
             *  corners); Newton it back down before it becomes a
             *  shared constrained vertex.  */
            for (int it = 0; it < 3; ++it)
            {
                const float h = 0.01f * sp;
                float xs[7], ys[7], zs[7];
                for (int q = 0; q < 7; ++q)
                {
                    xs[q] = mx;
                    ys[q] = my;
                    zs[q] = mz;
                }
                xs[1] += h;  xs[2] -= h;
                ys[3] += h;  ys[4] -= h;
                zs[5] += h;  zs[6] -= h;
                Region dummy = {};
                dummy.voxels = 7;
                dummy.X = xs;
                dummy.Y = ys;
                dummy.Z = zs;
                const float* rf = tape_eval_r(base, ctx, dummy);
                const float f = rf[0];
                const float h2 = 2 * h;
                const float gx = (rf[1] - rf[2]) / h2;
                const float gy = (rf[3] - rf[4]) / h2;
                const float gz = (rf[5] - rf[6]) / h2;
                const float g2 = gx*gx + gy*gy + gz*gz;
                if (!(g2 > 0) || !std::isfinite(f))
                    break;
                float step = -f / g2;
                const float cap = 0.5f * sp / sqrtf(g2);
                if (step > cap)
                    step = cap;
                if (step < -cap)
                    step = -cap;
                mx += gx * step;
                my += gy * step;
                mz += gz * step;
            }
            for (const size_t q : cl)
            {
                pt_of(ends[q]) = { mx, my, mz };
                used[q] = 1;
            }
        }
    }
    tape_ctx_free(ctx);

    /*  Append traced points as the constrained vertex set; the
     *  chains reference them by absolute surface index.  QEF
     *  features stay in the soup as plain surface points (the
     *  crease band drop retires the redundant ones).  */
    for (auto& poly : polys)
    {
        std::vector<uint32_t> chain;
        chain.reserve(poly.size());
        for (const DSurfPoint& q : poly)
        {
            chain.push_back(uint32_t(soup->surface.size()));
            soup->surface.push_back(q);
            ++soup->traced;
        }
        soup->tchains.push_back(std::move(chain));
    }
    soup->tclosed = closed;
    return true;
}

////////////////////////////////////////////////////////////////////////////
//  Stages B + C: Delaunay, refinement, extraction (CGAL)

#ifdef STIBIUM_HAS_CGAL

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Conforming_constrained_Delaunay_triangulation_3.h>
#include <CGAL/spatial_sort.h>
#include <CGAL/Spatial_sort_traits_adapter_3.h>
#include <CGAL/property_map.h>

namespace {

using K   = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb  = CGAL::Triangulation_vertex_base_with_info_3<int8_t, K>;
using Cb  = CGAL::Triangulation_cell_base_with_info_3<int8_t, K>;
using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
using DT  = CGAL::Delaunay_triangulation_3<K, Tds>;
using DPoint3 = DT::Point;

/*  The constrained path (MESH-NEXT, constrained-crease round):
 *  CGAL 6.2's conforming constrained Delaunay triangulation, with
 *  our info bases STACKED under its bases (both take their base
 *  class as the second template parameter and derive from it).
 *  The documented wrapper class only builds from a finished PLC;
 *  incremental insert() and insert_constrained_edge() live on the
 *  _impl class, which IS the triangulation (derives from
 *  Conforming_Delaunay_triangulation_3<T_3> deriving from T_3), so
 *  the whole DT code path ports onto it.  Header-declared public,
 *  version-pinned by the system CGAL; GPL like the rest of us.  */
using CVb  = CGAL::Conforming_constrained_Delaunay_triangulation_vertex_base_3<K, Vb>;
using CCb  = CGAL::Conforming_constrained_Delaunay_triangulation_cell_base_3<K, Cb>;
using CTds = CGAL::Triangulation_data_structure_3<CVb, CCb>;
/*  The impl's T_3 must be a Delaunay_triangulation_3 (it uses the
 *  conflict-region machinery), exactly as the wrapper builds its
 *  own: Delaunay over the CCDT-based Tds.  Delaunay under the hood
 *  also keeps nearest_vertex() available for the crowding guard.  */
using CTr  = CGAL::Delaunay_triangulation_3<K, CTds>;
using CCDT = CGAL::Conforming_constrained_Delaunay_triangulation_3_impl<CTr>;

/*  Vertex info: -1 inside, +1 outside, 0 surface.  */

std::vector<DSurfPoint> bisect_pending(const Deck* deck, TapeCtx* ctx,
                                       std::vector<PendingEdge>& edges,
                                       volatile int* halt)
{
    std::vector<DSurfPoint> out;
    const size_t n = edges.size();
    if (n == 0)
        return out;
    Tape* base = deck_base(deck);
    std::vector<float> xs(n), ys(n), zs(n), vals;
    for (int round = 0; round < BISECT_ROUNDS; ++round)
    {
        if (*halt)
            return out;
        for (size_t e = 0; e < n; ++e)
        {
            xs[e] = (edges[e].ax + edges[e].bx) * 0.5f;
            ys[e] = (edges[e].ay + edges[e].by) * 0.5f;
            zs[e] = (edges[e].az + edges[e].bz) * 0.5f;
        }
        eval_points_mt(deck, ctx, xs, ys, zs, vals);
        for (size_t e = 0; e < n; ++e)
        {
            PendingEdge& ed = edges[e];
            if (vals[e] < 0)
                { ed.ax = xs[e]; ed.ay = ys[e]; ed.az = zs[e]; }
            else
                { ed.bx = xs[e]; ed.by = ys[e]; ed.bz = zs[e]; }
        }
    }
    out.reserve(n);
    for (const PendingEdge& ed : edges)
        out.push_back({ (ed.ax + ed.bx) * 0.5f,
                        (ed.ay + ed.by) * 0.5f,
                        (ed.az + ed.bz) * 0.5f });
    return out;
}

/*  Stages B+C+repair, templated over the triangulation.  The plain
 *  DT path is the landed reference; CCDT_MODE additionally inserts
 *  the crease chains as constrained edges BEFORE refinement, making
 *  crease-crossing chords structurally impossible (the complete
 *  cure for sharp-point pinches; MESH-NEXT constrained-crease
 *  round).  The conforming machinery adds Steiner vertices ON the
 *  constrained segments whenever a constraint would otherwise be
 *  missing from the Delaunay structure; their info is indeterminate
 *  (with_info does not initialize), and their correct value is 0 -
 *  they sit on the crease, which IS the surface - so every insert
 *  batch is followed by a sweep marking them.  */
template <bool CCDT_MODE, class Tri>
bool mesh_impl(const Deck* deck, const DSoup& soup,
               volatile int* halt, DMesh* out,
               const std::vector<std::array<float, 3>>* quarantine
                       = nullptr)
{
    using TPoint = typename Tri::Point;
    PhaseTimer pt;
    prog_stage(4);
    using VH     = typename Tri::Vertex_handle;
    using CH     = typename Tri::Cell_handle;

    TapeCtx* ctx = tape_ctx_new(deck);
    Tri dt;

    /*  Vertex provenance, for the non-manifold debug dump
     *  (STIBIUM_DMESH_NM_DEBUG): 1 sample, 2 bisected surface,
     *  3 feature, 4 refinement insert, 5 repair insert.  Vertices
     *  the conforming machinery created carry no tag and identify
     *  themselves via ccdt_3_data().  */
    std::unordered_map<const void*, uint8_t> prov;

    const auto sweep_steiner = [&dt]() {
        if constexpr (CCDT_MODE)
            for (auto v = dt.finite_vertices_begin();
                 v != dt.finite_vertices_end(); ++v)
                if (v->ccdt_3_data().is_Steiner_vertex_on_edge())
                    v->info() = 0;
    };

    /*  Constrained-crease round, step 0 (before ANY insertion):
     *  extract the chains and let the oracle referee each segment -
     *  a chain segment is only constrained if it actually lies on
     *  the crease (|f|/|grad| at 1/4, 1/2, 3/4 within the wart
     *  tolerance).  The linker can shortcut where two creases cross
     *  (csg's sharp point: union seam meets the plane circle), and
     *  a wrong constraint is worse than any wart - it is FORCED
     *  into the triangulation.  Accepted segments also define the
     *  crease band used to drop redundant surface vertices and to
     *  keep repair out.  */
    std::vector<std::pair<uint32_t, uint32_t>> accepted;
    std::vector<std::array<float, 6>> cseg;
    std::vector<uint32_t> cseg_chain;   // owning chain per segment
    if constexpr (CCDT_MODE)
    {
        /*  Traced chains (exact, ordered, junction-sharing) when
         *  the tracer ran; the QEF chain extractor otherwise.  */
        DChains chains;
        if (!soup.tchains.empty())
        {
            chains.chains = soup.tchains;
            chains.closed.assign(soup.tclosed.begin(),
                                 soup.tclosed.end());
        }
        else
            chains = delaunay_chains(soup);
        std::vector<std::pair<uint32_t, uint32_t>> cand;
        std::vector<uint32_t> cand_chain;
        for (size_t c = 0; c < chains.chains.size(); ++c)
        {
            const auto& chain = chains.chains[c];
            for (size_t j = 0; j + 1 < chain.size(); ++j)
            {
                cand.push_back({ chain[j], chain[j + 1] });
                cand_chain.push_back(uint32_t(c));
            }
            if (chains.closed[c] && chain.size() > 2)
            {
                cand.push_back({ chain.back(), chain.front() });
                cand_chain.push_back(uint32_t(c));
            }
        }
        Tape* base = deck_base(deck);
        constexpr int NS = 3;
        std::vector<float> sx, sy, sz, sv, slen, shs;
        sx.reserve(cand.size() * NS * 7);
        for (const auto& [a, b] : cand)
        {
            const DSurfPoint& A = soup.surface[a];
            const DSurfPoint& B = soup.surface[b];
            const float ex = B.x - A.x, ey = B.y - A.y,
                        ez = B.z - A.z;
            const float len = sqrtf(ex*ex + ey*ey + ez*ez);
            slen.push_back(len);
            const float h = len > 0 ? len * 0.01f : 1e-5f;
            shs.push_back(h);
            for (int u = 0; u < NS; ++u)
            {
                const float t = float(u + 1) / (NS + 1);
                const float mx = A.x + t * ex;
                const float my = A.y + t * ey;
                const float mz = A.z + t * ez;
                /*  7-tap: value + central differences  */
                for (int q = 0; q < 7; ++q)
                {
                    sx.push_back(mx);
                    sy.push_back(my);
                    sz.push_back(mz);
                }
                const size_t o = sx.size() - 7;
                sx[o+1] += h;  sx[o+2] -= h;
                sy[o+3] += h;  sy[o+4] -= h;
                sz[o+5] += h;  sz[o+6] -= h;
            }
        }
        if (!sx.empty())
            eval_points_mt(deck, ctx, sx, sy, sz, sv);
        /*  Cross-pair duplicate-coverage dedupe (zeiss run,
         *  2026-07-15): a merged CSG tree repeats surfaces, so the
         *  SAME physical crease is traced by several min/max pairs
         *  with different vertex spacing.  Overlapping collinear
         *  constraints throw the conforming machinery into an
         *  encroachment cascade - Steiner points ping-pong toward
         *  the mutual corner until femtometer segments trip the
         *  through-vertex assertion (measured: Steiner triplet
         *  7e-15 apart at (-5, 3.75, 26); the coordinates are
         *  doubles, so they can only be CGAL constructions, not
         *  soup points).  A segment whose endpoints AND midpoint
         *  all lie on already-accepted coverage adds nothing: drop
         *  it.  0.2 sp keeps chain continuations (their far
         *  endpoint sticks 0.5 sp out) and junction spokes (they
         *  leave the corner perpendicular).  */
        uint64_t dup_dropped = 0, cross_dropped = 0;
        uint64_t quar_dropped = 0;
        /*  Interior-interior closest approach between the candidate
         *  and an accepted segment.  The three-point coverage probe
         *  is blind to a mid-segment crossing (endpoints and
         *  midpoint all sit far from the other segment; only the
         *  crossing point is near), and CCDT_3's contract forbids
         *  constraints meeting anywhere but shared endpoints - a
         *  shallow-angle crossing Steiner-cascades to femtometer
         *  segments and the through-vertex assertion.  Closest
         *  approach with BOTH parameters interior and distance
         *  under 0.05 sp = crossing or partial collinear overlap;
         *  shared junction corners approach at t ~ 0/1 and are
         *  exempt.  */
        /*  Round two of the crossing rules (zeiss again): the
         *  t-based corner exemption excused a REAL crossing at
         *  t = 0.046 ("Steiner point coincides with an existing
         *  vertex ... two input segments are too close").  What
         *  legalizes corner contact is SHARING a vertex, not
         *  parameter proximity to any endpoint:
         *   - segments sharing one endpoint (welded junction
         *     spokes) are judged by their FAR parts - midpoint or
         *     far endpoint hugging the other segment = shallow
         *     spoke, dropped;
         *   - segments sharing nothing get full closest-approach
         *     with NO exemption, except endpoint-to-interior
         *     contacts within the pre-split tolerance (5e-3 sp),
         *     which insertion legalizes into a shared vertex
         *     (T-junctions).  */
        const auto d2_point_seg = [](float px, float py, float pz,
                                     const std::array<float, 6>& s) {
            const float ax = s[0], ay = s[1], az = s[2];
            const float bx = s[3] - ax, by = s[4] - ay,
                        bz = s[5] - az;
            const float qx = px - ax, qy = py - ay, qz = pz - az;
            const float bb = bx*bx + by*by + bz*bz;
            float t = bb > 0 ? (qx*bx + qy*by + qz*bz) / bb : 0;
            t = t < 0 ? 0 : t > 1 ? 1 : t;
            const float dx = qx - t*bx, dy = qy - t*by,
                        dz = qz - t*bz;
            return dx*dx + dy*dy + dz*dz;
        };
        const auto crosses = [&](float ax, float ay, float az,
                                 float bx, float by, float bz,
                                 float r) {
            const float weld = 1e-3f * soup.spacing;
            const float weld2 = weld * weld;
            const float tj = 5e-3f * soup.spacing;
            const float tj2 = tj * tj;
            const float ux = bx - ax, uy = by - ay, uz = bz - az;
            const float uu = ux*ux + uy*uy + uz*uz;
            const auto d2pp = [](float x1, float y1, float z1,
                                 float x2, float y2, float z2) {
                const float dx = x1-x2, dy = y1-y2, dz = z1-z2;
                return dx*dx + dy*dy + dz*dz;
            };
            for (const auto& s : cseg)
            {
                const bool a0 = d2pp(ax, ay, az,
                                     s[0], s[1], s[2]) < weld2;
                const bool a1 = d2pp(ax, ay, az,
                                     s[3], s[4], s[5]) < weld2;
                const bool b0 = d2pp(bx, by, bz,
                                     s[0], s[1], s[2]) < weld2;
                const bool b1 = d2pp(bx, by, bz,
                                     s[3], s[4], s[5]) < weld2;
                if ((a0 || a1) && (b0 || b1))
                    continue;      // same span: dup check's turf
                if (a0 || a1 || b0 || b1)
                {
                    /*  One shared corner: judge the far parts  */
                    const float mx = (ax + bx) * 0.5f,
                                my = (ay + by) * 0.5f,
                                mz = (az + bz) * 0.5f;
                    const float fx = (a0 || a1) ? bx : ax,
                                fy = (a0 || a1) ? by : ay,
                                fz = (a0 || a1) ? bz : az;
                    if (d2_point_seg(mx, my, mz, s) < r * r ||
                        d2_point_seg(fx, fy, fz, s) < r * r)
                        return true;
                    continue;
                }
                /*  No sharing: full closest approach  */
                const float cx = s[0], cy = s[1], cz = s[2];
                const float vx = s[3] - cx, vy = s[4] - cy,
                            vz = s[5] - cz;
                const float vv = vx*vx + vy*vy + vz*vz;
                const float wx = ax - cx, wy = ay - cy,
                            wz = az - cz;
                const float uv = ux*vx + uy*vy + uz*vz;
                const float uw = ux*wx + uy*wy + uz*wz;
                const float vw = vx*wx + vy*wy + vz*wz;
                const float den = uu*vv - uv*uv;
                float t1, t2;
                if (den > 1e-12f * uu * vv)
                {
                    t1 = (uv*vw - vv*uw) / den;
                    t2 = (uu*vw - uv*uw) / den;
                }
                else
                {
                    t1 = 0.5f;
                    t2 = vv > 0 ? (vw + 0.5f*uv) / vv : 0;
                }
                t1 = t1 < 0 ? 0 : t1 > 1 ? 1 : t1;
                t2 = t2 < 0 ? 0 : t2 > 1 ? 1 : t2;
                const float px = ax + t1*ux - (cx + t2*vx);
                const float py = ay + t1*uy - (cy + t2*vy);
                const float pz = az + t1*uz - (cz + t2*vz);
                const float d2m = px*px + py*py + pz*pz;
                if (d2m >= r * r)
                    continue;
                /*  T-contact the pre-split legalizes: an ENDPOINT
                 *  of either segment landing on the other's
                 *  interior within its tolerance  */
                if (d2m < tj2 &&
                    (t1 <= 0.02f || t1 >= 0.98f ||
                     t2 <= 0.02f || t2 >= 0.98f))
                    continue;
                return true;
            }
            return false;
        };
        const auto covered = [&](float x, float y, float z,
                                 float r) {
            for (const auto& s : cseg)
            {
                const float ax = s[0], ay = s[1], az = s[2];
                const float bx = s[3] - ax, by = s[4] - ay,
                            bz = s[5] - az;
                const float px = x - ax, py = y - ay, pz = z - az;
                const float bb = bx*bx + by*by + bz*bz;
                float t = bb > 0
                        ? (px*bx + py*by + pz*bz) / bb : 0;
                t = t < 0 ? 0 : t > 1 ? 1 : t;
                const float dx = px - t*bx, dy = py - t*by,
                            dz = pz - t*bz;
                if (dx*dx + dy*dy + dz*dz < r * r)
                    return true;
            }
            return false;
        };
        for (size_t i = 0; i < cand.size(); ++i)
        {
            /*  Degenerate segments carry no law (belt to the
             *  march guard's braces): a float-noise-length
             *  constraint can only feed the conforming machinery
             *  pathologies.  */
            bool good = slen[i] > 1e-5f * soup.spacing;
            int judged_samples = 0;
            int kink_samples = 0;
            for (int u = 0; good && u < NS; ++u)
            {
                const size_t o = (i * NS + u) * 7;
                const float f = sv[o];
                const float h2 = 2 * shs[i];
                const float gx = (sv[o+1] - sv[o+2]) / h2;
                const float gy = (sv[o+3] - sv[o+4]) / h2;
                const float gz = (sv[o+5] - sv[o+6]) / h2;
                const float gl = sqrtf(gx*gx + gy*gy + gz*gz);
                /*  Relative tolerance with an absolute floor: QEF
                 *  endpoints carry ~1e-3 placement noise (measured
                 *  |f| <= 1.3e-3), so a purely relative bar starves
                 *  SHORT segments - csg's plane circle lost all 96
                 *  of its densest (0.5-0.9 cell) segments at
                 *  3.0-4.1% while true shortcuts sit far beyond
                 *  either bound.  */
                /*  "Couldn't judge" is NOT "guilty" (bino model,
                 *  2026-07-17: 198 of 292 rejections carried 0.0%
                 *  measured error and tipped the trust gate,
                 *  torching 2,000+ good constraints).  A non-finite
                 *  sample gives no verdict; a zero-gradient sample
                 *  with f ~ 0 sits exactly ON the kink (the central
                 *  difference cancels across a crease - which is
                 *  where these polylines live).  Only samples with
                 *  a real distance reading can convict.  */
                if (!std::isfinite(f))
                    continue;
                if (!(gl > 0))
                {
                    /*  Spacing-relative with the old absolute as
                     *  floor (correctness review #5: a bare 1e-6
                     *  is scale-dependent - 100x model, 100x
                     *  stricter).  The on-kink population reads
                     *  <= ~1e-3, so the zone is wide.  */
                    if (fabsf(f) > std::max(1e-6f,
                                1e-5f * soup.spacing))
                        good = false;
                    else
                        ++kink_samples;
                    continue;
                }
                ++judged_samples;
                if (fabsf(f) / gl > std::max(slen[i] * 0.03f,
                            soup.spacing * 0.05f))
                    good = false;
            }
            /*  A segment with NO judgeable sample could hide a
             *  shortcut through a NaN void - demand at least one
             *  real reading.  A segment whose EVERY sample sits
             *  exactly on the kink (finite, f ~ 0, cancelled
             *  gradient) is the opposite of a void: a crease
             *  lying IN a grid-aligned face reads this way along
             *  its whole length, and convicting it stripped the
             *  plinth junction of all law - ten rounds of
             *  air-chords at a circle the tracer had chained
             *  perfectly (2026-07-18, every reject 0.0% len).  */
            if (judged_samples == 0 && kink_samples == 0)
                good = false;
            if (good)
            {
                const DSurfPoint& A = soup.surface[cand[i].first];
                const DSurfPoint& B = soup.surface[cand[i].second];
                /*  0.2 sp grazed REAL neighbour rings (they start
                 *  at 0.22 sp in the additive joints); true
                 *  re-traces of the same crease agree to Newton
                 *  noise, so 0.05 sp still catches every genuine
                 *  duplicate.  */
                const float rr = 0.05f * soup.spacing;
                if (covered(A.x, A.y, A.z, rr) &&
                    covered(B.x, B.y, B.z, rr) &&
                    covered((A.x + B.x) * 0.5f,
                            (A.y + B.y) * 0.5f,
                            (A.z + B.z) * 0.5f, rr))
                {
                    ++dup_dropped;
                    continue;
                }
                if (crosses(A.x, A.y, A.z, B.x, B.y, B.z,
                            0.05f * soup.spacing))
                {
                    ++cross_dropped;
                    continue;
                }
                /*  Quarantined cascade sites (see the dispatcher's
                 *  retry loop): a conforming-cascade corner takes
                 *  its local constraints out of play so the other
                 *  99% survive.  */
                if (quarantine)
                {
                    const float qr = soup.spacing;
                    bool quar = false;
                    for (const auto& Q : *quarantine)
                    {
                        const float bx = B.x - A.x, by = B.y - A.y,
                                    bz = B.z - A.z;
                        const float px = Q[0] - A.x,
                                    py = Q[1] - A.y,
                                    pz = Q[2] - A.z;
                        const float bb = bx*bx + by*by + bz*bz;
                        float t = bb > 0
                                ? (px*bx + py*by + pz*bz) / bb : 0;
                        t = t < 0 ? 0 : t > 1 ? 1 : t;
                        const float dx = px - t*bx, dy = py - t*by,
                                    dz = pz - t*bz;
                        if (dx*dx + dy*dy + dz*dz < qr * qr)
                        {
                            quar = true;
                            break;
                        }
                    }
                    if (quar)
                    {
                        ++quar_dropped;
                        continue;
                    }
                }
                accepted.push_back(cand[i]);
                cseg.push_back({ A.x, A.y, A.z, B.x, B.y, B.z });
                cseg_chain.push_back(cand_chain[i]);
            }
            else if (getenv("STIBIUM_DMESH_SEG_DEBUG"))
            {
                const DSurfPoint& A = soup.surface[cand[i].first];
                const DSurfPoint& B = soup.surface[cand[i].second];
                float worst = 0;
                for (int u = 0; u < NS; ++u)
                {
                    const size_t o = (i * NS + u) * 7;
                    const float h2 = 2 * shs[i];
                    const float gx = (sv[o+1] - sv[o+2]) / h2;
                    const float gy = (sv[o+3] - sv[o+4]) / h2;
                    const float gz = (sv[o+5] - sv[o+6]) / h2;
                    const float gl =
                            sqrtf(gx*gx + gy*gy + gz*gz);
                    if (gl > 0 && std::isfinite(sv[o]))
                        worst = std::max(worst,
                                fabsf(sv[o]) / gl);
                }
                fprintf(stderr, "SEG reject (%.4f,%.4f,%.4f)-"
                        "(%.4f,%.4f,%.4f) len %.4f (%.2f sp) "
                        "worst %.5f (%.1f%% len)\n",
                        A.x, A.y, A.z, B.x, B.y, B.z, slen[i],
                        slen[i] / soup.spacing, worst,
                        100.f * worst / slen[i]);
            }
        }
        /*  The gate: constraints become law only when the chains
         *  earn it.  A high rejection rate means the extractor
         *  mislinked (csg: two creases crossing at the sharp point
         *  - 35% of its segments are shortcuts), and a half-trusted
         *  polyline plus the band drop built around it does more
         *  damage than no polyline (measured: 7 pinches vs 3).
         *  Fall back to unconstrained behavior for the model;
         *  junction-aware chain extraction is the designed next
         *  round.  */
        /*  Duplicate-coverage drops are not rejections - the crease
         *  is constrained, just by its first copy - so they leave
         *  the gate arithmetic entirely.  */
        /*  The gate is FALLBACK-ONLY (zeiss shallow round,
         *  2026-07-18).  It exists because the radius-graph
         *  extractor can MISLINK - a high rejection rate means
         *  even its accepted segments connect the wrong points
         *  (csg's 35%-shortcut era; half-trusted polylines
         *  measured 7 pinches vs 3).  Traced chains cannot
         *  mislink (ordered by the marcher), every accepted
         *  segment passed the oracle individually, and their
         *  rejects are local chord-sag over tight curvature -
         *  already dropped one by one.  Collective punishment on
         *  traced chains was measured harmful BOTH ways: the
         *  model-global nuke cost zeiss all ~15K constraints over
         *  1,784 local rejects (4 open, 965 nm), and per-chain
         *  conviction cost the bino 163 verified segments
         *  (14 open vs 0).  Judgment stays individual.  */
        const size_t judged = cand.size() -
                size_t(dup_dropped + cross_dropped + quar_dropped);
        if (soup.tchains.empty() && judged > 0 &&
            float(judged - accepted.size()) > 0.10f * float(judged))
        {
            accepted.clear();
            cseg.clear();
            cseg_chain.clear();
        }
        if ((dup_dropped || cross_dropped || quar_dropped) &&
            getenv("STIBIUM_DMESH_SEG_DEBUG"))
            fprintf(stderr, "constraints: %llu duplicate-coverage, "
                    "%llu crossing, %llu quarantined segments "
                    "dropped\n",
                    (unsigned long long)dup_dropped,
                    (unsigned long long)cross_dropped,
                    (unsigned long long)quar_dropped);
        /*  Targeted forensics: STIBIUM_DMESH_PROBE="x,y,z[,r]"
         *  dumps everything the constraint pipeline knows within
         *  r cells (default 1) of the point: accepted segments,
         *  QEF feature points (the tracer's seeds), and traced
         *  polyline points.  Segments-without-features = seedless
         *  crease; features-without-polylines = seeds refused.  */
        if (const char* pe = getenv("STIBIUM_DMESH_PROBE"))
        {
            float px2, py2, pz2, prr = 1;
            const int nsc = sscanf(pe, "%f,%f,%f,%f",
                                   &px2, &py2, &pz2, &prr);
            if (nsc >= 3)
            {
                const float pr2 = prr * prr *
                        soup.spacing * soup.spacing;
                for (const auto& s : cseg)
                {
                    float best = 1e30f;
                    for (int q = 0; q < 2; ++q)
                    {
                        const float dx = s[3*q] - px2,
                                    dy = s[3*q+1] - py2,
                                    dz = s[3*q+2] - pz2;
                        best = std::min(best,
                                dx*dx + dy*dy + dz*dz);
                    }
                    if (best < pr2)
                        fprintf(stderr, "PROBE seg (%.6f,%.6f,%.6f)"
                                "-(%.6f,%.6f,%.6f)\n",
                                s[0], s[1], s[2],
                                s[3], s[4], s[5]);
                }
                const size_t pf0 = soup.surface.size() -
                        size_t(soup.feature_points);
                uint64_t nfeat = 0;
                for (size_t i = pf0; i < soup.surface.size(); ++i)
                {
                    const DSurfPoint& p = soup.surface[i];
                    const float dx = p.x - px2, dy = p.y - py2,
                                dz = p.z - pz2;
                    if (dx*dx + dy*dy + dz*dz < pr2)
                    {
                        ++nfeat;
                        fprintf(stderr, "PROBE feat "
                                "(%.6f,%.6f,%.6f)\n",
                                p.x, p.y, p.z);
                    }
                }
                uint64_t ntr = 0;
                for (const auto& ch : soup.tchains)
                    for (const uint32_t idx : ch)
                    {
                        const DSurfPoint& p = soup.surface[idx];
                        const float dx = p.x - px2,
                                    dy = p.y - py2,
                                    dz = p.z - pz2;
                        if (dx*dx + dy*dy + dz*dz < pr2)
                            ++ntr;
                    }
                fprintf(stderr, "PROBE summary: %llu features, "
                        "%llu traced points in radius\n",
                        (unsigned long long)nfeat,
                        (unsigned long long)ntr);
            }
        }
    }
    pt.mark("segment referee");
    prog_stage(5);
    const auto near_crease = [&](float x, float y, float z,
                                 float r) {
        for (const auto& s : cseg)
        {
            const float ax = s[0], ay = s[1], az = s[2];
            const float bx = s[3] - ax, by = s[4] - ay,
                        bz = s[5] - az;
            const float px = x - ax, py = y - ay, pz = z - az;
            const float bb = bx*bx + by*by + bz*bz;
            float t = bb > 0 ? (px*bx + py*by + pz*bz) / bb : 0;
            t = t < 0 ? 0 : t > 1 ? 1 : t;
            const float dx = px - t*bx, dy = py - t*by,
                        dz = pz - t*bz;
            if (dx*dx + dy*dy + dz*dz < r * r)
                return true;
        }
        return false;
    };
    /*  Redundant-vertex drop radius: a surface vertex this close to
     *  the constrained crease pairs with a chain/Steiner vertex
     *  into sliver tets - the measured pinch pairs sat 0.15-0.25
     *  cells apart, all on the crease.  The crease polyline IS the
     *  surface there; the vertex adds nothing.  Sign witnesses
     *  (inside/outside samples) are never dropped.  */
    /*  Scaled by the dense round: the corridor exists along the
     *  constrained creases, which is exactly where the dense lattice
     *  fired - sliver pairing distance tracks the local pitch.  */
    /*  Auto-density makes the factor LOCAL: outside the drill-down
     *  bands the corridor stays at base width; inside, it tracks
     *  the band's pitch.  Quick-reject with the widest possible
     *  radius before paying for the box lookup.  */
    const float drop_max = 0.35f * soup.spacing / dense_factor();
    /*  Two-chain witness rule (the eyepiece-ring autopsy,
     *  2026-07-17): a vertex is redundant only when ONE crease
     *  owns it.  Between two constrained rings closer than two
     *  corridor widths the old test deleted every witness of the
     *  surface BETWEEN them, and extraction spanned the strip
     *  with chords through air (the double ring of thorns, fuzzed
     *  full-circle only when constraints are on - bino A/B).
     *  Junction corners still drop: there the two chains' closest
     *  points coincide at the shared corner, so the foot points
     *  sit within the radius of each other.  */
    const auto sole_owner = [&](float x, float y, float z,
                                float r) {
        float d1 = 1e30f, d2 = 1e30f;
        float f1x = 0, f1y = 0, f1z = 0;
        float f2x = 0, f2y = 0, f2z = 0;
        uint32_t c1 = UINT32_MAX;
        for (size_t si = 0; si < cseg.size(); ++si)
        {
            const auto& s = cseg[si];
            const float ax = s[0], ay = s[1], az = s[2];
            const float bx = s[3] - ax, by = s[4] - ay,
                        bz = s[5] - az;
            const float px = x - ax, py = y - ay, pz = z - az;
            const float bb = bx*bx + by*by + bz*bz;
            float t = bb > 0 ? (px*bx + py*by + pz*bz) / bb : 0;
            t = t < 0 ? 0 : t > 1 ? 1 : t;
            const float fx = ax + t*bx, fy = ay + t*by,
                        fz = az + t*bz;
            const float dx = x - fx, dy = y - fy, dz = z - fz;
            const float d = dx*dx + dy*dy + dz*dz;
            const uint32_t ch = cseg_chain[si];
            if (ch == c1)
            {
                if (d < d1)
                {
                    d1 = d;
                    f1x = fx; f1y = fy; f1z = fz;
                }
            }
            else if (d < d1)
            {
                d2 = d1;
                f2x = f1x; f2y = f1y; f2z = f1z;
                d1 = d;
                f1x = fx; f1y = fy; f1z = fz;
                c1 = ch;
            }
            else if (d < d2)
            {
                d2 = d;
                f2x = fx; f2y = fy; f2z = fz;
            }
        }
        if (!(d1 < r * r))
            return false;             // near no crease: keep
        if (!(d2 < r * r))
            return true;              // one owner: redundant, drop
        const float gx = f2x - f1x, gy = f2y - f1y, gz = f2z - f1z;
        /*  Feet together = junction corner (drop as before);
         *  feet apart = the point sits BETWEEN two creases and is
         *  the strip's only witness.  */
        return gx*gx + gy*gy + gz*gz < r * r;
    };
    const auto in_corridor = [&](float x, float y, float z) {
        if (!near_crease(x, y, z, drop_max))
            return false;
        const float eff = 0.35f * soup.spacing /
                std::max(dense_factor(),
                         box_dense_factor(soup, x, y, z));
        return sole_owner(x, y, z, std::min(eff, drop_max));
    };
    /*  Feature points referenced by an accepted segment are the
     *  chain itself and always enter; the rest of the QEF clump
     *  (curved creases run 26-51% duplicates, merged out of the
     *  chains) is redundant by the merge's own judgment, and every
     *  duplicate near a constrained segment ENCROACHES on it,
     *  forcing the machinery to place a Steiner at its projection
     *  - a guaranteed near-coincidence (the measured spheres
     *  pinch: feature vs Steiner, 0.15 cells).  */
    std::unordered_set<uint32_t> chain_used;
    for (const auto& [a, b] : accepted)
    {
        chain_used.insert(a);
        chain_used.insert(b);
    }

    /*  Stage B: triangulate everything (spatial-sort batch insert),
     *  then refine: any tet edge joining an inside vertex directly
     *  to an outside vertex gets a bisected surface point, until no
     *  such edge remains.  At convergence, inside and outside
     *  vertices are separated everywhere by surface vertices.  */
    std::vector<std::pair<TPoint, int8_t>> pts;
    pts.reserve(soup.samples.size());
    for (const DSample& s : soup.samples)
    {
        if (s.on_surface && !cseg.empty() &&
            in_corridor(s.x, s.y, s.z))
            continue;
        pts.push_back({ TPoint(s.x, s.y, s.z),
                        int8_t(s.on_surface ? 0
                               : s.inside ? -1 : 1) });
    }
    if constexpr (CCDT_MODE)
    {
        /*  No bulk constructor here: spatial-sort ourselves and
         *  insert with the previous cell as hint, mirroring the
         *  class's own insert_vertices_range (including the CORNER
         *  tag, which marks input vertices for the constraint
         *  bookkeeping).  */
        using PMap = CGAL::First_of_pair_property_map<
                std::pair<TPoint, int8_t>>;
        CGAL::spatial_sort(pts.begin(), pts.end(),
                CGAL::Spatial_sort_traits_adapter_3<K, PMap>());
        CH hint{};
        size_t wpn = 0;
        for (const auto& pr : pts)
        {
            if ((++wpn & 4095) == 0)
                prog_frac(float(wpn) / float(pts.size()));
            const auto before = dt.number_of_vertices();
            VH vh = dt.insert(pr.first, hint, false);
            hint = vh->cell();
            vh->ccdt_3_data().set_vertex_type(
                    CGAL::CDT_3_vertex_type::CORNER);
            if (dt.number_of_vertices() > before)
            {
                vh->info() = pr.second;
                prov.emplace(&*vh, 1);
            }
        }
    }
    else
        dt.insert(pts.begin(), pts.end());

    pt.mark("insert samples");
    prog_stage(6);
    /*  Surface points go in one by one with a coincidence guard: on
     *  grid-aligned geometry a bisected point can converge exactly
     *  onto a lattice sample, and overwriting that vertex's inside/
     *  outside witness with 'surface' corrupts extraction.  Handles
     *  are recorded because chain entries (indices into
     *  soup.surface) become constrained-edge endpoints.  */
    uint64_t witness_overwrites = 0;
    std::vector<VH> surf_vh(soup.surface.size(), VH());
    const size_t feat_base =
            soup.surface.size() - size_t(soup.feature_points);
    /*  Insert in SPATIAL order with a walk hint (perf round 2,
     *  2026-07-18: this loop was 51 s of the r2 bino wall - every
     *  insert AND its nearest_vertex guard did a cold full-tree
     *  point location; the witness path two screens up has used
     *  spatial_sort + hints since day one).  surf_vh stays in
     *  soup order via the index map - chain constraints depend on
     *  it.  The corridor drop touches no DT state and decides
     *  first.  */
    std::vector<std::pair<TPoint, size_t>> sorder;
    sorder.reserve(soup.surface.size());
    /*  Corridor-drop autopsy (STIBIUM_DMESH_DROP_PROBE=
     *  "x,y,z,r"): count survivors vs drops near a site - the
     *  witness-starvation instrument (dense-details class).  */
    static const char* dp_env = getenv("STIBIUM_DMESH_DROP_PROBE");
    float dpx = 0, dpy = 0, dpz = 0, dpr = 0;
    size_t dp_kept = 0, dp_drop = 0;
    if (dp_env)
        sscanf(dp_env, "%f,%f,%f,%f", &dpx, &dpy, &dpz, &dpr);
    for (size_t si = 0; si < soup.surface.size(); ++si)
    {
        const DSurfPoint& s = soup.surface[si];
        /*  Surface points inside the crease band are redundant
         *  with the constrained polyline and pair into pinch
         *  slivers - drop them, EXCEPT the chain members
         *  themselves.  */
        const bool probed = dp_env &&
                fabsf(s.x - dpx) < dpr &&
                fabsf(s.y - dpy) < dpr &&
                fabsf(s.z - dpz) < dpr;
        if (!cseg.empty() && !chain_used.count(uint32_t(si)) &&
            in_corridor(s.x, s.y, s.z))
        {
            if (probed)
                ++dp_drop;
            continue;
        }
        if (probed)
        {
            ++dp_kept;
            fprintf(stderr, "DPKEEP %.5g %.5g %.5g\n",
                    s.x, s.y, s.z);
        }
        sorder.push_back({ TPoint(s.x, s.y, s.z), si });
    }
    if (dp_env)
        fprintf(stderr, "DROP_PROBE (%.3f,%.3f,%.3f r%.2f): "
                "%zu kept, %zu corridor-dropped\n",
                dpx, dpy, dpz, dpr, dp_kept, dp_drop);
    {
        using PMap = CGAL::First_of_pair_property_map<
                std::pair<TPoint, size_t>>;
        CGAL::spatial_sort(sorder.begin(), sorder.end(),
                CGAL::Spatial_sort_traits_adapter_3<K, PMap>());
    }
    CH shint{};
    size_t spn = 0;
    /*  Hash-grid coincidence guard (perf round 3 rock 2a): the
     *  1e-3 sp near-coincidence check paid a full nearest_vertex
     *  walk per point.  A spatial hash of every vertex placed so
     *  far answers the same question in O(1); the DT walk is
     *  reserved for actual inserts.  Nearest-of-candidates keeps
     *  the semantics of nearest_vertex exactly.  */
    const double hw = 1e-3 * soup.spacing;
    std::unordered_map<uint64_t, std::vector<VH>> vhash;
    const auto hkey = [&](double x, double y, double z) {
        const int64_t cx = int64_t(std::floor(x / hw));
        const int64_t cy = int64_t(std::floor(y / hw));
        const int64_t cz = int64_t(std::floor(z / hw));
        return (uint64_t(cx & 0x1FFFFF) << 42) |
               (uint64_t(cy & 0x1FFFFF) << 21) |
                uint64_t(cz & 0x1FFFFF);
    };
    if constexpr (CCDT_MODE)
    {
        vhash.reserve(dt.number_of_vertices());
        for (auto vit = dt.finite_vertices_begin();
             vit != dt.finite_vertices_end(); ++vit)
            vhash[hkey(vit->point().x(), vit->point().y(),
                       vit->point().z())]
                    .push_back(vit);
    }
    for (const auto& [sp3, si] : sorder)
    {
        if ((++spn & 1023) == 0)
            prog_frac(float(spn) / float(sorder.size()));
        /*  Near-coincidence weld (zeiss run, 2026-07-15): traced
         *  corners from different chains Newton-converge onto the
         *  same machined corner within double noise (measured
         *  7e-15 apart at (-5, 3.75, 26)).  Two distinct DT
         *  vertices that close produce femtometer constraint
         *  segments (whose length UNDERFLOWS float) and
         *  constraints passing through a near-coincident third
         *  copy - the CGAL through-vertex assertion.  Within
         *  1e-3 sp everything is the same point: reuse the
         *  existing vertex, and the micro-segment collapses into
         *  the va == vb guard below.  */
        if constexpr (CCDT_MODE)
        {
            {
                VH nv{};
                double best = hw * hw;
                const double px2 = sp3.x(), py2 = sp3.y(),
                             pz2 = sp3.z();
                for (int dx = -1; dx <= 1; ++dx)
                 for (int dy = -1; dy <= 1; ++dy)
                  for (int dz = -1; dz <= 1; ++dz)
                  {
                    const auto it = vhash.find(hkey(
                            px2 + dx*hw, py2 + dy*hw,
                            pz2 + dz*hw));
                    if (it == vhash.end())
                        continue;
                    for (const VH& cand2 : it->second)
                    {
                        const double d2 = CGAL::squared_distance(
                                cand2->point(), sp3);
                        if (d2 < best)
                        {
                            best = d2;
                            nv = cand2;
                        }
                    }
                  }
                if (nv != VH())
                {
                    {
                        /*  Correctness review #3: this reuse path
                         *  overwrites an existing vertex's info
                         *  unconditionally - if nv is a +-1 sign
                         *  witness its vote dies silently.  The
                         *  siblings only stamp NEW vertices.
                         *  INSTRUMENT (fix waits on data): count
                         *  the actual witness kills.  */
                        if (nv->info() != 0)
                            ++witness_overwrites;
                        nv->info() = 0;
                        surf_vh[si] = nv;
                        shint = nv->cell();
                        continue;
                    }
                }
            }
        }
        const auto before = dt.number_of_vertices();
        VH vh;
        if constexpr (CCDT_MODE)
        {
            vh = dt.insert(sp3, shint, false);
            vh->ccdt_3_data().set_vertex_type(
                    CGAL::CDT_3_vertex_type::CORNER);
        }
        else
            vh = dt.insert(sp3);
        if (dt.number_of_vertices() > before)
        {
            vh->info() = 0;
            prov.emplace(&*vh, si < feat_base ? 2 : 3);
        }
        surf_vh[si] = vh;
        shint = vh->cell();
        if constexpr (CCDT_MODE)
            vhash[hkey(sp3.x(), sp3.y(), sp3.z())].push_back(vh);
    }
    if (witness_overwrites &&
        (getenv("STIBIUM_DMESH_TIME") ||
         getenv("STIBIUM_DMESH_CHIP_DEBUG")))
        fprintf(stderr, "WITNESS: %llu sign witnesses overwritten "
                "to surface by the reuse path (review #3)\n",
                (unsigned long long)witness_overwrites);

    pt.mark("insert points");
    prog_stage(7);
    /*  Insert the refereed chain segments as constrained edges,
     *  batched with one Delaunay restoration at the end (each
     *  insert_constrained_edge(…, true) would restore
     *  individually).
     *
     *  Through-vertex pre-split (zeiss run, 2026-07-15): mechanical
     *  models put lattice vertices EXACTLY on straight creases
     *  (grid-aligned geometry at rational coordinates), and a
     *  constraint passing exactly through a non-endpoint vertex is
     *  a CGAL assertion ("The constraint passes through a
     *  vertex!") - one such vertex out of 345K threw the whole
     *  model back to unconstrained DT.  A vertex on the crease
     *  line IS a crease point: split the constraint there and
     *  chain through it.  */
    uint64_t constrained = 0, through_splits = 0;
    if constexpr (CCDT_MODE)
    {
        /*  Bin-grid of every finite vertex, one cell per lattice
         *  pitch, so each segment only inspects its own tube.  */
        const float bin = soup.spacing;
        const auto cell_key = [&](float x, float y, float z) {
            const int64_t ix = int64_t(floorf(x / bin));
            const int64_t iy = int64_t(floorf(y / bin));
            const int64_t iz = int64_t(floorf(z / bin));
            return (uint64_t(ix & 0x1FFFFF) << 42) |
                   (uint64_t(iy & 0x1FFFFF) << 21) |
                    uint64_t(iz & 0x1FFFFF);
        };
        std::unordered_map<uint64_t, std::vector<VH>> vgrid;
        vgrid.reserve(size_t(dt.number_of_vertices()));
        for (auto v = dt.finite_vertices_begin();
             v != dt.finite_vertices_end(); ++v)
            vgrid[cell_key(float(v->point().x()),
                           float(v->point().y()),
                           float(v->point().z()))].push_back(v);
        /*  Radius sized for T-JUNCTIONS (the probe dump's verdict):
         *  a chain terminating ON another chain's crease stops via
         *  the tight corner bisection (5e-4 sp), so its endpoint
         *  sits within ~1e-3 sp of the through-chain's segment
         *  interior (measured 3.47e-4 sp at (-5, 3.75, 26)) - an
         *  endpoint-to-interior contact CCDT_3 forbids, ground into
         *  a femtometer Steiner cascade.  Splitting the through-
         *  segment at that endpoint turns the T into a legal shared
         *  vertex.  5e-3 sp covers the class with margin and stays
         *  far below the drop-band floor (0.35 sp) that keeps
         *  legitimate vertices off constraint interiors.  */
        const float eps = 5e-3f * soup.spacing;
        const float eps2 = eps * eps;
        for (const auto& [a, b] : accepted)
        {
            VH va = surf_vh[a], vb = surf_vh[b];
            if (va == VH() || vb == VH() || va == vb)
                continue;
            const float ax = float(va->point().x()),
                        ay = float(va->point().y()),
                        az = float(va->point().z());
            const float bx = float(vb->point().x()),
                        by = float(vb->point().y()),
                        bz = float(vb->point().z());
            const float ux = bx - ax, uy = by - ay, uz = bz - az;
            const float uu = ux*ux + uy*uy + uz*uz;
            std::vector<std::pair<float, VH>> mids;
            if (uu > 0)
            {
                const int ix0 = int(floorf(std::min(ax, bx)/bin)) - 1,
                          ix1 = int(floorf(std::max(ax, bx)/bin)) + 1,
                          iy0 = int(floorf(std::min(ay, by)/bin)) - 1,
                          iy1 = int(floorf(std::max(ay, by)/bin)) + 1,
                          iz0 = int(floorf(std::min(az, bz)/bin)) - 1,
                          iz1 = int(floorf(std::max(az, bz)/bin)) + 1;
                for (int ix = ix0; ix <= ix1; ++ix)
                for (int iy = iy0; iy <= iy1; ++iy)
                for (int iz = iz0; iz <= iz1; ++iz)
                {
                    const uint64_t k =
                            (uint64_t(int64_t(ix) & 0x1FFFFF) << 42) |
                            (uint64_t(int64_t(iy) & 0x1FFFFF) << 21) |
                             uint64_t(int64_t(iz) & 0x1FFFFF);
                    const auto gi = vgrid.find(k);
                    if (gi == vgrid.end())
                        continue;
                    for (const VH& v : gi->second)
                    {
                        if (v == va || v == vb)
                            continue;
                        const float px = float(v->point().x()) - ax,
                                    py = float(v->point().y()) - ay,
                                    pz = float(v->point().z()) - az;
                        const float t = (px*ux + py*uy + pz*uz) / uu;
                        if (t <= 1e-4f || t >= 1.f - 1e-4f)
                            continue;
                        const float dx = px - t*ux, dy = py - t*uy,
                                    dz = pz - t*uz;
                        if (dx*dx + dy*dy + dz*dz < eps2)
                            mids.push_back({ t, v });
                    }
                }
            }
            std::sort(mids.begin(), mids.end(),
                      [](const auto& l, const auto& r) {
                          return l.first < r.first;
                      });
            VH prev = va;
            for (const auto& [t, v] : mids)
            {
                (void)t;
                if (v == prev)
                    continue;
                v->info() = 0;   // on the crease line, by measure
                dt.insert_constrained_edge(prev, v, false);
                ++constrained;
                ++through_splits;
                prev = v;
            }
            if (prev != vb)
            {
                dt.insert_constrained_edge(prev, vb, false);
                ++constrained;
            }
        }
        if (through_splits && getenv("STIBIUM_DMESH_SEG_DEBUG"))
            fprintf(stderr, "constraints: %llu split at "
                    "through-vertices\n",
                    (unsigned long long)through_splits);
        dt.restore_Delaunay();
        sweep_steiner();
    }

    pt.mark("constraints+restore");
    prog_stage(8);
    constexpr int MAX_ROUNDS = 48;
    uint64_t inserted = 0;
    int round = 0;
    for (; round < MAX_ROUNDS; ++round)
    {
        prog_frac(std::min(0.95f, float(round) / 8.f));
        if (*halt)
            break;
        std::vector<PendingEdge> pending;
        static size_t seen_hint = 0;
        std::unordered_set<uint64_t> seen;
        seen.reserve(seen_hint + 64);
        for (auto e = dt.finite_edges_begin();
             e != dt.finite_edges_end(); ++e)
        {
            const auto v1 = e->first->vertex(e->second);
            const auto v2 = e->first->vertex(e->third);
            const int s1 = v1->info(), s2 = v2->info();
            if (s1 * s2 != -1)
                continue;   // only direct inside<->outside edges
            const auto& p1 = v1->point();
            const auto& p2 = v2->point();
            PendingEdge ed;
            if (s1 < 0)
                ed = { float(p1.x()), float(p1.y()), float(p1.z()),
                       float(p2.x()), float(p2.y()), float(p2.z()) };
            else
                ed = { float(p2.x()), float(p2.y()), float(p2.z()),
                       float(p1.x()), float(p1.y()), float(p1.z()) };
            const uint64_t h = coord_hash(ed.ax, ed.ay, ed.az) ^
                    (coord_hash(ed.bx, ed.by, ed.bz) * 31);
            if (seen.insert(h).second)
                pending.push_back(ed);
        }
        seen_hint = seen.size();
        if (pending.empty())
            break;
        auto fresh = bisect_pending(deck, ctx, pending, halt);

        /*  Slide-away pass (the pinch preventer at the source): a
         *  separator landing within a fraction of a cell of an
         *  existing surface vertex makes a sliver tet pair whose
         *  extraction is a 4-triangle edge - the measured mechanism
         *  behind every pinch, DT and CCDT alike (provenance dump
         *  2026-07-15: refine-insert vs chain/Steiner vertex,
         *  0.15-0.25 cells).  The separator itself is needed - so
         *  slide it along the chord away from the crowding vertex
         *  and Newton-project it back onto the surface.  It still
         *  breaks its inside/outside edge (stays in the conflict
         *  zone); the near-coincidence is gone.  */
        {
            /*  0 disables: measured 2026-07-15, sliding stalls on
             *  short band edges (slid point exits the conflict
             *  zone, the i/o edge survives, holes at cap).  Kept
             *  for the next design round.  */
            static const char* slide_env =
                    getenv("STIBIUM_DMESH_SLIDE");
            const float pinch_r = (slide_env ? atof(slide_env)
                    : 0.0f) * soup.spacing;
            std::vector<float> mx, my, mz, mh2, mc2;
            std::vector<size_t> midx;
            for (size_t fi = 0; fi < fresh.size(); ++fi)
            {
                const DSurfPoint& s = fresh[fi];
                const auto nv = dt.nearest_vertex(
                        TPoint(s.x, s.y, s.z));
                if (nv == VH())
                    continue;
                /*  Nearest SURFACE vertex - the nearest vertex
                 *  outright may be a sign witness masking a
                 *  surface vertex a hair further out.  */
                VH q = VH();
                float d = 0;
                std::vector<VH> ring{ nv };
                dt.finite_adjacent_vertices(
                        nv, std::back_inserter(ring));
                for (const VH& v : ring)
                {
                    if (v->info() != 0)
                        continue;
                    const float ddx =
                            s.x - float(v->point().x());
                    const float ddy =
                            s.y - float(v->point().y());
                    const float ddz =
                            s.z - float(v->point().z());
                    const float dd = sqrtf(ddx*ddx + ddy*ddy +
                                           ddz*ddz);
                    if (q == VH() || dd < d)
                    {
                        q = v;
                        d = dd;
                    }
                }
                if (q == VH() || d >= pinch_r)
                    continue;
                const float dx = s.x - float(q->point().x());
                const float dy = s.y - float(q->point().y());
                const float dz = s.z - float(q->point().z());
                /*  Away from the crowding vertex; if coincident,
                 *  along the bisected edge instead.  */
                float ux = dx, uy = dy, uz = dz, ul = d;
                if (!(ul > 1e-12f))
                {
                    ux = pending[fi].bx - pending[fi].ax;
                    uy = pending[fi].by - pending[fi].ay;
                    uz = pending[fi].bz - pending[fi].az;
                    ul = sqrtf(ux*ux + uy*uy + uz*uz);
                    if (!(ul > 0))
                        continue;
                }
                const float slide = 1.5f * pinch_r - d;
                mx.push_back(s.x + ux / ul * slide);
                my.push_back(s.y + uy / ul * slide);
                mz.push_back(s.z + uz / ul * slide);
                mh2.push_back(soup.spacing * 0.01f);
                mc2.push_back(soup.spacing);
                midx.push_back(fi);
            }
            if (!mx.empty())
            {
                project_points_impl(deck, ctx, mx, my, mz,
                                    mh2, mc2);
                for (size_t q = 0; q < midx.size(); ++q)
                    fresh[midx[q]] = { mx[q], my[q], mz[q] };
            }
        }

        uint64_t round_added = 0;
        for (const DSurfPoint& s : fresh)
        {
            /*  A bisected point can converge onto an EXISTING vertex
             *  (grid-aligned geometry: surfaces lying on lattice
             *  planes).  Never overwrite a sign witness's info -
             *  count vertices to detect coincidence.  */
            const auto before = dt.number_of_vertices();
            auto vh = dt.insert(TPoint(s.x, s.y, s.z));
            if (dt.number_of_vertices() > before)
            {
                vh->info() = 0;
                prov.emplace(&*vh, 4);
                ++inserted;
                ++round_added;
            }
        }
        sweep_steiner();
        if (!round_added)
            break;   // every insert coincided: no progress possible
    }
    /*  Which exit: converged (no pending in/out edges), stalled
     *  (inserts all coincided - the survivors get spanned by
     *  extraction), or round-capped.  A stall/cap leaves i/o edges
     *  in the DT and extraction chords them - chip factory.  */
    if (getenv("STIBIUM_DMESH_TIME"))
    {
        size_t leftover = 0;
        for (auto e = dt.finite_edges_begin();
             e != dt.finite_edges_end(); ++e)
            if (e->first->vertex(e->second)->info() *
                e->first->vertex(e->third)->info() == -1)
                ++leftover;
        fprintf(stderr, "REFINE: %d rounds, %llu inserted, "
                "%zu i/o edges leftover%s\n", round,
                (unsigned long long)inserted, leftover,
                round >= MAX_ROUNDS ? " (ROUND CAP)"
                : leftover ? " (STALLED)" : "");
    }
    pt.mark("refinement");
    prog_stage(9);

    /*  Stage C runs inside the repair loop: extract, then hunt WARTS
     *  - fold edges (adjacent-triangle normals disagreeing) whose
     *  midpoint is OFF the surface.  A fold on a true crease has its
     *  edge ON the surface and is left alone.  Wart midpoints are
     *  Newton-projected onto the surface and inserted, and the mesh
     *  re-extracted, until clean or capped (error-driven insertion,
     *  after Wang et al. 2025).  */
    /*  Repair-round budget (STIBIUM_DMESH_ROUNDS overrides): the
     *  bench minis converge in 3-4 rounds, but zeiss exits by CAP
     *  with work left on the table (333K inserts and climbing) -
     *  Nate's budget theory, and the v1-kindness finding (brute
     *  repair volume doubled as detail rescue), both want this
     *  dial.  */
    static const char* rounds_env = getenv("STIBIUM_DMESH_ROUNDS");
    const int MAX_REPAIR = rounds_env && atoi(rounds_env) > 0
            ? atoi(rounds_env) : 16;
    uint64_t repaired_total = 0;
    int repair_round = 0;
    std::vector<VH> xvh;   // extraction order -> vertex handle
    std::vector<std::pair<CH, CH>> tri_cells;  // per emitted tri
    /*  Residual chip edges of the FINAL extraction, stashed for the
     *  crease-snap pass (refilled every detection round; the
     *  vertex/triangle indices go stale the moment a re-extraction
     *  happens, so the MAX_REPAIR exit clears them).  */
    std::vector<uint32_t> snap_va, snap_vb, snap_t1, snap_t2;
    /*  Stall exit (perf round: repair was 64% of the export and
     *  zeiss's depth flatlines by round ~3 while rounds 4-15 pay
     *  full re-detection price for nothing - and the budget
     *  experiment proved extra rounds only convert churn into
     *  topology damage).  Two consecutive rounds without depth
     *  improvement = churn ahead, stop; the stash keeps the
     *  residual chips for the snap pass.  */
    float stall_prev = 1e30f;
    int stall_rounds = 0;
    std::vector<float> snap_d;   // chip depth at detection
    for (;; ++repair_round)
    {
    prog_frac(std::min(0.95f, float(repair_round) / 5.f));
    /*  Stage C: cell signs.  After convergence a finite cell cannot
     *  contain both signs (except in the crease band, where i/o
     *  edges are deliberately shadowed), so a non-surface vertex
     *  decides it; all-surface AND mixed cells ask the oracle
     *  (batched).  Infinite cells are outside by definition.
     *
     *  All-surface cells get DC semantics (STIBIUM_DMESH_DC=0
     *  reverts to centroid-only): at a concave crease the wedge
     *  corner lives inside a tet whose four vertices are all on
     *  the surface and whose CENTROID is in air - centroid-only
     *  classification clips the corner, and the clipped facet is
     *  exactly the measured chip (worst-depth plateau ~0.196 sp,
     *  invariant under 1x/2x/4x band density and both repair
     *  strategies: the repair keep-out and the chip target scale
     *  together, so insert-repair can never reach it).  Probe the
     *  four face barycenters too, and call the cell INSIDE when
     *  any probe reads materially negative: the chip facet's own
     *  barycenter is the deep-material witness.  Bias is one-
     *  sided by design - at a convex crease the mirror case fills
     *  an air notch (volume-adding, the visually-benign species,
     *  and the wart machinery already handles it).  */
    static const char* dc_env = getenv("STIBIUM_DMESH_DC");
    const bool dc_cells = !dc_env || atoi(dc_env) != 0;
    std::vector<CH> oracle_cells;
    std::vector<uint8_t> oracle_allsurf;
    std::vector<float> cxs, cys, czs, cvals;
    std::vector<float> ctol;   // per-cell material threshold
    for (auto c = dt.finite_cells_begin();
         c != dt.finite_cells_end(); ++c)
    {
        int8_t sign = 0;
        bool mixed = false;
        for (int i = 0; i < 4; ++i)
        {
            const int8_t s = c->vertex(i)->info();
            if (!s)
                continue;
            if (sign && s != sign)
            {
                mixed = true;
                break;
            }
            sign = s;
        }
        if (!sign || mixed)
        {
            double p[4][3];
            for (int i = 0; i < 4; ++i)
            {
                p[i][0] = c->vertex(i)->point().x();
                p[i][1] = c->vertex(i)->point().y();
                p[i][2] = c->vertex(i)->point().z();
            }
            const double x = p[0][0]+p[1][0]+p[2][0]+p[3][0],
                         y = p[0][1]+p[1][1]+p[2][1]+p[3][1],
                         z = p[0][2]+p[1][2]+p[2][2]+p[3][2];
            oracle_cells.push_back(c);
            /*  DC bias is suppressed INSIDE the crease band
             *  (plinth conviction, 2026-07-18): the bias was
             *  built the morning the tracer landed, for corner
             *  tets at UNCONSTRAINED creases.  Beside a
             *  constrained crease the corner is already
             *  articulated by law, and the barycenter probe
             *  swallows the junction's AIR wedge instead -
             *  facets roof from wall to face over the crease
             *  fence (88 tilted "teeth" on the plinth, 0 with
             *  the band suppression; chips unchanged - law+snap
             *  already cover the band).  The volume-ADDING
             *  species Nate convicted by eye days ago.
             *  STIBIUM_DMESH_DC_BAND, cells; 0 = old blanket
             *  bias.  */
            static const char* dcb_env =
                    getenv("STIBIUM_DMESH_DC_BAND");
            const float dcb = (dcb_env ? float(atof(dcb_env))
                                       : 0.75f) * soup.spacing;
            const bool allsurf = dc_cells && !mixed && !sign &&
                    !(dcb > 0 && !cseg.empty() &&
                      near_crease(float(x / 4), float(y / 4),
                                  float(z / 4), dcb));
            oracle_allsurf.push_back(allsurf ? 1 : 0);
            cxs.push_back(float(x / 4));
            cys.push_back(float(y / 4));
            czs.push_back(float(z / 4));
            if (allsurf)
            {
                /*  Four face barycenters (face i omits vertex i)  */
                for (int i = 0; i < 4; ++i)
                {
                    cxs.push_back(float((x - p[i][0]) / 3));
                    cys.push_back(float((y - p[i][1]) / 3));
                    czs.push_back(float((z - p[i][2]) / 3));
                }
                /*  Material bar: same 3%-of-edge species test the
                 *  chip detector uses, on the cell's longest edge
                 *  (near-tangent faces of legitimate air cells
                 *  read |f| ~ sagitta, far below it).  */
                double emax2 = 0;
                for (int i = 0; i < 4; ++i)
                    for (int j = i + 1; j < 4; ++j)
                    {
                        const double dx = p[i][0]-p[j][0],
                                     dy = p[i][1]-p[j][1],
                                     dz = p[i][2]-p[j][2];
                        const double e2 = dx*dx + dy*dy + dz*dz;
                        if (e2 > emax2)
                            emax2 = e2;
                    }
                ctol.push_back(0.03f * float(sqrt(emax2)));
            }
            else
                ctol.push_back(0);
            c->info() = 0;
        }
        else
            c->info() = sign;
    }
    if (!oracle_cells.empty())
    {
        eval_points_mt(deck, ctx, cxs, cys, czs, cvals);
        size_t at = 0;
        for (size_t i = 0; i < oracle_cells.size(); ++i)
        {
            if (oracle_allsurf[i])
            {
                float fmin = cvals[at];
                for (int q = 1; q < 5; ++q)
                    fmin = std::min(fmin, cvals[at + q]);
                oracle_cells[i]->info() =
                        fmin < -ctol[i] ? -1
                        : cvals[at] < 0 ? -1 : 1;
                at += 5;
            }
            else
            {
                oracle_cells[i]->info() = cvals[at] < 0 ? -1 : 1;
                ++at;
            }
        }
    }

    /*  Degenerate-sliver disenfranchisement (the z-fighting
     *  conviction, 2026-07-18): lattice-aligned flats breed
     *  ZERO-VOLUME tets lying flat in the surface plane - all
     *  four vertices on-surface, so the oracle's probes all read
     *  |f| ~ 0 and the cell's sign is rounding noise.  A noise
     *  sign disagreeing with a neighbor emits a facet through a
     *  cell with no interior: a second coincident layer in the
     *  plane (fightpix: ~3,600 fighting pixels on screws, base
     *  and THIN alike), with noise orientation (the apex dot is
     *  the same degenerate volume).  A cell with no volume has
     *  no geometric interior - ANY sign is consistent - so it
     *  adopts the majority sign of its finite neighbors, which
     *  cancels every facet through it.  Scale referee: flat
     *  slivers read |V|/emax^3 ~ 1e-6, the DC wedge tets the
     *  bias exists for read ~1e-2 - the 1e-4 bar splits the
     *  measured populations by two orders each side.  Stacked
     *  slivers settle by fixpoint (guarded).  STIBIUM_DMESH_
     *  DEGEN=0 disables.  */
    {
        static const char* dg_env = getenv("STIBIUM_DMESH_DEGEN");
        const float dg_bar = dg_env ? float(atof(dg_env)) : 1e-4f;
        /*  Pancake discriminator (fold-cure round 2, 2026-07-18):
         *  the REPAIR-MINTED pancake population (rel-vol 1e-3..
         *  1e-2, growing per round - Newton-projected inserts land
         *  exactly on flats and mint near-zero cells between old
         *  surface vertices) shares a volume decade with the DC
         *  wedge tets, so volume alone cannot disenfranchise it.
         *  The field can: a PANCAKE is thinner than the razor bar
         *  (height < STIBIUM_DMESH_PANCAKE * sp, default 0.02 -
         *  sub-print noise by definition) AND gradient-UNIFORM
         *  across its four corners (one flat patch).  A wedge
         *  holds a real corner: its corner gradients SPLIT - the
         *  very reason DC semantics exist (screws: 90 degrees,
         *  dot ~ 0, nowhere near the 0.9 bar).  Thin + uniform
         *  adopts neighbor signs like the float-degenerates; thin
         *  + split keeps its DC vote.  0 disables.  */
        static const char* pk_env = getenv("STIBIUM_DMESH_PANCAKE");
        const float pk_bar = (pk_env ? float(atof(pk_env)) : 0.02f)
                * soup.spacing;
        if (dg_bar > 0 || pk_bar > 0)
        {
            uint64_t vhist[10] = {};   // log10 rel-vol buckets
            std::vector<CH> degen;
            std::vector<CH> pk_cells;
            std::vector<std::array<VH, 4>> pk_verts;
            for (auto c = dt.finite_cells_begin();
                 c != dt.finite_cells_end(); ++c)
            {
                bool allsurf = true;
                for (int i = 0; i < 4 && allsurf; ++i)
                    allsurf = c->vertex(i)->info() == 0;
                if (!allsurf)
                    continue;
                double p[4][3];
                for (int i = 0; i < 4; ++i)
                {
                    p[i][0] = c->vertex(i)->point().x();
                    p[i][1] = c->vertex(i)->point().y();
                    p[i][2] = c->vertex(i)->point().z();
                }
                double emax2 = 0;
                for (int i = 0; i < 4; ++i)
                    for (int j = i + 1; j < 4; ++j)
                    {
                        const double dx = p[i][0]-p[j][0],
                                     dy = p[i][1]-p[j][1],
                                     dz = p[i][2]-p[j][2];
                        emax2 = std::max(emax2,
                                dx*dx + dy*dy + dz*dz);
                    }
                const double ax = p[1][0]-p[0][0],
                             ay = p[1][1]-p[0][1],
                             az = p[1][2]-p[0][2];
                const double bx = p[2][0]-p[0][0],
                             by = p[2][1]-p[0][1],
                             bz = p[2][2]-p[0][2];
                const double cx2 = p[3][0]-p[0][0],
                             cy2 = p[3][1]-p[0][1],
                             cz2 = p[3][2]-p[0][2];
                const double vol = fabs(
                        ax*(by*cz2 - bz*cy2) +
                        ay*(bz*cx2 - bx*cz2) +
                        az*(bx*cy2 - by*cx2)) / 6.0;
                const double e3 = emax2 * sqrt(emax2);
                if (!(e3 > 0))
                    continue;
                const double rel = vol / e3;
                int b = rel <= 0 ? 0
                        : int(9 + floor(log10(rel)));
                ++vhist[b < 0 ? 0 : b > 9 ? 9 : b];
                if (dg_bar > 0 && rel < dg_bar)
                {
                    degen.push_back(c);
                    continue;
                }
                if (!(pk_bar > 0))
                    continue;
                /*  Height over the largest face: h = 3V / Amax  */
                double amax = 0;
                for (int i = 0; i < 4; ++i)
                {
                    const double* q0 = p[(i + 1) & 3];
                    const double* q1 = p[(i + 2) & 3];
                    const double* q2 = p[(i + 3) & 3];
                    const double ux2 = q1[0]-q0[0],
                            uy2 = q1[1]-q0[1], uz2 = q1[2]-q0[2];
                    const double vx3 = q2[0]-q0[0],
                            vy3 = q2[1]-q0[1], vz3 = q2[2]-q0[2];
                    const double fx = uy2*vz3 - uz2*vy3;
                    const double fy = uz2*vx3 - ux2*vz3;
                    const double fz = ux2*vy3 - uy2*vx3;
                    amax = std::max(amax, 0.5 *
                            sqrt(fx*fx + fy*fy + fz*fz));
                }
                if (amax > 0 && 3.0 * vol / amax < pk_bar)
                {
                    pk_cells.push_back(c);
                    pk_verts.push_back({ c->vertex(0),
                            c->vertex(1), c->vertex(2),
                            c->vertex(3) });
                }
            }
            uint64_t pancakes = 0;
            if (!pk_cells.empty())
            {
                /*  One gradient per unique vertex (7-pt stencil,
                 *  batched); a cell is a pancake only when all
                 *  four corner gradients agree (min pairwise dot
                 *  > 0.9 ~ 26 degrees - a concave wedge reads
                 *  ~0).  Zero-length gradients keep the DC vote
                 *  (bad stencil is not evidence).  */
                std::unordered_map<const void*, uint32_t> gidx;
                std::vector<float> gx2, gy2, gz2, gv2;
                const float gh = 0.01f * soup.spacing;
                for (const auto& vv : pk_verts)
                    for (const VH& v : vv)
                        if (!gidx.count(&*v))
                        {
                            const uint32_t id =
                                    uint32_t(gidx.size());
                            gidx.emplace(&*v, id);
                            const float vx4 =
                                    float(v->point().x());
                            const float vy4 =
                                    float(v->point().y());
                            const float vz4 =
                                    float(v->point().z());
                            for (int o = 0; o < 6; ++o)
                            {
                                gx2.push_back(vx4 +
                                        (o==0?gh:o==1?-gh:0));
                                gy2.push_back(vy4 +
                                        (o==2?gh:o==3?-gh:0));
                                gz2.push_back(vz4 +
                                        (o==4?gh:o==5?-gh:0));
                            }
                        }
                eval_points_mt(deck, ctx, gx2, gy2, gz2,
                            gv2);
                std::vector<std::array<float, 3>> grad(
                        gidx.size());
                for (size_t i = 0; i < gidx.size(); ++i)
                {
                    const float dx2 = gv2[i*6] - gv2[i*6+1];
                    const float dy2 = gv2[i*6+2] - gv2[i*6+3];
                    const float dz2 = gv2[i*6+4] - gv2[i*6+5];
                    const float l = sqrtf(dx2*dx2 + dy2*dy2 +
                                          dz2*dz2);
                    grad[i] = l > 0
                            ? std::array<float, 3>{ dx2/l, dy2/l,
                                                    dz2/l }
                            : std::array<float, 3>{ 0, 0, 0 };
                }
                uint64_t dhist[12] = {};  // min-dot [-1,1] buckets
                for (size_t i = 0; i < pk_cells.size(); ++i)
                {
                    float mind = 1;
                    bool bad = false;
                    for (int a2 = 0; a2 < 4 && !bad; ++a2)
                    {
                        const auto& ga =
                                grad[gidx[&*pk_verts[i][a2]]];
                        if (ga[0] == 0 && ga[1] == 0 &&
                            ga[2] == 0)
                            bad = true;
                        for (int b2 = a2 + 1; b2 < 4; ++b2)
                        {
                            const auto& gb =
                                    grad[gidx[&*pk_verts[i][b2]]];
                            mind = std::min(mind,
                                    ga[0]*gb[0] + ga[1]*gb[1] +
                                    ga[2]*gb[2]);
                        }
                    }
                    if (!bad)
                    {
                        int b3 = int((mind + 1) * 6);
                        ++dhist[b3 < 0 ? 0 : b3 > 11 ? 11 : b3];
                    }
                    /*  The min-dot histogram (screws, 2026-07-18)
                     *  killed the uniform-gradient theory: the
                     *  thin population is ~0 uniform, ~1,050 at
                     *  dot 0 (crease-straddlers, GROWING per
                     *  repair round - the fighting-layer growth
                     *  signature) and ~100 at dot -1 (thin-wall
                     *  spanners).  The clean rule is HEIGHT
                     *  alone: a cell under 0.02 sp cannot hold a
                     *  product-bar feature whatever its
                     *  gradients say - its vote moves the
                     *  surface by less than its own sub-print
                     *  thickness.  The DC wedge cure lives at
                     *  0.03 sp+ with real cell heights - orders
                     *  away from the bar.  The hist stays as the
                     *  population instrument.  */
                    (void)mind;
                    degen.push_back(pk_cells[i]);
                    ++pancakes;
                }
                if (getenv("STIBIUM_DMESH_CHIP_DEBUG"))
                {
                    fprintf(stderr, "PANCAKE min-dot hist "
                            "[-1..1]: ");
                    for (int b3 = 0; b3 < 12; ++b3)
                        fprintf(stderr, "%llu ",
                                (unsigned long long)dhist[b3]);
                    fprintf(stderr, "\n");
                }
            }
            uint64_t adopted = 0;
            for (int round = 0; round < 8; ++round)
            {
                bool changed = false;
                for (const CH& c : degen)
                {
                    int votes = 0;
                    for (int i = 0; i < 4; ++i)
                    {
                        const CH n = c->neighbor(i);
                        votes += dt.is_infinite(n) ? 1
                                : int(n->info());
                    }
                    /*  Ties break INSIDE, deterministically: a
                     *  pancake between solid-below and air-above
                     *  reads 2 v 2, and stacked pancakes with
                     *  alternating noise signs emit a coincident
                     *  facet at EVERY flip (the doubled layers).
                     *  Any uniform choice makes a stack monotone
                     *  - one emission, at the stack's air side.  */
                    const int8_t want = votes > 0 ? 1 : -1;
                    if (c->info() != want)
                    {
                        c->info() = want;
                        changed = true;
                        ++adopted;
                    }
                }
                if (!changed)
                    break;
            }
            if (getenv("STIBIUM_DMESH_CHIP_DEBUG"))
            {
                fprintf(stderr, "DEGEN: %llu adopted of %zu degen "
                        "(%llu pancakes of %zu thin candidates); "
                        "all-surf rel-vol hist (<=1e-9..>=1): ",
                        (unsigned long long)adopted, degen.size(),
                        (unsigned long long)pancakes,
                        pk_cells.size());
                for (int b = 0; b < 10; ++b)
                    fprintf(stderr, "%llu ",
                            (unsigned long long)vhist[b]);
                fprintf(stderr, "\n");
            }
        }
    }

    /*  Extraction: facets whose three corners are surface vertices,
     *  between cells of opposite sign; oriented so the triangle is
     *  CCW seen from the outside cell.  Facets of a tet complex
     *  cannot self-intersect.  */
    std::unordered_map<const void*, uint32_t> vidx;
    out->verts.clear();
    out->tris.clear();
    xvh.clear();
    tri_cells.clear();
    const auto vertex_index = [&](VH v) -> uint32_t {
        const auto found = vidx.find(&*v);
        if (found != vidx.end())
            return found->second;
        const uint32_t id = uint32_t(out->verts.size() / 3);
        vidx.emplace(&*v, id);
        xvh.push_back(v);
        out->verts.push_back(float(v->point().x()));
        out->verts.push_back(float(v->point().y()));
        out->verts.push_back(float(v->point().z()));
        return id;
    };

    for (auto f = dt.finite_facets_begin();
         f != dt.finite_facets_end(); ++f)
    {
        const CH c = f->first;
        const int i = f->second;
        const CH n = c->neighbor(i);
        const int8_t cs = dt.is_infinite(c) ? 1 : c->info();
        const int8_t ns = dt.is_infinite(n) ? 1 : n->info();
        if (cs == ns || cs == 0 || ns == 0)
            continue;
        VH vs[3];
        for (int q = 0, w = 0; q < 4; ++q)
            if (q != i)
                vs[w++] = c->vertex(q);
        if (vs[0]->info() != 0 || vs[1]->info() != 0 ||
            vs[2]->info() != 0)
            continue;

        /*  Orient: normal toward the OUTSIDE cell's fourth vertex
         *  side.  c's fourth vertex is c->vertex(i); if c is the
         *  inside cell, the outside direction is away from it.  */
        const auto& p0 = vs[0]->point();
        const auto& p1 = vs[1]->point();
        const auto& p2 = vs[2]->point();
        const double ux = p1.x() - p0.x(), uy = p1.y() - p0.y(),
                     uz = p1.z() - p0.z();
        const double wx = p2.x() - p0.x(), wy = p2.y() - p0.y(),
                     wz = p2.z() - p0.z();
        const double nx = uy * wz - uz * wy;
        const double ny = uz * wx - ux * wz;
        const double nz = ux * wy - uy * wx;
        const auto& apex = c->vertex(i)->point();   // on cell c's side
        const double dx = apex.x() - p0.x(), dy = apex.y() - p0.y(),
                     dz = apex.z() - p0.z();
        const double d = nx * dx + ny * dy + nz * dz;
        /*  If c is outside, normal should point toward apex (d>0);
         *  if c is inside, away from it.  */
        const bool flip = (cs > 0) ? (d < 0) : (d > 0);
        const uint32_t a = vertex_index(vs[0]);
        const uint32_t b = vertex_index(flip ? vs[2] : vs[1]);
        const uint32_t cc = vertex_index(flip ? vs[1] : vs[2]);
        out->tris.push_back(a);
        out->tris.push_back(b);
        out->tris.push_back(cc);
        tri_cells.push_back({ c, n });
    }

    snap_va.clear();   // stale after this re-extraction
    snap_vb.clear();
    snap_t1.clear();
    snap_t2.clear();
    snap_d.clear();
    {
        char nb[32];
        snprintf(nb, sizeof(nb), "r%d signs+extract",
                 repair_round);
        pt.sub(nb);
    }
    if (repair_round >= MAX_REPAIR || *halt)
        break;

    /*  Wart hunt on the freshly extracted mesh.  Two strategies
     *  behind STIBIUM_DMESH_REPAIR:
     *   1 - fold-gated: only edges whose adjacent triangles
     *     disagree past ~78 degrees are checked.  Conservative;
     *     misses shallow-crease chords (on-surface warts remain).
     *   2 (default) - gate-free crease-seek: EVERY edge midpoint
     *     checked; wart edges split at the |f| peak along the
     *     chord (the branch switch = the crease), with a crowding
     *     guard (no inserts within a quarter edge of an existing
     *     vertex).  The guard is what makes it converge: without
     *     it, crease-crossing chords are self-similar and the loop
     *     churns into sliver pinches; with it, union-crease warts
     *     press down in ~5 rounds with zero topology damage.  */
    static const char* rep_env = getenv("STIBIUM_DMESH_REPAIR");
    const int repair_mode = rep_env ? atoi(rep_env) : 2;

    std::unordered_map<uint64_t, std::pair<uint32_t, uint32_t>> em;
    for (uint32_t t = 0; t < out->tris.size() / 3; ++t)
        for (int e = 0; e < 3; ++e)
        {
            uint64_t a = out->tris[3*t + e];
            uint64_t b = out->tris[3*t + (e + 1) % 3];
            if (a > b)
                std::swap(a, b);
            const uint64_t k = (a << 32) | b;
            auto it = em.find(k);
            if (it == em.end())
                em[k] = { t, UINT32_MAX };
            else if (it->second.second == UINT32_MAX)
                it->second.second = t;
        }
    const auto tri_normal = [&](uint32_t t, double n[3]) {
        const float* a = &out->verts[3 * out->tris[3*t]];
        const float* b = &out->verts[3 * out->tris[3*t + 1]];
        const float* c2 = &out->verts[3 * out->tris[3*t + 2]];
        const double ux = b[0]-a[0], uy = b[1]-a[1], uz = b[2]-a[2];
        const double wx = c2[0]-a[0], wy = c2[1]-a[1], wz = c2[2]-a[2];
        n[0] = uy*wz - uz*wy;
        n[1] = uz*wx - ux*wz;
        n[2] = ux*wy - uy*wx;
        const double l = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        if (l > 0) { n[0] /= l; n[1] /= l; n[2] /= l; }
    };
    std::vector<float> wx2, wy2, wz2, whs, wclamp;
    std::vector<float> eax, eay, eaz, ebx, eby, ebz;
    std::vector<uint32_t> wt1, wt2, wva, wvb;
    for (const auto& [k, pr] : em)
    {
        if (repair_mode < 2)
        {
            if (pr.second == UINT32_MAX)
                continue;
            double n1[3], n2[3];
            tri_normal(pr.first, n1);
            tri_normal(pr.second, n2);
            if (n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2] >= 0.2)
                continue;   // no fold here
        }
        const uint32_t va = uint32_t(k >> 32), vb = uint32_t(k);
        const float* A2 = &out->verts[3*va];
        const float* B2 = &out->verts[3*vb];
        const float ex = B2[0]-A2[0], ey = B2[1]-A2[1],
                    ez = B2[2]-A2[2];
        const float len = sqrtf(ex*ex + ey*ey + ez*ez);
        if (!(len > 0))
            continue;
        wx2.push_back((A2[0]+B2[0]) * 0.5f);
        wy2.push_back((A2[1]+B2[1]) * 0.5f);
        wz2.push_back((A2[2]+B2[2]) * 0.5f);
        whs.push_back(len * 0.01f);
        wclamp.push_back(len);
        eax.push_back(A2[0]);  eay.push_back(A2[1]);
        eaz.push_back(A2[2]);
        ebx.push_back(B2[0]);  eby.push_back(B2[1]);
        ebz.push_back(B2[2]);
        wt1.push_back(pr.first);
        wt2.push_back(pr.second);
        wva.push_back(va);
        wvb.push_back(vb);
    }
    if (wx2.empty())
        break;

    /*  Off-surface test: |f| / |grad| estimates distance; only
     *  midpoints beyond 5% of their edge length are warts.  */
    {
        Tape* base2 = deck_base(deck);
        std::vector<float> gxs(wx2.size()*7), gys(wx2.size()*7),
                gzs(wx2.size()*7), gv;
        for (size_t i = 0; i < wx2.size(); ++i)
        {
            for (int q = 0; q < 7; ++q)
            {
                gxs[i*7+q] = wx2[i];
                gys[i*7+q] = wy2[i];
                gzs[i*7+q] = wz2[i];
            }
            gxs[i*7+1] += whs[i];  gxs[i*7+2] -= whs[i];
            gys[i*7+3] += whs[i];  gys[i*7+4] -= whs[i];
            gzs[i*7+5] += whs[i];  gzs[i*7+6] -= whs[i];
        }
        eval_points_mt(deck, ctx, gxs, gys, gzs, gv);
        std::vector<size_t> cand;
        std::vector<uint8_t> cspec;   // 1 = chip (f < 0)
        std::vector<float> cdist;     // |f|/|grad| at detection
        uint64_t n_chip = 0, n_wart = 0;
        float worst_chip = 0;
        size_t worst_i = SIZE_MAX;
        for (size_t i = 0; i < wx2.size(); ++i)
        {
            const float f = gv[i*7];
            const float h2 = 2 * whs[i];
            const float gx = (gv[i*7+1]-gv[i*7+2]) / h2;
            const float gy = (gv[i*7+3]-gv[i*7+4]) / h2;
            const float gz = (gv[i*7+5]-gv[i*7+6]) / h2;
            const float gl = sqrtf(gx*gx + gy*gy + gz*gz);
            if (!(gl > 0) || !std::isfinite(f))
                continue;
            const float dist = fabsf(f) / gl;
            /*  Credibility gate (zeiss run, 2026-07-15): a chord
             *  with both endpoints on the surface cannot sag more
             *  than half its length; central differences at a
             *  wedge underestimate |grad| (they average the two
             *  face normals), inflating honest readings by
             *  1/cos(half-angle) - csg's real chips read
             *  0.53 x edge.  2 x edge admits every plausible
             *  wedge and still rejects gradient garbage (blend
             *  plateaus, scale nodes - the merged zeiss read
             *  576 sp on a 1.24 sp edge and the repair loop
             *  churned 300K inserts chasing it; repairs aimed by
             *  a garbage gradient insert garbage points).  */
            if (dist > 2.0f * wclamp[i])
                continue;
            if (dist > wclamp[i] * 0.03f)
            {
                cand.push_back(i);
                cdist.push_back(dist);
                /*  Sign of f at the sagging midpoint: f > 0 =
                 *  chord through EMPTY space (a wart), f < 0 =
                 *  chord through MATERIAL (a chip - the concave
                 *  species, diagnosed by eyeball 2026-07-15).  */
                cspec.push_back(f < 0 ? 1 : 0);
                (f < 0 ? n_chip : n_wart)++;
                if (f < 0 && dist > worst_chip)
                {
                    worst_chip = dist;
                    worst_i = i;
                }
            }
        }
        static const char* chip_env =
                getenv("STIBIUM_DMESH_CHIP_DEBUG");
        if (chip_env && (n_chip || n_wart))
            fprintf(stderr, "REPAIR round %d: %llu chips (f<0), "
                    "%llu warts (f>0), worst depth %.3f sp\n",
                    repair_round,
                    (unsigned long long)n_chip,
                    (unsigned long long)n_wart,
                    worst_chip / soup.spacing);
        if (chip_env && worst_i != SIZE_MAX)
            fprintf(stderr, "  worst chip at (%.4f, %.4f, %.4f) "
                    "edge %.3f sp  A(%.4f, %.4f, %.4f) "
                    "B(%.4f, %.4f, %.4f)\n",
                    wx2[worst_i], wy2[worst_i], wz2[worst_i],
                    wclamp[worst_i] / soup.spacing,
                    eax[worst_i], eay[worst_i], eaz[worst_i],
                    ebx[worst_i], eby[worst_i], ebz[worst_i]);
        for (size_t q = 0; q < cand.size(); ++q)
            if (cspec[q])
            {
                const size_t i = cand[q];
                snap_va.push_back(wva[i]);
                snap_vb.push_back(wvb[i]);
                snap_t1.push_back(wt1[i]);
                snap_t2.push_back(wt2[i]);
                snap_d.push_back(cdist[q]);
            }
        if (cand.empty())
            break;
        if (worst_chip >= stall_prev * 0.98f)
        {
            /*  Stall tolerance (STIBIUM_DMESH_STALL, rounds of
             *  non-improvement to forgive; default 1, old
             *  behavior = 2).  TIME=2 anatomy, r2 bino: depth
             *  plateaus by round 2 and rounds 3-4 cost ~20 s
             *  each for zero depth gain.  The round >= 2 floor
             *  protects the measured stall-then-improve at
             *  rounds 1 -> 2 (bino 0.339 -> 0.242).  */
            const char* st_env =
                    getenv("STIBIUM_DMESH_STALL");
            const int st_bar = st_env ? atoi(st_env) : 1;
            if (++stall_rounds >= st_bar && repair_round >= 2)
            {
                if (chip_env)
                    fprintf(stderr, "REPAIR: depth stalled at "
                            "%.3f sp, stopping after round %d\n",
                            worst_chip / soup.spacing,
                            repair_round);
                break;
            }
        }
        else
            stall_rounds = 0;
        stall_prev = std::min(stall_prev, worst_chip);

        /*  TBP refinement targets for chips (DelIso lineage; see
         *  doc/research/2026-07-15-pinch-manifoldness.md): a divot
         *  facet's dual Voronoi edge crosses the surface more than
         *  once, and the crossing point - the surface-Delaunay-
         *  ball center - is ON the surface and inside both
         *  adjacent circumballs, so inserting it destroys the
         *  facet.  A principled target, unlike the chord-derived
         *  midpoint.  Dense-sample each chip's adjacent facets'
         *  dual segments, linear-interp the crossings (projection
         *  Newton-polishes afterward), take the crossing nearest
         *  the sagging midpoint.  */
        std::vector<float> sdbx(cand.size()), sdby(cand.size()),
                sdbz(cand.size());
        std::vector<uint8_t> sdb_ok(cand.size(), 0);
        /*  Opt-in (STIBIUM_DMESH_SDB=1): measured 2026-07-15 on
         *  csg at 96^3, the depth trajectory (worst chip 0.500 ->
         *  0.196 sp over 4 rounds) is IDENTICAL to plain midpoint
         *  repair - the plateau is set by the guards and the
         *  sampling density, not the target choice.  Kept as a
         *  knob for revisiting at other densities.  */
        static const char* sdb_env = getenv("STIBIUM_DMESH_SDB");
        if (repair_mode >= 2 && sdb_env && atoi(sdb_env) != 0)
        {
            constexpr int DS = 33;
            struct Seg
            {
                size_t q;
                float ax, ay, az, bx, by, bz;
                size_t o;
            };
            std::vector<Seg> segs;
            std::vector<float> px, py, pz, pv;
            for (size_t q = 0; q < cand.size(); ++q)
            {
                if (!cspec[q])
                    continue;
                const size_t i = cand[q];
                for (const uint32_t t : { wt1[i], wt2[i] })
                {
                    if (t == UINT32_MAX || t >= tri_cells.size())
                        continue;
                    const CH tc = tri_cells[t].first;
                    const CH tn = tri_cells[t].second;
                    if (dt.is_infinite(tc) || dt.is_infinite(tn))
                        continue;
                    /*  Kernel circumcenter (our cell bases carry
                     *  info, not the Delaunay cell base's cached
                     *  circumcenter, so dt.dual() won't compile).  */
                    const auto d1 = CGAL::circumcenter(
                            tc->vertex(0)->point(),
                            tc->vertex(1)->point(),
                            tc->vertex(2)->point(),
                            tc->vertex(3)->point());
                    const auto d2 = CGAL::circumcenter(
                            tn->vertex(0)->point(),
                            tn->vertex(1)->point(),
                            tn->vertex(2)->point(),
                            tn->vertex(3)->point());
                    Seg s;
                    s.q = q;
                    s.o = px.size();
                    s.ax = float(d1.x());
                    s.ay = float(d1.y());
                    s.az = float(d1.z());
                    s.bx = float(d2.x());
                    s.by = float(d2.y());
                    s.bz = float(d2.z());
                    /*  Degenerate tets throw circumcenters far
                     *  away; a dual segment much longer than the
                     *  local scale is useless.  */
                    const float sx2 = s.bx - s.ax,
                                sy2 = s.by - s.ay,
                                sz2 = s.bz - s.az;
                    if (sx2*sx2 + sy2*sy2 + sz2*sz2 >
                        25.f * wclamp[i] * wclamp[i])
                        continue;
                    for (int u = 0; u < DS; ++u)
                    {
                        const float tt = float(u) / (DS - 1);
                        px.push_back(s.ax + tt * sx2);
                        py.push_back(s.ay + tt * sy2);
                        pz.push_back(s.az + tt * sz2);
                    }
                    segs.push_back(s);
                }
            }
            if (!px.empty())
            {
                eval_points_mt(deck, ctx, px, py, pz, pv);
                std::vector<float> bestd(cand.size(), 1e30f);
                for (const Seg& s : segs)
                {
                    const size_t i = cand[s.q];
                    for (int u = 0; u + 1 < DS; ++u)
                    {
                        const float f0 = pv[s.o + u];
                        const float f1 = pv[s.o + u + 1];
                        if (!std::isfinite(f0) ||
                            !std::isfinite(f1) ||
                            (f0 < 0) == (f1 < 0))
                            continue;
                        const float tt =
                                (float(u) + f0 / (f0 - f1)) /
                                (DS - 1);
                        const float cx2 =
                                s.ax + tt * (s.bx - s.ax);
                        const float cy2 =
                                s.ay + tt * (s.by - s.ay);
                        const float cz2 =
                                s.az + tt * (s.bz - s.az);
                        const float dx = cx2 - wx2[i],
                                    dy = cy2 - wy2[i],
                                    dz = cz2 - wz2[i];
                        const float d2m = dx*dx + dy*dy + dz*dz;
                        if (d2m < bestd[s.q] &&
                            d2m < 4.f * wclamp[i] * wclamp[i])
                        {
                            bestd[s.q] = d2m;
                            sdbx[s.q] = cx2;
                            sdby[s.q] = cy2;
                            sdbz[s.q] = cz2;
                            sdb_ok[s.q] = 1;
                        }
                    }
                }
            }
        }

        std::vector<float> kx, ky, kz, kh, kc;
        std::vector<uint8_t> kspec;
        if (repair_mode >= 2)
        {
            /*  Warts (chords through EMPTY space) crease-seek: a
             *  midpoint split of a crease-crossing chord is
             *  self-similar and never converges; split at the |f|
             *  peak along the chord instead - the branch switch,
             *  i.e. the crease.  Chips (chords through MATERIAL)
             *  do NOT seek the crease: their |f| peak IS the
             *  corner, which the constraint already owns - a
             *  crease-seek target there is exactly what the
             *  keep-out must block (measured: 87 of csg's 90
             *  candidates blocked, chips shipped).  Chips take
             *  their midpoint, which Newton-projects OUTWARD onto
             *  the nearest face, safely off the polyline.  */
            constexpr int NSAMP = 15;
            std::vector<float> sx(cand.size() * NSAMP),
                    sy(cand.size() * NSAMP),
                    sz(cand.size() * NSAMP), sv;
            for (size_t q = 0; q < cand.size(); ++q)
            {
                const size_t i = cand[q];
                for (int u = 0; u < NSAMP; ++u)
                {
                    const float t = float(u + 1) / (NSAMP + 1);
                    sx[q*NSAMP + u] = eax[i] + t * (ebx[i] - eax[i]);
                    sy[q*NSAMP + u] = eay[i] + t * (eby[i] - eay[i]);
                    sz[q*NSAMP + u] = eaz[i] + t * (ebz[i] - eaz[i]);
                }
            }
            eval_points_mt(deck, ctx, sx, sy, sz, sv);
            for (size_t q = 0; q < cand.size(); ++q)
            {
                const size_t i = cand[q];
                if (cspec[q])
                {
                    if (sdb_ok[q])
                    {
                        kx.push_back(sdbx[q]);
                        ky.push_back(sdby[q]);
                        kz.push_back(sdbz[q]);
                    }
                    else
                    {
                        kx.push_back(wx2[i]);
                        ky.push_back(wy2[i]);
                        kz.push_back(wz2[i]);
                    }
                }
                else
                {
                    int best = 0;
                    float bf = -1;
                    for (int u = 0; u < NSAMP; ++u)
                    {
                        const float av = fabsf(sv[q*NSAMP + u]);
                        if (std::isfinite(av) && av > bf)
                        {
                            bf = av;
                            best = u;
                        }
                    }
                    kx.push_back(sx[q*NSAMP + best]);
                    ky.push_back(sy[q*NSAMP + best]);
                    kz.push_back(sz[q*NSAMP + best]);
                }
                kh.push_back(whs[i]);
                kc.push_back(wclamp[i]);
                kspec.push_back(cspec[q] ? (sdb_ok[q] ? 2 : 1) : 0);
            }
        }
        else
        {
            for (size_t q = 0; q < cand.size(); ++q)
            {
                const size_t i = cand[q];
                kx.push_back(wx2[i]);
                ky.push_back(wy2[i]);
                kz.push_back(wz2[i]);
                kh.push_back(whs[i]);
                kc.push_back(wclamp[i]);
                kspec.push_back(cspec[q] ? (sdb_ok[q] ? 2 : 1) : 0);
            }
        }
        project_points_impl(deck, ctx, kx, ky, kz, kh, kc);
        /*  Post-projection oracle (plinth tooth autopsy,
         *  2026-07-18): Newton projection OSCILLATES at concave
         *  creases (the gradient flips across the kink) and a
         *  stalled projection parks in the corner's AIR wedge -
         *  measured on the plinth junction: repair vertices at
         *  quarter-cell heights 0.1-0.7 mm outside the wall,
         *  minting the tilted "foot teeth" (88 at r1, gone at
         *  r2 - they scale with the lattice, not the model).
         *  Same question stage-A's phantom oracle asks of QEF
         *  minimizers, one stage downstream: is the projected
         *  point ON the surface?  A miss is dropped, not
         *  inserted - the chip it wanted to press is bounded
         *  (3%-of-edge species) while an air vertex is
         *  unbounded damage.  */
        std::vector<uint8_t> offsurf(kx.size(), 0);
        if (!kx.empty())
        {
            const float vh2 = 0.01f * soup.spacing;
            std::vector<float> vx(kx.size() * 7), vy(kx.size() * 7),
                    vz(kx.size() * 7), vv;
            for (size_t i = 0; i < kx.size(); ++i)
            {
                for (int q = 0; q < 7; ++q)
                {
                    vx[i*7 + q] = kx[i];
                    vy[i*7 + q] = ky[i];
                    vz[i*7 + q] = kz[i];
                }
                vx[i*7 + 1] += vh2;   vx[i*7 + 2] -= vh2;
                vy[i*7 + 3] += vh2;   vy[i*7 + 4] -= vh2;
                vz[i*7 + 5] += vh2;   vz[i*7 + 6] -= vh2;
            }
            eval_points_mt(deck, ctx, vx, vy, vz, vv);
            for (size_t i = 0; i < kx.size(); ++i)
            {
                const float f = vv[i*7];
                const float gx = vv[i*7 + 1] - vv[i*7 + 2];
                const float gy = vv[i*7 + 3] - vv[i*7 + 4];
                const float gz = vv[i*7 + 5] - vv[i*7 + 6];
                const float gl = sqrtf(gx*gx + gy*gy + gz*gz) /
                        (2 * vh2);
                if (!std::isfinite(f) ||
                    (gl > 1e-12f &&
                     fabsf(f) / gl > 0.05f * soup.spacing))
                    offsurf[i] = 1;
            }
        }
        uint64_t added = 0;
        uint64_t blk_keepout = 0, blk_crowd = 0, blk_air = 0;
        for (size_t i = 0; i < kx.size(); ++i)
        {
            if (offsurf[i])
            {
                ++blk_air;
                continue;
            }
            /*  Repair keep-out: the constrained crease repairs
             *  itself; only the smooth field is ours to press (an
             *  insert on a constrained edge destroys it and the
             *  re-conforming Steiner point lands a sliver away -
             *  the pinch factory, measured 2026-07-15).  Chip
             *  repairs land ON a face, not the crease, so they
             *  only need the sliver radius - 0.75 would starve
             *  them entirely (measured).  */
            /*  Chip keep-out swept 2026-07-15: 0.75 starves chips
             *  entirely (87 blocked, 0 fixed); 0.15 escalates
             *  (near-crease inserts spawn new sliver chords faster
             *  than they fix - 17 -> 31 over rounds, pinch splits
             *  appear); 0.35 presses monotonically 87 -> 17 with
             *  zero topology damage.  A 0.10 shield + 0.3 sp
             *  crease-vertex-only crowding floor was also tried:
             *  36 more repairs land but detection PLATEAUS at ~20
             *  (fresh shallow chips replace pressed ones - the
             *  relative tolerance is self-similar and the guards
             *  rightly stop the descent) - reverted, 0.35 ends
             *  better.  The residue is the greedy-repair
             *  asymptote; the structural cure is extraction-level
             *  (MESH-NEXT).  */
            /*  Species 2 (SDB centers) skip the segment keep-out:
             *  they are principled insert locations, gated by the
             *  crowding guard and refereed by the counters.  */
            /*  Density round 2026-07-15: the sliver radius the
             *  keep-out protects scales with the LOCAL pitch, not
             *  the base pitch - chips live in the crease band,
             *  where the dense lattice divides it.  Keeping the
             *  RATIO at the measured-good 0.35 x pitch preserves
             *  the sweep's verdict at the band's scale (fixed
             *  0.35 sp starves the band entirely: 281/281 blocked,
             *  depth stuck at 0.375 sp; 0.6 x chord-length proxy
             *  ditto - chords ARE one local pitch, so it resolved
             *  to 0.30 sp and still blocked everything).  */
            /*  STIBIUM_DMESH_KEEPOUT scales the repair keep-out.
             *  DEFAULT 0 (off) since the eyepiece-ring autopsy
             *  (2026-07-17): the keep-out was starving groove
             *  repairs beside every constrained crease (the ring
             *  "thorn crown"; zeiss worst depth 0.573 -> 0.212 sp
             *  and nm 168 -> 124 with it off, bench better
             *  everywhere) while its anti-sliver duty is already
             *  covered by the crowding guard below.  =1 restores
             *  the old radii.  */
            static const char* ko_env =
                    getenv("STIBIUM_DMESH_KEEPOUT");
            const float ko_scale = ko_env ? float(atof(ko_env))
                                          : 0.f;
            const float ko = ko_scale * (kspec[i] == 2 ? 0.0f
                    : kspec[i]
                        ? 0.35f * soup.spacing /
                          std::max(dense_factor(),
                                   box_dense_factor(soup, kx[i],
                                                    ky[i], kz[i]))
                        : 0.75f * soup.spacing);
            if (ko > 0 && !cseg.empty() &&
                near_crease(kx[i], ky[i], kz[i], ko))
            {
                ++blk_keepout;
                if (chip_env && repair_round >= 3)
                    fprintf(stderr, "BLOCKED %s target "
                            "(%.4f, %.4f, %.4f) edge %.3f sp\n",
                            kspec[i] ? "chip" : "wart",
                            kx[i], ky[i], kz[i],
                            kc[i] / soup.spacing);
                continue;
            }
            const TPoint pt(kx[i], ky[i], kz[i]);
            if (repair_mode >= 2)
            {
                /*  Crowding guard: an insert within a quarter of
                 *  its edge of an existing vertex makes the sliver
                 *  pinches - skip it (the pinch preventer).  */
                const auto nv = dt.nearest_vertex(pt);
                if (nv != VH())
                {
                    const double d2 = CGAL::squared_distance(
                            nv->point(), pt);
                    /*  Law-blindness dial (plinth autopsy,
                     *  2026-07-18): the guard counted constraint
                     *  and Steiner vertices as crowding and
                     *  vetoed every repair beside the junction it
                     *  had just been given law for (blocked
                     *  26->61 across rounds while detected chips
                     *  ROSE).  When the nearest vertex is LAW,
                     *  its crowding radius scales by
                     *  STIBIUM_DMESH_CROWD_LAW; default 1 =
                     *  unchanged (the fossil at the keep-out
                     *  comment warns a naive exemption plateaued
                     *  in the sparse-law era - re-refereed now on
                     *  the plinth before any default moves).  */
                    static const char* cl_env =
                            getenv("STIBIUM_DMESH_CROWD_LAW");
                    const double cls = cl_env ? atof(cl_env) : 1.0;
                    double r = 0.25 * kc[i];
                    if (cls != 1.0)
                    {
                        const auto pit = prov.find(&*nv);
                        bool law = pit != prov.end() &&
                                pit->second == 1;
                        if constexpr (CCDT_MODE)
                            if (!law && nv->ccdt_3_data()
                                        .is_Steiner_vertex_on_edge())
                                law = true;
                        if (law)
                            r *= cls;
                    }
                    if (d2 < r * r)
                    {
                        ++blk_crowd;
                        continue;
                    }
                }
            }
            const auto before = dt.number_of_vertices();
            auto vh = dt.insert(pt);
            if (dt.number_of_vertices() > before)
            {
                vh->info() = 0;
                prov.emplace(&*vh, 5);
                ++added;
            }
        }
        if (chip_env)
            fprintf(stderr, "REPAIR round %d: inserted %llu, "
                    "keep-out blocked %llu, crowding blocked "
                    "%llu, off-surface dropped %llu\n",
                    repair_round,
                    (unsigned long long)added,
                    (unsigned long long)blk_keepout,
                    (unsigned long long)blk_crowd,
                    (unsigned long long)blk_air);
        {
            char nb[32];
            snprintf(nb, sizeof(nb), "r%d detect+insert",
                     repair_round);
            pt.sub(nb);
        }
        if (!added)
            break;
        repaired_total += added;
        sweep_steiner();
    }
    }   // repair loop
    pt.mark("extract+repair");
    prog_stage(10);
    /*  Stage dumps (STIBIUM_DMESH_STAGES=prefix): the formation
     *  film - drop each stage into a slicer and watch where a
     *  defect is born instead of inferring it from deltas.  */
    static const char* stages_env = getenv("STIBIUM_DMESH_STAGES");
    if (stages_env)
        write_stl_raw((std::string(stages_env) +
                       "_1_extract_repair.stl").c_str(),
                      out->verts, out->tris);

    /*  Manifold pass (Manifold-DC after Schaefer/Ju/Warren 2007;
     *  see doc/research/2026-07-15-pinch-manifoldness.md): our
     *  extraction is restricted-Delaunay class, so a separator or
     *  tangency can pinch the surface into a 4-triangle edge.  The
     *  tet complex knows the truth: circulating the cells around a
     *  pinch edge, the two facets bounding the same INSIDE run are
     *  one wedge of the solid - one sheet.  Pair them, then split
     *  every vertex whose facet fan is disconnected into one output
     *  vertex per sheet (coincident copies, zero geometric change).
     *  The result is 2-manifold at every vertex by construction.
     *  STIBIUM_DMESH_MANIFOLD=0 disables.  */
    {
        static const char* man_env = getenv("STIBIUM_DMESH_MANIFOLD");
        const bool manifold_on = man_env ? atoi(man_env) != 0 : true;
        const size_t ntri = out->tris.size() / 3;
        if (manifold_on && ntri)
        {
            /*  Edge -> incident triangle list  */
            std::unordered_map<uint64_t, std::vector<uint32_t>> etris;
            for (uint32_t t = 0; t < ntri; ++t)
                for (int e = 0; e < 3; ++e)
                {
                    uint64_t a = out->tris[3*t + e];
                    uint64_t b = out->tris[3*t + (e + 1) % 3];
                    if (a > b)
                        std::swap(a, b);
                    etris[(a << 32) | b].push_back(t);
                }

            /*  Triangle lookup by sorted vertex triple.  EXACT
             *  key (correctness review finding 2): the old
             *  21-bit-per-index packing silently aliased above
             *  2,097,151 vertices - multi-million-tri exports
             *  already ship.  Struct key with equality resolves
             *  hash collisions correctly at any size.  */
            struct TriKeyHash
            {
                size_t operator()(
                        const std::array<uint32_t, 3>& k) const
                {
                    uint64_t h = 0xcbf29ce484222325ull;
                    for (const uint32_t w : k)
                        h = (h ^ w) * 0x100000001b3ull;
                    return size_t(h);
                }
            };
            std::unordered_map<std::array<uint32_t, 3>, uint32_t,
                               TriKeyHash> trikey;
            const auto tri_key = [](uint32_t x, uint32_t y,
                                    uint32_t z) {
                if (x > y) std::swap(x, y);
                if (y > z) std::swap(y, z);
                if (x > y) std::swap(x, y);
                return std::array<uint32_t, 3>{ x, y, z };
            };
            for (uint32_t t = 0; t < ntri; ++t)
                trikey[tri_key(out->tris[3*t], out->tris[3*t + 1],
                               out->tris[3*t + 2])] = t;

            /*  Reverse map: vertex handle -> output index  */
            std::unordered_map<const void*, uint32_t> vrev;
            for (uint32_t i = 0; i < xvh.size(); ++i)
                vrev[&*xvh[i]] = i;

            /*  Sheet pairing at pinch edges: edge key ->
             *  (tri -> sheet id).  Unrecoverable rings leave the
             *  edge unpaired = conservative split.  */
            std::unordered_map<uint64_t,
                    std::unordered_map<uint32_t, int>> sheet;
            for (const auto& [k, ts] : etris)
            {
                if (ts.size() <= 2)
                    continue;
                const uint32_t ia = uint32_t(k >> 32);
                const uint32_t ib = uint32_t(k);
                const VH vp = xvh[ia], vq = xvh[ib];
                typename Tri::Cell_handle ec;
                int ei, ej;
                if (!dt.is_edge(vp, vq, ec, ei, ej))
                    continue;
                /*  Ring of cells around the edge, with signs  */
                auto circ = dt.incident_cells(ec, ei, ej);
                const auto done = circ;
                std::vector<CH> ring;
                do
                    ring.push_back(CH(circ));
                while (++circ != done);
                const auto cell_sign = [&](CH c) -> int {
                    return dt.is_infinite(c) ? 1 : int(c->info());
                };
                /*  Transitions between consecutive ring cells;
                 *  each is an extracted facet {vp, vq, w}.  */
                struct Trans { uint32_t tri; int before; };
                std::vector<Trans> trans;
                bool ok = true;
                for (size_t r = 0; r < ring.size() && ok; ++r)
                {
                    const CH c1 = ring[r];
                    const CH c2 = ring[(r + 1) % ring.size()];
                    const int s1 = cell_sign(c1);
                    const int s2 = cell_sign(c2);
                    if (s1 == s2 || s1 == 0 || s2 == 0)
                        continue;
                    /*  Shared facet's third vertex  */
                    int ni = -1;
                    for (int q = 0; q < 4; ++q)
                        if (c1->neighbor(q) == c2)
                            ni = q;
                    if (ni < 0) { ok = false; break; }
                    VH w = VH();
                    for (int q = 0; q < 4; ++q)
                    {
                        if (q == ni)
                            continue;
                        VH v = c1->vertex(q);
                        if (v != vp && v != vq)
                            w = v;
                    }
                    const auto wf = w == VH() ? vrev.end()
                                              : vrev.find(&*w);
                    if (wf == vrev.end()) { ok = false; break; }
                    const auto tf = trikey.find(
                            tri_key(ia, ib, wf->second));
                    if (tf == trikey.end()) { ok = false; break; }
                    trans.push_back({ tf->second, s1 });
                }
                if (!ok || trans.size() != ts.size())
                    continue;
                /*  Pair each out->in transition with the next
                 *  transition (the inside run between them is one
                 *  wedge of the solid).  */
                auto& sh = sheet[k];
                int sid = 0;
                for (size_t r = 0; r < trans.size(); ++r)
                    if (trans[r].before > 0)
                    {
                        sh[trans[r].tri] = sid;
                        sh[trans[(r + 1) % trans.size()].tri] = sid;
                        ++sid;
                    }
            }

            /*  Geometric fallback pairing (prov-5 pinch round,
             *  2026-07-17: 62% of residual nm vertices were
             *  repair-loop inserts whose tet-ring walk failed,
             *  landing in the conservative keep-the-pinch path).
             *  Two sheets crossing an edge are each locally
             *  planar, so greedy max-dot normal matching recovers
             *  the pairing with no DT walk at all.  The map is
             *  shared by both endpoint fans, so duplication stays
             *  consistent and cannot tear a hole.  */
            static const char* fb_env =
                    getenv("STIBIUM_DMESH_PAIR_FALLBACK");
            const bool fb_on = fb_env && atoi(fb_env) != 0;
            for (const auto& [k, ts] : etris)
            {
                if (!fb_on)
                    break;
                if (ts.size() <= 2 || (ts.size() & 1))
                    continue;
                /*  The walk loop default-constructs sheet[k] even
                 *  when its pairing produced nothing - absent OR
                 *  empty means unpaired.  */
                const auto sf0 = sheet.find(k);
                if (sf0 != sheet.end() && !sf0->second.empty())
                    continue;
                std::vector<std::array<float, 3>> nrm(ts.size());
                for (size_t q = 0; q < ts.size(); ++q)
                {
                    const uint32_t t = ts[q];
                    const float* A = &out->verts[3 *
                            out->tris[3*t]];
                    const float* B = &out->verts[3 *
                            out->tris[3*t + 1]];
                    const float* C = &out->verts[3 *
                            out->tris[3*t + 2]];
                    const float ux = B[0]-A[0], uy = B[1]-A[1],
                                uz = B[2]-A[2];
                    const float vx = C[0]-A[0], vy = C[1]-A[1],
                                vz = C[2]-A[2];
                    float nx = uy*vz - uz*vy;
                    float ny = uz*vx - ux*vz;
                    float nz = ux*vy - uy*vx;
                    const float l = sqrtf(nx*nx + ny*ny + nz*nz);
                    if (l > 0)
                    {
                        nx /= l;
                        ny /= l;
                        nz /= l;
                    }
                    nrm[q] = { nx, ny, nz };
                }
                std::vector<char> used(ts.size(), 0);
                auto& sh = sheet[k];
                int sid = 0;
                bool bad = false;
                for (size_t q = 0; q < ts.size() && !bad; ++q)
                {
                    if (used[q])
                        continue;
                    used[q] = 1;
                    float bestd = -2.f;
                    size_t bj = SIZE_MAX;
                    for (size_t j = q + 1; j < ts.size(); ++j)
                    {
                        if (used[j])
                            continue;
                        const float d =
                                nrm[q][0]*nrm[j][0] +
                                nrm[q][1]*nrm[j][1] +
                                nrm[q][2]*nrm[j][2];
                        if (d > bestd)
                        {
                            bestd = d;
                            bj = j;
                        }
                    }
                    if (bj == SIZE_MAX)
                    {
                        bad = true;
                        break;
                    }
                    used[bj] = 1;
                    sh[ts[q]] = sid;
                    sh[ts[bj]] = sid;
                    ++sid;
                }
                if (bad)
                    sheet.erase(k);
            }
            if (getenv("STIBIUM_DMESH_NM_DEBUG"))
            {
                size_t nm_e = 0, walk = 0, fall = 0, odd = 0;
                for (const auto& [k, ts] : etris)
                {
                    if (ts.size() <= 2)
                        continue;
                    ++nm_e;
                    if (ts.size() & 1)
                        ++odd;
                    const auto sf0 = sheet.find(k);
                    if (sf0 != sheet.end() &&
                        !sf0->second.empty())
                        ++walk;
                }
                (void)fall;
                fprintf(stderr, "MANIFOLD: %zu pinch edges, %zu "
                        "paired, %zu odd-count\n", nm_e, walk,
                        odd);
            }

            /*  Fan split: per vertex, union incident triangles
             *  across manifold edges (and across paired pinch
             *  edges); each extra component gets a coincident
             *  duplicate vertex.  Adjacency is read from a
             *  SNAPSHOT - rewrites must not perturb later
             *  vertices' edge keys.  */
            const std::vector<uint32_t> tris0 = out->tris;
            std::vector<std::vector<uint32_t>> vtris(
                    out->verts.size() / 3);
            for (uint32_t t = 0; t < ntri; ++t)
                for (int e = 0; e < 3; ++e)
                    vtris[tris0[3*t + e]].push_back(t);

            uint64_t split = 0;
            for (uint32_t v = 0; v < vtris.size(); ++v)
            {
                const auto& ts = vtris[v];
                if (ts.size() < 2)
                    continue;
                /*  Local union-find  */
                std::unordered_map<uint32_t, uint32_t> id;
                std::vector<uint32_t> parent(ts.size());
                for (uint32_t q = 0; q < ts.size(); ++q)
                {
                    id[ts[q]] = q;
                    parent[q] = q;
                }
                std::function<uint32_t(uint32_t)> find =
                        [&](uint32_t x) -> uint32_t {
                    while (parent[x] != x)
                        x = parent[x] = parent[parent[x]];
                    return x;
                };
                for (const uint32_t t : ts)
                    for (int e = 0; e < 3; ++e)
                    {
                        uint64_t a = tris0[3*t + e];
                        uint64_t b = tris0[3*t + (e + 1) % 3];
                        if (a != v && b != v)
                            continue;
                        if (a > b)
                            std::swap(a, b);
                        const uint64_t ek = (a << 32) | b;
                        const auto& ets = etris[ek];
                        if (ets.size() == 2)
                        {
                            parent[find(id[ets[0]])] =
                                    find(id[ets[1]]);
                        }
                        else if (ets.size() > 2)
                        {
                            const auto sf = sheet.find(ek);
                            if (sf == sheet.end())
                            {
                                /*  Pairing unrecoverable: keep the
                                 *  pinch (union everything) - an
                                 *  honest non-manifold edge beats a
                                 *  torn hole.  */
                                for (size_t q = 1; q < ets.size();
                                     ++q)
                                    parent[find(id[ets[q]])] =
                                            find(id[ets[0]]);
                                continue;
                            }
                            /*  Union same-sheet pairs; any tri the
                             *  pairing missed joins sheet 0.  */
                            std::unordered_map<int, uint32_t> first;
                            for (const uint32_t et : ets)
                            {
                                const auto pf = sf->second.find(et);
                                const int sid =
                                        pf == sf->second.end()
                                        ? 0 : pf->second;
                                const auto ff = first.find(sid);
                                if (ff == first.end())
                                    first[sid] = et;
                                else
                                    parent[find(id[et])] =
                                            find(id[ff->second]);
                            }
                        }
                    }
                /*  Components -> duplicates  */
                std::unordered_map<uint32_t, uint32_t> comp_vert;
                bool any = false;
                for (uint32_t q = 0; q < ts.size(); ++q)
                {
                    const uint32_t root = find(q);
                    auto cf = comp_vert.find(root);
                    if (cf == comp_vert.end())
                    {
                        const uint32_t nv = comp_vert.empty()
                                ? v
                                : uint32_t(out->verts.size() / 3);
                        if (nv != v)
                        {
                            out->verts.push_back(
                                    out->verts[3*v]);
                            out->verts.push_back(
                                    out->verts[3*v + 1]);
                            out->verts.push_back(
                                    out->verts[3*v + 2]);
                            ++split;
                            any = true;
                        }
                        comp_vert[root] = nv;
                        cf = comp_vert.find(root);
                    }
                    if (cf->second != v)
                    {
                        const uint32_t t = ts[q];
                        for (int e = 0; e < 3; ++e)
                            if (out->tris[3*t + e] == v)
                                out->tris[3*t + e] = cf->second;
                    }
                }
                (void)any;
            }
            out->split_verts = split;
        }
    }

    pt.mark("manifold pass");
    prog_stage(11);
    if (stages_env)
        write_stl_raw((std::string(stages_env) +
                       "_2_manifold.stl").c_str(),
                      out->verts, out->tris);
    /*  Crease-snap pass (STIBIUM_DMESH_SNAP=0 disables): the
     *  residual chips are chords of MOSTLY-AIR tets whose bottom
     *  facet dips under a concave crease - a facet-level defect no
     *  cell classification or insert-repair can reach (the chip
     *  target and the repair keep-out both scale with local pitch,
     *  so the blocked residue is density-invariant: plateau
     *  0.196/0.198/0.192 sp at 1x/2x/4x band density, both repair
     *  strategies).  But the traced polylines know EXACTLY where
     *  the surface corner is: split each residual chip edge's two
     *  triangles at the nearest traced-crease point, tenting the
     *  chord up onto the crease.  Pure output surgery - no CGAL
     *  insert, no Steiner, no keep-out, manifoldness preserved by
     *  construction (the tent's new edges each get exactly two
     *  triangles, the apex fan is closed).  */
    {
        /*  FRESH reads (correctness r2 #1): the app dialog
         *  setenvs these per export; a static here silently
         *  reuses the FIRST export's values.  */
        const char* snap_env = getenv("STIBIUM_DMESH_SNAP");
        const bool snap_on = !snap_env || atoi(snap_env) != 0;
        uint64_t snapped = 0, snap_skipped = 0;
        uint64_t skip_far = 0, skip_hits = 0, snap_surf = 0;
        uint64_t skip_damage = 0, skip_churn = 0;
        uint64_t snap_valley = 0, churn_saved = 0;
        /*  Instrument: per-chip outcome + depth (sp units) to a
         *  file - the skip-class anatomy (which population is
         *  actually deep) drives the valley-snap calibration.  */
        FILE* sdump = nullptr;
        if (const char* sd_env = getenv("STIBIUM_DMESH_SNAP_DUMP"))
            sdump = fopen(sd_env, "w");
        /*  Self-feeding stash (zeiss run, 2026-07-15): when the
         *  repair loop exhausts MAX_REPAIR (real models churn), the
         *  per-round stash is stale and was cleared - but the chips
         *  are still there and the polylines still know where the
         *  creases are.  Sweep the FINAL mesh once and build the
         *  snap list fresh; the wave lookup re-finds triangles by
         *  position anyway.  */
        static const char* rep_env2 = getenv("STIBIUM_DMESH_REPAIR");
        const int rmode2 = rep_env2 ? atoi(rep_env2) : 2;
        /*  Runs as the classic stash-sweep fallback when the
         *  per-round list is empty, AND in warts-only mode when
         *  contact chains exist (the per-round list records
         *  chips; air-chords over contact seams are warts and
         *  need their own sweep).  */
        const bool warts_only = !snap_va.empty();
        uint64_t wart_admitted = 0;
        if (snap_on &&
            (snap_va.empty() || !soup.contact_chains.empty()) &&
            (!soup.tchains.empty() ||
             !soup.contact_chains.empty()) &&
            rmode2 >= 2 && !*halt)
        {
            std::unordered_map<uint64_t,
                    std::pair<uint32_t, uint32_t>> fem;
            for (uint32_t t = 0;
                 t < uint32_t(out->tris.size() / 3); ++t)
                for (int e = 0; e < 3; ++e)
                {
                    uint64_t a = out->tris[3*t + e];
                    uint64_t b = out->tris[3*t + (e + 1) % 3];
                    if (a > b)
                        std::swap(a, b);
                    const uint64_t k = (a << 32) | b;
                    auto it = fem.find(k);
                    if (it == fem.end())
                        fem[k] = { t, UINT32_MAX };
                    else if (it->second.second == UINT32_MAX)
                        it->second.second = t;
                }
            std::vector<uint64_t> keys;
            std::vector<float> qxs, qys, qzs, qhs, qlen;
            for (const auto& [k, pr] : fem)
            {
                if (pr.second == UINT32_MAX)
                    continue;
                const float* A2 =
                        &out->verts[3 * uint32_t(k >> 32)];
                const float* B2 = &out->verts[3 * uint32_t(k)];
                const float ex = B2[0]-A2[0], ey = B2[1]-A2[1],
                            ez = B2[2]-A2[2];
                const float len = sqrtf(ex*ex + ey*ey + ez*ez);
                if (!(len > 0))
                    continue;
                keys.push_back(k);
                qxs.push_back((A2[0]+B2[0]) * 0.5f);
                qys.push_back((A2[1]+B2[1]) * 0.5f);
                qzs.push_back((A2[2]+B2[2]) * 0.5f);
                qhs.push_back(len * 0.01f);
                qlen.push_back(len);
            }
            std::vector<float> gxs(keys.size()*7),
                    gys(keys.size()*7), gzs(keys.size()*7), gv;
            for (size_t i = 0; i < keys.size(); ++i)
            {
                for (int q = 0; q < 7; ++q)
                {
                    gxs[i*7+q] = qxs[i];
                    gys[i*7+q] = qys[i];
                    gzs[i*7+q] = qzs[i];
                }
                gxs[i*7+1] += qhs[i];  gxs[i*7+2] -= qhs[i];
                gys[i*7+3] += qhs[i];  gys[i*7+4] -= qhs[i];
                gzs[i*7+5] += qhs[i];  gzs[i*7+6] -= qhs[i];
            }
            if (!keys.empty())
                eval_points_mt(deck, ctx, gxs, gys, gzs,
                            gv);
            /*  #4a wart admission: air-chords at contact loci
             *  read f >= 0 at their midpoints (chord spans AIR
             *  over the tangent seam) and were always filtered
             *  out of snap - repair's territory, except repair
             *  has no law there either.  Admit a WART edge only
             *  when a contact chain runs within one chord length
             *  of its midpoint; everywhere else the chips-only
             *  gate stands.  */
            std::vector<std::array<float, 6>> csegs;
            for (const auto& ch : soup.contact_chains)
                for (size_t i = 0; i + 1 < ch.size(); ++i)
                    csegs.push_back({ ch[i][0], ch[i][1],
                                      ch[i][2], ch[i+1][0],
                                      ch[i+1][1], ch[i+1][2] });
            const auto near_contact = [&](float mx, float my,
                                          float mz, float rr) {
                const float r2 = rr * rr;
                for (const auto& s : csegs)
                {
                    const float bx2 = s[3]-s[0], by2 = s[4]-s[1],
                                bz2 = s[5]-s[2];
                    const float bb = bx2*bx2 + by2*by2 + bz2*bz2;
                    float t2 = bb > 0
                            ? ((mx-s[0])*bx2 + (my-s[1])*by2 +
                               (mz-s[2])*bz2) / bb
                            : 0;
                    t2 = t2 < 0 ? 0 : t2 > 1 ? 1 : t2;
                    const float dx = mx - (s[0] + t2*bx2),
                                dy = my - (s[1] + t2*by2),
                                dz = mz - (s[2] + t2*bz2);
                    if (dx*dx + dy*dy + dz*dz <= r2)
                        return true;
                }
                return false;
            };
            for (size_t i = 0; i < keys.size(); ++i)
            {
                const float f = gv[i*7];
                const float h2 = 2 * qhs[i];
                const float gx = (gv[i*7+1]-gv[i*7+2]) / h2;
                const float gy = (gv[i*7+3]-gv[i*7+4]) / h2;
                const float gz = (gv[i*7+5]-gv[i*7+6]) / h2;
                const float gl = sqrtf(gx*gx + gy*gy + gz*gz);
                if (!(gl > 0) || !std::isfinite(f))
                    continue;
                if (f < 0 && warts_only)
                    continue;   /* chips already stashed per-round */
                if (f >= 0 &&
                    (csegs.empty() ||
                     !near_contact(qxs[i], qys[i], qzs[i],
                                   qlen[i])))
                    continue;
                const float dist = fabsf(f) / gl;
                if (dist <= qlen[i] * 0.03f ||
                    dist > 2.0f * qlen[i])
                    continue;
                snap_va.push_back(uint32_t(keys[i] >> 32));
                snap_vb.push_back(uint32_t(keys[i]));
                snap_t1.push_back(fem[keys[i]].first);
                snap_t2.push_back(fem[keys[i]].second);
                snap_d.push_back(dist);
                if (f >= 0)
                    ++wart_admitted;
            }
            if (getenv("STIBIUM_DMESH_CHIP_DEBUG"))
                fprintf(stderr, "SWEEP: %zu edges scanned, %zu "
                        "contact segs, %llu warts admitted "
                        "(warts_only=%d)\n",
                        keys.size(), csegs.size(),
                        (unsigned long long)wart_admitted,
                        int(warts_only));
            pt.sub("snap: stash sweep");
        }
        if (snap_on && !snap_va.empty() &&
            (!soup.tchains.empty() ||
             !soup.contact_chains.empty()))
        {
            struct Seg { float ax, ay, az, bx, by, bz; };
            std::vector<Seg> segs;
            for (size_t c = 0; c < soup.tchains.size(); ++c)
            {
                const auto& ch = soup.tchains[c];
                for (size_t i = 0; i + 1 < ch.size(); ++i)
                {
                    const DSurfPoint& a = soup.surface[ch[i]];
                    const DSurfPoint& b = soup.surface[ch[i + 1]];
                    segs.push_back({ a.x, a.y, a.z, b.x, b.y, b.z });
                }
                if (soup.tclosed[c] && ch.size() >= 3)
                {
                    const DSurfPoint& a = soup.surface[ch.back()];
                    const DSurfPoint& b = soup.surface[ch.front()];
                    segs.push_back({ a.x, a.y, a.z, b.x, b.y, b.z });
                }
            }
            /*  #4a: contact chains join the snap law - same
             *  attribution cap, same referees, no constraints.  */
            for (const auto& ch : soup.contact_chains)
                for (size_t i = 0; i + 1 < ch.size(); ++i)
                    segs.push_back({ ch[i][0], ch[i][1],
                                     ch[i][2], ch[i+1][0],
                                     ch[i+1][1], ch[i+1][2] });
            const auto pos_eq = [&](uint32_t v, const float* p) {
                return out->verts[3*v] == p[0] &&
                       out->verts[3*v + 1] == p[1] &&
                       out->verts[3*v + 2] == p[2];
            };
            const auto vpos_hash = [&](uint32_t v) {
                return coord_hash(out->verts[3*v],
                                  out->verts[3*v + 1],
                                  out->verts[3*v + 2]);
            };
            /*  Consecutive chips along a crease share triangles, so
             *  the stashed triangle indices go stale as soon as a
             *  neighbour is snapped.  Snap in WAVES: each wave
             *  re-finds every pending edge by POSITION (stable
             *  under both the manifold pass's coincident-duplicate
             *  rewrites and earlier snaps), and skips triangles
             *  already touched this wave.  */
            std::vector<uint8_t> done(snap_va.size(), 0);
            std::vector<std::array<float, 3>> far_pts;
            pt.sub("snap: seg build");
            bool progress = true;
            while (progress)
            {
                prog_frac(snap_va.empty() ? 0.f :
                        float(snapped + snap_skipped) /
                        float(snap_va.size()));
                progress = false;
                std::unordered_map<uint64_t,
                        std::vector<std::pair<uint32_t, int>>> pem;
                for (uint32_t t = 0;
                     t < uint32_t(out->tris.size() / 3); ++t)
                    for (int e = 0; e < 3; ++e)
                    {
                        const uint64_t h =
                                vpos_hash(out->tris[3*t + e]) ^
                                vpos_hash(out->tris[3*t + (e+1)%3]);
                        pem[h].push_back({ t, e });
                    }
                std::unordered_set<uint32_t> touched;
                for (size_t s = 0; s < snap_va.size(); ++s)
                {
                    if (done[s])
                        continue;
                    bool tent_curved = false;
                    const float A2[3] = {
                            out->verts[3*snap_va[s]],
                            out->verts[3*snap_va[s] + 1],
                            out->verts[3*snap_va[s] + 2] };
                    const float B2[3] = {
                            out->verts[3*snap_vb[s]],
                            out->verts[3*snap_vb[s] + 1],
                            out->verts[3*snap_vb[s] + 2] };
                    const float mx = (A2[0] + B2[0]) * 0.5f,
                                my = (A2[1] + B2[1]) * 0.5f,
                                mz = (A2[2] + B2[2]) * 0.5f;
                    /*  Nearest point on the traced polylines  */
                    float best = 1e30f, px = 0, py = 0, pz = 0;
                    for (const Seg& sg : segs)
                    {
                        const float bx2 = sg.bx - sg.ax,
                                    by2 = sg.by - sg.ay,
                                    bz2 = sg.bz - sg.az;
                        const float qx = mx - sg.ax,
                                    qy = my - sg.ay,
                                    qz = mz - sg.az;
                        const float bb = bx2*bx2 + by2*by2 + bz2*bz2;
                        float t2 = bb > 0
                                ? (qx*bx2 + qy*by2 + qz*bz2) / bb
                                : 0;
                        t2 = t2 < 0 ? 0 : t2 > 1 ? 1 : t2;
                        const float cx2 = sg.ax + t2*bx2,
                                    cy2 = sg.ay + t2*by2,
                                    cz2 = sg.az + t2*bz2;
                        const float dx = mx - cx2, dy = my - cy2,
                                    dz = mz - cz2;
                        const float d2 = dx*dx + dy*dy + dz*dz;
                        if (d2 < best)
                        {
                            best = d2;
                            px = cx2; py = cy2; pz = cz2;
                        }
                    }
                    /*  Attribution gate: the tent apex must be
                     *  CLOSE relative to the chip's own depth - a
                     *  wedge corner sits within ~1.6x the |f|/grad
                     *  reading (the central-difference inflation),
                     *  so a polyline further than 2.5x the depth
                     *  belongs to some OTHER feature, and tenting
                     *  to it builds a ridge taller than the divot
                     *  (measured on zeiss: worst depth ROSE 0.975
                     *  -> 1.258 sp under unconditional 1 sp
                     *  attribution).  The absolute 1-cell cap
                     *  stays.  */
                    /*  Floor at the drop-band radius: within
                     *  0.35 sp the crease owns the surface by the
                     *  pipeline's own corridor definition, so
                     *  attribution there is never wrong (csg's
                     *  shallow chips sit at ~0.2 sp with ~0.06 sp
                     *  depth and must still tent).  */
                    /*  ALL radii at the LOCAL pitch (screws fold
                     *  autopsy, 2026-07-18): inside a level-2
                     *  band the global floor builds tents three
                     *  sizes too wide for quarter-cell chips -
                     *  overlapping tents FOLD (4,621 rough edges
                     *  with 180.0-degree dihedrals; snap-off
                     *  reads 2,855 at 161).  The corridor made
                     *  this move on 2026-07-16; snap missed it.  */
                    const float lsp = soup.spacing /
                            std::max(dense_factor(),
                                     box_dense_factor(soup, mx,
                                                      my, mz));
                    const float acap = std::min(lsp,
                            std::max(2.5f * snap_d[s],
                                     0.35f * lsp));
                    if (best > acap * acap)
                    {
                        /*  No polyline owns this chip: a divot in
                         *  SMOOTH surface (untraced tangent blends
                         *  - the zeiss FAR clusters; the seam
                         *  system {f_A = f_B, f = 0} that would
                         *  trace their centerline has a DOUBLE
                         *  root at G1 contact and cannot be
                         *  marched).  The apex that can never
                         *  build a ridge is the surface itself:
                         *  Newton-project the midpoint onto f = 0
                         *  and tent there (STIBIUM_DMESH_SNAP_SURF
                         *  = 0 reverts to skipping).  The same
                         *  depth-proportional acceptance applies -
                         *  a projection that travels further than
                         *  the divot is deep is a bad gradient,
                         *  not a cure.  */
                        static const char* ss_env =
                                getenv("STIBIUM_DMESH_SNAP_SURF");
                        /*  Depth floor: tenting a near-zero divot
                         *  onto the surface is pure churn (the
                         *  surface is already there), and on
                         *  repair-churned bands it self-feeds
                         *  through the stash sweep (engrave_0_5:
                         *  7,795 tents of 0.000-depth noise, +9K
                         *  tris).  0.01 sp splits the measured
                         *  populations: the churn noise reads
                         *  exactly 0.000, the real zeiss FAR
                         *  divots 0.01-0.59 (a 0.05 floor left
                         *  5.2K of them unfixed).  */
                        bool surf_ok = (!ss_env ||
                                        atoi(ss_env) != 0) &&
                                snap_d[s] >= 0.01f * soup.spacing;
                        if (surf_ok)
                        {
                            std::vector<float>
                                    sx{ mx }, sy{ my }, sz{ mz },
                                    hs2{ 0.01f * soup.spacing },
                                    cl2{ 0.75f * soup.spacing };
                            project_points_impl(deck, ctx, sx, sy,
                                                sz, hs2, cl2);
                            const float dx2 = sx[0] - mx,
                                        dy2 = sy[0] - my,
                                        dz2 = sz[0] - mz;
                            const float d2 = dx2*dx2 + dy2*dy2 +
                                             dz2*dz2;
                            surf_ok = std::isfinite(d2) &&
                                      d2 <= acap * acap;
                            if (surf_ok)
                            {
                                px = sx[0];
                                py = sy[0];
                                pz = sz[0];
                                ++snap_surf;
                            }
                        }
                        /*  Surface projection overshot the cap
                         *  (concave-throat gradients under-read -
                         *  the historical 993 FAR skips).  Last
                         *  chance: project onto the CURVATURE
                         *  VALLEY of the blend (quality P1) and
                         *  judge it by the POLYLINE attribution
                         *  cap - a valley line is a feature line,
                         *  same semantics as a traced crease.
                         *  STIBIUM_DMESH_SNAP_VALLEY=0 reverts to
                         *  skipping.  */
                        if (!surf_ok)
                        {
                            static const char* vy_env = getenv(
                                    "STIBIUM_DMESH_SNAP_VALLEY");
                            bool val_ok = (!vy_env ||
                                           atoi(vy_env) != 0) &&
                                    snap_d[s] >=
                                        0.01f * soup.spacing;
                            if (val_ok)
                            {
                                float vx2 = mx, vy2 = my,
                                      vz2 = mz;
                                val_ok = valley_project(deck, ctx,
                                        soup.spacing, lsp,
                                        vx2, vy2, vz2);
                                if (val_ok)
                                {
                                    const float dxv = vx2 - mx,
                                                dyv = vy2 - my,
                                                dzv = vz2 - mz;
                                    val_ok = dxv*dxv + dyv*dyv +
                                             dzv*dzv <=
                                             acap * acap;
                                }
                                if (val_ok)
                                {
                                    px = vx2;
                                    py = vy2;
                                    pz = vz2;
                                    ++snap_valley;
                                }
                            }
                            if (!val_ok)
                            {
                                done[s] = 1;
                                ++snap_skipped;
                                ++skip_far;
                                far_pts.push_back({ mx, my, mz });
                                if (sdump)
                                    fprintf(sdump, "far %.4f %.4f "
                                            "%.4f %.4f\n", mx, my,
                                            mz, snap_d[s] /
                                                soup.spacing);
                                continue;
                            }
                        }
                    }
                    /*  Locate the edge's two triangles by position.
                     *  Coincident sheet duplicates (pinch splits)
                     *  alias the position key - demand EXACTLY two
                     *  real matches or leave the edge alone.  */
                    const uint64_t h = coord_hash(A2[0], A2[1],
                                                  A2[2]) ^
                                       coord_hash(B2[0], B2[1],
                                                  B2[2]);
                    const auto pf = pem.find(h);
                    std::pair<uint32_t, int> hit[2];
                    int nhit = 0;
                    bool stale = false;
                    if (pf != pem.end())
                        for (const auto& [t, e] : pf->second)
                        {
                            /*  An entry whose triangle a snap
                             *  already rewrote this wave no longer
                             *  matches by position - the edge moved
                             *  into a tent triangle pem can't see
                             *  until the next wave.  Defer, don't
                             *  conclude.  */
                            if (touched.count(t))
                            {
                                stale = true;
                                continue;
                            }
                            const uint32_t i0 = out->tris[3*t + e];
                            const uint32_t i1 =
                                    out->tris[3*t + (e+1) % 3];
                            const bool fwd = pos_eq(i0, A2) &&
                                             pos_eq(i1, B2);
                            const bool rev = pos_eq(i0, B2) &&
                                             pos_eq(i1, A2);
                            if (!fwd && !rev)
                                continue;
                            if (nhit < 2)
                                hit[nhit] = { t, e };
                            ++nhit;
                        }
                    if (nhit != 2)
                    {
                        if (stale)
                            continue;     // retry next wave
                        done[s] = 1;
                        ++snap_skipped;
                        ++skip_hits;
                        if (getenv("STIBIUM_DMESH_CHIP_DEBUG"))
                            fprintf(stderr, "  SNAP skip nhit=%d "
                                    "at (%.4f, %.4f, %.4f)\n",
                                    nhit, mx, my, mz);
                        continue;
                    }
                    if (touched.count(hit[0].first) ||
                        touched.count(hit[1].first))
                        continue;         // next wave
                    /*  Tent damage referee (2026-07-17 night, the
                     *  thinning autopsy): a tent whose apex
                     *  attributes to the wrong feature at a step
                     *  corner BUILDS a divot deeper than the chip
                     *  it cures (bino: one tent minted 0.309 sp
                     *  from a sub-0.242 site; the 2.5x attribution
                     *  gate cannot see it - the wrong crease was
                     *  close enough).  Same law as everything that
                     *  survived here: measure the damage.  Probe
                     *  the four new edges' midpoints; if any reads
                     *  a feature deeper than the chip being cured
                     *  (plus a noise floor), refuse the tent - the
                     *  chip stays at its known bounded depth.  */
                    {
                        const uint32_t w0 =
                                out->tris[3*hit[0].first +
                                          (hit[0].second+2) % 3];
                        const uint32_t w1 =
                                out->tris[3*hit[1].first +
                                          (hit[1].second+2) % 3];
                        const float* C2 = &out->verts[3*w0];
                        const float* D2 = &out->verts[3*w1];
                        /*  Churn gate (Nate's eyes, 2026-07-17
                         *  night: "razor-thin z-fighting flat-
                         *  surface scars", THIN screws): a tent
                         *  whose apex barely leaves the plane of
                         *  the triangles it splits adds no
                         *  geometry - just a sliver bump one
                         *  depth-quantum above a flat, flickering
                         *  in the render.  The FIELD cannot see
                         *  this class (|f| ~ 0 either way); the
                         *  plane can.  Refuse tents under 0.02 sp
                         *  of tent height over BOTH split
                         *  triangles.  */
                        const auto plane_h = [&](const float* P,
                                                 const float* Q,
                                                 const float* R) {
                            const float ux = Q[0]-P[0],
                                    uy = Q[1]-P[1], uz = Q[2]-P[2];
                            const float vx2 = R[0]-P[0],
                                    vy2 = R[1]-P[1],
                                    vz2 = R[2]-P[2];
                            float nx = uy*vz2 - uz*vy2;
                            float ny = uz*vx2 - ux*vz2;
                            float nz = ux*vy2 - uy*vx2;
                            const float nl = sqrtf(nx*nx + ny*ny +
                                                   nz*nz);
                            if (!(nl > 0))
                                return 0.f;
                            return fabsf(((px-P[0])*nx +
                                          (py-P[1])*ny +
                                          (pz-P[2])*nz) / nl);
                        };
                        const float th = std::max(
                                plane_h(A2, B2, C2),
                                plane_h(A2, B2, D2));
                        if (th < 0.02f * soup.spacing)
                        {
                            /*  Curvature carve-out (quality
                             *  P2.2): the churn bar exists for
                             *  FLATS (sliver bumps z-fight in the
                             *  render), but the refused
                             *  population bunches at 0.012-0.019
                             *  sp right under the bar - on a
                             *  curved blend that tent height IS
                             *  the chord sagitta (L^2 kappa / 8),
                             *  real bend geometry.  Allow the
                             *  tent where the measured curvature
                             *  accounts for the apex height; a
                             *  flat reads kappa ~ 0 and still
                             *  refuses.  0.005 sp absolute floor
                             *  keeps zero-depth churn out.
                             *  STIBIUM_DMESH_CHURN_CURVE=0
                             *  disables.  */
                            bool curved_ok = false;
                            static const char* cc_env = getenv(
                                    "STIBIUM_DMESH_CHURN_CURVE");
                            /*  CALIBRATION HUNT OPEN (Nate's
                             *  eyes on curve-r1: a few fresh
                             *  chips along small steps, a few of
                             *  the same class cleaned - net
                             *  visually neutral, metrics way up).
                             *  Two separators REFUTED trying to
                             *  exclude the step class: kappa
                             *  ceiling 2/lsp (blend fillets are
                             *  SUB-CELL radius and read as sharp
                             *  as steps; saves 875 -> 13, worst
                             *  back to 0.170) and law-routing
                             *  (screws' step-base saves are
                             *  UNTRACED-routed 183/186, while the
                             *  deepest bino cure - worst 0.170 ->
                             *  0.093 at the (-2.7, 37.1, 61.7)
                             *  mystery site - is LAW-routed).
                             *  Next move needs eye-referee
                             *  coordinates of an actual fresh
                             *  step chip.  */
                            /*  DEPTH floor 0.02 sp (the step-row
                             *  autopsy, Nate's coordinate
                             *  (-5.88, 48.79, 70.53)): the
                             *  carve-out was saving tents on
                             *  0.009 sp chips along step bases -
                             *  3-micron divots nobody can see,
                             *  cured by lone off-plane slivers
                             *  everybody can.  The 0.093 sp win
                             *  came from the deep tail (0.05-
                             *  0.17 sp) alone.  Depth separates
                             *  what three geometric separators
                             *  (kappa ceiling, law routing,
                             *  distance) could not.  */
                            if ((!cc_env || atoi(cc_env) != 0) &&
                                snap_d[s] >= 0.02f * soup.spacing &&
                                th >= 0.005f * soup.spacing)
                            {
                                KOut kc;
                                if (curvature_probe(deck, ctx,
                                        0.1f * soup.spacing,
                                        px, py, pz, &kc))
                                {
                                    const float ex2 = B2[0]-A2[0],
                                            ey2 = B2[1]-A2[1],
                                            ez2 = B2[2]-A2[2];
                                    const float el2 =
                                            ex2*ex2 + ey2*ey2 +
                                            ez2*ez2;
                                    const float kmag = std::max(
                                            fabsf(kc.k),
                                            fabsf(kc.kmax));
                                    curved_ok = th <=
                                            0.25f * kmag * el2;
                                }
                            }
                            if (curved_ok)
                            {
                                ++churn_saved;
                                tent_curved = true;
                            }
                            else
                            {
                                done[s] = 1;
                                ++snap_skipped;
                                ++skip_churn;
                                if (sdump)
                                    fprintf(sdump, "churn %.4f "
                                            "%.4f %.4f %.4f\n",
                                            mx, my, mz,
                                            snap_d[s] /
                                                soup.spacing);
                                continue;
                            }
                        }
                        const float ref[4][3] = {
                            { (A2[0]+px)*0.5f, (A2[1]+py)*0.5f,
                              (A2[2]+pz)*0.5f },
                            { (B2[0]+px)*0.5f, (B2[1]+py)*0.5f,
                              (B2[2]+pz)*0.5f },
                            { (C2[0]+px)*0.5f, (C2[1]+py)*0.5f,
                              (C2[2]+pz)*0.5f },
                            { (D2[0]+px)*0.5f, (D2[1]+py)*0.5f,
                              (D2[2]+pz)*0.5f } };
                        const float rh = 0.01f * soup.spacing;
                        std::vector<float> rx(28), ry(28), rz(28),
                                rv;
                        for (int q = 0; q < 4; ++q)
                            for (int o = 0; o < 7; ++o)
                            {
                                rx[q*7+o] = ref[q][0];
                                ry[q*7+o] = ref[q][1];
                                rz[q*7+o] = ref[q][2];
                            }
                        for (int q = 0; q < 4; ++q)
                        {
                            rx[q*7+1] += rh;  rx[q*7+2] -= rh;
                            ry[q*7+3] += rh;  ry[q*7+4] -= rh;
                            rz[q*7+5] += rh;  rz[q*7+6] -= rh;
                        }
                        eval_points(deck_base(deck), ctx, rx, ry,
                                    rz, rv);
                        const float bar = std::max(snap_d[s],
                                0.02f * soup.spacing);
                        bool damage = false;
                        for (int q = 0; q < 4 && !damage; ++q)
                        {
                            const float f = rv[q*7];
                            const float gx = (rv[q*7+1]-rv[q*7+2]);
                            const float gy = (rv[q*7+3]-rv[q*7+4]);
                            const float gz = (rv[q*7+5]-rv[q*7+6]);
                            const float gl = sqrtf(gx*gx + gy*gy +
                                    gz*gz) / (2 * rh);
                            if (std::isfinite(f) && gl > 0 &&
                                fabsf(f) / gl > bar)
                                damage = true;
                        }
                        if (damage)
                        {
                            done[s] = 1;
                            ++snap_skipped;
                            ++skip_damage;
                            if (sdump)
                                fprintf(sdump, "damage %.4f %.4f "
                                        "%.4f %.4f\n", mx, my, mz,
                                        snap_d[s] / soup.spacing);
                            continue;
                        }
                    }
                    const uint32_t pi =
                            uint32_t(out->verts.size() / 3);
                    out->verts.push_back(px);
                    out->verts.push_back(py);
                    out->verts.push_back(pz);
                    for (int q = 0; q < 2; ++q)
                    {
                        const uint32_t t = hit[q].first;
                        const int e = hit[q].second;
                        const uint32_t v0 = out->tris[3*t + e];
                        const uint32_t v1 =
                                out->tris[3*t + (e+1) % 3];
                        const uint32_t v2 =
                                out->tris[3*t + (e+2) % 3];
                        out->tris[3*t] = v0;
                        out->tris[3*t + 1] = pi;
                        out->tris[3*t + 2] = v2;
                        out->tris.push_back(pi);
                        out->tris.push_back(v1);
                        out->tris.push_back(v2);
                        touched.insert(t);
                    }
                    done[s] = 1;
                    ++snapped;
                    if (sdump)
                        fprintf(sdump, "%s %.4f %.4f %.4f "
                                "%.4f\n",
                                tent_curved ? "tent-curved"
                                            : "tent",
                                mx, my, mz,
                                snap_d[s] / soup.spacing);
                    progress = true;
                }
            }
            if (sdump)
                fclose(sdump);
            if (getenv("STIBIUM_DMESH_CHIP_DEBUG"))
            {
                fprintf(stderr, "SNAP: %llu chip edges tented "
                        "(%llu onto the surface, %llu onto "
                        "curvature valleys), %llu skipped "
                        "(%llu far, %llu hits, %llu damage, "
                        "%llu churn; %llu churn-refusals saved "
                        "by curvature)\n",
                        (unsigned long long)snapped,
                        (unsigned long long)snap_surf,
                        (unsigned long long)snap_valley,
                        (unsigned long long)snap_skipped,
                        (unsigned long long)skip_far,
                        (unsigned long long)skip_hits,
                        (unsigned long long)skip_damage,
                        (unsigned long long)skip_churn,
                        (unsigned long long)churn_saved);
                if (g_vly_try)
                    fprintf(stderr, "VALLEY: %llu attempts, %llu "
                            "proj-fail, %llu floor-reject, best "
                            "kmin*sp %.4f\n",
                            (unsigned long long)g_vly_try,
                            (unsigned long long)g_vly_proj,
                            (unsigned long long)g_vly_floor,
                            g_vly_kmin);
                /*  Where do the UNTRACED chips live?  5-cell bins,
                 *  top clusters - probe targets for the seed-
                 *  coverage hunt.  */
                std::unordered_map<uint64_t, uint32_t> bins;
                const float bw = 5.f * soup.spacing;
                for (const auto& P : far_pts)
                {
                    const int64_t bx2 = int64_t(floorf(P[0] / bw));
                    const int64_t by2 = int64_t(floorf(P[1] / bw));
                    const int64_t bz2 = int64_t(floorf(P[2] / bw));
                    ++bins[(uint64_t(bx2 & 0x1FFFFF) << 42) |
                           (uint64_t(by2 & 0x1FFFFF) << 21) |
                            uint64_t(bz2 & 0x1FFFFF)];
                }
                std::vector<std::pair<uint32_t, uint64_t>> top;
                for (const auto& [k, n2] : bins)
                    top.push_back({ n2, k });
                std::sort(top.rbegin(), top.rend());
                for (size_t i = 0; i < top.size() && i < 6; ++i)
                {
                    const auto sx17 = [](uint64_t v) {
                        int64_t s = int64_t(v & 0x1FFFFF);
                        if (s & 0x100000)
                            s -= 0x200000;
                        return s;
                    };
                    const uint64_t k = top[i].second;
                    fprintf(stderr, "  FAR cluster: %u chips near "
                            "(%.1f, %.1f, %.1f)\n", top[i].first,
                            (sx17(k >> 42) + 0.5f) * bw,
                            (sx17(k >> 21) + 0.5f) * bw,
                            (sx17(k) + 0.5f) * bw);
                }
            }
        }
        out->snapped = snapped;
        pt.sub("snap: waves+referees");
    }
    pt.mark("snap pass");
    prog_stage(12);
    if (pt.level >= 2)
        fprintf(stderr, "TIME2   %-24s %8.2f s  (%llu pts, "
                "%llu calls)\n", "EVAL total so far",
                g_eval.secs, (unsigned long long)g_eval.pts,
                (unsigned long long)g_eval.calls);

    /*  End-of-pipeline chip referee (diagnostic): the repair-loop
     *  metric measures the raw extraction, so the snap pass needs
     *  its own honest number - a full sweep of the FINAL mesh.  */
    if (getenv("STIBIUM_DMESH_CHIP_DEBUG"))
    {
        std::unordered_set<uint64_t> eseen;
        std::vector<float> fxs, fys, fzs, fhs, flen, fv;
        for (size_t t = 0; t < out->tris.size(); t += 3)
            for (int e = 0; e < 3; ++e)
            {
                uint64_t a = out->tris[t + e];
                uint64_t b = out->tris[t + (e + 1) % 3];
                if (a > b)
                    std::swap(a, b);
                if (!eseen.insert((a << 32) | b).second)
                    continue;
                const float* A2 = &out->verts[3*a];
                const float* B2 = &out->verts[3*b];
                const float ex = B2[0]-A2[0], ey = B2[1]-A2[1],
                            ez = B2[2]-A2[2];
                const float len =
                        sqrtf(ex*ex + ey*ey + ez*ez);
                if (!(len > 0))
                    continue;
                fxs.push_back((A2[0]+B2[0]) * 0.5f);
                fys.push_back((A2[1]+B2[1]) * 0.5f);
                fzs.push_back((A2[2]+B2[2]) * 0.5f);
                fhs.push_back(len * 0.01f);
                flen.push_back(len);
            }
        std::vector<float> gxs(fxs.size()*7), gys(fxs.size()*7),
                gzs(fxs.size()*7);
        for (size_t i = 0; i < fxs.size(); ++i)
        {
            for (int q = 0; q < 7; ++q)
            {
                gxs[i*7+q] = fxs[i];
                gys[i*7+q] = fys[i];
                gzs[i*7+q] = fzs[i];
            }
            gxs[i*7+1] += fhs[i];  gxs[i*7+2] -= fhs[i];
            gys[i*7+3] += fhs[i];  gys[i*7+4] -= fhs[i];
            gzs[i*7+5] += fhs[i];  gzs[i*7+6] -= fhs[i];
        }
        eval_points_mt(deck, ctx, gxs, gys, gzs, fv);
        float worst = 0;
        float wpos[3] = { 0, 0, 0 };
        uint64_t nfinal = 0;
        /*  Site dump (STIBIUM_DMESH_CHIP_DUMP=path): every
         *  defective edge midpoint, BOTH signs - chips (edge
         *  under the surface, f < 0) AND air-chords (edge in
         *  air, f > 0) - "x y z signed_depth_sp len_sp" rows.
         *  The localization instrument for the densely-packed-
         *  details class (defects cluster where features crowd;
         *  the dump is what turns "abounds" into coordinates).  */
        FILE* cdump = nullptr;
        if (const char* cd = getenv("STIBIUM_DMESH_CHIP_DUMP"))
            cdump = fopen(cd, "w");
        for (size_t i = 0; i < fxs.size(); ++i)
        {
            const float f = fv[i*7];
            const float h2 = 2 * fhs[i];
            const float gx = (fv[i*7+1]-fv[i*7+2]) / h2;
            const float gy = (fv[i*7+3]-fv[i*7+4]) / h2;
            const float gz = (fv[i*7+5]-fv[i*7+6]) / h2;
            const float gl = sqrtf(gx*gx + gy*gy + gz*gz);
            if (!(gl > 0) || !std::isfinite(f))
                continue;
            const float dist = fabsf(f) / gl;
            if (dist > 2.0f * flen[i])
                continue;   // gradient not credible (see repair gate)
            if (dist <= flen[i] * 0.03f)
                continue;
            if (cdump)
                fprintf(cdump, "%.6g %.6g %.6g %.4f %.3f\n",
                        fxs[i], fys[i], fzs[i],
                        (f < 0 ? -dist : dist) / soup.spacing,
                        flen[i] / soup.spacing);
            if (f >= 0)
                continue;
            ++nfinal;
            if (dist > worst)
            {
                worst = dist;
                wpos[0] = fxs[i];
                wpos[1] = fys[i];
                wpos[2] = fzs[i];
            }
        }
        if (cdump)
            fclose(cdump);
        fprintf(stderr, "FINAL mesh: %llu chip edges, worst depth "
                "%.3f sp at (%.4f, %.4f, %.4f)\n",
                (unsigned long long)nfinal, worst / soup.spacing,
                wpos[0], wpos[1], wpos[2]);

        dmesh_face_sweep(deck, out->verts, out->tris,
                         soup.spacing,
                         getenv("STIBIUM_DMESH_FACE_DUMP"));
    }

    /*  Quality accounting: closed 2-manifold means every edge is
     *  shared by exactly two triangles.  OPEN edges are counted
     *  on GEOMETRIC (welded) ids - a hole is a geometric fact,
     *  and index-space seams between manifold-pass split sheets
     *  are not holes.  Non-manifoldness stays counted on INDEX
     *  ids: index-manifold structure is exactly what the split
     *  vertices provide, and re-welding them would un-measure
     *  the pass's own work.  */
    std::unordered_map<uint64_t, uint32_t> edge_count;
    for (size_t t = 0; t < out->tris.size(); t += 3)
        for (int e = 0; e < 3; ++e)
        {
            uint64_t a = out->tris[t + e];
            uint64_t b = out->tris[t + (e + 1) % 3];
            if (a > b)
                std::swap(a, b);
            ++edge_count[(a << 32) | b];
        }
    out->nonmanifold_edges = 0;
    for (const auto& [k, n2] : edge_count)
    {
        (void)k;
        if (n2 > 2)
            ++out->nonmanifold_edges;
    }
    {
        const std::vector<uint32_t> wid = weld_ids(out->verts);
        std::unordered_map<uint64_t, uint32_t> gec;
        for (size_t t = 0; t < out->tris.size(); t += 3)
            for (int e = 0; e < 3; ++e)
            {
                uint64_t a = wid[out->tris[t + e]];
                uint64_t b = wid[out->tris[t + (e + 1) % 3]];
                if (a == b)
                    continue;   // degenerate under weld
                if (a > b)
                    std::swap(a, b);
                ++gec[(a << 32) | b];
            }
        out->open_edges = 0;
        for (const auto& [k, n2] : gec)
        {
            (void)k;
            if (n2 < 2)
                ++out->open_edges;
        }
    }
    /*  Pinch-site anatomy (sheet-separation campaign, opening
     *  instrument): at every non-manifold edge midpoint, march the
     *  surface normal both ways and measure how thick the SOLID
     *  web (inward) and the AIR gap (outward) are.  The cure's
     *  shape depends on which species dominates and at what
     *  thickness relative to the local pitch.  */
    if (getenv("STIBIUM_DMESH_NM_DEBUG"))
    {
        std::vector<std::array<float, 3>> mids;
        for (const auto& [k, n2] : edge_count)
        {
            if (n2 <= 2)
                continue;
            const uint32_t va = uint32_t(k >> 32);
            const uint32_t vb = uint32_t(k);
            mids.push_back({
                0.5f * (out->verts[3*va] + out->verts[3*vb]),
                0.5f * (out->verts[3*va+1] + out->verts[3*vb+1]),
                0.5f * (out->verts[3*va+2] + out->verts[3*vb+2]) });
        }
        if (!mids.empty())
        {
            const float sp = soup.spacing;
            const float h = 0.01f * sp;
            const int NT = 40;               // 40 x 0.05 sp = 2 sp
            /*  Per site: 7 gradient taps + 2*NT march points  */
            std::vector<float> qx, qy, qz, qv;
            for (const auto& m : mids)
                for (int q = 0; q < 7; ++q)
                {
                    qx.push_back(m[0] + (q==1?h:q==2?-h:0));
                    qy.push_back(m[1] + (q==3?h:q==4?-h:0));
                    qz.push_back(m[2] + (q==5?h:q==6?-h:0));
                }
            eval_points_mt(deck, ctx, qx, qy, qz, qv);
            std::vector<std::array<float, 3>> nrm(mids.size());
            for (size_t i = 0; i < mids.size(); ++i)
            {
                float gx = qv[i*7+1] - qv[i*7+2];
                float gy = qv[i*7+3] - qv[i*7+4];
                float gz = qv[i*7+5] - qv[i*7+6];
                const float l = sqrtf(gx*gx + gy*gy + gz*gz);
                nrm[i] = l > 0
                    ? std::array<float,3>{ gx/l, gy/l, gz/l }
                    : std::array<float,3>{ 0, 0, 0 };
            }
            qx.clear(); qy.clear(); qz.clear();
            for (size_t i = 0; i < mids.size(); ++i)
                for (int s2 = 0; s2 < 2; ++s2)
                    for (int j = 1; j <= NT; ++j)
                    {
                        const float t = (s2 ? 1.f : -1.f) *
                                0.05f * sp * j;
                        qx.push_back(mids[i][0] + t*nrm[i][0]);
                        qy.push_back(mids[i][1] + t*nrm[i][1]);
                        qz.push_back(mids[i][2] + t*nrm[i][2]);
                    }
            eval_points_mt(deck, ctx, qx, qy, qz, qv);
            int hist_solid[5] = {}, hist_air[5] = {}, none = 0;
            const auto bucket = [](float t) {
                return t < 0.25f ? 0 : t < 0.5f ? 1 : t < 1.f ? 2
                     : t < 2.f ? 3 : 4;
            };
            for (size_t i = 0; i < mids.size(); ++i)
            {
                float solid = 1e9f, air = 1e9f;
                for (int j = 0; j < NT; ++j)     // inward (-n)
                    if (qv[i*2*NT + j] >= 0)
                    {
                        solid = 0.05f * (j + 1);
                        break;
                    }
                for (int j = 0; j < NT; ++j)     // outward (+n)
                    if (qv[i*2*NT + NT + j] <= 0)
                    {
                        air = 0.05f * (j + 1);
                        break;
                    }
                if (solid > 100 && air > 100)
                    ++none;
                else if (solid <= air)
                    ++hist_solid[bucket(solid)];
                else
                    ++hist_air[bucket(air)];
                fprintf(stderr, "NM site (%.4f, %.4f, %.4f) "
                        "solid %.2f sp, air %.2f sp\n",
                        mids[i][0], mids[i][1], mids[i][2],
                        solid > 100 ? -1.f : solid,
                        air > 100 ? -1.f : air);
            }
            fprintf(stderr, "NM anatomy (%zu sites, thinnest "
                    "species; buckets <.25/<.5/<1/<2/none sp):\n"
                    "  solid web: %d %d %d %d\n"
                    "  air gap:   %d %d %d %d\n"
                    "  thick both ways: %d\n",
                    mids.size(),
                    hist_solid[0], hist_solid[1], hist_solid[2],
                    hist_solid[3],
                    hist_air[0], hist_air[1], hist_air[2],
                    hist_air[3], none);
        }
    }
    if (getenv("STIBIUM_DMESH_NM_DEBUG"))
        for (const auto& [k, n2] : edge_count)
        {
            if (n2 == 2)
                continue;
            for (const uint32_t vi : { uint32_t(k >> 32),
                                       uint32_t(k) })
            {
                if (vi >= xvh.size())
                {
                    fprintf(stderr, "NM edge (count %u): "
                            "(%.4f, %.4f, %.4f) split duplicate\n",
                            unsigned(n2), out->verts[3*vi],
                            out->verts[3*vi + 1],
                            out->verts[3*vi + 2]);
                    continue;
                }
                const VH v = xvh[vi];
                const auto it = prov.find(&*v);
                int steiner = -1, ncstr = -1;
                if constexpr (CCDT_MODE)
                {
                    steiner = v->ccdt_3_data()
                            .is_Steiner_vertex_on_edge();
                    ncstr = int(v->ccdt_3_data()
                            .number_of_incident_constraints());
                }
                fprintf(stderr, "NM edge (count %u): "
                        "(%.4f, %.4f, %.4f) prov %d steiner %d "
                        "ncstr %d\n", unsigned(n2),
                        out->verts[3*vi], out->verts[3*vi + 1],
                        out->verts[3*vi + 2],
                        it == prov.end() ? 0 : int(it->second),
                        steiner, ncstr);
            }
        }
    out->iterations = uint64_t(round);
    out->inserted = inserted;
    out->repaired = repaired_total;
    out->repair_rounds = uint64_t(repair_round);
    out->constrained = constrained;
    if constexpr (CCDT_MODE)
        for (auto v = dt.finite_vertices_begin();
             v != dt.finite_vertices_end(); ++v)
            if (v->ccdt_3_data().is_Steiner_vertex_on_edge())
                ++out->steiner;

    tape_ctx_free(ctx);
    return !*halt;
}

}  // namespace

/*  Flat-fan decimation (the cost lever, opened 2026-07-18 from
 *  the screws loop: 91.8% of a dense-band export's triangles are
 *  strictly axis-aligned flats - eighth-cell tessellation of
 *  planes that three vertices could carry).  A vertex whose
 *  ENTIRE incident fan lies in one plane adds nothing: collapse
 *  it into a neighbour.  Creases, law, and every curved surface
 *  survive geometrically - their fans are never coplanar - so
 *  the pass needs no provenance.  In-plane orientation check
 *  guards against fold-overs; boundary-of-plane vertices keep
 *  themselves (their fans contain off-plane triangles).
 *  STIBIUM_DMESH_DECIMATE=0 disables.  */
/*  Edge-collapse link condition (correctness review finding 1 +
 *  Nate's z-fighting sighting, 2026-07-17): collapsing (drop ->
 *  keep) is topologically safe only when their one-rings share
 *  EXACTLY the two vertices opposite the shared edge(s).  A
 *  third shared neighbour means the collapse mints a doubled
 *  coplanar facet - both normals point up, the orientation
 *  guard is blind to it, and it renders as z-fighting.  */
static bool link_ok(const std::vector<uint32_t>& T,
                    const std::vector<std::vector<uint32_t>>& inc,
                    uint32_t drop, uint32_t keep)
{
    std::unordered_set<uint32_t> ring;
    for (const uint32_t t : inc[drop])
        for (int e = 0; e < 3; ++e)
        {
            const uint32_t w = T[t + e];
            if (w != drop && w != keep)
                ring.insert(w);
        }
    int shared = 0;
    std::unordered_set<uint32_t> seen;
    for (const uint32_t t : inc[keep])
        for (int e = 0; e < 3; ++e)
        {
            const uint32_t w = T[t + e];
            if (w != drop && w != keep && ring.count(w) &&
                seen.insert(w).second)
                ++shared;
        }
    return shared <= 2;
}

static void decimate_flats(const Deck* deck, DMesh* out)
{
    /*  Winding referee (MARATHON 2 Tier A cure 1, autopsy: the
     *  150 reversed facets on bino are decimation-minted;
     *  DECIMATE=0 reads 6).  Candidates are SELECTED under the
     *  existing guards, then a batched field-gradient referee
     *  vetoes any collapse that would leave a surviving facet
     *  wound against the field at its post-collapse centroid.
     *  Vetoed candidates forfeit their independence claim for
     *  the pass (retried next pass).  */
    TapeCtx* wctx = deck ? tape_ctx_new(deck) : nullptr;
    float mesh_edge = 0;
    for (size_t t = 0; t + 3 <= out->tris.size() &&
                       t < 64 * 3; t += 3)
    {
        const uint32_t a = out->tris[t], b = out->tris[t+1];
        const float dx = out->verts[3*a] - out->verts[3*b];
        const float dy = out->verts[3*a+1] - out->verts[3*b+1];
        const float dz = out->verts[3*a+2] - out->verts[3*b+2];
        mesh_edge = std::max(mesh_edge,
                sqrtf(dx*dx + dy*dy + dz*dz));
    }
    const size_t nv = out->verts.size() / 3;
    if (!nv)
        return;
    std::vector<uint32_t>& T = out->tris;
    auto vp = [&](uint32_t v, int q) {
        return out->verts[3 * v + q];
    };
    uint64_t dk_veto = 0, dk_fail = 0;
    float dk_kmax_seen = 0;
    for (int pass = 0; pass < 30; ++pass)
    {
        /*  Incident-triangle lists  */
        std::vector<std::vector<uint32_t>> inc(nv);
        for (size_t t = 0; t < T.size(); t += 3)
            for (int e = 0; e < 3; ++e)
                inc[T[t + e]].push_back(uint32_t(t));
        std::vector<uint32_t> remap(nv);
        for (size_t v = 0; v < nv; ++v)
            remap[v] = uint32_t(v);
        std::vector<uint8_t> touched(nv, 0);
        size_t collapsed = 0;
        std::vector<uint32_t> cand_v, cand_t;
        for (size_t v = 0; v < nv; ++v)
        {
            const auto& fan = inc[v];
            if (fan.size() < 3 || touched[v])
                continue;
            /*  Fan plane from the first triangle  */
            float n0[3] = { 0, 0, 0 };
            bool flat = true;
            float scale = 0;
            for (const uint32_t t : fan)
            {
                const uint32_t a = T[t], b = T[t+1], c2 = T[t+2];
                const float ux = vp(b,0)-vp(a,0),
                            uy = vp(b,1)-vp(a,1),
                            uz = vp(b,2)-vp(a,2);
                const float wx = vp(c2,0)-vp(a,0),
                            wy = vp(c2,1)-vp(a,1),
                            wz = vp(c2,2)-vp(a,2);
                float nx = uy*wz - uz*wy, ny = uz*wx - ux*wz,
                      nz = ux*wy - uy*wx;
                const float nl = sqrtf(nx*nx + ny*ny + nz*nz);
                if (!(nl > 0))
                    continue;
                nx /= nl; ny /= nl; nz /= nl;
                if (!(n0[0] || n0[1] || n0[2]))
                {
                    n0[0] = nx; n0[1] = ny; n0[2] = nz;
                    scale = sqrtf(nl);
                    continue;
                }
                if (nx*n0[0] + ny*n0[1] + nz*n0[2] < 1.f - 1e-6f)
                {
                    flat = false;
                    break;
                }
            }
            if (!flat || !(n0[0] || n0[1] || n0[2]))
                continue;
            /*  Collapse v into its nearest fan neighbour, with an
             *  in-plane orientation guard: every surviving fan
             *  triangle must keep its normal direction.  */
            uint32_t tgt = UINT32_MAX;
            float bd = 1e30f;
            for (const uint32_t t : fan)
                for (int e = 0; e < 3; ++e)
                {
                    const uint32_t w = T[t + e];
                    if (w == v || touched[w])
                        continue;
                    const float dx = vp(w,0)-vp(uint32_t(v),0),
                                dy = vp(w,1)-vp(uint32_t(v),1),
                                dz = vp(w,2)-vp(uint32_t(v),2);
                    const float d = dx*dx + dy*dy + dz*dz;
                    if (d < bd)
                    {
                        bd = d;
                        tgt = w;
                    }
                }
            if (tgt == UINT32_MAX)
                continue;
            bool ok = true;
            for (const uint32_t t : fan)
            {
                uint32_t a = T[t], b = T[t+1], c2 = T[t+2];
                if (a == v) a = tgt;
                if (b == v) b = tgt;
                if (c2 == v) c2 = tgt;
                if (a == b || b == c2 || a == c2)
                    continue;        // degenerate: will be culled
                const float ux = vp(b,0)-vp(a,0),
                            uy = vp(b,1)-vp(a,1),
                            uz = vp(b,2)-vp(a,2);
                const float wx = vp(c2,0)-vp(a,0),
                            wy = vp(c2,1)-vp(a,1),
                            wz = vp(c2,2)-vp(a,2);
                const float nx = uy*wz - uz*wy,
                            ny = uz*wx - ux*wz,
                            nz = ux*wy - uy*wx;
                const float dot = nx*n0[0] + ny*n0[1] + nz*n0[2];
                /*  Not just flips: NEAR-degenerate survivors are
                 *  razor slivers whose normals flicker - Nate's
                 *  "gently poked with a razorblade" hairline cuts
                 *  across the flats (autod31, hundreds).  A
                 *  surviving triangle keeps a real fraction of
                 *  the reference area or the collapse is
                 *  refused.  */
                if (dot <= 0.01f * scale * scale)
                {
                    ok = false;      // flip, zero, or razor sliver
                    break;
                }
            }
            if (!ok || !link_ok(T, inc, uint32_t(v), tgt))
                continue;
            cand_v.push_back(uint32_t(v));
            cand_t.push_back(tgt);
            touched[v] = 1;
            touched[tgt] = 1;
            /*  Claim the whole fan (fold conviction, 2026-07-18):
             *  collapses are decided on ORIGINAL positions but
             *  applied at pass end, so a triangle with two
             *  vertices collapsing in one pass is double-remapped
             *  into a configuration NEITHER guard evaluated - the
             *  in-plane fold mint behind every z-fighting scar
             *  (fightpix: decimate_flats solo 522, weld solo 0).
             *  Independent-set rule: every vertex of every fan
             *  triangle becomes untouchable this pass; no triangle
             *  is rewritten twice between checks.  The 30-pass
             *  loop recovers the collapse rate.  */
            for (const uint32_t t : fan)
                for (int e = 0; e < 3; ++e)
                    touched[T[t + e]] = 1;
        }
        /*  Field referee over all candidates' surviving facets,
         *  one batch.  */
        /*  DEFAULT OFF (2026-07-18): the veto was built on the
         *  autopsy's decimation attribution, which dissolved -
         *  the reversed count is a prusa pinch-walk artifact +
         *  weld's handful.  Measured cost on screws: ~20 blocked
         *  collapses -> 68 prusa opens (manifold yes -> no) for
         *  zero measured benefit.  STIBIUM_DMESH_DECVETO=1
         *  re-arms it for experiments.  */
        static const char* dv_env = getenv("STIBIUM_DMESH_DECVETO");
        if (dv_env && atoi(dv_env) != 0 && wctx && !cand_v.empty())
        {
            std::vector<float> rx, ry, rz, rv2;
            std::vector<size_t> rcand;
            std::vector<float> rnrm;
            const float gh = std::max(1e-5f,
                    0.01f * mesh_edge);
            for (size_t ci = 0; ci < cand_v.size(); ++ci)
            {
                const uint32_t v2 = cand_v[ci],
                               tg = cand_t[ci];
                for (const uint32_t t : inc[v2])
                {
                    uint32_t q[3] = { T[t], T[t+1], T[t+2] };
                    for (int e = 0; e < 3; ++e)
                        if (q[e] == v2)
                            q[e] = tg;
                    if (q[0] == q[1] || q[1] == q[2] ||
                        q[0] == q[2])
                        continue;
                    const float ux = vp(q[1],0)-vp(q[0],0),
                            uy = vp(q[1],1)-vp(q[0],1),
                            uz = vp(q[1],2)-vp(q[0],2);
                    const float wx = vp(q[2],0)-vp(q[0],0),
                            wy = vp(q[2],1)-vp(q[0],1),
                            wz = vp(q[2],2)-vp(q[0],2);
                    const float cx = (vp(q[0],0)+vp(q[1],0)+
                                      vp(q[2],0)) / 3;
                    const float cy = (vp(q[0],1)+vp(q[1],1)+
                                      vp(q[2],1)) / 3;
                    const float cz = (vp(q[0],2)+vp(q[1],2)+
                                      vp(q[2],2)) / 3;
                    rcand.push_back(ci);
                    rnrm.push_back(uy*wz - uz*wy);
                    rnrm.push_back(uz*wx - ux*wz);
                    rnrm.push_back(ux*wy - uy*wx);
                    for (int o = 0; o < 6; ++o)
                    {
                        rx.push_back(cx + (o==0?gh:o==1?-gh:0));
                        ry.push_back(cy + (o==2?gh:o==3?-gh:0));
                        rz.push_back(cz + (o==4?gh:o==5?-gh:0));
                    }
                }
            }
            if (!rx.empty())
            {
                eval_points_mt(deck, wctx, rx, ry, rz, rv2);
                std::vector<uint8_t> veto(cand_v.size(), 0);
                for (size_t i = 0; i < rcand.size(); ++i)
                {
                    const float gx = rv2[i*6] - rv2[i*6+1];
                    const float gy = rv2[i*6+2] - rv2[i*6+3];
                    const float gz = rv2[i*6+4] - rv2[i*6+5];
                    /*  DECISIVE opposition only: the candidate
                    *  pool lives beside razors whose survivor
                    *  normals are near-degenerate - a knife-edge
                    *  <= 0 test vetoes on sign noise (measured:
                    *  7.6K blocked welds, razors resurrected).  */
                    const float dt2 = rnrm[i*3]*gx +
                            rnrm[i*3+1]*gy + rnrm[i*3+2]*gz;
                    const float nl2 = sqrtf(rnrm[i*3]*rnrm[i*3] +
                            rnrm[i*3+1]*rnrm[i*3+1] +
                            rnrm[i*3+2]*rnrm[i*3+2]);
                    const float gl2 = sqrtf(gx*gx + gy*gy +
                                            gz*gz);
                    if (dt2 < -0.1f * nl2 * gl2)
                        veto[rcand[i]] = 1;
                }
                size_t kept = 0;
                for (size_t ci = 0; ci < cand_v.size(); ++ci)
                    if (!veto[ci])
                    {
                        cand_v[kept] = cand_v[ci];
                        cand_t[kept] = cand_t[ci];
                        ++kept;
                    }
                cand_v.resize(kept);
                cand_t.resize(kept);
            }
        }
        /*  CURVATURE veto (cylinder quilting, 2026-07-18: the
         *  wireframe splotches - big coplanar plates tiling
         *  cylinder walls - are legal merges by the bit-tight
         *  flatness test, because extraction facets the cylinder
         *  into exactly-planar patches.  Normals cannot see the
         *  underlying curvature; the FIELD can.  Batched Hessian
         *  stencil at each candidate: refuse the collapse when
         *  the merged fan's chord sagitta (kappa * span^2 / 8)
         *  exceeds the bar - a flat reads kappa ~ 0 and merges
         *  freely, a cylinder wall caps its plate size at the
         *  radius the bar allows.  STIBIUM_DMESH_DECCURVE sets
         *  the bar in mesh-edge units (default 0.05; 0
         *  disables).  */
        static const char* dk_env =
                getenv("STIBIUM_DMESH_DECCURVE");
        const float dk_bar = (dk_env ? float(atof(dk_env))
                                     : 0.05f) * mesh_edge;
        if (dk_bar > 0 && wctx && !cand_v.empty())
        {
            const float kh = 0.1f * mesh_edge;
            std::vector<float> kx(cand_v.size() * 19),
                    ky(cand_v.size() * 19),
                    kz(cand_v.size() * 19), kv;
            std::vector<float> span2(cand_v.size(), 0.f);
            for (size_t ci = 0; ci < cand_v.size(); ++ci)
            {
                const uint32_t v2 = cand_v[ci];
                for (int q = 0; q < 19; ++q)
                {
                    kx[ci*19+q] = vp(v2,0) + K_SO[q][0]*kh;
                    ky[ci*19+q] = vp(v2,1) + K_SO[q][1]*kh;
                    kz[ci*19+q] = vp(v2,2) + K_SO[q][2]*kh;
                }
                /*  Post-collapse span: max pairwise distance
                 *  among the fan's vertices.  */
                std::vector<uint32_t> fv;
                for (const uint32_t t : inc[v2])
                    for (int e = 0; e < 3; ++e)
                        fv.push_back(T[t + e]);
                std::sort(fv.begin(), fv.end());
                fv.erase(std::unique(fv.begin(), fv.end()),
                         fv.end());
                float s2 = 0;
                for (size_t i = 0; i < fv.size(); ++i)
                    for (size_t j = i + 1; j < fv.size(); ++j)
                    {
                        const float dx = vp(fv[i],0)-vp(fv[j],0),
                                dy = vp(fv[i],1)-vp(fv[j],1),
                                dz = vp(fv[i],2)-vp(fv[j],2);
                        s2 = std::max(s2,
                                dx*dx + dy*dy + dz*dz);
                    }
                span2[ci] = s2;
            }
            eval_points_mt(deck, wctx, kx, ky, kz, kv);
            size_t kept = 0;
            for (size_t ci = 0; ci < cand_v.size(); ++ci)
            {
                KOut kd;
                bool keep_c = true;
                if (kout_from_vals(&kv[ci*19], kh, &kd))
                {
                    const float kmag = std::max(fabsf(kd.k),
                                                fabsf(kd.kmax));
                    dk_kmax_seen = std::max(dk_kmax_seen, kmag);
                    keep_c = 0.125f * kmag * span2[ci] <= dk_bar;
                }
                else
                    ++dk_fail;
                if (keep_c)
                {
                    cand_v[kept] = cand_v[ci];
                    cand_t[kept] = cand_t[ci];
                    ++kept;
                }
                else
                    ++dk_veto;
            }
            cand_v.resize(kept);
            cand_t.resize(kept);
        }
        for (size_t ci = 0; ci < cand_v.size(); ++ci)
        {
            remap[cand_v[ci]] = cand_t[ci];
            ++collapsed;
        }
        cand_v.clear();
        cand_t.clear();
        if (!collapsed)
            break;
        /*  Apply remaps, cull degenerates  */
        size_t w2 = 0;
        for (size_t t = 0; t < T.size(); t += 3)
        {
            uint32_t a = remap[T[t]], b = remap[T[t+1]],
                     c2 = remap[T[t+2]];
            if (a == b || b == c2 || a == c2)
                continue;
            T[w2] = a; T[w2+1] = b; T[w2+2] = c2;
            w2 += 3;
        }
        T.resize(w2);
    }
    if ((dk_veto || dk_fail) && (getenv("STIBIUM_DMESH_TIME") ||
                                 getenv("STIBIUM_DMESH_CHIP_DEBUG")))
        fprintf(stderr, "DECCURVE: %llu vetoed, %llu probe-fail, "
                "kmax %.4f, edge %.4f\n",
                (unsigned long long)dk_veto,
                (unsigned long long)dk_fail,
                dk_kmax_seen, mesh_edge);
    if (wctx)
        tape_ctx_free(wctx);
}

/*  Sliver weld (the razor-gash class, 2026-07-17): the raw
 *  extraction mints degenerate slivers along crease-adjacent
 *  flats (screws census: 834 pre-decimation; ~90% have a
 *  shortest edge under 0.15 sp) - invisible while buried among
 *  small coplanar neighbours, but decimation leaves the
 *  survivors lying alone across big clean planes where they
 *  render as hairline cuts (Nate's "gently poked with a
 *  razorblade").  Collapse the short edge of any razor triangle
 *  (height < 0.02 sp), moving the FLATTER endpoint into the
 *  more-featured one so creases and law never move; the motion
 *  is sub-bar.  Per-triangle orientation guard against
 *  fold-overs.  Runs before flat decimation.  */
static void weld_slivers(const Deck* deck, DMesh* out, float sp)
{
    /*  Field referee for collapses (attribution 2026-07-18:
     *  WELDSLIV=0 reads 4 prusa-reversed vs 150 with weld on -
     *  THIS pass is the reversed-facet mint; its own-normal
     *  orientation guard judges razor fans with garbage
     *  references, the flip_slivers lesson).  Same batched
     *  post-collapse gradient veto as decimate_flats.  */
    TapeCtx* wwctx = deck ? tape_ctx_new(deck) : nullptr;
    const size_t nv = out->verts.size() / 3;
    if (!nv || !(sp > 0))
        return;
    std::vector<uint32_t>& T = out->tris;
    auto vp = [&](uint32_t v, int q) {
        return out->verts[3 * v + q];
    };
    const float hbar = 0.02f * sp;
    const float ebar = 0.15f * sp;
    /*  24 passes (was 10): the independent-set fan-claim defers
     *  same-pass neighbour welds; the queue needs more rounds to
     *  drain (measured: +20 surviving structures on screws at 10
     *  passes cost prusa manifold=yes).  Converges early via the
     *  !welded break.  */
    for (int pass = 0; pass < 24; ++pass)
    {
        std::vector<std::vector<uint32_t>> inc(nv);
        for (size_t t = 0; t < T.size(); t += 3)
            for (int e = 0; e < 3; ++e)
                inc[T[t + e]].push_back(uint32_t(t));
        /*  Fan planarity per vertex: min pairwise normal dot
         *  across incident triangles (1 = perfectly flat).  */
        const auto tnormal = [&](uint32_t t, float n[3]) {
            const uint32_t a = T[t], b = T[t+1], c2 = T[t+2];
            const float ux = vp(b,0)-vp(a,0), uy = vp(b,1)-vp(a,1),
                        uz = vp(b,2)-vp(a,2);
            const float wx = vp(c2,0)-vp(a,0),
                        wy = vp(c2,1)-vp(a,1),
                        wz = vp(c2,2)-vp(a,2);
            n[0] = uy*wz - uz*wy;
            n[1] = uz*wx - ux*wz;
            n[2] = ux*wy - uy*wx;
            const float l = sqrtf(n[0]*n[0] + n[1]*n[1] +
                                  n[2]*n[2]);
            if (l > 0)
            {
                n[0] /= l; n[1] /= l; n[2] /= l;
            }
            return l > 0;
        };
        const auto planarity = [&](uint32_t v) {
            float ref[3];
            bool have = false;
            float mind = 1.f;
            for (const uint32_t t : inc[v])
            {
                float n[3];
                if (!tnormal(t, n))
                    continue;
                if (!have)
                {
                    ref[0] = n[0]; ref[1] = n[1]; ref[2] = n[2];
                    have = true;
                    continue;
                }
                mind = std::min(mind, n[0]*ref[0] + n[1]*ref[1] +
                                      n[2]*ref[2]);
            }
            return have ? mind : -1.f;
        };
        std::vector<uint32_t> remap(nv);
        for (size_t v = 0; v < nv; ++v)
            remap[v] = uint32_t(v);
        std::vector<uint8_t> touched(nv, 0);
        size_t welded = 0;
        std::vector<uint32_t> wc_v, wc_t;
        for (size_t t = 0; t < T.size(); t += 3)
        {
            const uint32_t vs[3] = { T[t], T[t+1], T[t+2] };
            if (touched[vs[0]] || touched[vs[1]] || touched[vs[2]])
                continue;
            float el[3];
            float emin = 1e30f, emax = 0;
            int emini = 0;
            for (int e = 0; e < 3; ++e)
            {
                const uint32_t a = vs[e], b = vs[(e+1)%3];
                const float dx = vp(b,0)-vp(a,0),
                            dy = vp(b,1)-vp(a,1),
                            dz = vp(b,2)-vp(a,2);
                el[e] = sqrtf(dx*dx + dy*dy + dz*dz);
                if (el[e] < emin)
                {
                    emin = el[e];
                    emini = e;
                }
                emax = std::max(emax, el[e]);
            }
            if (!(emax > 0) || emin > ebar)
                continue;
            float ux = vp(vs[1],0)-vp(vs[0],0),
                  uy = vp(vs[1],1)-vp(vs[0],1),
                  uz = vp(vs[1],2)-vp(vs[0],2);
            float wx = vp(vs[2],0)-vp(vs[0],0),
                  wy = vp(vs[2],1)-vp(vs[0],1),
                  wz = vp(vs[2],2)-vp(vs[0],2);
            const float cxn = uy*wz - uz*wy, cyn = uz*wx - ux*wz,
                        czn = ux*wy - uy*wx;
            const float area2 = sqrtf(cxn*cxn + cyn*cyn + czn*czn);
            if (area2 / emax > 2 * hbar)
                continue;            // not a razor
            const uint32_t a = vs[emini], b = vs[(emini+1)%3];
            /*  Flatter endpoint moves; ties move b into a.  */
            const float pa = planarity(a), pb = planarity(b);
            const uint32_t keep = pa < pb ? a : b;
            const uint32_t drop = pa < pb ? b : a;
            /*  Orientation guard: every surviving triangle of
             *  drop's fan keeps its own normal direction.  */
            bool ok = true;
            for (const uint32_t t2 : inc[drop])
            {
                uint32_t q[3] = { T[t2], T[t2+1], T[t2+2] };
                float n0[3];
                if (!tnormal(t2, n0))
                    continue;
                for (int e = 0; e < 3; ++e)
                    if (q[e] == drop)
                        q[e] = keep;
                if (q[0] == q[1] || q[1] == q[2] || q[0] == q[2])
                    continue;
                const float ax2 = vp(q[1],0)-vp(q[0],0),
                            ay2 = vp(q[1],1)-vp(q[0],1),
                            az2 = vp(q[1],2)-vp(q[0],2);
                const float bx2 = vp(q[2],0)-vp(q[0],0),
                            by2 = vp(q[2],1)-vp(q[0],1),
                            bz2 = vp(q[2],2)-vp(q[0],2);
                const float nx2 = ay2*bz2 - az2*by2,
                            ny2 = az2*bx2 - ax2*bz2,
                            nz2 = ax2*by2 - ay2*bx2;
                if (nx2*n0[0] + ny2*n0[1] + nz2*n0[2] <= 0)
                {
                    ok = false;
                    break;
                }
            }
            if (!ok || !link_ok(T, inc, drop, keep))
                continue;
            wc_v.push_back(drop);
            wc_t.push_back(keep);
            touched[drop] = 1;
            touched[keep] = 1;
            /*  Claim the whole fan (correctness r2 #3: same
             *  decide-on-original/apply-at-pass-end double-remap
             *  hazard decimate_flats was convicted of - weld read
             *  clean empirically, but structural beats lucky).  */
            for (const uint32_t t2 : inc[drop])
                for (int e = 0; e < 3; ++e)
                    touched[T[t2 + e]] = 1;
        }
        /*  Field veto REVERTED for weld (2026-07-18: it blocked
         *  7.6K welds and TRIPLED the razor census 627 -> 1,744
         *  - resurrecting the class weld exists to kill.  The
         *  reversed-facet cure moved to repair_winding's
         *  local-majority flip instead.)  */
        if (false && wwctx && !wc_v.empty())
        {
            std::vector<float> rx, ry, rz, rv2, rnrm;
            std::vector<size_t> rcand;
            const float gh = std::max(1e-5f, 0.05f * sp);
            for (size_t ci = 0; ci < wc_v.size(); ++ci)
            {
                const uint32_t v2 = wc_v[ci], tg = wc_t[ci];
                for (const uint32_t t : inc[v2])
                {
                    uint32_t q[3] = { T[t], T[t+1], T[t+2] };
                    for (int e = 0; e < 3; ++e)
                        if (q[e] == v2)
                            q[e] = tg;
                    if (q[0] == q[1] || q[1] == q[2] ||
                        q[0] == q[2])
                        continue;
                    const float ux = vp(q[1],0)-vp(q[0],0),
                            uy = vp(q[1],1)-vp(q[0],1),
                            uz = vp(q[1],2)-vp(q[0],2);
                    const float wx = vp(q[2],0)-vp(q[0],0),
                            wy = vp(q[2],1)-vp(q[0],1),
                            wz = vp(q[2],2)-vp(q[0],2);
                    rcand.push_back(ci);
                    rnrm.push_back(uy*wz - uz*wy);
                    rnrm.push_back(uz*wx - ux*wz);
                    rnrm.push_back(ux*wy - uy*wx);
                    const float cx = (vp(q[0],0)+vp(q[1],0)+
                                      vp(q[2],0)) / 3;
                    const float cy = (vp(q[0],1)+vp(q[1],1)+
                                      vp(q[2],1)) / 3;
                    const float cz = (vp(q[0],2)+vp(q[1],2)+
                                      vp(q[2],2)) / 3;
                    for (int o = 0; o < 6; ++o)
                    {
                        rx.push_back(cx + (o==0?gh:o==1?-gh:0));
                        ry.push_back(cy + (o==2?gh:o==3?-gh:0));
                        rz.push_back(cz + (o==4?gh:o==5?-gh:0));
                    }
                }
            }
            if (!rx.empty())
            {
                eval_points_mt(deck, wwctx, rx, ry, rz, rv2);
                std::vector<uint8_t> veto(wc_v.size(), 0);
                for (size_t i = 0; i < rcand.size(); ++i)
                {
                    const float gx = rv2[i*6] - rv2[i*6+1];
                    const float gy = rv2[i*6+2] - rv2[i*6+3];
                    const float gz = rv2[i*6+4] - rv2[i*6+5];
                    /*  DECISIVE opposition only: the candidate
                    *  pool lives beside razors whose survivor
                    *  normals are near-degenerate - a knife-edge
                    *  <= 0 test vetoes on sign noise (measured:
                    *  7.6K blocked welds, razors resurrected).  */
                    const float dt2 = rnrm[i*3]*gx +
                            rnrm[i*3+1]*gy + rnrm[i*3+2]*gz;
                    const float nl2 = sqrtf(rnrm[i*3]*rnrm[i*3] +
                            rnrm[i*3+1]*rnrm[i*3+1] +
                            rnrm[i*3+2]*rnrm[i*3+2]);
                    const float gl2 = sqrtf(gx*gx + gy*gy +
                                            gz*gz);
                    if (dt2 < -0.1f * nl2 * gl2)
                        veto[rcand[i]] = 1;
                }
                size_t kept = 0;
                for (size_t ci = 0; ci < wc_v.size(); ++ci)
                    if (!veto[ci])
                    {
                        wc_v[kept] = wc_v[ci];
                        wc_t[kept] = wc_t[ci];
                        ++kept;
                    }
                wc_v.resize(kept);
                wc_t.resize(kept);
            }
        }
        for (size_t ci = 0; ci < wc_v.size(); ++ci)
        {
            remap[wc_v[ci]] = wc_t[ci];
            ++welded;
        }
        wc_v.clear();
        wc_t.clear();
        if (!welded)
            break;
        size_t w2 = 0;
        for (size_t t = 0; t < T.size(); t += 3)
        {
            const uint32_t a = remap[T[t]], b = remap[T[t+1]],
                           c2 = remap[T[t+2]];
            if (a == b || b == c2 || a == c2)
                continue;
            T[w2] = a; T[w2+1] = b; T[w2+2] = c2;
            w2 += 3;
        }
        T.resize(w2);
    }
    if (wwctx)
        tape_ctx_free(wwctx);
}

/*  Coplanar sliver flip (the razor class, round 2, 2026-07-17):
 *  the weld handles short-edge razors, but the LONG-THIN needles
 *  beside constraint lines survive it - their endpoints are far
 *  apart, and collapsing would move real geometry.  Nate's read:
 *  razors radiate through QEM as flat-face scratches, most
 *  visible exactly there.  The cure that moves NOTHING: when a
 *  razor triangle and its neighbour across the long edge are
 *  coplanar, flip the diagonal - the sliver quad re-splits into
 *  two healthy triangles.  Flip accepted only when the minimum
 *  triangle height IMPROVES and both new triangles keep the
 *  plane's orientation (convexity falls out of that test).  */
static void flip_slivers(DMesh* out, float sp)
{
    if (!(sp > 0))
        return;
    std::vector<uint32_t>& T = out->tris;
    auto vp = [&](uint32_t v, int q) {
        return out->verts[3 * v + q];
    };
    const float hbar = 0.05f * sp;
    const auto tri_h = [&](uint32_t a, uint32_t b, uint32_t c2,
                           float n[3]) {
        const float ux = vp(b,0)-vp(a,0), uy = vp(b,1)-vp(a,1),
                    uz = vp(b,2)-vp(a,2);
        const float wx = vp(c2,0)-vp(a,0), wy = vp(c2,1)-vp(a,1),
                    wz = vp(c2,2)-vp(a,2);
        n[0] = uy*wz - uz*wy;
        n[1] = uz*wx - ux*wz;
        n[2] = ux*wy - uy*wx;
        const float a2 = sqrtf(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        float emax = 0;
        const uint32_t vs[3] = { a, b, c2 };
        for (int e = 0; e < 3; ++e)
        {
            const uint32_t p = vs[e], q = vs[(e+1)%3];
            const float dx = vp(q,0)-vp(p,0),
                        dy = vp(q,1)-vp(p,1),
                        dz = vp(q,2)-vp(p,2);
            emax = std::max(emax,
                    sqrtf(dx*dx + dy*dy + dz*dz));
        }
        return emax > 0 ? a2 / emax : 0.f;
    };
    for (int pass = 0; pass < 8; ++pass)
    {
        std::unordered_map<uint64_t, std::vector<uint32_t>> et;
        for (size_t t = 0; t < T.size(); t += 3)
            for (int e = 0; e < 3; ++e)
            {
                uint64_t a = T[t + e], b = T[t + (e+1)%3];
                if (a > b)
                    std::swap(a, b);
                et[(a << 32) | b].push_back(uint32_t(t));
            }
        std::vector<uint8_t> dirty(T.size() / 3, 0);
        size_t flips = 0;
        for (size_t t = 0; t < T.size(); t += 3)
        {
            if (dirty[t / 3])
                continue;
            float n1[3];
            const float h1 = tri_h(T[t], T[t+1], T[t+2], n1);
            if (!(h1 > 0) || h1 > hbar)
                continue;
            /*  Longest edge of the razor  */
            int le = 0;
            float lel = 0;
            for (int e = 0; e < 3; ++e)
            {
                const uint32_t p = T[t + e],
                               q = T[t + (e+1)%3];
                const float dx = vp(q,0)-vp(p,0),
                            dy = vp(q,1)-vp(p,1),
                            dz = vp(q,2)-vp(p,2);
                const float l = dx*dx + dy*dy + dz*dz;
                if (l > lel)
                {
                    lel = l;
                    le = e;
                }
            }
            const uint32_t a = T[t + le], b = T[t + (le+1)%3],
                           c2 = T[t + (le+2)%3];
            uint64_t ka = a, kb = b;
            if (ka > kb)
                std::swap(ka, kb);
            const auto it = et.find((ka << 32) | kb);
            if (it == et.end() || it->second.size() != 2)
                continue;
            const uint32_t t2 = it->second[0] == uint32_t(t)
                    ? it->second[1] : it->second[0];
            if (dirty[t2 / 3])
                continue;
            uint32_t d = UINT32_MAX;
            for (int e = 0; e < 3; ++e)
            {
                const uint32_t w = T[t2 + e];
                if (w != a && w != b)
                    d = w;
            }
            if (d == UINT32_MAX || d == c2)
                continue;
            float n2[3];
            const float h2 = tri_h(T[t2], T[t2+1], T[t2+2], n2);
            /*  Coplanarity against the NEIGHBOUR's plane - the
             *  razor's own normal is numerically unstable
             *  precisely because it is a razor (first cut of
             *  this pass tested razor-normal dot and flipped
             *  NOTHING).  The razor's apex must lie in the
             *  neighbour's plane.  */
            const float l2 = sqrtf(n2[0]*n2[0] + n2[1]*n2[1] +
                                   n2[2]*n2[2]);
            if (!(l2 > 0))
                continue;
            const float px = vp(c2,0) - vp(d,0),
                        py = vp(c2,1) - vp(d,1),
                        pz = vp(c2,2) - vp(d,2);
            const float pdist = fabsf(px*n2[0] + py*n2[1] +
                                      pz*n2[2]) / l2;
            if (pdist > 1e-4f * sp)
                continue;
            /*  Flipped pair: (c2, d, a) and (d, c2, b) - both
             *  must keep the NEIGHBOUR plane's orientation and
             *  beat the old minimum height (convexity falls out
             *  of the orientation test).  */
            float na[3], nb[3];
            const float ha = tri_h(c2, d, a, na);
            const float hb = tri_h(d, c2, b, nb);
            if (!(ha > 0) || !(hb > 0) ||
                std::min(ha, hb) <= std::min(h1, h2) ||
                na[0]*n2[0] + na[1]*n2[1] + na[2]*n2[2] <= 0 ||
                nb[0]*n2[0] + nb[1]*n2[1] + nb[2]*n2[2] <= 0)
                continue;
            T[t] = c2;   T[t+1] = d;    T[t+2] = a;
            T[t2] = d;   T[t2+1] = c2;  T[t2+2] = b;
            dirty[t / 3] = 1;
            dirty[t2 / 3] = 1;
            ++flips;
        }
        if (!flips)
            break;
    }
}

/*  Post-mutation quality recount (correctness review finding 1,
 *  2026-07-17): decimation was the only mesh-mutating stage
 *  running AFTER the quality accounting - its damage shipped
 *  uncounted (the razor-sliver gashes reached Nate's slicer with
 *  a clean ledger row).  Same semantics as the in-pipeline
 *  accounting: opens on welded GEOMETRIC ids, nm on index ids.  */
/*  Slicer-hygiene seal (2026-07-18, the "0 auto-repairs" bar):
 *  manifold-pass twins start position-identical, then a later
 *  vertex-moving pass (sliver weld, snap tent) shifts ONE copy by
 *  ~1e-4 sp, and every slicer's exact matcher reads the seam as
 *  open edges (prusa: bino 528, kinea 1).  Merge vertices within
 *  STIBIUM_DMESH_SEAL * sp (default 3e-4 - a tenth of a MICRON at
 *  product scale, three orders under the 0.1 mm feature bar) onto
 *  the first-seen position, and cull the zero-area slivers the
 *  merge exposes.  0 disables.  */
static void seal_seams(DMesh* out, float sp)
{
    static const char* env = getenv("STIBIUM_DMESH_SEAL");
    const float eps = (env ? float(atof(env)) : 3e-4f) * sp;
    if (!(eps > 0) || out->verts.empty())
        return;
    const size_t nv = out->verts.size() / 3;
    const auto cell = [&](float v) {
        return int64_t(floorf(v / eps));
    };
    const auto key = [&](int64_t x, int64_t y, int64_t z) {
        return (uint64_t(x & 0x1FFFFF) << 42) |
               (uint64_t(y & 0x1FFFFF) << 21) |
                uint64_t(z & 0x1FFFFF);
    };
    std::unordered_map<uint64_t, std::vector<uint32_t>> grid;
    grid.reserve(nv);
    std::vector<uint32_t> remap(nv);
    uint64_t sealed = 0;
    for (size_t v = 0; v < nv; ++v)
    {
        const float x = out->verts[3*v], y = out->verts[3*v + 1],
                    z = out->verts[3*v + 2];
        const int64_t cx = cell(x), cy = cell(y), cz = cell(z);
        uint32_t canon = UINT32_MAX;
        for (int dx = -1; dx <= 1 && canon == UINT32_MAX; ++dx)
            for (int dy = -1; dy <= 1 && canon == UINT32_MAX; ++dy)
                for (int dz = -1; dz <= 1; ++dz)
                {
                    const auto it = grid.find(key(cx+dx, cy+dy,
                                                  cz+dz));
                    if (it == grid.end())
                        continue;
                    for (const uint32_t w : it->second)
                    {
                        const float ex = out->verts[3*w] - x;
                        const float ey = out->verts[3*w+1] - y;
                        const float ez = out->verts[3*w+2] - z;
                        if (ex*ex + ey*ey + ez*ez <= eps * eps)
                        {
                            canon = w;
                            break;
                        }
                    }
                    if (canon != UINT32_MAX)
                        break;
                }
        if (canon == UINT32_MAX)
        {
            grid[key(cx, cy, cz)].push_back(uint32_t(v));
            remap[v] = uint32_t(v);
        }
        else
        {
            remap[v] = canon;
            if (out->verts[3*canon] != x ||
                out->verts[3*canon+1] != y ||
                out->verts[3*canon+2] != z)
                ++sealed;
        }
    }
    std::vector<uint32_t>& T = out->tris;
    /*  First pass: remap + cull point/edge-degenerates.  */
    size_t w2 = 0;
    for (size_t t = 0; t < T.size(); t += 3)
    {
        const uint32_t a = remap[T[t]], b = remap[T[t+1]],
                       c = remap[T[t+2]];
        if (a == b || b == c || a == c)
            continue;
        T[w2] = a; T[w2+1] = b; T[w2+2] = c;
        w2 += 3;
    }
    T.resize(w2);
    /*  Second pass: annihilate opposite-facing twins.  Sealing a
     *  pinch seam brings two sheets into exact contact; their
     *  facing triangles become coincident zero-thickness walls
     *  (same vertex triple, opposite winding) that slicers count
     *  as reversed facets.  Cancel them in pairs, like the soap
     *  film they are.  Same-winding duplicates collapse to one.  */
    std::unordered_map<uint64_t, int> par;
    par.reserve(T.size() / 3);
    const auto tkey = [](uint32_t a, uint32_t b, uint32_t c) {
        uint32_t v[3] = { a, b, c };
        std::sort(v, v + 3);
        uint64_t h = 0xcbf29ce484222325ull;
        for (const uint32_t w : v)
            h = (h ^ w) * 0x100000001b3ull;
        return h;
    };
    const auto parity = [](uint32_t a, uint32_t b, uint32_t c) {
        /*  +1 / -1 by whether (a,b,c) is an even permutation of
         *  the sorted triple.  */
        int p = 1;
        if (a > b) { std::swap(a, b); p = -p; }
        if (b > c) { std::swap(b, c); p = -p; }
        if (a > b) { std::swap(a, b); p = -p; }
        return p;
    };
    for (size_t t = 0; t < T.size(); t += 3)
        par[tkey(T[t], T[t+1], T[t+2])] +=
                parity(T[t], T[t+1], T[t+2]);
    uint64_t annihilated = 0;
    std::unordered_map<uint64_t, int> left = par;
    size_t w3 = 0;
    for (size_t t = 0; t < T.size(); t += 3)
    {
        const uint64_t k = tkey(T[t], T[t+1], T[t+2]);
        const int p = parity(T[t], T[t+1], T[t+2]);
        int& net = par[k];
        int& budget = left[k];
        /*  Keep only |net| copies of the majority orientation;
         *  everything else cancels.  */
        if (net != 0 && (p > 0) == (net > 0) && budget != 0)
        {
            T[w3] = T[t]; T[w3+1] = T[t+1]; T[w3+2] = T[t+2];
            w3 += 3;
            budget -= p;
        }
        else
            ++annihilated;
    }
    T.resize(w3);
    if ((sealed || annihilated) &&
        (getenv("STIBIUM_DMESH_TIME") ||
         getenv("STIBIUM_DMESH_CHIP_DEBUG")))
        fprintf(stderr, "SEAL: %llu near-coincident vertices "
                "merged, %llu twin facets annihilated\n",
                (unsigned long long)sealed,
                (unsigned long long)annihilated);
}

/*  Winding repair (2026-07-18, same bar): a facet whose normal
 *  opposes the field gradient at its centroid is wound backwards
 *  BY DEFINITION (the fold-war's gradient-vs-normal detector,
 *  landed as a cure).  Slicers count these as facets_reversed
 *  (screws: 228).  Batch-eval centroid gradients, flip where
 *  the disagreement is decisive (dot < -0.1 of |n||g|); leave
 *  ambiguous razor noise alone.  STIBIUM_DMESH_WINDING=0
 *  disables.  */
static void repair_winding(const Deck* deck, DMesh* out)
{
    static const char* env = getenv("STIBIUM_DMESH_WINDING");
    if (env && atoi(env) == 0)
        return;
    const size_t nt = out->tris.size() / 3;
    if (!nt)
        return;
    TapeCtx* ctx = tape_ctx_new(deck);
    const std::vector<uint32_t>& T = out->tris;
    const std::vector<float>& V = out->verts;
    /*  Stencil pitch from the mesh scale.  */
    float ex0 = 0;
    for (size_t t = 0; t < std::min<size_t>(nt, 64); ++t)
    {
        const float dx = V[3*T[3*t]] - V[3*T[3*t+1]];
        const float dy = V[3*T[3*t]+1] - V[3*T[3*t+1]+1];
        const float dz = V[3*T[3*t]+2] - V[3*T[3*t+1]+2];
        ex0 = std::max(ex0, sqrtf(dx*dx + dy*dy + dz*dz));
    }
    const float h = std::max(1e-5f, 0.01f * ex0);
    std::vector<float> xs(nt * 6), ys(nt * 6), zs(nt * 6), vals;
    for (size_t t = 0; t < nt; ++t)
    {
        const uint32_t a = T[3*t], b = T[3*t+1], c = T[3*t+2];
        const float cx = (V[3*a] + V[3*b] + V[3*c]) / 3;
        const float cy = (V[3*a+1] + V[3*b+1] + V[3*c+1]) / 3;
        const float cz = (V[3*a+2] + V[3*b+2] + V[3*c+2]) / 3;
        for (int o = 0; o < 6; ++o)
        {
            xs[t*6+o] = cx + (o==0 ? h : o==1 ? -h : 0);
            ys[t*6+o] = cy + (o==2 ? h : o==3 ? -h : 0);
            zs[t*6+o] = cz + (o==4 ? h : o==5 ? -h : 0);
        }
    }
    eval_points_mt(deck, ctx, xs, ys, zs, vals);
    tape_ctx_free(ctx);

    /*  Winding is a COHERENCE property (prusa counts facets whose
     *  traversal disagrees with their neighbors around shared
     *  edges), so per-facet field tests cannot fix it - measured:
     *  145 field-flips left 191 incoherent.  The classic cure:
     *  flood-fill orientation across every clean 2-facet edge
     *  (shared edges must be traversed in opposite directions),
     *  then ONE field vote per connected component sets the
     *  global sign (sum of per-facet normal-dot-gradient,
     *  area-weighted - individual razor centroids are noisy, the
     *  component sum is not).  Pinch edges (> 2 facets) are
     *  barriers, not links.  */
    const auto ekey = [](uint32_t a, uint32_t b) {
        return (uint64_t(std::min(a, b)) << 32) | std::max(a, b);
    };
    std::unordered_map<uint64_t, std::array<uint32_t, 2>> einc;
    std::unordered_map<uint64_t, uint32_t> ecount;
    einc.reserve(nt * 3 / 2);
    for (uint32_t t = 0; t < nt; ++t)
        for (int e = 0; e < 3; ++e)
        {
            const uint64_t k = ekey(T[3*t + e],
                                    T[3*t + (e + 1) % 3]);
            auto& n2 = ecount[k];
            if (n2 < 2)
            {
                auto& pr = einc[k];
                pr[n2] = t;
            }
            ++n2;
        }
    std::vector<uint8_t> flip(nt, 0), seen(nt, 0);
    std::vector<uint32_t> stack;
    uint64_t flipped = 0, comps = 0;
    const auto traverses = [&](uint32_t t, uint32_t a,
                               uint32_t b) {
        /*  Does triangle t (with current flip state) traverse the
         *  directed edge a->b?  */
        for (int e = 0; e < 3; ++e)
        {
            uint32_t p = T[3*t + e], q = T[3*t + (e + 1) % 3];
            if (flip[t])
                std::swap(p, q);
            if (p == a && q == b)
                return true;
        }
        return false;
    };
    for (uint32_t t0 = 0; t0 < nt; ++t0)
    {
        if (seen[t0])
            continue;
        ++comps;
        std::vector<uint32_t> comp;
        stack.assign(1, t0);
        seen[t0] = 1;
        while (!stack.empty())
        {
            const uint32_t t = stack.back();
            stack.pop_back();
            comp.push_back(t);
            for (int e = 0; e < 3; ++e)
            {
                const uint32_t a = T[3*t + e],
                               b = T[3*t + (e + 1) % 3];
                const uint64_t k = ekey(a, b);
                if (ecount[k] != 2)
                    continue;       // boundary or pinch: barrier
                const auto& pr = einc[k];
                const uint32_t u = pr[0] == t ? pr[1] : pr[0];
                if (seen[u])
                    continue;
                /*  t's EFFECTIVE direction (traverses applies
                 *  flip[t]); u evaluated unflipped.  Coherence
                 *  wants opposite directions, so same direction
                 *  means u flips.  */
                const bool t_ab = traverses(t, a, b);
                const bool u_ab = traverses(u, a, b);
                flip[u] = uint8_t(t_ab == u_ab);
                seen[u] = 1;
                stack.push_back(u);
            }
        }
        /*  Field vote for the component's global sign - but only
         *  overrule extraction (globally consistent by
         *  construction) on DECISIVE evidence: sizeable component,
         *  strong majority.  Razor-fan fragments bounded by pinch
         *  edges vote near zero, and flipping them minted the very
         *  backwards edges this pass exists to cure (bino: 632
         *  ambiguous flips -> +12 prusa opens.  Measured, of
         *  course).  */
        double vote = 0, mag = 0;
        for (const uint32_t t : comp)
        {
            const float gx = vals[t*6] - vals[t*6+1];
            const float gy = vals[t*6+2] - vals[t*6+3];
            const float gz = vals[t*6+4] - vals[t*6+5];
            const uint32_t a = T[3*t], b = T[3*t+1], c = T[3*t+2];
            const float ux = V[3*b]-V[3*a],
                        uy = V[3*b+1]-V[3*a+1],
                        uz = V[3*b+2]-V[3*a+2];
            const float wx = V[3*c]-V[3*a],
                        wy = V[3*c+1]-V[3*a+1],
                        wz = V[3*c+2]-V[3*a+2];
            float nx = uy*wz - uz*wy;
            float ny = uz*wx - ux*wz;
            float nz = ux*wy - uy*wx;
            if (flip[t])
            {
                nx = -nx; ny = -ny; nz = -nz;
            }
            const double c2 = double(nx*gx + ny*gy + nz*gz);
            vote += c2;
            mag += std::abs(c2);
        }
        const bool decisive = comp.size() >= 8 &&
                mag > 0 && std::abs(vote) > 0.5 * mag;
        if (!decisive)
            continue;   // leave exactly as extracted
        const uint8_t invert = vote < 0 ? 1 : 0;
        for (const uint32_t t : comp)
            if (flip[t] != invert)
            {
                std::swap(out->tris[3*t+1], out->tris[3*t+2]);
                ++flipped;
            }
    }
    /*  Local-majority singleton flip (2026-07-18): weld-minted
     *  reversed facets are scattered singletons whose every
     *  2-incident edge is traversed the SAME direction as its
     *  partner (prusa's facets_reversed).  A facet unanimously
     *  opposed by its partners is the odd one out - flip it.
     *  Strictly reduces incoherent edges; a few passes settle.  */
    uint64_t maj_flips = 0;
    for (int pass = 0; pass < 6; ++pass)
    {
        uint64_t pf = 0;
        for (uint32_t t = 0; t < nt; ++t)
        {
            int bad = 0, clean = 0;
            for (int e = 0; e < 3; ++e)
            {
                const uint32_t a = T[3*t + e],
                               b = T[3*t + (e + 1) % 3];
                const uint64_t k = ekey(a, b);
                if (ecount[k] != 2)
                    continue;
                const auto& pr = einc[k];
                const uint32_t u = pr[0] == t ? pr[1] : pr[0];
                ++clean;
                if (traverses(t, a, b) == traverses(u, a, b))
                    ++bad;
            }
            if (clean >= 2 && bad == clean)
            {
                flip[t] ^= 1;
                std::swap(out->tris[3*t+1], out->tris[3*t+2]);
                ++pf;
            }
        }
        maj_flips += pf;
        if (!pf)
            break;
    }
    if ((flipped || maj_flips) &&
        (getenv("STIBIUM_DMESH_TIME") ||
         getenv("STIBIUM_DMESH_CHIP_DEBUG")))
        fprintf(stderr, "WINDING: %llu component flips, %llu "
                "majority flips\n",
                (unsigned long long)flipped,
                (unsigned long long)maj_flips);
}

/*  FACE-DEVIATION sweep (Nate's metric audit, 2026-07-20): the
 *  edge-midpoint worst-depth metric cannot see the class he
 *  cares about - "points placed on the surface, but the
 *  triangle formed through air or through the volume" has
 *  on-surface edges and an off-surface INTERIOR - and it kept
 *  improving with no visual change.  This probes every triangle
 *  CENTROID against the tape: deviation |f|/|grad f| (sign
 *  says dips-into-material vs roofs-through-air), reported as
 *  AREA-WEIGHTED statistics, because the eye integrates area,
 *  not maxima (calibrated: screws configs Nate ranked visually
 *  separate 25x on %-area-over-bar and not at all on
 *  worst-depth).  Shared by the export tail and the --facedev
 *  CLI verb (offline referee for PAST meshes).  dump_path:
 *  "x y z signed_dev_sp area_sp2" rows for heatmaps.  */
void dmesh_face_sweep(const Deck* deck,
                      const std::vector<float>& verts,
                      const std::vector<uint32_t>& tris,
                      float spacing, const char* dump_path)
{
    TapeCtx* ctx = tape_ctx_new(deck);
    const size_t ntri = tris.size() / 3;
    std::vector<float> cxs(ntri * 7), cys(ntri * 7),
            czs(ntri * 7), harr(ntri), aarr(ntri), cv;
    for (size_t t = 0; t < ntri; ++t)
    {
        const float* A2 = &verts[3*tris[3*t]];
        const float* B2 = &verts[3*tris[3*t + 1]];
        const float* C2 = &verts[3*tris[3*t + 2]];
        const float gx = (A2[0] + B2[0] + C2[0]) / 3;
        const float gy = (A2[1] + B2[1] + C2[1]) / 3;
        const float gz = (A2[2] + B2[2] + C2[2]) / 3;
        const float ux = B2[0]-A2[0], uy = B2[1]-A2[1],
                uz = B2[2]-A2[2];
        const float vx2 = C2[0]-A2[0], vy2 = C2[1]-A2[1],
                vz2 = C2[2]-A2[2];
        const float nx2 = uy*vz2 - uz*vy2,
                ny2 = uz*vx2 - ux*vz2,
                nz2 = ux*vy2 - uy*vx2;
        aarr[t] = 0.5f * sqrtf(nx2*nx2 + ny2*ny2 + nz2*nz2);
        /*  Probe h rides the triangle's own scale.  */
        const float el = sqrtf(std::max(
                ux*ux + uy*uy + uz*uz,
                vx2*vx2 + vy2*vy2 + vz2*vz2));
        harr[t] = std::max(el * 0.01f, 1e-4f * spacing);
        for (int q = 0; q < 7; ++q)
        {
            cxs[t*7+q] = gx;
            cys[t*7+q] = gy;
            czs[t*7+q] = gz;
        }
        cxs[t*7+1] += harr[t];  cxs[t*7+2] -= harr[t];
        cys[t*7+3] += harr[t];  cys[t*7+4] -= harr[t];
        czs[t*7+5] += harr[t];  czs[t*7+6] -= harr[t];
    }
    eval_points_mt(deck, ctx, cxs, cys, czs, cv);
    FILE* fdump = dump_path ? fopen(dump_path, "w") : nullptr;
    double a_tot = 0, a_dev = 0, wsum = 0;
    float fworst = 0;
    float fwpos[3] = { 0, 0, 0 };
    uint64_t nface = 0;
    const float bar = 0.05f * spacing;
    for (size_t t = 0; t < ntri; ++t)
    {
        const float f = cv[t*7];
        const float h2 = 2 * harr[t];
        const float gx2 = (cv[t*7+1]-cv[t*7+2]) / h2;
        const float gy2 = (cv[t*7+3]-cv[t*7+4]) / h2;
        const float gz2 = (cv[t*7+5]-cv[t*7+6]) / h2;
        const float gl = sqrtf(gx2*gx2 + gy2*gy2 + gz2*gz2);
        if (!(gl > 0) || !std::isfinite(f) || !(aarr[t] > 0))
            continue;
        const float dist = fabsf(f) / gl;
        /*  Credibility gate, as in the edge sweep: a first-
         *  order distance is not trustworthy past ~2 local
         *  edge lengths.  */
        if (dist > 200.f * harr[t])
            continue;
        a_tot += aarr[t];
        wsum += double(dist) * aarr[t];
        if (dist > bar)
        {
            ++nface;
            a_dev += aarr[t];
        }
        if (dist > fworst)
        {
            fworst = dist;
            fwpos[0] = cxs[t*7];
            fwpos[1] = cys[t*7];
            fwpos[2] = czs[t*7];
        }
        if (fdump && dist > 0.01f * spacing)
            fprintf(fdump, "%.6g %.6g %.6g %.4f %.5f\n",
                    cxs[t*7], cys[t*7], czs[t*7],
                    (f < 0 ? -dist : dist) / spacing,
                    aarr[t] / (spacing * spacing));
    }
    if (fdump)
        fclose(fdump);
    fprintf(stderr, "FACE: %llu faces over 0.05 sp (%.3f%% of "
            "area), worst %.3f sp at (%.4f, %.4f, %.4f), "
            "area-weighted mean %.4f sp\n",
            (unsigned long long)nface,
            a_tot > 0 ? 100.0 * a_dev / a_tot : 0.0,
            fworst / spacing, fwpos[0], fwpos[1], fwpos[2],
            a_tot > 0 ? wsum / a_tot / spacing : 0.0);
    tape_ctx_free(ctx);
}

/*  Edge-pinch split, SEAM-CLOSURE-CONSISTENT (the "0
 *  auto-repairs" finish line; design per doc/reviews/
 *  2026-07-18-prusa-internals.md sec 4.2): residual prusa opens
 *  are REAL pinch seams (4-/6-incident edges, net winding 0)
 *  where thin features self-contact.  The naive split (nm edges
 *  as pure barriers in vertex-fan connectivity) un-shares each
 *  seam without re-closing the sheets: the two faces that should
 *  stay partners across a seam edge land in different vertex
 *  copies, so every seam half-edge becomes a boundary (bino:
 *  274 nm -> 2,202 open, measured).  Cure: PAIR the incident
 *  faces at each nm edge - radial sort around the edge axis,
 *  then couple each b->a-wound face with the next a->b-wound
 *  face counterclockwise (the two faces bounding each MATERIAL
 *  wedge; winding is trustworthy - repair_winding runs first) -
 *  and let those pairs COUNT as fan connectivity alongside clean
 *  (2-incident) edges.  Each sheet then re-closes by
 *  construction: paired faces share both endpoint copies, so
 *  every split edge copy keeps exactly 2 opposite-wound faces.
 *  Copies are COINCIDENT (zero vertex motion - geometric welded
 *  counters unchanged, watertight law intact); the split lives
 *  at the index level, which 3MF preserves (STL welds it away).
 *  Runs from the EXPORT path, strictly AFTER optional QEM:
 *  meshopt tears coincident copies into real boundary holes if
 *  the split runs first, and its collapses mint fresh pinches
 *  of their own (bino 274 -> 325 nm, measured) - last position
 *  cures both.  STIBIUM_DMESH_PINCHSPLIT=0 disables.  */
uint64_t dmesh_split_pinches(std::vector<float>& V,
                             std::vector<uint32_t>& T)
{
    /*  DEFAULT ON (2026-07-18, Nate's call after the referee
     *  matrix: prusa manifold=yes on screws/bino/bino+QEM, zero
     *  vertex motion, geometric welded counters unchanged).  */
    const char* env = getenv("STIBIUM_DMESH_PINCHSPLIT");
    if (env && atoi(env) == 0)
        return 0;
    const auto ekey = [](uint32_t a, uint32_t b) {
        return (uint64_t(std::min(a, b)) << 32) | std::max(a, b);
    };
    uint64_t total_split = 0;
    uint32_t first_mid = UINT32_MAX;
    for (int round = 0; round < 8; ++round)
    {
        const size_t nv = V.size() / 3;
        /*  All reads this round come from a FROZEN snapshot; the
         *  live array only receives index rewrites.  Each slot is
         *  owned by its original vertex, so it is written at most
         *  once per round and every lookup key stays valid (the
         *  old read-the-live-array pattern made later vertices
         *  see already-rewritten ids as unknown edges ->
         *  spurious barriers).  */
        const std::vector<uint32_t> T0 = T;
        std::unordered_map<uint64_t, uint32_t> ecount;
        ecount.reserve(T0.size());
        for (size_t t = 0; t < T0.size(); t += 3)
            for (int e = 0; e < 3; ++e)
                ++ecount[ekey(T0[t + e], T0[t + (e + 1) % 3])];
        /*  Face lists for nm (>2-incident) edges, then the
         *  material-wedge pairing per edge.  */
        std::unordered_map<uint64_t, std::vector<uint32_t>> nmf;
        for (uint32_t t = 0; t < T0.size(); t += 3)
            for (int e = 0; e < 3; ++e)
            {
                const uint32_t a = T0[t + e],
                               b = T0[t + (e + 1) % 3];
                if (a == b)
                    continue;
                const uint64_t k = ekey(a, b);
                if (ecount[k] > 2)
                    nmf[k].push_back(t);
            }
        std::unordered_map<uint64_t,
                std::vector<std::pair<uint32_t, uint32_t>>> nmpair;
        for (const auto& [k, faces] : nmf)
        {
            const uint32_t A = uint32_t(k >> 32),
                           B = uint32_t(k & 0xffffffffu);
            const float* pa = &V[3 * A];
            const float* pb = &V[3 * B];
            float dx = pb[0]-pa[0], dy = pb[1]-pa[1],
                  dz = pb[2]-pa[2];
            const float dl = sqrtf(dx*dx + dy*dy + dz*dz);
            if (dl <= 0)
                continue;
            dx /= dl; dy /= dl; dz /= dl;
            /*  Radial basis (e1, e2) perpendicular to the edge;
             *  CCW angle = right-hand rotation about a->b.  */
            float e1x, e1y, e1z;
            if (fabsf(dx) <= fabsf(dy) && fabsf(dx) <= fabsf(dz))
                { e1x = 0; e1y = -dz; e1z = dy; }
            else if (fabsf(dy) <= fabsf(dz))
                { e1x = dz; e1y = 0; e1z = -dx; }
            else
                { e1x = -dy; e1y = dx; e1z = 0; }
            const float e1l = sqrtf(e1x*e1x + e1y*e1y + e1z*e1z);
            e1x /= e1l; e1y /= e1l; e1z /= e1l;
            const float e2x = dy*e1z - dz*e1y,
                        e2y = dz*e1x - dx*e1z,
                        e2z = dx*e1y - dy*e1x;
            struct RF { float ang; uint32_t t; int dir; };
            std::vector<RF> rf;
            rf.reserve(faces.size());
            for (const uint32_t t : faces)
            {
                int dir = 0;
                uint32_t w = 0;
                for (int e = 0; e < 3; ++e)
                {
                    const uint32_t p = T0[t + e],
                                   q = T0[t + (e + 1) % 3];
                    if (p == A && q == B)
                        { dir = 1; w = T0[t + (e + 2) % 3]; }
                    else if (p == B && q == A)
                        { dir = -1; w = T0[t + (e + 2) % 3]; }
                }
                if (!dir)
                    continue;
                const float rx = V[3*w]   - pa[0],
                            ry = V[3*w+1] - pa[1],
                            rz = V[3*w+2] - pa[2];
                rf.push_back({ atan2f(rx*e2x + ry*e2y + rz*e2z,
                                      rx*e1x + ry*e1y + rz*e1z),
                               t, dir });
            }
            std::sort(rf.begin(), rf.end(),
                    [](const RF& x, const RF& y)
                    { return x.ang < y.ang; });
            /*  A dir- face has material on its CCW side, a dir+
             *  face on its clockwise side: each material wedge
             *  runs from a dir- face to the next dir+ face CCW.
             *  ZERO-width wedges are skipped: at a doubled-sheet
             *  kiss two walls' faces sit at angularly COINCIDENT
             *  positions, so sort order between them is float
             *  noise - a noise pairing couples wall A to wall B
             *  (a material wedge with no volume) and the fan
             *  union then glues both sheets, blocking the split.
             *  Greedy nearest-CCW past the eps keeps degenerate
             *  (odd-count) fans safe; leftovers stay barriers.  */
            std::vector<char> used(rf.size(), 0);
            auto& prs = nmpair[k];
            for (size_t s = 0; s < rf.size(); ++s)
            {
                if (used[s] || rf[s].dir != -1)
                    continue;
                for (size_t j = 1; j < rf.size(); ++j)
                {
                    const size_t c = (s + j) % rf.size();
                    if (used[c] || rf[c].dir != 1)
                        continue;
                    float dw = rf[c].ang - rf[s].ang;
                    if (dw < 0)
                        dw += 2.f * float(M_PI);
                    if (dw < 1e-3f)
                        continue;
                    prs.push_back({ rf[s].t, rf[c].t });
                    used[s] = used[c] = 1;
                    break;
                }
            }
        }
        std::vector<std::vector<uint32_t>> inc(nv);
        for (uint32_t t = 0; t < T0.size(); t += 3)
            for (int e = 0; e < 3; ++e)
                inc[T0[t + e]].push_back(t);
        uint64_t split = 0;
        for (size_t v = 0; v < nv; ++v)
        {
            const auto& fan = inc[v];
            if (fan.size() < 2)
                continue;
            /*  Union facets sharing a CLEAN (2-incident) edge at
             *  v, plus PAIRED facets across nm edges at v.  */
            int ncomp = 0;
            std::vector<int> parent(fan.size());
            for (size_t i = 0; i < fan.size(); ++i)
                parent[i] = int(i);
            std::function<int(int)> find = [&](int x) {
                while (parent[x] != x)
                    x = parent[x] = parent[parent[x]];
                return x;
            };
            std::unordered_map<uint64_t,
                    std::vector<size_t>> at_v;
            std::unordered_map<uint32_t, size_t> fanidx;
            std::vector<uint64_t> nm_at_v;
            for (size_t i = 0; i < fan.size(); ++i)
            {
                const uint32_t t = fan[i];
                fanidx[t] = i;
                for (int e = 0; e < 3; ++e)
                {
                    const uint32_t a = T0[t + e],
                                   b = T0[t + (e + 1) % 3];
                    if (a != v && b != v)
                        continue;
                    const uint64_t k = ekey(a, b);
                    if (ecount[k] == 2)
                        at_v[k].push_back(i);
                    else if (ecount[k] > 2)
                        nm_at_v.push_back(k);
                }
            }
            for (const auto& [k, lst] : at_v)
                for (size_t j = 1; j < lst.size(); ++j)
                    parent[find(int(lst[0]))] =
                            find(int(lst[j]));
            for (const uint64_t k : nm_at_v)
            {
                const auto it = nmpair.find(k);
                if (it == nmpair.end())
                    continue;
                for (const auto& pr : it->second)
                {
                    const auto ia = fanidx.find(pr.first);
                    const auto ib = fanidx.find(pr.second);
                    if (ia != fanidx.end() && ib != fanidx.end())
                        parent[find(int(ia->second))] =
                                find(int(ib->second));
                }
            }
            std::unordered_map<int, int> roots;
            for (size_t i = 0; i < fan.size(); ++i)
            {
                const int r = find(int(i));
                if (!roots.count(r))
                    roots[r] = ncomp++;
            }
            if (ncomp <= 1)
                continue;
            /*  Component 0 keeps v; each further component gets
             *  a coincident copy (zero vertex motion).  */
            std::vector<uint32_t> copy_of(ncomp, uint32_t(v));
            for (int c2 = 1; c2 < ncomp; ++c2)
            {
                copy_of[c2] = uint32_t(V.size() / 3);
                V.push_back(V[3*v]);
                V.push_back(V[3*v+1]);
                V.push_back(V[3*v+2]);
                ++split;
            }
            for (size_t i = 0; i < fan.size(); ++i)
            {
                const int c2 = roots[find(int(i))];
                if (c2 == 0)
                    continue;
                const uint32_t t = fan[i];
                for (int e = 0; e < 3; ++e)
                    if (T0[t + e] == v)
                        T[t + e] = copy_of[c2];
            }
        }
        total_split += split;
        /*  Stuck seams (bino: 157/274, measured): the two sheets
         *  reconnect through BOTH endpoint fans - the seam
         *  terminates where the surfaces merge, so no vertex
         *  partition separates the edge without tearing clean
         *  edges.  But the slicer manifold test is per-EDGE, not
         *  per-vertex: subdivide the stuck edge at its MIDPOINT
         *  (collinear - zero geometric deviation) and the new
         *  vertex, whose fan is entangled with nothing but the
         *  seam itself, splits cleanly by the pairing next
         *  round.  Endpoints stay fused (bowtie vertices are
         *  legal manifold-test topology).  Midpoint edges never
         *  re-subdivide; a face carrying TWO nm edges defers the
         *  second to the next stall round (its corner moved when
         *  the first was cut); anything stuck at the round cap
         *  is honest residue.  */
        if (!split)
        {
            uint64_t cut = 0;
            for (const auto& [k, faces] : nmf)
            {
                const uint32_t A = uint32_t(k >> 32),
                               B = uint32_t(k & 0xffffffffu);
                if (A >= first_mid || B >= first_mid)
                    continue;
                bool intact = true;
                for (const uint32_t t : faces)
                {
                    bool ha = false, hb = false;
                    for (int e = 0; e < 3; ++e)
                    {
                        ha |= T[t + e] == A;
                        hb |= T[t + e] == B;
                    }
                    intact &= ha && hb;
                }
                if (!intact)
                    continue;
                const uint32_t M = uint32_t(V.size() / 3);
                first_mid = std::min(first_mid, M);
                V.push_back((V[3*A] + V[3*B]) * 0.5f);
                V.push_back((V[3*A+1] + V[3*B+1]) * 0.5f);
                V.push_back((V[3*A+2] + V[3*B+2]) * 0.5f);
                for (const uint32_t t : faces)
                {
                    /*  (.., A, B, ..) -> in-place copy with B
                     *  := M, plus appended copy with A := M -
                     *  pure index substitution, winding kept;
                     *  the (M, w) diagonal is shared opposite-
                     *  wound by the two halves.  */
                    uint32_t half[3];
                    for (int e = 0; e < 3; ++e)
                    {
                        half[e] = T[t + e] == A ? M : T[t + e];
                        if (T[t + e] == B)
                            T[t + e] = M;
                    }
                    T.insert(T.end(), half, half + 3);
                }
                ++cut;
            }
            if (!cut)
                break;
        }
    }
    if (total_split && (getenv("STIBIUM_DMESH_TIME") ||
                        getenv("STIBIUM_DMESH_CHIP_DEBUG")))
        fprintf(stderr, "PINCHSPLIT: %llu sheet copies minted at "
                "pinch seams\n", (unsigned long long)total_split);
    return total_split;
}

static void recount_quality(DMesh* out)
{
    std::unordered_map<uint64_t, uint32_t> ec;
    for (size_t t = 0; t < out->tris.size(); t += 3)
        for (int e = 0; e < 3; ++e)
        {
            uint64_t a = out->tris[t + e];
            uint64_t b = out->tris[t + (e + 1) % 3];
            if (a > b)
                std::swap(a, b);
            ++ec[(a << 32) | b];
        }
    out->nonmanifold_edges = 0;
    for (const auto& [k, n2] : ec)
    {
        (void)k;
        if (n2 > 2)
            ++out->nonmanifold_edges;
    }
    const std::vector<uint32_t> wid = weld_ids(out->verts);
    std::unordered_map<uint64_t, uint32_t> gec;
    for (size_t t = 0; t < out->tris.size(); t += 3)
        for (int e = 0; e < 3; ++e)
        {
            uint64_t a = wid[out->tris[t + e]];
            uint64_t b = wid[out->tris[t + (e + 1) % 3]];
            if (a == b)
                continue;
            if (a > b)
                std::swap(a, b);
            ++gec[(a << 32) | b];
        }
    out->open_edges = 0;
    for (const auto& [k, n2] : gec)
    {
        (void)k;
        if (n2 < 2)
            ++out->open_edges;
    }
}

/*  Observability pack (Nate's ask, 2026-07-17): binary STL
 *  writers so every intermediate product can be dropped into a
 *  normal slicer next to the export.  */
static void write_stl_raw(const char* path,
                          const std::vector<float>& V,
                          const std::vector<uint32_t>& T)
{
    FILE* f = fopen(path, "wb");
    if (!f)
        return;
    char hdr[80] = "stibium stage dump";
    fwrite(hdr, 1, 80, f);
    const uint32_t nt = uint32_t(T.size() / 3);
    fwrite(&nt, 4, 1, f);
    for (uint32_t t = 0; t < nt; ++t)
    {
        float rec[12] = { 0, 0, 0 };
        for (int e = 0; e < 3; ++e)
            for (int q = 0; q < 3; ++q)
                rec[3 + 3*e + q] = V[3 * T[3*t + e] + q];
        const float ux = rec[6]-rec[3], uy = rec[7]-rec[4],
                    uz = rec[8]-rec[5];
        const float vx = rec[9]-rec[3], vy = rec[10]-rec[4],
                    vz = rec[11]-rec[5];
        rec[0] = uy*vz - uz*vy;
        rec[1] = uz*vx - ux*vz;
        rec[2] = ux*vy - uy*vx;
        fwrite(rec, 4, 12, f);
        const uint16_t attr = 0;
        fwrite(&attr, 2, 1, f);
    }
    fclose(f);
}

/*  Traced polylines as thin triangular tubes - visible in any
 *  STL viewer next to the mesh (STIBIUM_DMESH_DUMP_CHAINS).  */
static void emit_tube(std::vector<float>& V,
                      std::vector<uint32_t>& T,
                      const float A[3], const float B[3],
                      float rr)
{
    float d[3] = { B[0]-A[0], B[1]-A[1], B[2]-A[2] };
    const float l = sqrtf(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
    if (!(l > 0))
        return;
    for (int q = 0; q < 3; ++q)
        d[q] /= l;
    float u[3] = { fabsf(d[0]) < 0.9f ? 1.f : 0.f,
                   fabsf(d[0]) < 0.9f ? 0.f : 1.f, 0 };
    float w[3] = { d[1]*u[2]-d[2]*u[1],
                   d[2]*u[0]-d[0]*u[2],
                   d[0]*u[1]-d[1]*u[0] };
    const float wl = sqrtf(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
    for (int q = 0; q < 3; ++q)
        w[q] /= wl > 0 ? wl : 1;
    float v2[3] = { d[1]*w[2]-d[2]*w[1],
                    d[2]*w[0]-d[0]*w[2],
                    d[0]*w[1]-d[1]*w[0] };
    const uint32_t base = uint32_t(V.size() / 3);
    for (int end = 0; end < 2; ++end)
        for (int k = 0; k < 3; ++k)
        {
            const float ang = float(k) * 2.0944f;
            const float cx2 = cosf(ang), sx2 = sinf(ang);
            const float* P = end ? B : A;
            V.push_back(P[0] + rr * (cx2*w[0] + sx2*v2[0]));
            V.push_back(P[1] + rr * (cx2*w[1] + sx2*v2[1]));
            V.push_back(P[2] + rr * (cx2*w[2] + sx2*v2[2]));
        }
    for (int k = 0; k < 3; ++k)
    {
        const uint32_t a0 = base + k;
        const uint32_t a1 = base + (k + 1) % 3;
        const uint32_t b0 = base + 3 + k;
        const uint32_t b1 = base + 3 + (k + 1) % 3;
        T.insert(T.end(), { a0, b0, a1 });
        T.insert(T.end(), { a1, b0, b1 });
    }
    T.insert(T.end(), { base, base + 1, base + 2 });
    T.insert(T.end(), { base + 3, base + 5, base + 4 });
}

static void dump_chains_stl(const char* path, const DSoup& soup)
{
    std::vector<float> V;
    std::vector<uint32_t> T;
    const float rr = 0.06f * soup.spacing;
    for (size_t c = 0; c < soup.tchains.size(); ++c)
    {
        const auto& ch = soup.tchains[c];
        const size_t nseg = ch.size() - (soup.tclosed[c] ? 0 : 1);
        for (size_t j = 0; j < nseg && ch.size() >= 2; ++j)
        {
            const DSurfPoint& A = soup.surface[ch[j]];
            const DSurfPoint& B =
                    soup.surface[ch[(j + 1) % ch.size()]];
            const float pa[3] = { A.x, A.y, A.z };
            const float pb[3] = { B.x, B.y, B.z };
            emit_tube(V, T, pa, pb, rr);
        }
    }
    write_stl_raw(path, V, T);
    fprintf(stderr, "CHAINS dump: %zu chains, %zu tris -> %s\n",
            soup.tchains.size(), T.size() / 3, path);
}

bool delaunay_mesh(const Deck* deck, Region r, volatile int* halt,
                   DMesh* out)
{
    /*  Optimism with rollback (2026-07-17): no a-priori gate can
     *  cleanly split thin walls from concave grooves at leaf-corner
     *  resolution (the wall-gap estimator's populations overlap in
     *  0.2-1.5 mm - measured).  So the level-2 cores run
     *  optimistically, and if the extraction reports OPEN edges,
     *  the guilty leaves (open-edge midpoints inside level-2 dense
     *  boxes) are demoted and stage A re-runs.  Damage is measured,
     *  not predicted; models that produce no holes never pay.  */
    std::unordered_map<uint64_t, int> demote;
    std::vector<std::array<float, 3>> noweld;
    std::unordered_map<uint64_t, int> promote;
    /*  Ship the LEAST-DAMAGED attempt (2026-07-18): the loop used
     *  to ship whatever the final attempt produced, and a rollback
     *  that backfired shipped 12 opens when attempt 2 had 4.
     *  Damage is measured - use the measurement at the exit.  */
    DMesh best;
    bool have_best = false;
    const auto better = [](const DMesh& a, const DMesh& b) {
        return a.open_edges != b.open_edges
                ? a.open_edges < b.open_edges
                : a.nonmanifold_edges < b.nonmanifold_edges;
    };
    /*  Trace carry (perf round 3, rock 1a): the crease polylines
     *  are geometric objects of the FIELD - invariant to density
     *  promotion/demotion - yet every strips/retreat re-run
     *  re-traced from scratch (~18 s on bino, the log's continue
     *  swallowed it).  Trace once; later attempts re-append the
     *  carried chain POINTS into the fresh soup (indices are
     *  soup-relative and must be rebuilt).  REFEREE VERDICT
     *  (2026-07-18): REFUTED as a default - the re-trace on the
     *  PROMOTED soup mints MORE law (denser seeds complete more
     *  chains; constrained 12,042 -> 11,504 with carry, a -4.5%
     *  law loss), even though worst depth read BETTER (0.097 vs
     *  0.170 - a lead worth a fresh-eyes look: is some law
     *  HURTING depth?).  OPT-IN via STIBIUM_DMESH_TRACE_CARRY=1
     *  for experiments; the 18 s re-trace is the price of law.  */
    std::vector<std::vector<std::array<float, 3>>> carried;
    std::vector<uint8_t> carried_closed;
    bool have_carry = false;
    for (int attempt = 0; attempt < 4; ++attempt)
    {
        PhaseTimer pt;
        DSoup soup = delaunay_sample(deck, r, halt,
                demote.empty() ? nullptr : &demote,
                noweld.empty() ? nullptr : &noweld,
                promote.empty() ? nullptr : &promote);
        pt.mark("sample+bisect+QEF");
        if (*halt)
            return false;
        static const char* tr_env = getenv("STIBIUM_DMESH_TRACE");
        static const char* tc_env =
                getenv("STIBIUM_DMESH_TRACE_CARRY");
        const bool carry_ok = tc_env && atoi(tc_env) != 0;
        if (have_carry && carry_ok)
        {
            for (size_t c2 = 0; c2 < carried.size(); ++c2)
            {
                std::vector<uint32_t> chain;
                chain.reserve(carried[c2].size());
                for (const auto& p : carried[c2])
                {
                    chain.push_back(
                            uint32_t(soup.surface.size()));
                    soup.surface.push_back({ p[0], p[1], p[2] });
                    ++soup.traced;
                }
                soup.tchains.push_back(std::move(chain));
                soup.tclosed.push_back(carried_closed[c2]);
            }
        }
        else if (!tr_env || atoi(tr_env) != 0)
        {
            prog_stage(3);
            delaunay_trace(deck, r, &soup, halt);
            if (carry_ok && !soup.tchains.empty())
            {
                carried.clear();
                carried_closed.assign(soup.tclosed.begin(),
                                      soup.tclosed.end());
                for (const auto& ch : soup.tchains)
                {
                    std::vector<std::array<float, 3>> pts2;
                    pts2.reserve(ch.size());
                    for (const uint32_t i : ch)
                        pts2.push_back({
                                soup.surface[i].x,
                                soup.surface[i].y,
                                soup.surface[i].z });
                    carried.push_back(std::move(pts2));
                }
                have_carry = true;
            }
        }
        /*  #4a round 4: tseed step-tracing - fillet-boundary
         *  seams marched from the shallow channel, after crease
         *  law has claimed its cells.  Snap targets only.  */
        trace_step_seams(deck, &soup, halt);
        /*  #4a VERTEX RAIL: the part that actually kills seam-
         *  crossing air-chords (measured on the off-axis rim:
         *  tangency warts read second-order field depth, below
         *  every snap floor - snap cannot own this class).
         *  Every traced contact point is an exact f = 0 point;
         *  injected as PLAIN surface vertices - never
         *  constraints, never tracer seeds (we are past
         *  delaunay_trace here) - they give Delaunay a rail of
         *  sites along the seam, and chords across it lose to
         *  the rail on the empty-circumball test.
         *  STIBIUM_DMESH_CONTACT_POINTS=0 keeps the rail off.  */
        {
            /*  OPT-IN (first referee round: off-axis near-site
             *  skinny 1.2 -> 2.0% - the rail-to-lattice tris
             *  are thin by construction; seam-crossing verdict
             *  awaits the eye referee).  */
            const char* cp_env =
                    getenv("STIBIUM_DMESH_CONTACT_POINTS");
            if (cp_env && atoi(cp_env) != 0 &&
                !soup.contact_chains.empty())
            {
                size_t injected = 0;
                for (const auto& ch : soup.contact_chains)
                    for (const auto& p : ch)
                    {
                        soup.surface.push_back(
                                { p[0], p[1], p[2] });
                        ++injected;
                    }
                if (getenv("STIBIUM_DMESH_TIME") ||
                    getenv("STIBIUM_DMESH_CHIP_DEBUG"))
                    fprintf(stderr, "CONTACT: %zu rail vertices "
                            "injected\n", injected);
            }
        }
        if (const char* cp = getenv("STIBIUM_DMESH_DUMP_CHAINS"))
            dump_chains_stl(cp, soup);
        /*  #4a viz: contact chains as slightly fatter tubes
         *  (distinguishable beside the crease tubes).  */
        if (const char* cc = getenv("STIBIUM_DMESH_DUMP_CONTACT"))
        {
            std::vector<float> V;
            std::vector<uint32_t> T;
            for (const auto& ch : soup.contact_chains)
                for (size_t i = 0; i + 1 < ch.size(); ++i)
                    emit_tube(V, T, ch[i].data(),
                              ch[i + 1].data(),
                              0.1f * soup.spacing);
            write_stl_raw(cc, V, T);
        }
        /*  Close-ring strip cores (the additive-joint diagnosis):
         *  chain pairs closer than ~3 dense pitches bound a strip
         *  of REAL sub-articulation geometry (two parts' rims at
         *  slightly different heights - 6-8 rings in a 1.5 mm band
         *  at the bino collar).  The strip between them has no
         *  witnesses and chords ring-to-ring: the visible
         *  interference.  Promote exactly the strip leaves to
         *  level 3 and re-run stage A - tens of leaves, not the
         *  blanket-level-3 disaster; the retreat loop demotes any
         *  that tear.  */
        /*  DEFAULT ON (2026-07-17, Nate's verdict): bino v15
         *  (dedupe + strips) is "by FAR the cleanest
         *  representation we've had" - the 588 nm that benched
         *  this to opt-in are the sub-visual pinch species, and
         *  count is a lying metric in BOTH directions.  */
        static const char* strips_env =
                getenv("STIBIUM_DMESH_STRIPS");
        if ((!strips_env || atoi(strips_env) != 0) &&
            attempt == 0 && promote.empty() && autodense() &&
            !soup.tchains.empty())
        {
            const float strip_r = 0.75f * soup.spacing;
            /*  Level from the MEASURED rim gap (A/B 2026-07-18:
             *  blanket level 3 = clean joints at 3.04M tris,
             *  strips off = ratty joints at 979K - the dial
             *  between them is the gap the detector already
             *  measures).  A strip articulates with ~2 lattice
             *  rows between its rims: gaps of half a cell and up
             *  get quarter-cell pitch (level 2); only tighter
             *  pairs pay for eighth-cells.  STIBIUM_DMESH_
             *  STRIP_GAP sets the threshold (cells); DEFAULT 0 =
             *  blanket level 3, because the dial was measured
             *  nearly moot: 895/1077 zeiss strip gaps sit under a
             *  QUARTER cell (they are the near-tangent assembly
             *  contacts of the pinch census) - bar 0.25 refunds
             *  only 9% of the tris.  The 2.5x is the honest price
             *  of sub-quarter-cell fits; the tri-count lever is
             *  decimation, not strip starvation.  */
            static const char* sg2_env =
                    getenv("STIBIUM_DMESH_STRIP_GAP");
            const float gap2 = (sg2_env ? float(atof(sg2_env))
                                        : 0.f) * soup.spacing;
            std::vector<std::array<float, 4>> pts;  // x,y,z,chain
            for (size_t c = 0; c < soup.tchains.size(); ++c)
                for (const uint32_t ix : soup.tchains[c])
                {
                    const DSurfPoint& p = soup.surface[ix];
                    pts.push_back({ p.x, p.y, p.z, float(c) });
                }
            for (size_t i = 0; i < pts.size(); ++i)
                for (size_t j = i + 1; j < pts.size(); ++j)
                {
                    if (pts[i][3] == pts[j][3])
                        continue;
                    const float dx = pts[i][0] - pts[j][0];
                    const float dy = pts[i][1] - pts[j][1];
                    const float dz = pts[i][2] - pts[j][2];
                    const float d2 = dx*dx + dy*dy + dz*dz;
                    if (d2 >= strip_r * strip_r)
                        continue;
                    const int lvl = d2 >= gap2 * gap2 ? 2 : 3;
                    const float mx = 0.5f * (pts[i][0] +
                                             pts[j][0]);
                    const float my = 0.5f * (pts[i][1] +
                                             pts[j][1]);
                    const float mz = 0.5f * (pts[i][2] +
                                             pts[j][2]);
                    for (const auto& db : soup.dense_boxes)
                        if (mx >= db.lo[0] && mx <= db.hi[0] &&
                            my >= db.lo[1] && my <= db.hi[1] &&
                            mz >= db.lo[2] && mz <= db.hi[2])
                        {
                            int& lv = promote[db.key];
                            lv = std::max(lv, lvl);
                        }
                }
            if (!promote.empty())
            {
                if (getenv("STIBIUM_DMESH_TIME") ||
                    getenv("STIBIUM_DMESH_CHIP_DEBUG"))
                {
                    size_t p2 = 0, p3 = 0;
                    for (const auto& [k, l] : promote)
                        ++(l >= 3 ? p3 : p2);
                    fprintf(stderr, "STRIPS: %zu close-ring strip "
                            "leaves promoted (%zu @2, %zu @3, "
                            "gap bar %.2f sp)\n",
                            promote.size(), p2, p3,
                            gap2 / soup.spacing);
                }
                continue;   // re-run stage A before paying for B+C
            }
        }
        pt.mark("crease tracer");
        if (*halt)
            return false;
        *out = DMesh();
        const bool ok = delaunay_mesh_soup(deck, soup, halt, out);
        pt.mark("mesh (B+C total)");
        if (!ok || out->open_edges == 0 || *halt)
        {
            const char* de =
                    getenv("STIBIUM_DMESH_DECIMATE");
            if (ok && !*halt && (!de || atoi(de) != 0))
            {
                /*  Individual gates for fold bisection (the
                 *  z-fight hunt): the trio shares the DECIMATE
                 *  master switch; each stage can be soloed off.  */
                const char* ws_env =
                        getenv("STIBIUM_DMESH_WELDSLIV");
                const char* fs_env =
                        getenv("STIBIUM_DMESH_FLIPSLIV");
                const char* df_env =
                        getenv("STIBIUM_DMESH_DECFLATS");
                if (!ws_env || atoi(ws_env) != 0)
                    weld_slivers(deck, out, soup.spacing);
                pt.sub("weld_slivers");
                if (!fs_env || atoi(fs_env) != 0)
                    flip_slivers(out, soup.spacing);
                pt.sub("flip_slivers");
                if (!df_env || atoi(df_env) != 0)
                    decimate_flats(deck, out);
                pt.sub("decimate_flats");
                seal_seams(out, soup.spacing);
                repair_winding(deck, out);
                /*  Pinch split moved to the EXPORT tail (post-
                 *  QEM): dmesh_split_pinches.  */
                pt.sub("seal+winding");
                recount_quality(out);
                pt.sub("recount_quality");
                pt.mark("fix stages");
                prog_stage(13);
                if (pt.level >= 2)
                    fprintf(stderr, "TIME2   %-24s %8.2f s  "
                            "(%llu pts, %llu calls)\n",
                            "EVAL grand total", g_eval.secs,
                            (unsigned long long)g_eval.pts,
                            (unsigned long long)g_eval.calls);
            }
            return ok;
        }
        if (!have_best || better(*out, best))
        {
            best = *out;
            have_best = true;
        }

        /*  Map open edges back to their level-2 leaves - on
         *  GEOMETRIC ids, matching the referee: demotion may only
         *  be driven by real holes, never by split-sheet seams
         *  (the phantom-hole vandalism of zeiss autod22).  */
        const std::vector<uint32_t> wid = weld_ids(out->verts);
        std::unordered_map<uint64_t, uint32_t> ec;
        for (size_t t = 0; t < out->tris.size(); t += 3)
            for (int e = 0; e < 3; ++e)
            {
                uint64_t a = wid[out->tris[t + e]];
                uint64_t b = wid[out->tris[t + (e + 1) % 3]];
                if (a == b)
                    continue;
                if (a > b)
                    std::swap(a, b);
                ++ec[(a << 32) | b];
            }
        size_t newly = 0, newly_weld = 0;
        /*  One conviction per leaf per ATTEMPT (a leaf with many
         *  open edges is one failure, not many rollback steps).  */
        std::unordered_set<uint64_t> hit;
        for (const auto& [k, n2] : ec)
        {
            if (n2 >= 2)
                continue;
            const uint32_t va = uint32_t(k >> 32), vb = uint32_t(k);
            const float mx = 0.5f * (out->verts[3*va] +
                                     out->verts[3*vb]);
            const float my = 0.5f * (out->verts[3*va + 1] +
                                     out->verts[3*vb + 1]);
            const float mz = 0.5f * (out->verts[3*va + 2] +
                                     out->verts[3*vb + 2]);
            /*  Half-cell inflation: seam holes sit ON the box
             *  boundary and read epsilon-outside it (measured:
             *  6 opens, 0 demotions, retreat stuck).  */
            const float be = 0.5f * soup.spacing;
            bool claimed = false;
            for (const auto& db : soup.dense_boxes)
            {
                if (db.level < 2 ||
                    mx < db.lo[0] - be || mx > db.hi[0] + be ||
                    my < db.lo[1] - be || my > db.hi[1] + be ||
                    mz < db.lo[2] - be || mz > db.hi[2] + be)
                    continue;
                claimed = true;
                /*  A repeat conviction of a still-dense leaf is
                 *  progress too (it steps down one more level) -
                 *  but a leaf already at flood cannot retreat
                 *  further and must not count as new work.  */
                if (hit.insert(db.key).second &&
                    demote[db.key] < db.level)
                    ++newly;
            }
            /*  Unclaimed holes near welded samples are weld-seam
             *  damage: ban welding there and re-run.  */
            if (!claimed && soup.welded)
            {
                noweld.push_back({ mx, my, mz });
                ++newly_weld;
            }
        }
        for (const uint64_t k : hit)
            ++demote[k];
        if (getenv("STIBIUM_DMESH_TIME") ||
            getenv("STIBIUM_DMESH_CHIP_DEBUG"))
            fprintf(stderr, "AUTOD retreat %d: %llu open edges, "
                    "%zu level-2 leaves demoted, %zu no-weld "
                    "zones\n", attempt + 1,
                    (unsigned long long)out->open_edges, newly,
                    newly_weld);
        /*  Holes neither core- nor weld-attributable cannot be
         *  fixed by retreating.  */
        if (newly == 0 && newly_weld == 0)
            break;
    }
    if (have_best && better(best, *out))
        *out = std::move(best);
    {
        const char* de = getenv("STIBIUM_DMESH_DECIMATE");
        if (!de || atoi(de) != 0)
        {
            PhaseTimer fixpt;
            const char* ws_env =
                    getenv("STIBIUM_DMESH_WELDSLIV");
            const char* fs_env =
                    getenv("STIBIUM_DMESH_FLIPSLIV");
            const char* df_env =
                    getenv("STIBIUM_DMESH_DECFLATS");
            if (!ws_env || atoi(ws_env) != 0)
                weld_slivers(deck, out, r.X[1] - r.X[0]);
            fixpt.sub("weld_slivers");
            if (!fs_env || atoi(fs_env) != 0)
                flip_slivers(out, r.X[1] - r.X[0]);
            fixpt.sub("flip_slivers");
            if (!df_env || atoi(df_env) != 0)
                decimate_flats(deck, out);
            fixpt.sub("decimate_flats");
            seal_seams(out, r.X[1] - r.X[0]);
            repair_winding(deck, out);
            /*  Pinch split moved to the EXPORT tail (post-QEM):
             *  dmesh_split_pinches.  */
            fixpt.sub("seal+winding");
            recount_quality(out);
            fixpt.sub("recount_quality");
            fixpt.mark("fix stages");
            prog_stage(13);
        }
    }
    return true;
}

bool delaunay_mesh_soup(const Deck* deck, const DSoup& soup,
                        volatile int* halt, DMesh* out)
{
    static const char* env = getenv("STIBIUM_DMESH_CCDT");
    const bool use_ccdt = env ? atoi(env) != 0 : true;
    if (use_ccdt)
    {
        /*  Self-healing constraint retries (zeiss run, 2026-07-15):
         *  real models hand the conforming machinery pathologies
         *  one corner at a time (T-junctions, tangential chain
         *  meetings), and each cascade site used to cost EVERY
         *  constraint on the model.  The CCDT code prints the
         *  offending vertex to std::cerr before asserting: capture
         *  it, quarantine the site's local constraints, rebuild.
         *  Up to 8 sites; an unparseable failure degrades to the
         *  unconstrained path exactly as before.  */
        std::vector<std::array<float, 3>> quarantine;
        for (int attempt = 0; attempt < 8; ++attempt)
        {
            std::stringstream capture;
            std::streambuf* old_buf =
                    std::cerr.rdbuf(capture.rdbuf());
            bool ok = false, threw = false;
            std::string what;
            try
            {
                ok = mesh_impl<true, CCDT>(deck, soup, halt, out,
                        quarantine.empty() ? nullptr : &quarantine);
            }
            catch (const std::exception& e)
            {
                threw = true;
                what = e.what();
            }
            std::cerr.rdbuf(old_buf);
            const std::string text = capture.str();
            if (!text.empty())
                fputs(text.c_str(), stderr);
            if (!threw)
                return ok;
            double qx, qy, qz;
            bool parsed = false;
            const size_t vp = text.rfind("vertex #");
            if (vp != std::string::npos)
            {
                const size_t eq = text.find('=', vp);
                if (eq != std::string::npos &&
                    sscanf(text.c_str() + eq + 1, "%lf %lf %lf",
                           &qx, &qy, &qz) == 3)
                    parsed = true;
            }
            if (!parsed)
            {
                /*  Second dialect: "insert_Steiner_point_on_
                 *  subconstraint: Steiner point coincides with an
                 *  existing vertex\n  -> Steiner point: x y z"
                 *  (thrown via what(), which our capture also
                 *  holds when the CCDT prints it).  */
                for (const std::string& hay : { text, what })
                {
                    const size_t sp2 = hay.rfind("Steiner point:");
                    if (sp2 != std::string::npos &&
                        sscanf(hay.c_str() + sp2 + 14,
                               "%lf %lf %lf", &qx, &qy, &qz) == 3)
                    {
                        parsed = true;
                        break;
                    }
                }
            }
            if (!parsed)
            {
                fprintf(stderr, "delaunay: constrained path failed "
                        "(%s), no cascade site parsed; falling "
                        "back unconstrained\n", what.c_str());
                break;
            }
            quarantine.push_back({ float(qx), float(qy),
                                   float(qz) });
            fprintf(stderr, "delaunay: conforming cascade at "
                    "(%.4f, %.4f, %.4f); quarantining and "
                    "retrying (%zu site%s)\n", qx, qy, qz,
                    quarantine.size(),
                    quarantine.size() == 1 ? "" : "s");
            if (attempt == 7)
                fprintf(stderr, "delaunay: quarantine retries "
                        "exhausted; falling back unconstrained\n");
        }
    }
    return mesh_impl<false, DT>(deck, soup, halt, out);
}

#else  // !STIBIUM_HAS_CGAL

bool delaunay_mesh(const Deck*, Region, volatile int*, DMesh*)
{
    return false;
}

#endif

/*
 *  Adaptive-Delaunay meshing, stages A-C (doc/MESH-NEXT.md).
 *  Blueprint: Keeter, "Please Steal My Meshing Algorithm Idea"
 *  (2026-07-03), annotated with this kernel's machinery: the octree
 *  descent runs on STANDARD-pushed tapes (pointwise-exact by the
 *  fuzzer's standard), interval culls hand the far field over as
 *  proven-sign box corners, and bisection batches through eval_r.
 */

#include <cmath>
#include <cstring>
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

    void add_sample(float x, float y, float z, bool inside)
    {
        if (seen_samples.insert(coord_hash(x, y, z)).second)
            soup.samples.push_back({ x, y, z, inside });
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

/*  Evaluate an arbitrary point list in MIN_VOLUME-sized batches on
 *  the given tape (same dummy-Region pattern as get_normals).  */
void eval_points(const Tape* tape, TapeCtx* ctx,
                 const std::vector<float>& xs,
                 const std::vector<float>& ys,
                 const std::vector<float>& zs,
                 std::vector<float>& out)
{
    const size_t n = xs.size();
    out.resize(n);
    for (size_t at = 0; at < n; at += MIN_VOLUME)
    {
        const unsigned count =
                unsigned(n - at < MIN_VOLUME ? n - at : MIN_VOLUME);
        Region dummy = {};
        dummy.voxels = count;
        dummy.X = const_cast<float*>(xs.data() + at);
        dummy.Y = const_cast<float*>(ys.data() + at);
        dummy.Z = const_cast<float*>(zs.data() + at);
        const float* v = tape_eval_r(tape, ctx, dummy);
        memcpy(out.data() + at, v, count * sizeof(float));
    }
}

/*  Leaf block: evaluate the corner lattice, tag signs, and queue
 *  every sign-change lattice edge for bisection.  */
void sample_block(Collector& c, const Region& r, Tape* tape)
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

    const auto idx = [&](uint32_t i, uint32_t j, uint32_t k) {
        return (size_t(k) * cy + j) * cx + i;
    };
    const auto inside = [&](size_t p) { return vals[p] < 0; };

    bool any_in = false, any_out = false;
    for (size_t p = 0; p < n; ++p)
        (inside(p) ? any_in : any_out) = true;

    for (uint32_t k = 0; k < cz; ++k)
        for (uint32_t j = 0; j < cy; ++j)
            for (uint32_t i = 0; i < cx; ++i)
            {
                const size_t p = idx(i, j, k);
                c.add_sample(xs[p], ys[p], zs[p], inside(p));
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
        ++c.soup.hidden_candidates;   // interval said maybe; samples say no

    if (c.spacing == 0)
        c.spacing = r.X[1] - r.X[0];

    /*  Candidate feature cells: any voxel with >= 3 crossing edges.
     *  Judged later, once bisected positions and normals exist.  */
    for (uint32_t k = 0; k < nk; ++k)
        for (uint32_t j = 0; j < nj; ++j)
            for (uint32_t i = 0; i < ni; ++i)
            {
                FeatCell fc;
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

/*  Octree descent, mirroring the production mesher's structure:
 *  scalar interval eval + STANDARD push per level, cull on decided
 *  sign, leaf blocks small enough that their corner lattice fits a
 *  few eval_r batches.  */
constexpr uint64_t LEAF_VOXELS = 64;

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
        sample_block(c, r, sub);
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

void bisect_edges(Collector& c)
{
    const size_t n = c.edges.size();
    if (n == 0)
        return;
    Tape* base = deck_base(c.deck);
    std::vector<float> xs(n), ys(n), zs(n), vals;
    for (int round = 0; round < BISECT_ROUNDS; ++round)
    {
        if (*c.halt)
            return;
        for (size_t e = 0; e < n; ++e)
        {
            xs[e] = (c.edges[e].ax + c.edges[e].bx) * 0.5f;
            ys[e] = (c.edges[e].ay + c.edges[e].by) * 0.5f;
            zs[e] = (c.edges[e].az + c.edges[e].bz) * 0.5f;
        }
        eval_points(base, c.ctx, xs, ys, zs, vals);
        for (size_t e = 0; e < n; ++e)
        {
            PendingEdge& ed = c.edges[e];
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
    c.soup.surface.reserve(n);
    for (const PendingEdge& ed : c.edges)
        c.soup.surface.push_back({ (ed.ax + ed.bx) * 0.5f,
                                   (ed.ay + ed.by) * 0.5f,
                                   (ed.az + ed.bz) * 0.5f });
}

}  // namespace

namespace {

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
        if (min_dot > SPREAD_DOT)
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

        /*  Reject escapes (with a small margin): a QEF that leaves
         *  its cell is extrapolating noise.  */
        const float margin = c.spacing * 0.25f;
        bool ok = true;
        for (int q = 0; q < 3; ++q)
            if (x(q) < fc.lo[q] - margin || x(q) > fc.hi[q] + margin)
                ok = false;
        if (!ok)
            continue;

        c.soup.surface.push_back({ x(0), x(1), x(2) });
        ++added;
        for (uint8_t q = 0; q < fc.n; ++q)
            suppress[fc.pts[q]] = 1;
    }
    c.soup.feature_points = added;

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

DSoup delaunay_sample(const Deck* deck, Region r, volatile int* halt)
{
    Collector c;
    c.deck = deck;
    c.ctx = tape_ctx_new(deck);
    c.halt = halt;
    descend(c, r, deck_base(deck));
    bisect_edges(c);
    feature_points(c);
    tape_ctx_free(c.ctx);
    return std::move(c.soup);
}

////////////////////////////////////////////////////////////////////////////
//  Stages B + C: Delaunay, refinement, extraction (CGAL)

#ifdef STIBIUM_HAS_CGAL

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/Triangulation_cell_base_with_info_3.h>

namespace {

using K   = CGAL::Exact_predicates_inexact_constructions_kernel;
using Vb  = CGAL::Triangulation_vertex_base_with_info_3<int8_t, K>;
using Cb  = CGAL::Triangulation_cell_base_with_info_3<int8_t, K>;
using Tds = CGAL::Triangulation_data_structure_3<Vb, Cb>;
using DT  = CGAL::Delaunay_triangulation_3<K, Tds>;
using DPoint3 = DT::Point;

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
        eval_points(base, ctx, xs, ys, zs, vals);
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

}  // namespace

bool delaunay_mesh(const Deck* deck, Region r, volatile int* halt,
                   DMesh* out)
{
    DSoup soup = delaunay_sample(deck, r, halt);
    if (*halt)
        return false;
    return delaunay_mesh_soup(deck, soup, halt, out);
}

bool delaunay_mesh_soup(const Deck* deck, const DSoup& soup,
                        volatile int* halt, DMesh* out)
{
    TapeCtx* ctx = tape_ctx_new(deck);

    /*  Stage B: triangulate everything (spatial-sort batch insert),
     *  then refine: any tet edge joining an inside vertex directly
     *  to an outside vertex gets a bisected surface point, until no
     *  such edge remains.  At convergence, inside and outside
     *  vertices are separated everywhere by surface vertices.  */
    std::vector<std::pair<DPoint3, int8_t>> pts;
    pts.reserve(soup.samples.size());
    for (const DSample& s : soup.samples)
        pts.push_back({ DPoint3(s.x, s.y, s.z),
                        int8_t(s.inside ? -1 : 1) });
    DT dt(pts.begin(), pts.end());

    /*  Surface points go in one by one with a coincidence guard: on
     *  grid-aligned geometry a bisected point can converge exactly
     *  onto a lattice sample, and overwriting that vertex's inside/
     *  outside witness with 'surface' corrupts extraction.  */
    for (const DSurfPoint& s : soup.surface)
    {
        const auto before = dt.number_of_vertices();
        auto vh = dt.insert(DPoint3(s.x, s.y, s.z));
        if (dt.number_of_vertices() > before)
            vh->info() = 0;
    }

    constexpr int MAX_ROUNDS = 48;
    uint64_t inserted = 0;
    int round = 0;
    for (; round < MAX_ROUNDS; ++round)
    {
        if (*halt)
            break;
        std::vector<PendingEdge> pending;
        std::unordered_set<uint64_t> seen;
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
        if (pending.empty())
            break;
        const auto fresh = bisect_pending(deck, ctx, pending, halt);
        for (const DSurfPoint& s : fresh)
        {
            /*  A bisected point can converge onto an EXISTING vertex
             *  (grid-aligned geometry: surfaces lying on lattice
             *  planes).  Never overwrite a sign witness's info -
             *  count vertices to detect coincidence.  */
            const auto before = dt.number_of_vertices();
            auto vh = dt.insert(DPoint3(s.x, s.y, s.z));
            if (dt.number_of_vertices() > before)
            {
                vh->info() = 0;
                ++inserted;
            }
        }
    }

    /*  Stage C: cell signs.  After convergence a finite cell cannot
     *  contain both signs, so any non-surface vertex decides it;
     *  all-surface cells ask the oracle at their centroid (batched).
     *  Infinite cells are outside by definition.  */
    std::vector<DT::Cell_handle> oracle_cells;
    std::vector<float> cxs, cys, czs, cvals;
    for (auto c = dt.finite_cells_begin();
         c != dt.finite_cells_end(); ++c)
    {
        int8_t sign = 0;
        for (int i = 0; i < 4 && !sign; ++i)
            sign = c->vertex(i)->info();
        if (!sign)
        {
            double x = 0, y = 0, z = 0;
            for (int i = 0; i < 4; ++i)
            {
                x += c->vertex(i)->point().x();
                y += c->vertex(i)->point().y();
                z += c->vertex(i)->point().z();
            }
            oracle_cells.push_back(c);
            cxs.push_back(float(x / 4));
            cys.push_back(float(y / 4));
            czs.push_back(float(z / 4));
            c->info() = 0;
        }
        else
            c->info() = sign;
    }
    if (!oracle_cells.empty())
    {
        eval_points(deck_base(deck), ctx, cxs, cys, czs, cvals);
        for (size_t i = 0; i < oracle_cells.size(); ++i)
            oracle_cells[i]->info() = cvals[i] < 0 ? -1 : 1;
    }

    /*  Extraction: facets whose three corners are surface vertices,
     *  between cells of opposite sign; oriented so the triangle is
     *  CCW seen from the outside cell.  Facets of a tet complex
     *  cannot self-intersect.  */
    std::unordered_map<const void*, uint32_t> vidx;
    out->verts.clear();
    out->tris.clear();
    const auto vertex_index = [&](DT::Vertex_handle v) -> uint32_t {
        const auto found = vidx.find(&*v);
        if (found != vidx.end())
            return found->second;
        const uint32_t id = uint32_t(out->verts.size() / 3);
        vidx.emplace(&*v, id);
        out->verts.push_back(float(v->point().x()));
        out->verts.push_back(float(v->point().y()));
        out->verts.push_back(float(v->point().z()));
        return id;
    };

    for (auto f = dt.finite_facets_begin();
         f != dt.finite_facets_end(); ++f)
    {
        const DT::Cell_handle c = f->first;
        const int i = f->second;
        const DT::Cell_handle n = c->neighbor(i);
        const int8_t cs = dt.is_infinite(c) ? 1 : c->info();
        const int8_t ns = dt.is_infinite(n) ? 1 : n->info();
        if (cs == ns || cs == 0 || ns == 0)
            continue;
        DT::Vertex_handle vs[3];
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
    }

    /*  Quality accounting: closed 2-manifold means every edge is
     *  shared by exactly two triangles.  */
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
    out->open_edges = 0;
    out->nonmanifold_edges = 0;
    for (const auto& [k, n2] : edge_count)
    {
        (void)k;
        if (n2 < 2)
            ++out->open_edges;
        else if (n2 > 2)
            ++out->nonmanifold_edges;
    }
    out->iterations = uint64_t(round);
    out->inserted = inserted;

    tape_ctx_free(ctx);
    return !*halt;
}

#else  // !STIBIUM_HAS_CGAL

bool delaunay_mesh(const Deck*, Region, volatile int*, DMesh*)
{
    return false;
}

#endif

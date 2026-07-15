/*
 *  Adaptive-Delaunay meshing, stages A-C (doc/MESH-NEXT.md).
 *  Blueprint: Keeter, "Please Steal My Meshing Algorithm Idea"
 *  (2026-07-03), annotated with this kernel's machinery: the octree
 *  descent runs on STANDARD-pushed tapes (pointwise-exact by the
 *  fuzzer's standard), interval culls hand the far field over as
 *  proven-sign box corners, and bisection batches through eval_r.
 */

#include <array>
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
                c.add_sample(xs[p], ys[p], zs[p], inside(p),
                             vals[p] == 0);
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
        eval_points(base, ctx, xs, ys, zs, vals);
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
            if (linked(cand, uint32_t(cur)) || true)
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

DSoup delaunay_sample(const Deck* deck, Region r, volatile int* halt)
{
    Collector c;
    c.deck = deck;
    c.ctx = tape_ctx_new(deck);
    c.halt = halt;
    descend(c, r, deck_base(deck));
    bisect_edges(c);
    feature_points(c);
    c.soup.spacing = c.spacing;
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
               volatile int* halt, DMesh* out)
{
    using TPoint = typename Tri::Point;
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
    if constexpr (CCDT_MODE)
    {
        const DChains chains = delaunay_chains(soup);
        std::vector<std::pair<uint32_t, uint32_t>> cand;
        for (size_t c = 0; c < chains.chains.size(); ++c)
        {
            const auto& chain = chains.chains[c];
            for (size_t j = 0; j + 1 < chain.size(); ++j)
                cand.push_back({ chain[j], chain[j + 1] });
            if (chains.closed[c] && chain.size() > 2)
                cand.push_back({ chain.back(), chain.front() });
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
            eval_points(base, ctx, sx, sy, sz, sv);
        for (size_t i = 0; i < cand.size(); ++i)
        {
            bool good = slen[i] > 0;
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
                if (!std::isfinite(f) || !(gl > 0) ||
                    fabsf(f) / gl > std::max(slen[i] * 0.03f,
                            soup.spacing * 0.05f))
                    good = false;
            }
            if (good)
            {
                accepted.push_back(cand[i]);
                const DSurfPoint& A = soup.surface[cand[i].first];
                const DSurfPoint& B = soup.surface[cand[i].second];
                cseg.push_back({ A.x, A.y, A.z, B.x, B.y, B.z });
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
        if (!cand.empty() &&
            float(cand.size() - accepted.size()) >
                    0.10f * float(cand.size()))
        {
            accepted.clear();
            cseg.clear();
        }
    }
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
    const float drop_r = 0.35f * soup.spacing;
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
            near_crease(s.x, s.y, s.z, drop_r))
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
        for (const auto& pr : pts)
        {
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

    /*  Surface points go in one by one with a coincidence guard: on
     *  grid-aligned geometry a bisected point can converge exactly
     *  onto a lattice sample, and overwriting that vertex's inside/
     *  outside witness with 'surface' corrupts extraction.  Handles
     *  are recorded because chain entries (indices into
     *  soup.surface) become constrained-edge endpoints.  */
    std::vector<VH> surf_vh;
    surf_vh.reserve(soup.surface.size());
    const size_t feat_base =
            soup.surface.size() - size_t(soup.feature_points);
    for (size_t si = 0; si < soup.surface.size(); ++si)
    {
        const DSurfPoint& s = soup.surface[si];
        /*  Surface points inside the crease band are redundant
         *  with the constrained polyline and pair into pinch
         *  slivers - drop them, EXCEPT the chain members
         *  themselves.  */
        if (!cseg.empty() && !chain_used.count(uint32_t(si)) &&
            near_crease(s.x, s.y, s.z, drop_r))
        {
            surf_vh.push_back(VH());
            continue;
        }
        const auto before = dt.number_of_vertices();
        VH vh;
        if constexpr (CCDT_MODE)
        {
            vh = dt.insert(TPoint(s.x, s.y, s.z), CH{}, false);
            vh->ccdt_3_data().set_vertex_type(
                    CGAL::CDT_3_vertex_type::CORNER);
        }
        else
            vh = dt.insert(TPoint(s.x, s.y, s.z));
        if (dt.number_of_vertices() > before)
        {
            vh->info() = 0;
            prov.emplace(&*vh, si < feat_base ? 2 : 3);
        }
        surf_vh.push_back(vh);
    }

    /*  The constrained-crease round: chain segments become
     *  constrained edges, batched with one Delaunay restoration at
     *  the end (each insert_constrained_edge(…, true) would restore
     *  individually).  */
    /*  Insert the refereed chain segments as constrained edges,
     *  batched with one Delaunay restoration at the end.  */
    uint64_t constrained = 0;
    if constexpr (CCDT_MODE)
    {
        for (const auto& [a, b] : accepted)
        {
            VH va = surf_vh[a], vb = surf_vh[b];
            if (va == VH() || vb == VH() || va == vb)
                continue;
            dt.insert_constrained_edge(va, vb, false);
            ++constrained;
        }
        dt.restore_Delaunay();
        sweep_steiner();
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

    /*  Stage C runs inside the repair loop: extract, then hunt WARTS
     *  - fold edges (adjacent-triangle normals disagreeing) whose
     *  midpoint is OFF the surface.  A fold on a true crease has its
     *  edge ON the surface and is left alone.  Wart midpoints are
     *  Newton-projected onto the surface and inserted, and the mesh
     *  re-extracted, until clean or capped (error-driven insertion,
     *  after Wang et al. 2025).  */
    constexpr int MAX_REPAIR = 16;
    uint64_t repaired_total = 0;
    int repair_round = 0;
    std::vector<VH> xvh;   // extraction order -> vertex handle
    for (;; ++repair_round)
    {
    /*  Stage C: cell signs.  After convergence a finite cell cannot
     *  contain both signs (except in the crease band, where i/o
     *  edges are deliberately shadowed), so a non-surface vertex
     *  decides it; all-surface AND mixed cells ask the oracle at
     *  their centroid (batched).  Infinite cells are outside by
     *  definition.  */
    std::vector<CH> oracle_cells;
    std::vector<float> cxs, cys, czs, cvals;
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
    xvh.clear();
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
        eval_points(base2, ctx, gxs, gys, gzs, gv);
        std::vector<size_t> cand;
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
            if (dist > wclamp[i] * 0.03f)
                cand.push_back(i);
        }
        if (cand.empty())
            break;

        std::vector<float> kx, ky, kz, kh, kc;
        if (repair_mode >= 2)
        {
            /*  Crease-seek: a midpoint split of a crease-crossing
             *  chord is self-similar and never converges; split at
             *  the |f| peak along the chord instead - the branch
             *  switch, i.e. the crease.  */
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
            eval_points(base2, ctx, sx, sy, sz, sv);
            for (size_t q = 0; q < cand.size(); ++q)
            {
                const size_t i = cand[q];
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
                kh.push_back(whs[i]);
                kc.push_back(wclamp[i]);
            }
        }
        else
        {
            for (const size_t i : cand)
            {
                kx.push_back(wx2[i]);
                ky.push_back(wy2[i]);
                kz.push_back(wz2[i]);
                kh.push_back(whs[i]);
                kc.push_back(wclamp[i]);
            }
        }
        project_points_impl(deck, ctx, kx, ky, kz, kh, kc);
        uint64_t added = 0;
        for (size_t i = 0; i < kx.size(); ++i)
        {
            /*  Repair keep-out: the constrained crease repairs
             *  itself; only the smooth field is ours to press (an
             *  insert on a constrained edge destroys it and the
             *  re-conforming Steiner point lands a sliver away -
             *  the pinch factory, measured 2026-07-15).  */
            if (!cseg.empty() &&
                near_crease(kx[i], ky[i], kz[i],
                            0.75f * soup.spacing))
                continue;
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
                    const double r = 0.25 * kc[i];
                    if (d2 < r * r)
                        continue;
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
        if (!added)
            break;
        repaired_total += added;
        sweep_steiner();
    }
    }   // repair loop

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
    if (getenv("STIBIUM_DMESH_NM_DEBUG"))
        for (const auto& [k, n2] : edge_count)
        {
            if (n2 == 2)
                continue;
            for (const uint32_t vi : { uint32_t(k >> 32),
                                       uint32_t(k) })
            {
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
    static const char* env = getenv("STIBIUM_DMESH_CCDT");
    const bool use_ccdt = env ? atoi(env) != 0 : true;
    if (use_ccdt)
        return mesh_impl<true, CCDT>(deck, soup, halt, out);
    return mesh_impl<false, DT>(deck, soup, halt, out);
}

#else  // !STIBIUM_HAS_CGAL

bool delaunay_mesh(const Deck*, Region, volatile int*, DMesh*)
{
    return false;
}

#endif

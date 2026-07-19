/*
 *  Referees for the adaptive-Delaunay meshing campaign
 *  (doc/MESH-NEXT.md).  Stage gates:
 *
 *  "[.dmesh]"   stage A - every bisected surface point re-evaluates
 *               to ~zero on the tape, signs and counts sane.
 *  "[.dmeshBC]" stages B+C - triangulate, refine, extract; then the
 *               gauntlet v1: closed (no open edges), manifold edge
 *               count, signed volume vs the analytic value.
 */

#include <chrono>
#include <functional>
#include <unordered_set>
#include <cmath>
#include <cstdint>
#include <vector>

#include "catch/catch.hpp"

#include "fab/tree/triangulate/delaunay.h"
#include "fab/tree/triangulate.h"
#include "fab/formats/stl.h"
#include "fab/tree/tape.h"
#include "fab/tree/tree.h"
#include "fab/tree/parser.h"
#include "fab/util/region.h"

namespace {

struct DeckRegion
{
    Deck* deck = nullptr;
    Region r = {};

    DeckRegion(const char* expr, uint32_t n, float half)
    {
        MathTree* tree = parse(expr);
        REQUIRE(tree != nullptr);
        deck = deck_from_tree(tree);
        free_tree(tree);
        r.ni = n;
        r.nj = n;
        r.nk = n;
        r.voxels = uint64_t(n) * n * n;
        build_arrays(&r, -half, -half, -half, half, half, half);
    }
    ~DeckRegion()
    {
        free_arrays(&r);
        deck_free(deck);
    }
};

/*  Connected components over shared vertices: floating slivers are
 *  CLOSED, so they pass open/non-manifold edge checks - only a
 *  component count exposes them (spotted by eyeball on the cube's
 *  chamfered edges).  */
uint32_t mesh_components(const DMesh& m)
{
    const uint32_t nv = uint32_t(m.verts.size() / 3);
    std::vector<uint32_t> parent(nv);
    for (uint32_t v = 0; v < nv; ++v)
        parent[v] = v;
    std::function<uint32_t(uint32_t)> find =
            [&](uint32_t v) -> uint32_t {
        while (parent[v] != v)
        {
            parent[v] = parent[parent[v]];
            v = parent[v];
        }
        return v;
    };
    for (size_t t = 0; t < m.tris.size(); t += 3)
    {
        const uint32_t a = find(m.tris[t]);
        const uint32_t b = find(m.tris[t + 1]);
        const uint32_t c = find(m.tris[t + 2]);
        parent[b] = a;
        parent[c] = a;
    }
    std::unordered_set<uint32_t> roots;
    std::unordered_set<uint32_t> used;
    for (size_t t = 0; t < m.tris.size(); ++t)
        used.insert(m.tris[t]);
    for (uint32_t v : used)
        roots.insert(find(v));
    return uint32_t(roots.size());
}

/*  Wart detector: sharp folds between adjacent triangles.  A crease
 *  is a fold on purpose; a WART is a fold where two faces nearly
 *  double back (dihedral far past any real crease).  Count pairs
 *  whose normals oppose (dot < -0.2 ~ fold sharper than 100
 *  degrees) - Nate counts ~20 by eye on the spheres model.  */
uint32_t mesh_warts(const DMesh& m, float fold_dot = -0.2f)
{
    std::unordered_map<uint64_t, std::pair<uint32_t, uint32_t>> em;
    const auto tri_normal = [&](uint32_t t, double n[3]) {
        const float* a = &m.verts[3 * m.tris[3*t]];
        const float* b = &m.verts[3 * m.tris[3*t + 1]];
        const float* c = &m.verts[3 * m.tris[3*t + 2]];
        const double ux = b[0]-a[0], uy = b[1]-a[1], uz = b[2]-a[2];
        const double wx = c[0]-a[0], wy = c[1]-a[1], wz = c[2]-a[2];
        n[0] = uy*wz - uz*wy;
        n[1] = uz*wx - ux*wz;
        n[2] = ux*wy - uy*wx;
        const double l = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
        if (l > 0) { n[0] /= l; n[1] /= l; n[2] /= l; }
    };
    for (uint32_t t = 0; t < m.tris.size() / 3; ++t)
        for (int e = 0; e < 3; ++e)
        {
            uint64_t a = m.tris[3*t + e];
            uint64_t b = m.tris[3*t + (e + 1) % 3];
            if (a > b) std::swap(a, b);
            const uint64_t k = (a << 32) | b;
            auto it = em.find(k);
            if (it == em.end())
                em[k] = { t, UINT32_MAX };
            else if (it->second.second == UINT32_MAX)
                it->second.second = t;
        }
    uint32_t warts = 0;
    for (const auto& [k, pr] : em)
    {
        (void)k;
        if (pr.second == UINT32_MAX)
            continue;
        double n1[3], n2[3];
        tri_normal(pr.first, n1);
        tri_normal(pr.second, n2);
        const double d = n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2];
        if (d < fold_dot)
            ++warts;
    }
    return warts;
}

double mesh_volume(const DMesh& m)
{
    /*  Divergence theorem over the closed surface: V = sum of
     *  signed tet volumes (origin, v0, v1, v2), positive for
     *  outward orientation.  */
    double vol = 0;
    for (size_t t = 0; t < m.tris.size(); t += 3)
    {
        const float* a = &m.verts[3 * m.tris[t]];
        const float* b = &m.verts[3 * m.tris[t + 1]];
        const float* c = &m.verts[3 * m.tris[t + 2]];
        const double cx = double(b[1]) * c[2] - double(b[2]) * c[1];
        const double cy = double(b[2]) * c[0] - double(b[0]) * c[2];
        const double cz = double(b[0]) * c[1] - double(b[1]) * c[0];
        vol += (a[0] * cx + a[1] * cy + a[2] * cz) / 6.0;
    }
    return vol;
}

}  // namespace

TEST_CASE("Delaunay stage A: surface points sit on the surface",
          "[.dmesh]")
{
    const char* MODELS[] = {
        "-r++qXqYqZf1",                                  // sphere
        "aa-f-0.6X-Xf0.6aa-f-0.6Y-Yf0.6a-f-0.6Z-Zf0.6",  // cube
        "i-r++qXqYqZf1-r++q-Xf0.5q-Yf0.25q-Zf0.1f0.8",   // sphere union
        "ai-r++qXqYqZf1-r++q-Xf0.5qYqZf0.7n-Zf0.2",      // nested CSG
    };
    for (const char* m : MODELS)
    {
        CAPTURE(m);
        DeckRegion d(m, 48, 1.1f);
        volatile int halt = 0;
        const DSoup soup = delaunay_sample(d.deck, d.r, &halt);

        WARN("model " << m << ": " << soup.samples.size()
             << " samples, " << soup.surface.size() << " surface pts, "
             << soup.leaf_blocks << " leaf blocks, "
             << soup.culled_empty << "/" << soup.culled_full
             << " culled E/F, " << soup.hidden_candidates
             << " hidden candidates, " << soup.dense_blocks
             << " dense");

        REQUIRE(soup.surface.size() > 100);
        bool any_in = false, any_out = false;
        for (const DSample& s : soup.samples)
            (s.inside ? any_in : any_out) = true;
        REQUIRE(any_in);
        REQUIRE(any_out);

        /*  Bisected points sit ON the surface; QEF feature points
         *  (the tail of the array) sit on crease intersections and
         *  carry O(solver) error - judged to a looser bound.  */
        TapeCtx* ctx = tape_ctx_new(d.deck);
        const size_t n_bisected =
                soup.surface.size() - soup.feature_points;
        double worst = 0, worst_feat = 0;
        for (size_t i = 0; i < soup.surface.size(); ++i)
        {
            const DSurfPoint& p = soup.surface[i];
            const float v =
                    tape_eval_f(deck_base(d.deck), ctx, p.x, p.y, p.z);
            double& w = i < n_bisected ? worst : worst_feat;
            if (fabs(v) > w)
                w = fabs(v);
        }
        tape_ctx_free(ctx);
        WARN("worst |f|: bisected " << worst << ", feature "
             << worst_feat << " (" << soup.feature_points
             << " feature pts)");
        CHECK(worst < 1e-3);
        if (soup.feature_points)
            CHECK(worst_feat < 0.05);
    }

    /*  Crease-band density smoke (the STIBIUM_DMESH_DENSE knob is
     *  off by default - keep the machinery referee'd): dense leaves
     *  fire exactly where the pushed tape keeps a live min/max
     *  choice, so the sphere has none and csg has many, and the
     *  band's midpoint lattice must grow the point soup.  */
    {
        DeckRegion d(MODELS[3], 48, 1.1f);
        volatile int halt = 0;
        /*  Isolate the DENSE-knob machinery from the chainless-
         *  curvature trigger: at bench pitch every curved
         *  surface reads several degrees per cell and the
         *  trigger fires legitimately, but this smoke test
         *  referees the knob, not the trigger.  */
        setenv("STIBIUM_DMESH_CURVEBAR", "0", 1);
        const DSoup base = delaunay_sample(d.deck, d.r, &halt);
        setenv("STIBIUM_DMESH_DENSE", "1", 1);
        const DSoup dense = delaunay_sample(d.deck, d.r, &halt);
        DeckRegion ds(MODELS[0], 48, 1.1f);
        const DSoup sphere = delaunay_sample(ds.deck, ds.r, &halt);
        unsetenv("STIBIUM_DMESH_DENSE");
        unsetenv("STIBIUM_DMESH_CURVEBAR");
        WARN("dense smoke: csg " << base.dense_blocks << " -> "
             << dense.dense_blocks << " dense blocks, "
             << base.samples.size() << " -> " << dense.samples.size()
             << " samples; sphere " << sphere.dense_blocks
             << " dense blocks");
        CHECK(base.dense_blocks == 0);
        CHECK(dense.dense_blocks > 0);
        CHECK(dense.samples.size() > base.samples.size());
        CHECK(dense.spacing == base.spacing);   // base pitch kept
        CHECK(sphere.dense_blocks == 0);        // no min/max clauses
    }
}

TEST_CASE("Delaunay refinement rebuilds a stripped surface",
          "[.dmeshBC]")
{
    /*  The mechanism test for Keeter step 5: hand stage B a soup
     *  with EVERY surface point removed.  All separation must then
     *  come from the refinement loop bisecting direct inside<->
     *  outside tet edges - and the result must still be closed and
     *  volumetrically sane.  */
    DeckRegion d("-r++qXqYqZf1", 24, 1.6f);
    volatile int halt = 0;
    DSoup soup = delaunay_sample(d.deck, d.r, &halt);
    const size_t had = soup.surface.size();
    soup.surface.clear();

    DMesh m;
    if (!delaunay_mesh_soup(d.deck, soup, &halt, &m))
    {
        WARN("delaunay_mesh unavailable (built without CGAL?)");
        return;
    }
    WARN("stripped sphere: dropped " << had << " surface pts; "
         << m.iterations << " refine rounds rebuilt " << m.inserted
         << " points; " << m.tris.size() / 3 << " tris, "
         << m.open_edges << " open edges, "
         << m.nonmanifold_edges << " non-manifold");
    CHECK(m.iterations > 0);
    CHECK(m.inserted > 100);
    CHECK(m.open_edges == 0);
    const double vol = mesh_volume(m);
    const double rel_err = fabs(vol - 4.18879) / 4.18879;
    WARN("stripped sphere volume " << vol << " rel err " << rel_err);
    CHECK(rel_err < 0.05);
}

TEST_CASE("Delaunay stages B+C: closed mesh, correct volume",
          "[.dmeshBC]")
{
    struct Case
    {
        const char* expr;
        double volume;   // analytic
    };
    const Case CASES[] = {
        /*  sphere r=1: 4/3 pi  */
        { "-r++qXqYqZf1", 4.18879 },
        /*  cube [-0.6,0.6]^3  */
        { "aa-f-0.6X-Xf0.6aa-f-0.6Y-Yf0.6a-f-0.6Z-Zf0.6", 1.728 },
        /*  union of two spheres r=1 and r=0.8 at (0.5,0.25,0.1):
            overlapping - no closed form here; volume checked > max
            of either and < sum  */
        { "i-r++qXqYqZf1-r++q-Xf0.5q-Yf0.25q-Zf0.1f0.8", -1 },
    };
    for (const Case& tc : CASES)
    {
        CAPTURE(tc.expr);
        DeckRegion d(tc.expr, 48, 1.6f);
        volatile int halt = 0;
        DMesh m;
        const bool ok = delaunay_mesh(d.deck, d.r, &halt, &m);
        if (!ok)
        {
            WARN("delaunay_mesh unavailable (built without CGAL?)");
            return;
        }

        const double vol = mesh_volume(m);
        WARN("model " << tc.expr << ": "
             << m.verts.size() / 3 << " verts, "
             << m.tris.size() / 3 << " tris, "
             << m.iterations << " refine rounds, "
             << m.inserted << " inserted, "
             << m.open_edges << " open edges, "
             << m.nonmanifold_edges << " non-manifold edges, "
             << "volume " << vol);

        REQUIRE(m.tris.size() > 100 * 3);
        CHECK(m.open_edges == 0);
        CHECK(m.iterations < 48);
        if (tc.volume > 0)
        {
            const double rel_err = fabs(vol - tc.volume) / tc.volume;
            CHECK(rel_err < 0.02);
        }
        else
        {
            const double lo = 4.18879;            // sphere r=1
            const double hi = 4.18879 + 2.14466;  // sum of both
            CHECK(vol > lo);
            CHECK(vol < hi);
        }
    }
}

TEST_CASE("Delaunay vs Manifold DC: head to head", "[.dmeshVS]")
{
    /*  Same models, same regions: runtime, triangle count, and the
     *  structural properties each can promise.  STLs land in the
     *  working directory for human eyeballs.  */
    struct Case
    {
        const char* name;
        const char* expr;
    };
    const Case CASES[] = {
        { "sphere",   "-r++qXqYqZf1" },
        { "cube",     "aa-f-0.6X-Xf0.6aa-f-0.6Y-Yf0.6a-f-0.6Z-Zf0.6" },
        { "csg",      "ai-r++qXqYqZf1-r++q-Xf0.5qYqZf0.7n-Zf0.2" },
        { "spheres",  "i-r++qXqYqZf1-r++q-Xf0.5q-Yf0.25q-Zf0.1f0.8" },
    };
    using clk = std::chrono::steady_clock;
    for (const Case& tc : CASES)
    {
        CAPTURE(tc.name);
        DeckRegion d(tc.expr, 64, 1.6f);
        volatile int halt = 0;

        auto t0 = clk::now();
        DMesh m;
        if (!delaunay_mesh(d.deck, d.r, &halt, &m))
        {
            WARN("delaunay unavailable");
            return;
        }
        auto t1 = clk::now();

        /*  The incumbent (plain Manifold DC, detect_features off),
         *  through its own public entry.  */
        MathTree* tree = parse(tc.expr);
        REQUIRE(tree != nullptr);
        float* dc_verts = nullptr;
        unsigned dc_count = 0;
        auto t2 = clk::now();
        triangulate(tree, d.r, false, &halt, &dc_verts, &dc_count);
        auto t3 = clk::now();
        free_tree(tree);

        const auto ms = [](clk::time_point a, clk::time_point b) {
            return std::chrono::duration<double, std::milli>(b - a)
                    .count();
        };
        const uint32_t comps = mesh_components(m);
        WARN("folds(dot<-0.2): " << mesh_warts(m)
             << "; repaired " << m.repaired << " in "
             << m.repair_rounds << " rounds");
        WARN(tc.name << ": delaunay " << m.tris.size() / 3
             << " tris in " << ms(t0, t1) << " ms ("
             << m.open_edges << " open, " << m.nonmanifold_edges
             << " non-manifold, " << comps
             << " components) | manifold-DC " << dc_count / 9
             << " tris in " << ms(t2, t3) << " ms");

        char fname[128];
        snprintf(fname, sizeof(fname), "dmesh_%s.stl", tc.name);
        save_stl_indexed(m.verts.data(), m.tris.data(),
                         uint32_t(m.tris.size() / 3), fname);
        snprintf(fname, sizeof(fname), "dcref_%s.stl", tc.name);
        save_stl(dc_verts, dc_count, fname);   // count = FLOATS
        free(dc_verts);
    }
}

TEST_CASE("Thin feature below the lattice: the stage-D case",
          "[.dmeshVS]")
{
    /*  A plate 0.01 thick in x - far below the 24^3 lattice pitch
     *  (~0.13).  Point sampling CANNOT see it, for us or for
     *  Manifold DC; the difference is that our interval oracle
     *  KNOWS it might be there: the hidden-candidate counter is the
     *  stage-D drill-down trigger, measured here from day one.  */
    /*  (Centered off-lattice on purpose: a plate centered at 0 gets
     *  sampled by the lattice plane at exactly x = 0.)  */
    DeckRegion d("aa-f0.032X-Xf0.042aa-f-0.6Y-Yf0.6a-f-0.6Z-Zf0.6",
                 24, 1.6f);
    volatile int halt = 0;
    const DSoup soup = delaunay_sample(d.deck, d.r, &halt);
    WARN("thin plate at coarse lattice: " << soup.surface.size()
         << " surface pts, " << soup.hidden_candidates
         << " hidden-feature candidate blocks of "
         << soup.leaf_blocks);
    /*  Day one this asserted surface.empty(): point sampling
     *  cannot see the plate, only the interval bound can, and
     *  the candidate counter was the unexploited trigger.  The
     *  certify arc (2026-07-19) closed the loop: candidates are
     *  interval-certified (proven opposite-sign sub-box) and
     *  drill-down RESOLVES the plate - the once-invisible
     *  feature now arrives as real surface.  */
    CHECK(!soup.surface.empty());
}

TEST_CASE("Delaunay showcase STLs at higher resolution", "[.dmeshSTL]")
{
    /*  One-off generation for human eyeballs - resolution via
     *  STIBIUM_DMESH_N (default 96).  */
    uint32_t n = 96;
    if (const char* env = getenv("STIBIUM_DMESH_N"))
        if (atoi(env) > 0)
            n = uint32_t(atoi(env));
    const struct { const char* name; const char* expr; } CASES[] = {
        { "sphere",  "-r++qXqYqZf1" },
        { "cube",    "aa-f-0.6135X-Xf0.6135aa-f-0.6135Y-Yf0.6135a-f-0.6135Z-Zf0.6135" },
        { "cube_aligned",
                     "aa-f-0.6X-Xf0.6aa-f-0.6Y-Yf0.6a-f-0.6Z-Zf0.6" },
        { "csg",     "ai-r++qXqYqZf1-r++q-Xf0.5qYqZf0.7n-Zf0.2" },
        { "spheres", "i-r++qXqYqZf1-r++q-Xf0.5q-Yf0.25q-Zf0.1f0.8" },
    };
    for (const auto& tc : CASES)
    {
        DeckRegion d(tc.expr, n, 1.6f);
        volatile int halt = 0;
        DMesh m;
        if (!delaunay_mesh(d.deck, d.r, &halt, &m))
            return;
        char fname[128];
        snprintf(fname, sizeof(fname), "dmesh_hi_%s.stl", tc.name);
        save_stl_indexed(m.verts.data(), m.tris.data(),
                         uint32_t(m.tris.size() / 3), fname);
        WARN(tc.name << " @" << n << ": " << m.tris.size() / 3
             << " tris, " << m.open_edges << " open, "
             << m.nonmanifold_edges << " non-manifold, "
             << mesh_components(m) << " components; folds sharp/"
             << "medium/gentle: " << mesh_warts(m, -0.2f) << "/"
             << mesh_warts(m, 0.2f) << "/" << mesh_warts(m, 0.5f)
             << "; repaired " << m.repaired << " in "
             << m.repair_rounds << " rounds; constrained "
             << m.constrained << " edges, " << m.steiner
             << " steiner, " << m.split_verts
             << " split, " << m.snapped << " snapped -> " << fname);
    }
}

TEST_CASE("Crease tracer: exact curves with junctions", "[.dtrace]")
{
    /*  The SSI predictor-corrector tracer judged against models
     *  whose crease topology is known a priori, and against the
     *  oracle: every traced point must BE on the surface.  */
    struct Case
    {
        const char* name;
        const char* expr;
        int min_open;
        int min_closed;
    };
    const Case CASES[] = {
        /*  de-aligned cube: 4 pillar edges (one min/max pair) +
         *  2 closed face loops (the outer pair) = all 12 edges.
         *  (Under STIBIUM_DMESH_DENSE the loops arrive instead as
         *  corner-shared open arcs - 12 edge arcs, same wireframe,
         *  same constraint graph.)  */
        { "cube", "aa-f-0.6135X-Xf0.6135aa-f-0.6135Y-Yf0.6135a-f-0.6135Z-Zf0.6135",
          4, 2 },
        /*  union: one closed seam loop, no junctions  */
        { "union", "i-r++qXqYqZf1-r++q-Xf0.5q-Yf0.25q-Zf0.1f0.8",
          0, 1 },
        /*  csg: trimmed seam arc + the closed cut-face boundary
         *  (traced through its two branch-switch kinks)  */
        { "csg", "ai-r++qXqYqZf1-r++q-Xf0.5qYqZf0.7n-Zf0.2",
          1, 1 },
        /*  tent prism (abs round, 2026-07-15): the ridge
         *  {x = 0, z = 0.8} is an ABS crease (the surface kinks
         *  inside abs(x) with the clause active on both sides) -
         *  invisible to the min/max pair machinery, traced by the
         *  abs generator.  (The two 135-degree base creases sit
         *  below the QEF feature threshold, so the pair tracer
         *  has no seeds there - a seed-coverage limit, not an abs
         *  regression.)  */
        { "tent", "a-+bXZf0.8nZ",
          1, 0 },
        /*  extruded cylinder: rims are pair creases (2 closed
         *  loops); the abs(z) clause kinks at z = 0 where the
         *  CIRCLE branch owns the surface - the kink-activity
         *  gate must refuse that phantom equator.  */
        { "extrude", "a-r+qXqYf1-bZf0.6",
          0, 2 },
    };
    for (const Case& tc : CASES)
    {
        CAPTURE(tc.name);
        DeckRegion d(tc.expr, 64, 1.6f);
        volatile int halt = 0;
        DSoup soup = delaunay_sample(d.deck, d.r, &halt);
        REQUIRE(delaunay_trace(d.deck, d.r, &soup, &halt));

        int open = 0, closed = 0;
        size_t pts = 0;
        for (size_t c = 0; c < soup.tchains.size(); ++c)
        {
            (soup.tclosed[c] ? closed : open)++;
            pts += soup.tchains[c].size();
        }
        WARN(tc.name << ": " << soup.tchains.size()
             << " traced curves (" << open << " open, " << closed
             << " closed), " << pts << " points");
        CHECK(open >= tc.min_open);
        CHECK(closed >= tc.min_closed);

        /*  Every traced point on the surface, per the oracle  */
        TapeCtx* ctx = tape_ctx_new(d.deck);
        float worst = 0;
        for (const auto& chain : soup.tchains)
            for (const uint32_t idx : chain)
            {
                const DSurfPoint& p = soup.surface[idx];
                worst = std::max(worst, std::fabs(tape_eval_f(
                        deck_base(d.deck), ctx, p.x, p.y, p.z)));
            }
        tape_ctx_free(ctx);
        WARN(tc.name << ": worst traced |f| = " << worst);
        /*  Smooth-crease points converge to float exactness (the
         *  union measures 1e-7); the bound is set by points
         *  STRADDLING a branch-switch kink (csg's cut loop crosses
         *  the seam), where central differences mix the two
         *  branches' gradients and Newton converges loosely -
         *  measured 7.8e-4, still under the ~1.3e-3 QEF feature
         *  noise the pipeline already tolerates.  */
        CHECK(worst < 1.5e-3f);
    }
}

TEST_CASE("Crease chains: known topology recovered", "[.dchain]")
{
    /*  The chain extractor judged against models whose crease
     *  topology is known a priori.  */
    struct Case
    {
        const char* name;
        const char* expr;
        int min_loops;      // closed chains expected (at least)
        int min_open;       // open chains expected (at least)
    };
    const Case CASES[] = {
        /*  de-aligned cube: 12 edges between 8 corner junctions  */
        { "cube", "aa-f-0.6135X-Xf0.6135aa-f-0.6135Y-Yf0.6135a-f-0.6135Z-Zf0.6135",
          0, 8 },
        /*  two-sphere union: one closed intersection circle  */
        { "union", "i-r++qXqYqZf1-r++q-Xf0.5q-Yf0.25q-Zf0.1f0.8",
          1, 0 },
        /*  csg: sphere cut by plane - at least one closed crease  */
        { "csg", "ai-r++qXqYqZf1-r++q-Xf0.5qYqZf0.7n-Zf0.2",
          1, 0 },
    };
    /*  The fallback extractor's contract is BASE-density features:
     *  the crease-band dense lattice (2026-07-15) halves feature
     *  spacing and its covariance junction-split dissolves (cube
     *  corners sail through).  Production is protected by the
     *  oracle trust gate (>10% rejected segments -> DT semantics);
     *  dense-input robustness is queued in MESH-NEXT.  */
    setenv("STIBIUM_DMESH_DENSE", "0", 1);
    /*  Same contract for stage-D: the phantom-oracle flags (curved
     *  creases report their plane-fit sagitta) would otherwise
     *  drill csg's crease band dense.  */
    setenv("STIBIUM_DMESH_AUTODENSE", "0", 1);
    for (const Case& tc : CASES)
    {
        CAPTURE(tc.name);
        DeckRegion d(tc.expr, 64, 1.6f);
        volatile int halt = 0;
        const DSoup soup = delaunay_sample(d.deck, d.r, &halt);
        const DChains ch = delaunay_chains(soup);

        int loops = 0, open = 0;
        size_t chained = 0, longest = 0;
        for (size_t c = 0; c < ch.chains.size(); ++c)
        {
            (ch.closed[c] ? loops : open)++;
            chained += ch.chains[c].size();
            longest = std::max(longest, ch.chains[c].size());
        }
        WARN(tc.name << ": " << soup.feature_points << " features ("
             << ch.reps << " after merge) -> " << loops
             << " loops + " << open << " open chains ("
             << ch.stray << " stray, longest " << longest
             << ", chained " << chained << ")");
        CHECK(loops >= tc.min_loops);
        CHECK(open >= tc.min_open);
        /*  The chains should account for most REPRESENTATIVES
         *  (curved creases clump duplicates; merge collapses them) */
        if (ch.reps > 0)
        {
            const double frac = double(chained) / double(ch.reps);
            CHECK(frac > 0.7);
        }
    }
    unsetenv("STIBIUM_DMESH_DENSE");
    unsetenv("STIBIUM_DMESH_AUTODENSE");
}

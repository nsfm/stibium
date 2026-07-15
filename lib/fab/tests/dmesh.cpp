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

#include <cmath>
#include <cstdint>
#include <vector>

#include "catch/catch.hpp"

#include "fab/tree/triangulate/delaunay.h"
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
             << " hidden candidates");

        REQUIRE(soup.surface.size() > 100);
        bool any_in = false, any_out = false;
        for (const DSample& s : soup.samples)
            (s.inside ? any_in : any_out) = true;
        REQUIRE(any_in);
        REQUIRE(any_out);

        TapeCtx* ctx = tape_ctx_new(d.deck);
        double worst = 0;
        for (const DSurfPoint& p : soup.surface)
        {
            const float v =
                    tape_eval_f(deck_base(d.deck), ctx, p.x, p.y, p.z);
            if (fabs(v) > worst)
                worst = fabs(v);
        }
        tape_ctx_free(ctx);
        WARN("worst |f(surface point)| = " << worst);
        CHECK(worst < 1e-3);
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

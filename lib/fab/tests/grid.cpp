#include <Python.h>
#include <catch/catch.hpp>

#include <cmath>
#include <cstdio>
#include <vector>

#include "fab/tree/tree.h"
#include "fab/tree/parser.h"
#include "fab/tree/eval.h"
#include "fab/tree/grid.h"
#include "fab/tree/math/math_g.h"
#include "fab/tree/node/node.h"
#include "fab/util/region.h"

/*
 *  Engine tests for OP_GRID: a signed-distance grid sampled from the
 *  analytic sphere |p| - 1, evaluated through the ordinary tree
 *  machinery and compared against the analytic field.
 */

namespace {

/*  Registers an N^3 grid of the unit sphere's SDF over [-2, 2]^3.
 *  Returns the grid id. */
uint32_t register_sphere_grid(unsigned n, const std::string& key = "")
{
    std::vector<float> d(size_t(n) * n * n);
    const float lo = -2, hi = 2;
    const float step = (hi - lo) / (n - 1);
    size_t idx = 0;
    for (unsigned k = 0; k < n; ++k) {
        const float z = lo + k * step;
        for (unsigned j = 0; j < n; ++j) {
            const float y = lo + j * step;
            for (unsigned i = 0; i < n; ++i) {
                const float x = lo + i * step;
                d[idx++] = sqrtf(x*x + y*y + z*z) - 1;
            }
        }
    }
    return grid_register(std::move(d), n, n, n,
                         lo, lo, lo, hi, hi, hi, key);
}

}  // namespace

TEST_CASE("Grid registry", "[grid]")
{
    SECTION("Registration and lookup")
    {
        const uint32_t id = register_sphere_grid(9);
        REQUIRE(id != 0);
        MeshGrid* g = grid_lookup(id);
        REQUIRE(g != nullptr);
        REQUIRE(grid_id(g) == id);
        MeshGrid* missing = grid_lookup(999999);
        REQUIRE(missing == nullptr);
    }

    SECTION("Key-based reuse")
    {
        const uint32_t a = register_sphere_grid(9, "grid-test-reuse");
        REQUIRE(grid_find("grid-test-reuse") == a);
        REQUIRE(grid_find("grid-test-missing") == 0);
    }

    SECTION("Bad dimensions rejected")
    {
        std::vector<float> tiny(4, 0.f);
        const uint32_t id = grid_register(std::move(tiny), 2, 2, 1,
                                          0, 0, 0, 1, 1, 1, "");
        REQUIRE(id == 0);
    }

    SECTION("Trim frees unreferenced grids only")
    {
        grid_registry_trim();    // clear leftovers from other sections

        const uint32_t loose = register_sphere_grid(9);
        const uint32_t held = register_sphere_grid(9);
        char math[32];
        snprintf(math, sizeof(math), "g%u", held);
        MathTree* t = parse(math);
        REQUIRE(t != nullptr);

        grid_registry_trim();
        MeshGrid* loose_g = grid_lookup(loose);
        MeshGrid* held_g = grid_lookup(held);
        REQUIRE(loose_g == nullptr);
        REQUIRE(held_g != nullptr);

        free_tree(t);
        grid_registry_trim();
        held_g = grid_lookup(held);
        REQUIRE(held_g == nullptr);
    }
}

TEST_CASE("Grid parsing", "[grid]")
{
    SECTION("g<id> parses to a working tree")
    {
        const uint32_t id = register_sphere_grid(33);
        char math[32];
        snprintf(math, sizeof(math), "g%u", id);
        MathTree* t = parse(math);
        REQUIRE(t != nullptr);
        REQUIRE(t->num_levels == 2);    // X/Y/Z leaves + the grid node
        free_tree(t);
        grid_registry_trim();
    }

    SECTION("Unknown id fails the parse")
    {
        MathTree* t = parse("g999999");
        REQUIRE(t == nullptr);
    }

    SECTION("Duplicate grid nodes dedupe; distinct grids don't")
    {
        const uint32_t a = register_sphere_grid(9);
        const uint32_t b = register_sphere_grid(9);
        char math[64];

        // Same grid twice: one node, so min(g, g) folds to one leaf
        snprintf(math, sizeof(math), "ig%ug%u", a, a);
        MathTree* t = parse(math);
        REQUIRE(t != nullptr);
        const unsigned nodes_same = t->active[1];
        REQUIRE(nodes_same == 1);
        free_tree(t);

        // Different grids: two grid nodes under the min
        snprintf(math, sizeof(math), "ig%ug%u", a, b);
        t = parse(math);
        REQUIRE(t != nullptr);
        const unsigned nodes_diff = t->active[1];
        REQUIRE(nodes_diff == 2);
        free_tree(t);
        grid_registry_trim();
    }
}

TEST_CASE("Grid evaluation matches the analytic field", "[grid]")
{
    const uint32_t id = register_sphere_grid(65);
    char math[32];
    snprintf(math, sizeof(math), "g%u", id);
    MathTree* t = parse(math);
    REQUIRE(t != nullptr);

    SECTION("Point evaluation, inside the sampled box")
    {
        // Trilinear error on a smooth SDF at h = 4/64 stays well
        // under a hundredth here; the curvature term is h^2/8.
        float worst = 0;
        for (float x = -1.9f; x < 2; x += 0.37f) {
            for (float y = -1.9f; y < 2; y += 0.41f) {
                for (float z = -1.9f; z < 2; z += 0.43f) {
                    const float expected = sqrtf(x*x + y*y + z*z) - 1;
                    const float got = eval_f(t, x, y, z);
                    const float err = fabsf(got - expected);
                    if (err > worst)    worst = err;
                }
            }
        }
        CAPTURE(worst);
        REQUIRE(worst < 0.01f);
    }

    SECTION("Point evaluation, outside the sampled box")
    {
        // Beyond the box the field is boundary sample + distance to
        // the box: still positive and growing, so CSG against other
        // shapes can't hallucinate material out there.
        const float near_face = eval_f(t, 2.5f, 0, 0);
        const float far_corner = eval_f(t, 5, 5, 5);
        REQUIRE(near_face > 1.0f);
        REQUIRE(far_corner > near_face);
    }

    SECTION("Interval evaluation is conservative and useful")
    {
        // A box straddling the surface must contain 0
        Interval X = {0.5f, 1.5f}, Y = {-0.25f, 0.25f}, Z = {-0.25f, 0.25f};
        Interval out = eval_i(t, X, Y, Z);
        REQUIRE(out.lower <= 0.f);
        REQUIRE(out.upper >= 0.f);

        // A box fully inside the sphere must be all-negative
        X = {-0.2f, 0.2f};  Y = {-0.2f, 0.2f};  Z = {-0.2f, 0.2f};
        out = eval_i(t, X, Y, Z);
        REQUIRE(out.upper < 0.f);

        // A box fully outside (but inside the sampled bounds)
        X = {1.5f, 1.9f};  Y = {1.5f, 1.9f};  Z = {1.5f, 1.9f};
        out = eval_i(t, X, Y, Z);
        REQUIRE(out.lower > 0.f);

        // A box fully outside the sampled bounds
        X = {3.f, 4.f};  Y = {0.f, 1.f};  Z = {0.f, 1.f};
        out = eval_i(t, X, Y, Z);
        REQUIRE(out.lower > 0.f);
    }

    SECTION("Interval bounds actually contain point samples")
    {
        // Sliding boxes: every point value inside the box must lie
        // within the interval result (the conservative contract that
        // all the pruning machinery relies on).
        for (int trial = 0; trial < 20; ++trial) {
            const float cx = -1.8f + 0.19f * trial;
            const float cy = -1.5f + 0.16f * trial;
            const float cz = 1.7f - 0.17f * trial;
            Interval X = {cx, cx + 0.6f};
            Interval Y = {cy, cy + 0.5f};
            Interval Z = {cz - 0.4f, cz};
            const Interval out = eval_i(t, X, Y, Z);
            for (float fx = 0; fx <= 1; fx += 0.5f) {
                for (float fy = 0; fy <= 1; fy += 0.5f) {
                    for (float fz = 0; fz <= 1; fz += 0.5f) {
                        const float v = eval_f(t,
                                X.lower + fx * (X.upper - X.lower),
                                Y.lower + fy * (Y.upper - Y.lower),
                                Z.lower + fz * (Z.upper - Z.lower));
                        const bool contained =
                                v >= out.lower - 1e-5f &&
                                v <= out.upper + 1e-5f;
                        CAPTURE(trial);
                        CAPTURE(v);
                        REQUIRE(contained);
                    }
                }
            }
        }
    }

    free_tree(t);
    grid_registry_trim();
}

TEST_CASE("Grid nodes compose with transforms", "[grid]")
{
    const uint32_t id = register_sphere_grid(65);
    char math[64];

    SECTION("Translation via map")
    {
        // m -Xf0.5 _ _ g<id>: X' = X - 0.5, so the sphere's surface
        // sits at x = 1.5 along +x
        snprintf(math, sizeof(math), "m-Xf0.5__g%u", id);
        MathTree* t = parse(math);
        REQUIRE(t != nullptr);

        const float at_new_surface = eval_f(t, 1.5f, 0, 0);
        const float at_old_surface = eval_f(t, 1.0f, 0, 0);
        CAPTURE(at_new_surface);
        CAPTURE(at_old_surface);
        REQUIRE(fabsf(at_new_surface) < 0.01f);
        REQUIRE(fabsf(at_old_surface + 0.5f) < 0.01f);
        free_tree(t);
    }

    SECTION("Scale via map")
    {
        // X' = X * 0.5: sphere surface lands at x = 2
        snprintf(math, sizeof(math), "m*Xf0.5__g%u", id);
        MathTree* t = parse(math);
        REQUIRE(t != nullptr);
        const float v = eval_f(t, 2.0f, 0, 0);
        CAPTURE(v);
        REQUIRE(fabsf(v) < 0.01f);
        free_tree(t);
    }

    grid_registry_trim();
}

TEST_CASE("Grid trees clone and share the payload", "[grid]")
{
    const uint32_t id = register_sphere_grid(33);
    char math[32];
    snprintf(math, sizeof(math), "g%u", id);
    MathTree* t = parse(math);
    REQUIRE(t != nullptr);

    MathTree* c = clone_tree(t);
    REQUIRE(c != nullptr);

    // Same values through both trees
    const float a = eval_f(t, 0.3f, -0.2f, 0.7f);
    const float b = eval_f(c, 0.3f, -0.2f, 0.7f);
    REQUIRE(a == b);

    // The grid survives freeing either tree first
    free_tree(t);
    const float after = eval_f(c, 0.3f, -0.2f, 0.7f);
    REQUIRE(after == b);
    free_tree(c);

    grid_registry_trim();
    MeshGrid* gone = grid_lookup(id);
    REQUIRE(gone == nullptr);
}

TEST_CASE("Grid gradient evaluation", "[grid]")
{
    const uint32_t id = register_sphere_grid(65);
    char math[32];
    snprintf(math, sizeof(math), "g%u", id);
    MathTree* t = parse(math);
    REQUIRE(t != nullptr);

    // Gradient of a sphere SDF is the unit radial direction; check
    // through the region-based gradient evaluator.
    Region r = {};
    r.ni = r.nj = r.nk = 1;
    r.voxels = 1;

    float xs[1] = {0.8f}, ys[1] = {0.6f}, zs[1] = {0.f};
    uint16_t Ls[1] = {0};
    r.X = xs;  r.Y = ys;  r.Z = zs;  r.L = Ls;

    derivative* d = eval_g(t, r);
    REQUIRE(d != nullptr);
    CAPTURE(d[0].dx);
    CAPTURE(d[0].dy);
    CAPTURE(d[0].dz);

    // |p| = 1.0 -> value ~0, gradient ~(0.8, 0.6, 0)
    REQUIRE(fabsf(d[0].v) < 0.01f);
    REQUIRE(fabsf(d[0].dx - 0.8f) < 0.05f);
    REQUIRE(fabsf(d[0].dy - 0.6f) < 0.05f);
    REQUIRE(fabsf(d[0].dz) < 0.05f);

    free_tree(t);
    grid_registry_trim();
}

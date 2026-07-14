#include <Python.h>
#include <cstring>
#include <catch/catch.hpp>

#include "fab/fab.h"
#include "fab/tree/tree.h"
#include "fab/tree/parser.h"
#include "fab/tree/eval.h"
#include "fab/tree/tape.h"
#include "fab/tree/node/node.h"
#include "fab/util/switches.h"

/*  Expressions chosen to cover every eval path the renderer and
 *  mesher exercise: min/max CSG (pruning targets), trig, powers,
 *  and plain arithmetic.  (Old-style prefix strings; see parser.c.) */
static const char* CORPUS[] = {
    "-r++qXqYqZf1",                                  // sphere
    "i-r++qXqYqZf1-r++q-Xf0.5q-Yf0.25q-Zf0.1f0.8",   // union of spheres
    "a-r++qXqYqZf1nY",                               // sphere ∩ half-space
    "+s*Xf3.0c*Yf2.0",                               // wavy sheet
    "ai-r++qXqYqZf1-r++q-Xf0.5qYqZf0.7n-Zf0.2",      // nested min/max
    "*+Xf2.0+Yf2.0",                                 // dedup'd arithmetic
};

/*  min(sphere at origin r=1, sphere at x=10 r=0.5): interval eval
 *  over a region near the origin decides the min entirely in favor
 *  of the first sphere, so pushes have something real to prune.  */
static const char* PRUNABLE = "i-r++qXqYqZf1-r++q-Xf10qYqZf0.5";

TEST_CASE("Tape point evaluation matches eval_f", "[tape]")
{
    for (const char* expr : CORPUS)
    {
        MathTree* tree = parse(expr);
        REQUIRE(tree != nullptr);

        Deck* deck = deck_from_tree(tree);
        TapeCtx* ctx = tape_ctx_new(deck);
        Tape* tape = deck_base(deck);

        for (float x = -2; x <= 2; x += 0.8f)
        for (float y = -2; y <= 2; y += 0.8f)
        for (float z = -2; z <= 2; z += 0.8f)
        {
            const float a = eval_f(tree, x, y, z);
            const float b = tape_eval_f(tape, ctx, x, y, z);
            CAPTURE(expr);
            CAPTURE(x);
            CAPTURE(y);
            CAPTURE(z);
            REQUIRE(a == b);
        }

        tape_ctx_free(ctx);
        deck_free(deck);
        free_tree(tree);
    }
}

TEST_CASE("Tape interval evaluation matches eval_i", "[tape]")
{
    const Interval boxes[][3] = {
        { {-1, 1}, {-1, 1}, {-1, 1} },
        { {0.25f, 2}, {-0.5f, 0.5f}, {-2, -0.25f} },
        { {-0.1f, 0.1f}, {-0.1f, 0.1f}, {-0.1f, 0.1f} },
        { {5, 6}, {5, 6}, {5, 6} },
    };

    for (const char* expr : CORPUS)
    {
        MathTree* tree = parse(expr);
        REQUIRE(tree != nullptr);

        Deck* deck = deck_from_tree(tree);
        TapeCtx* ctx = tape_ctx_new(deck);
        Tape* tape = deck_base(deck);

        for (const auto& box : boxes)
        {
            const Interval a = eval_i(tree, box[0], box[1], box[2]);
            const Interval b = tape_eval_i(tape, ctx,
                                           box[0], box[1], box[2]);
            CAPTURE(expr);
            CAPTURE(box[0].lower);
            CAPTURE(box[0].upper);
            REQUIRE(a.lower == b.lower);
            REQUIRE(a.upper == b.upper);
        }

        tape_ctx_free(ctx);
        deck_free(deck);
        free_tree(tree);
    }
}

/*  Fills xs/ys/zs with a flattened grid of n^3 sample points  */
static unsigned sample_grid(float* xs, float* ys, float* zs, int n)
{
    unsigned q = 0;
    for (int k = 0; k < n; ++k)
    for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
    {
        xs[q] = -1.5f + 3.0f * i / (n - 1);
        ys[q] = -1.5f + 3.0f * j / (n - 1);
        zs[q] = -1.5f + 3.0f * k / (n - 1);
        q++;
    }
    return q;
}

TEST_CASE("Tape region evaluation matches eval_r", "[tape]")
{
    float xs[MIN_VOLUME], ys[MIN_VOLUME], zs[MIN_VOLUME];
    const unsigned count = sample_grid(xs, ys, zs, 5);
    REQUIRE(count < MIN_VOLUME);

    Region dummy = {};
    dummy.X = xs;
    dummy.Y = ys;
    dummy.Z = zs;
    dummy.voxels = count;

    for (const char* expr : CORPUS)
    {
        MathTree* tree = parse(expr);
        REQUIRE(tree != nullptr);

        Deck* deck = deck_from_tree(tree);
        TapeCtx* ctx = tape_ctx_new(deck);

        const float* a = eval_r(tree, dummy);
        const float* b = tape_eval_r(deck_base(deck), ctx, dummy);
        CAPTURE(expr);
        REQUIRE(memcmp(a, b, count * sizeof(float)) == 0);

        tape_ctx_free(ctx);
        deck_free(deck);
        free_tree(tree);
    }
}

TEST_CASE("Tape gradient evaluation matches eval_g", "[tape]")
{
    float xs[MIN_VOLUME], ys[MIN_VOLUME], zs[MIN_VOLUME];
    const unsigned count = sample_grid(xs, ys, zs, 3);
    REQUIRE(count <= MIN_VOLUME / 4);

    Region dummy = {};
    dummy.X = xs;
    dummy.Y = ys;
    dummy.Z = zs;
    dummy.voxels = count;

    for (const char* expr : CORPUS)
    {
        MathTree* tree = parse(expr);
        REQUIRE(tree != nullptr);

        // eval_g expects constants in derivative layout (see shaded8)
        for (unsigned i = 0; i < tree->num_constants; ++i)
            fill_results_g(tree->constants[i],
                           tree->constants[i]->results.f);

        Deck* deck = deck_from_tree(tree);
        TapeCtx* ctx = tape_ctx_new(deck);

        const derivative* a = eval_g(tree, dummy);
        const derivative* b = tape_eval_g(deck_base(deck), ctx, dummy);
        CAPTURE(expr);
        REQUIRE(memcmp(a, b, count * sizeof(derivative)) == 0);

        tape_ctx_free(ctx);
        deck_free(deck);
        free_tree(tree);
    }
}

TEST_CASE("Standard push prunes decided min/max and preserves values",
          "[tape]")
{
    MathTree* tree = parse(PRUNABLE);
    REQUIRE(tree != nullptr);

    Deck* deck = deck_from_tree(tree);
    TapeCtx* ctx = tape_ctx_new(deck);
    Tape* base = deck_base(deck);

    const Interval X = {-0.4f, 0.4f}, Y = {-0.4f, 0.4f}, Z = {-0.4f, 0.4f};
    tape_eval_i(base, ctx, X, Y, Z);
    Tape* pushed = tape_push(base, ctx, X, Y, Z, TAPE_PUSH_STANDARD);

    REQUIRE(pushed != base);
    REQUIRE(tape_length(pushed) < tape_length(base));
    // The surviving branch has no min/max left
    REQUIRE(tape_is_terminal(pushed));

    // Pruned evaluation must agree with the base tape and the
    // unpruned tree at points inside the region (pruning is exact)
    for (float x = -0.35f; x <= 0.35f; x += 0.17f)
    for (float y = -0.35f; y <= 0.35f; y += 0.17f)
    for (float z = -0.35f; z <= 0.35f; z += 0.17f)
    {
        const float t = tape_eval_f(pushed, ctx, x, y, z);
        const float b = tape_eval_f(base, ctx, x, y, z);
        const float d = eval_f(tree, x, y, z);
        CAPTURE(x);
        CAPTURE(y);
        CAPTURE(z);
        REQUIRE(t == b);
        REQUIRE(t == d);
    }

    tape_release(pushed);

    tape_ctx_free(ctx);
    deck_free(deck);
    free_tree(tree);
}

TEST_CASE("Terminal and unchanged pushes return the same tape", "[tape]")
{
    // No min/max at all: base tape is terminal
    MathTree* tree = parse("-r++qXqYqZf1");
    REQUIRE(tree != nullptr);
    Deck* deck = deck_from_tree(tree);
    TapeCtx* ctx = tape_ctx_new(deck);
    Tape* base = deck_base(deck);

    REQUIRE(tape_is_terminal(base));
    const Interval I = {-1, 1};
    tape_eval_i(base, ctx, I, I, I);
    Tape* pushed = tape_push(base, ctx, I, I, I, TAPE_PUSH_STANDARD);
    REQUIRE(pushed == base);
    tape_release(pushed);

    tape_ctx_free(ctx);
    deck_free(deck);
    free_tree(tree);

    // min/max present but undecidable over the region: same tape back
    tree = parse("iXY");
    REQUIRE(tree != nullptr);
    deck = deck_from_tree(tree);
    ctx = tape_ctx_new(deck);
    base = deck_base(deck);

    tape_eval_i(base, ctx, I, I, I);
    pushed = tape_push(base, ctx, I, I, I, TAPE_PUSH_STANDARD);
    REQUIRE(pushed == base);
    tape_release(pushed);

    tape_ctx_free(ctx);
    deck_free(deck);
    free_tree(tree);
}

TEST_CASE("Binary push collapses sign-determined regions", "[tape]")
{
    MathTree* tree = parse(PRUNABLE);
    REQUIRE(tree != nullptr);
    Deck* deck = deck_from_tree(tree);
    TapeCtx* ctx = tape_ctx_new(deck);
    Tape* base = deck_base(deck);

    // Entirely inside the origin sphere: the root's sign is fixed,
    // so in boolean context the whole tape collapses to one constant
    const Interval X = {-0.3f, 0.3f}, Y = {-0.3f, 0.3f}, Z = {-0.3f, 0.3f};
    const Interval before = tape_eval_i(base, ctx, X, Y, Z);
    REQUIRE(before.upper < 0);

    Tape* pushed = tape_push(base, ctx, X, Y, Z, TAPE_PUSH_BINARY);
    REQUIRE(pushed != base);
    REQUIRE(tape_length(pushed) == 1);

    // Signs (all a raster consumer reads) must be preserved
    for (float x = -0.25f; x <= 0.25f; x += 0.1f)
    for (float y = -0.25f; y <= 0.25f; y += 0.1f)
    {
        const float t = tape_eval_f(pushed, ctx, x, y, 0.1f);
        const float b = tape_eval_f(base, ctx, x, y, 0.1f);
        CAPTURE(x);
        CAPTURE(y);
        REQUIRE((t < 0) == (b < 0));
    }

    tape_release(pushed);
    tape_ctx_free(ctx);
    deck_free(deck);
    free_tree(tree);
}

TEST_CASE("Nested pushes and base walk-up", "[tape]")
{
    MathTree* tree = parse(PRUNABLE);
    REQUIRE(tree != nullptr);
    Deck* deck = deck_from_tree(tree);
    TapeCtx* ctx = tape_ctx_new(deck);
    Tape* base = deck_base(deck);

    const Interval A = {-0.4f, 0.4f};
    tape_eval_i(base, ctx, A, A, A);
    Tape* ta = tape_push(base, ctx, A, A, A, TAPE_PUSH_STANDARD);
    REQUIRE(ta != base);

    const Interval B = {-0.1f, 0.1f};
    tape_eval_i(ta, ctx, B, B, B);
    Tape* tb = tape_push(ta, ctx, B, B, B, TAPE_PUSH_STANDARD);

    // Inside A: stays on a tape whose recorded bounds contain it
    Tape* found = tape_base_for_point(tb, 0.3f, 0.3f, 0.3f);
    REQUIRE((found == ta || found == tb));

    // Far outside everything: walks all the way to the base
    found = tape_base_for_point(tb, 5, 5, 5);
    REQUIRE(found == base);

    // Values through the whole stack agree at a point inside B
    const float va = tape_eval_f(base, ctx, 0.05f, 0.05f, 0.05f);
    const float vb = tape_eval_f(ta, ctx, 0.05f, 0.05f, 0.05f);
    const float vc = tape_eval_f(tb, ctx, 0.05f, 0.05f, 0.05f);
    REQUIRE(va == vb);
    REQUIRE(va == vc);

    tape_release(tb);
    tape_release(ta);
    tape_ctx_free(ctx);
    deck_free(deck);
    free_tree(tree);
}

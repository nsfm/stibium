#include <catch/catch.hpp>

#include <cmath>

#include "fab/tree/analytics.h"
#include "fab/tree/parser.h"
#include "fab/tree/tree.h"

static FieldStats analyze(const char* math, bool flat=false)
{
    MathTree* tree = parse(math);
    REQUIRE(tree != nullptr);
    volatile int halt = 0;
    FieldStats stats;
    REQUIRE(analyze_field(tree, -1, -1, flat ? 0 : -1,
                          1, 1, flat ? 0 : 1,
                          100, flat, -1, &halt, &stats));
    free_tree(tree);
    return stats;
}

TEST_CASE("Analytics: sphere volume and center")
{
    const auto s = analyze("-r++qXqYqZf0.8");
    REQUIRE(s.volume == Approx(4 * M_PI / 3 * 0.8 * 0.8 * 0.8)
            .epsilon(0.01));
    REQUIRE(std::fabs(s.com[0]) < 0.01);
    REQUIRE(std::fabs(s.com[1]) < 0.01);
    REQUIRE(std::fabs(s.com[2]) < 0.01);
    REQUIRE(std::fabs(s.tight[0] + 0.8) < 0.02);
    REQUIRE(std::fabs(s.tight[3] - 0.8) < 0.02);
}

TEST_CASE("Analytics: offset shape moves the center of mass")
{
    // Sphere at (0.3, 0, 0), radius 0.5
    const auto s = analyze("-r++q-Xf0.3qYqZf0.5");
    REQUIRE(std::fabs(s.com[0] - 0.3) < 0.01);
    REQUIRE(std::fabs(s.com[1]) < 0.01);
    REQUIRE(s.volume == Approx(4 * M_PI / 3 * 0.125).epsilon(0.015));
}

TEST_CASE("Analytics: cube volume and tight bounds")
{
    const auto s = analyze("aa-f-0.6X-Xf0.6aa-f-0.6Y-Yf0.6a-f-0.6Z-Zf0.6");
    REQUIRE(s.volume == Approx(1.2 * 1.2 * 1.2).epsilon(0.01));
    for (int i = 0; i < 3; ++i)
    {
        REQUIRE(std::fabs(s.tight[i] + 0.6) < 0.02);
        REQUIRE(std::fabs(s.tight[i + 3] - 0.6) < 0.02);
    }
}

TEST_CASE("Analytics: flat shapes report area")
{
    const auto s = analyze("-r+qXqYf0.8", true);
    REQUIRE(s.volume == Approx(M_PI * 0.8 * 0.8).epsilon(0.01));
    REQUIRE(s.com[2] == 0.0);
}

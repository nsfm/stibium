#include <Python.h>
#include <catch/catch.hpp>

#include "fab/fab.h"
#include "fab/tree/tree.h"
#include "fab/tree/parser.h"

TEST_CASE("Basic parsing")
{
    MathTree* t;

    SECTION("Parsing 'X'")
    {
        t = parse("X");
        REQUIRE(t != nullptr);
        REQUIRE(t->num_levels == 1);
        REQUIRE(t->num_constants == 0);
        free(t);
    }

    SECTION("Parsing '+Xf1.0'")
    {
        t = parse("+Xf1.0");
        REQUIRE(t != nullptr);
        REQUIRE(t->num_constants == 1);
        REQUIRE(t->num_levels == 2);
        free(t);
    }

    SECTION("Duplicate pruning")
    {
        t = parse("*+Xf2.0+Yf2.0");
        REQUIRE(t != nullptr);
        REQUIRE(t->num_constants == 1);
        REQUIRE(t->num_levels == 3);
        free(t);
    }

    SECTION("Invalid parse")
    {
        t = parse("f");
        REQUIRE(t == nullptr);
        free(t);
    }

    SECTION("Parsing an empty string")
    {
        t = parse("");
        REQUIRE(t == nullptr);
        free(t);
    }

    SECTION("Incomplete expression")
    {
        t = parse("+X");
        REQUIRE(t == nullptr);
        free(t);
    }

    SECTION("Long expression")
    {
        std::string s;
        for (int i=0; i < 99; ++i)
            s += "i";
        for (int i=0; i < 100; ++i)
            s += "X";

        t = parse(s.c_str());
        REQUIRE(t != nullptr);
    }
}

#include "fab/tree/eval.h"

TEST_CASE("mod and floor opcodes")
{
    MathTree* t;

    SECTION("Prefix mod: MXf3.0")
    {
        t = parse("MXf3.0");
        REQUIRE(t != nullptr);
        REQUIRE(eval_f(t, 7.5, 0, 0) == Approx(1.5));
        REQUIRE(eval_f(t, -0.5, 0, 0) == Approx(2.5));
        free(t);
    }

    SECTION("Prefix floor: FX")
    {
        t = parse("FX");
        REQUIRE(t != nullptr);
        REQUIRE(eval_f(t, 2.75, 0, 0) == Approx(2.0));
        REQUIRE(eval_f(t, -0.25, 0, 0) == Approx(-1.0));
        free(t);
    }

    SECTION("Infix mod and floor")
    {
        t = parse("=mod(X, 3.0) + floor(Y);");
        REQUIRE(t != nullptr);
        REQUIRE(eval_f(t, 7.5, 1.9, 0) == Approx(2.5));
        free(t);
    }

    SECTION("Interval mod within one period")
    {
        t = parse("MXf3.0");
        REQUIRE(t != nullptr);
        Interval x = {0.5, 1.0}, y = {0, 0}, z = {0, 0};
        Interval out = eval_i(t, x, y, z);
        REQUIRE(out.lower == Approx(0.5));
        REQUIRE(out.upper == Approx(1.0));
        free(t);
    }

    SECTION("Interval mod across a period boundary")
    {
        t = parse("MXf3.0");
        REQUIRE(t != nullptr);
        Interval x = {2.5, 3.5}, y = {0, 0}, z = {0, 0};
        Interval out = eval_i(t, x, y, z);
        REQUIRE(out.lower <= 0.0f);
        REQUIRE(out.upper >= 3.0f - 1e-6f);
        free(t);
    }
}

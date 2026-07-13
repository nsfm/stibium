#include <catch/catch.hpp>

#include <cmath>
#include <fstream>
#include <string>

#include "fab/formats/svg.h"
#include "fab/tree/contour.h"
#include "fab/tree/parser.h"
#include "fab/tree/tree.h"

namespace {

std::vector<ContourPath> trace(const char* math, bool detect,
                               uint32_t n = 100)
{
    MathTree* tree = parse(math);
    REQUIRE(tree != nullptr);
    volatile int halt = 0;
    std::vector<ContourPath> paths;
    contour_field(tree, -1, -1, 1, 1, n, n, 0, detect, &halt, paths);
    free_tree(tree);
    return paths;
}

// Shoelace formula: positive for counterclockwise loops
double signed_area(const ContourPath& p)
{
    double area = 0;
    for (size_t i = 0; i < p.size(); ++i)
    {
        const auto& a = p[i];
        const auto& b = p[(i + 1) % p.size()];
        area += double(a[0]) * b[1] - double(b[0]) * a[1];
    }
    return area / 2;
}

double perimeter(const ContourPath& p)
{
    double len = 0;
    for (size_t i = 0; i < p.size(); ++i)
    {
        const auto& a = p[i];
        const auto& b = p[(i + 1) % p.size()];
        len += std::hypot(double(b[0]) - a[0], double(b[1]) - a[1]);
    }
    return len;
}

bool has_point_near(const std::vector<ContourPath>& paths,
                    float x, float y, float tol)
{
    for (const auto& p : paths)
        for (const auto& pt : p)
            if (std::hypot(pt[0] - x, pt[1] - y) < tol)
                return true;
    return false;
}

// 2D test shapes (prefix math strings)
const char* CIRCLE = "-r+qXqYf0.8";                       // radius 0.8
const char* SQUARE = "aa-f-0.6X-Xf0.6a-f-0.6Y-Yf0.6";     // [-0.6,0.6]^2
const char* ANNULUS = "a-r+qXqYf0.8-f0.4r+qXqY";          // r in [0.4,0.8]

}  // namespace

TEST_CASE("Contour: circle")
{
    auto paths = trace(CIRCLE, false);
    REQUIRE(paths.size() == 1);
    REQUIRE(signed_area(paths[0]) ==
            Approx(M_PI * 0.8 * 0.8).epsilon(0.01));
    REQUIRE(perimeter(paths[0]) == Approx(2 * M_PI * 0.8).epsilon(0.01));
    // Outer boundary winds counterclockwise (positive area)
    REQUIRE(signed_area(paths[0]) > 0);
}

TEST_CASE("Contour: square with sharp corners")
{
    auto plain = trace(SQUARE, false);
    REQUIRE(plain.size() == 1);
    REQUIRE(signed_area(plain[0]) == Approx(1.2 * 1.2).epsilon(0.01));

    auto sharp = trace(SQUARE, true);
    REQUIRE(sharp.size() == 1);
    REQUIRE(signed_area(sharp[0]) == Approx(1.2 * 1.2).epsilon(0.002));

    // Feature detection must recover the true corners
    for (float sx : {-0.6f, 0.6f})
        for (float sy : {-0.6f, 0.6f})
        {
            CAPTURE(sx);
            CAPTURE(sy);
            REQUIRE(has_point_near(sharp, sx, sy, 0.005f));
        }
}

TEST_CASE("Contour: annulus hole orientation")
{
    auto paths = trace(ANNULUS, false);
    REQUIRE(paths.size() == 2);

    // One CCW outer loop, one CW hole
    const double a0 = signed_area(paths[0]);
    const double a1 = signed_area(paths[1]);
    const double outer = std::max(a0, a1);
    const double hole = std::min(a0, a1);
    REQUIRE(outer == Approx(M_PI * 0.8 * 0.8).epsilon(0.01));
    REQUIRE(hole == Approx(-M_PI * 0.4 * 0.4).epsilon(0.01));

    // Net filled area
    const double net = outer + hole;
    REQUIRE(net == Approx(M_PI * (0.8 * 0.8 - 0.4 * 0.4)).epsilon(0.01));
}

TEST_CASE("Contour: clipped shape closes along the bounds")
{
    // Circle bigger than the sample window: the contour must still be
    // a closed loop, traced along the window border.
    auto paths = trace("-r+qXqYf1.5", false);
    REQUIRE(paths.size() == 1);
    // Filled area is the whole window (to within a cell)
    REQUIRE(signed_area(paths[0]) == Approx(4.0).epsilon(0.02));
}

TEST_CASE("SVG writer: valid document")
{
    auto paths = trace(ANNULUS, true);
    REQUIRE(save_svg(paths, -1, -1, 1, 1, "test_out.svg"));

    std::ifstream f("test_out.svg");
    std::string data((std::istreambuf_iterator<char>(f)),
                      std::istreambuf_iterator<char>());
    REQUIRE(data.find("<svg") != std::string::npos);
    REQUIRE(data.find("viewBox=\"-1 -1 2 2\"") != std::string::npos);
    REQUIRE(data.find("fill-rule=\"evenodd\"") != std::string::npos);
    // Two loops -> two closed subpaths
    REQUIRE(std::count(data.begin(), data.end(), 'M') >= 2);
    REQUIRE(std::count(data.begin(), data.end(), 'Z') == 2);
    std::remove("test_out.svg");
}

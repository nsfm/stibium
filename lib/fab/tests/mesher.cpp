#include <catch/catch.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iterator>
#include <map>
#include <string>
#include <vector>

#include "fab/formats/stl.h"
#include "fab/tree/parser.h"
#include "fab/tree/tree.h"
#include "fab/tree/triangulate.h"
#include "fab/util/region.h"

namespace {

// Flat triangle soup: 9 floats per triangle (xyz * 3 corners).
using Soup = std::vector<float>;

Soup mesh_shape(const char* math,
                float xmin, float ymin, float zmin,
                float xmax, float ymax, float zmax,
                float res, bool detect_features)
{
    MathTree* tree = parse(math);
    REQUIRE(tree != nullptr);

    Region r = {};
    r.ni = uint32_t((xmax - xmin) * res);
    r.nj = uint32_t((ymax - ymin) * res);
    r.nk = uint32_t((zmax - zmin) * res);
    r.voxels = uint64_t(r.ni) * r.nj * r.nk;
    build_arrays(&r, xmin, ymin, zmin, xmax, ymax, zmax);

    volatile int halt = 0;
    float* verts = nullptr;
    unsigned count = 0;
    triangulate(tree, r, detect_features, &halt, &verts, &count);

    Soup s(verts, verts + count);
    free(verts);
    free_arrays(&r);
    free_tree(tree);
    return s;
}

size_t tri_count(const Soup& s) { return s.size() / 9; }

// Signed volume via the divergence theorem (positive for outward winding).
double signed_volume(const Soup& s)
{
    double vol = 0;
    for (size_t i = 0; i + 9 <= s.size(); i += 9)
    {
        const double ax = s[i],   ay = s[i+1], az = s[i+2];
        const double bx = s[i+3], by = s[i+4], bz = s[i+5];
        const double cx = s[i+6], cy = s[i+7], cz = s[i+8];
        vol += ax * (by * cz - bz * cy)
             + ay * (bz * cx - bx * cz)
             + az * (bx * cy - by * cx);
    }
    return vol / 6;
}

double surface_area(const Soup& s)
{
    double area = 0;
    for (size_t i = 0; i + 9 <= s.size(); i += 9)
    {
        const double ux = s[i+3] - s[i],  uy = s[i+4] - s[i+1],
                     uz = s[i+5] - s[i+2];
        const double vx = s[i+6] - s[i],  vy = s[i+7] - s[i+1],
                     vz = s[i+8] - s[i+2];
        const double nx = uy * vz - uz * vy;
        const double ny = uz * vx - ux * vz;
        const double nz = ux * vy - uy * vx;
        area += std::sqrt(nx*nx + ny*ny + nz*nz) / 2;
    }
    return area;
}

// Number of directed edges without a matching reverse edge.
// Zero for a closed, consistently-wound mesh.
size_t unmatched_edges(const Soup& s)
{
    typedef std::array<float, 6> Edge;
    std::map<Edge, long> count;
    for (size_t i = 0; i + 9 <= s.size(); i += 9)
    {
        const float* t = &s[i];
        for (int e = 0; e < 3; ++e)
        {
            const float* u = t + 3 * e;
            const float* v = t + 3 * ((e + 1) % 3);
            count[{u[0], u[1], u[2], v[0], v[1], v[2]}]++;
            count[{v[0], v[1], v[2], u[0], u[1], u[2]}]--;
        }
    }
    size_t bad = 0;
    for (const auto& kv : count)
        if (kv.second != 0)
            bad++;
    return bad;
}

bool has_vertex_near(const Soup& s, float x, float y, float z, float tol)
{
    for (size_t i = 0; i + 3 <= s.size(); i += 3)
    {
        const float dx = s[i] - x, dy = s[i+1] - y, dz = s[i+2] - z;
        if (std::sqrt(dx*dx + dy*dy + dz*dz) < tol)
            return true;
    }
    return false;
}

// Canonical text form: each triangle rotated so its lexicographically
// smallest corner comes first (winding preserved), then all triangles
// sorted.  Coordinates printed as hex floats, so the dump is exact and
// insensitive to triangle emission order.
std::string canonical_dump(const Soup& s)
{
    std::vector<std::string> tris;
    tris.reserve(tri_count(s));
    char buf[256];
    for (size_t i = 0; i + 9 <= s.size(); i += 9)
    {
        std::array<std::array<float, 3>, 3> v =
            {{ {s[i],   s[i+1], s[i+2]},
               {s[i+3], s[i+4], s[i+5]},
               {s[i+6], s[i+7], s[i+8]} }};
        int first = 0;
        for (int j = 1; j < 3; ++j)
            if (v[j] < v[first])
                first = j;
        snprintf(buf, sizeof(buf), "%a %a %a  %a %a %a  %a %a %a\n",
                 v[first][0], v[first][1], v[first][2],
                 v[(first+1)%3][0], v[(first+1)%3][1], v[(first+1)%3][2],
                 v[(first+2)%3][0], v[(first+2)%3][1], v[(first+2)%3][2]);
        tris.push_back(buf);
    }
    std::sort(tris.begin(), tris.end());
    std::string out;
    for (const auto& t : tris)
        out += t;
    return out;
}

void write_dump(const std::string& name, const Soup& s)
{
    std::ofstream f(name);
    f << canonical_dump(s);
}

// Test shapes (prefix math strings, same encoding fab/shapes.py emits):
//   sphere, radius 1, centered at origin
const char* SPHERE = "-r++qXqYqZf1";
//   axis-aligned cube spanning [-0.6, 0.6]^3 (sharp edges + corners)
const char* CUBE = "aa-f-0.6X-Xf0.6aa-f-0.6Y-Yf0.6a-f-0.6Z-Zf0.6";

}  // namespace

TEST_CASE("Mesher: sphere invariants")
{
    for (bool detect : {false, true})
    {
        CAPTURE(detect);
        auto s = mesh_shape(SPHERE, -1.2f, -1.2f, -1.2f, 1.2f, 1.2f, 1.2f,
                            20, detect);
        REQUIRE(tri_count(s) > 1000);
        REQUIRE(signed_volume(s) == Approx(4 * M_PI / 3).epsilon(0.03));
        REQUIRE(surface_area(s) == Approx(4 * M_PI).epsilon(0.03));
        REQUIRE(unmatched_edges(s) == 0);
    }
}

TEST_CASE("Mesher: cube invariants")
{
    for (bool detect : {false, true})
    {
        CAPTURE(detect);
        auto s = mesh_shape(CUBE, -1, -1, -1, 1, 1, 1, 20, detect);
        REQUIRE(tri_count(s) > 100);
        REQUIRE(signed_volume(s) == Approx(1.2 * 1.2 * 1.2).epsilon(0.03));
        REQUIRE(surface_area(s) == Approx(6 * 1.2 * 1.2).epsilon(0.03));
        REQUIRE(unmatched_edges(s) == 0);
    }
}

TEST_CASE("Mesher: feature detection reconstructs sharp corners")
{
    // With detect_features on, the Kobbelt algorithm should place
    // vertices at (or very near) the cube's true corners; plain
    // marching tets cannot (it clips them at voxel-edge crossings).
    auto s = mesh_shape(CUBE, -1, -1, -1, 1, 1, 1, 20, true);
    for (float sx : {-0.6f, 0.6f})
        for (float sy : {-0.6f, 0.6f})
            for (float sz : {-0.6f, 0.6f})
            {
                CAPTURE(sx);
                CAPTURE(sy);
                CAPTURE(sz);
                REQUIRE(has_vertex_near(s, sx, sy, sz, 0.02f));
            }
}

TEST_CASE("Mesher: indexed mesh matches soup output")
{
    for (bool detect : {false, true})
    {
        CAPTURE(detect);

        MathTree* tree = parse(CUBE);
        REQUIRE(tree != nullptr);
        Region r = {};
        r.ni = r.nj = r.nk = 40;
        r.voxels = uint64_t(r.ni) * r.nj * r.nk;
        build_arrays(&r, -1, -1, -1, 1, 1, 1);
        volatile int halt = 0;

        std::vector<float> verts;
        std::vector<uint32_t> indices;
        triangulate_indexed(tree, r, detect, &halt, verts, indices);
        free_arrays(&r);
        free_tree(tree);

        auto soup = mesh_shape(CUBE, -1, -1, -1, 1, 1, 1, 20, detect);

        // Same triangles, in the same order, with identical coordinates.
        const size_t expanded = indices.size() * 3;
        REQUIRE(expanded == soup.size());
        for (size_t i = 0; i < indices.size(); ++i)
        {
            REQUIRE(indices[i] < verts.size() / 3);
            for (int j = 0; j < 3; ++j)
                REQUIRE(verts[indices[i]*3 + j] == soup[i*3 + j]);
        }

        // Every vertex-table entry is referenced (compact output).
        std::vector<bool> seen(verts.size() / 3, false);
        for (auto i : indices)
            seen[i] = true;
        REQUIRE(std::find(seen.begin(), seen.end(), false) == seen.end());

        // And the two STL writers emit byte-identical files.
        save_stl(soup.data(), soup.size(), "stl_soup.stl");
        save_stl_indexed(verts.data(), indices.data(), indices.size() / 3,
                         "stl_indexed.stl");
        std::ifstream a("stl_soup.stl", std::ios::binary);
        std::ifstream b("stl_indexed.stl", std::ios::binary);
        std::string sa((std::istreambuf_iterator<char>(a)),
                        std::istreambuf_iterator<char>());
        std::string sb((std::istreambuf_iterator<char>(b)),
                        std::istreambuf_iterator<char>());
        REQUIRE(sa.size() == sb.size());
        REQUIRE(sa == sb);
        std::remove("stl_soup.stl");
        std::remove("stl_indexed.stl");
    }
}

// Golden dumps: canonical mesh output written to the working directory,
// for exact before/after comparison across mesher refactors.
// Run with:  SbFabTest "[golden]"
TEST_CASE("Mesher: golden dumps", "[golden]")
{
    write_dump("golden_sphere_plain.txt",
               mesh_shape(SPHERE, -1.2f, -1.2f, -1.2f, 1.2f, 1.2f, 1.2f,
                          25, false));
    write_dump("golden_sphere_detect.txt",
               mesh_shape(SPHERE, -1.2f, -1.2f, -1.2f, 1.2f, 1.2f, 1.2f,
                          25, true));
    write_dump("golden_cube_plain.txt",
               mesh_shape(CUBE, -1, -1, -1, 1, 1, 1, 25, false));
    write_dump("golden_cube_detect.txt",
               mesh_shape(CUBE, -1, -1, -1, 1, 1, 1, 25, true));
    REQUIRE(true);
}

// Heavy benchmark, hidden from the default run.  Mesh a gyroid shell:
// lots of surface area, so triangle count scales hard with resolution.
// Run with:  /usr/bin/time -v SbFabTest "[bench]"  (peak RSS + wall time)
TEST_CASE("Mesher: gyroid benchmark", "[.bench]")
{
    const char* gyroid =
        "=cos(X)*sin(Y) + cos(Y)*sin(Z) + cos(Z)*sin(X) - 0.2;";
    const float res = std::getenv("MESHER_BENCH_RES")
        ? std::atof(std::getenv("MESHER_BENCH_RES")) : 6.0f;
    const bool detect = !std::getenv("MESHER_BENCH_PLAIN");

    if (std::getenv("MESHER_BENCH_INDEXED"))
    {
        // The indexed path, as used by mesh export.
        MathTree* tree = parse(gyroid);
        REQUIRE(tree != nullptr);
        Region r = {};
        r.ni = r.nj = r.nk = uint32_t(16 * res);
        r.voxels = uint64_t(r.ni) * r.nj * r.nk;
        build_arrays(&r, -8, -8, -8, 8, 8, 8);
        volatile int halt = 0;
        std::vector<float> verts;
        std::vector<uint32_t> indices;
        triangulate_indexed(tree, r, detect, &halt, verts, indices);
        free_arrays(&r);
        free_tree(tree);
        printf("gyroid @ res %g, detect=%d: %zu triangles, %zu verts "
               "(indexed)\n", res, detect, indices.size() / 3,
               verts.size() / 3);
        REQUIRE(indices.size() > 0);
    }
    else
    {
        auto s = mesh_shape(gyroid, -8, -8, -8, 8, 8, 8, res, detect);
        printf("gyroid @ res %g, detect=%d: %zu triangles\n",
               res, detect, tri_count(s));
        REQUIRE(tri_count(s) > 0);
    }
}

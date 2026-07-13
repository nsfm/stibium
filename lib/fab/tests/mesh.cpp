#include <catch/catch.hpp>

#include <array>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <random>
#include <string>
#include <vector>

#include "fab/mesh/trimesh.h"
#include "fab/mesh/mesh_query.h"

using namespace fab_mesh;

namespace {

// Flat triangle soup: 9 floats per triangle (xyz * 3 corners).
using Soup = std::vector<float>;

std::string g_tmpdir;

std::string tmp_path(const char* name)
{
    if (g_tmpdir.empty())
    {
        char tmpl[] = "/tmp/sbmesh_XXXXXX";
        char* d = mkdtemp(tmpl);
        g_tmpdir = d ? d : "/tmp";
    }
    return g_tmpdir + "/" + name;
}

struct V3 { float x, y, z; };

V3 tri_normal(const float* t)
{
    float ux = t[3]-t[0], uy = t[4]-t[1], uz = t[5]-t[2];
    float vx = t[6]-t[0], vy = t[7]-t[1], vz = t[8]-t[2];
    float nx = uy*vz - uz*vy, ny = uz*vx - ux*vz, nz = ux*vy - uy*vx;
    float l = std::sqrt(nx*nx + ny*ny + nz*nz);
    if (l == 0) l = 1;
    return {nx/l, ny/l, nz/l};
}

// --- geometry -----------------------------------------------------------

// Outward-wound cube spanning [-s, s]^3, 12 triangles.
Soup make_cube(float s)
{
    const float c[8][3] = {
        {-s,-s,-s}, { s,-s,-s}, { s, s,-s}, {-s, s,-s},
        {-s,-s, s}, { s,-s, s}, { s, s, s}, {-s, s, s}
    };
    const int quads[6][4] = {
        {1,2,6,5}, {0,4,7,3}, {3,7,6,2}, {0,1,5,4}, {4,5,6,7}, {0,3,2,1}
    };
    Soup s_out;
    auto push = [&](int i){ s_out.push_back(c[i][0]);
                            s_out.push_back(c[i][1]);
                            s_out.push_back(c[i][2]); };
    for (auto& q : quads)
    {
        push(q[0]); push(q[1]); push(q[2]);
        push(q[0]); push(q[2]); push(q[3]);
    }
    return s_out;
}

// Icosphere: subdivided icosahedron, all vertices on radius r.
Soup make_icosphere(int level, float r)
{
    const float t = (1.0f + std::sqrt(5.0f)) / 2.0f;
    std::vector<std::array<float,3>> vs = {
        {-1, t, 0}, {1, t, 0}, {-1,-t, 0}, {1,-t, 0},
        {0,-1, t}, {0, 1, t}, {0,-1,-t}, {0, 1,-t},
        {t, 0,-1}, {t, 0, 1}, {-t, 0,-1}, {-t, 0, 1}
    };
    std::vector<std::array<int,3>> fs = {
        {0,11,5},{0,5,1},{0,1,7},{0,7,10},{0,10,11},
        {1,5,9},{5,11,4},{11,10,2},{10,7,6},{7,1,8},
        {3,9,4},{3,4,2},{3,2,6},{3,6,8},{3,8,9},
        {4,9,5},{2,4,11},{6,2,10},{8,6,7},{9,8,1}
    };

    for (int it = 0; it < level; ++it)
    {
        std::map<std::pair<int,int>,int> mid;
        auto midpoint = [&](int a, int b) -> int {
            auto key = std::minmax(a, b);
            auto f = mid.find(key);
            if (f != mid.end()) return f->second;
            std::array<float,3> m = {
                (vs[a][0]+vs[b][0])/2, (vs[a][1]+vs[b][1])/2,
                (vs[a][2]+vs[b][2])/2 };
            int id = (int)vs.size();
            vs.push_back(m);
            mid[key] = id;
            return id;
        };
        std::vector<std::array<int,3>> nf;
        for (auto& f : fs)
        {
            int a = midpoint(f[0], f[1]);
            int b = midpoint(f[1], f[2]);
            int c = midpoint(f[2], f[0]);
            nf.push_back({f[0], a, c});
            nf.push_back({f[1], b, a});
            nf.push_back({f[2], c, b});
            nf.push_back({a, b, c});
        }
        fs.swap(nf);
    }

    for (auto& v : vs)
    {
        float l = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
        v[0] = v[0]/l*r; v[1] = v[1]/l*r; v[2] = v[2]/l*r;
    }

    Soup s;
    for (auto& f : fs)
        for (int k = 0; k < 3; ++k)
        {
            s.push_back(vs[f[k]][0]);
            s.push_back(vs[f[k]][1]);
            s.push_back(vs[f[k]][2]);
        }
    return s;
}

// --- STL writers --------------------------------------------------------

void write_binary_stl(const std::string& path, const Soup& soup,
                      const std::string& header)
{
    std::ofstream f(path, std::ios::binary);
    char hdr[80];
    std::memset(hdr, 0, sizeof(hdr));
    std::memcpy(hdr, header.data(), std::min(header.size(), sizeof(hdr)));
    f.write(hdr, 80);
    uint32_t count = (uint32_t)(soup.size() / 9);
    f.write((const char*)&count, 4);
    for (size_t i = 0; i + 9 <= soup.size(); i += 9)
    {
        V3 n = tri_normal(&soup[i]);
        f.write((const char*)&n, 12);
        f.write((const char*)&soup[i], 36);
        uint16_t attr = 0;
        f.write((const char*)&attr, 2);
    }
}

void write_ascii_stl(const std::string& path, const Soup& soup,
                     const std::string& name)
{
    std::ofstream f(path);
    f << "solid " << name << "\n";
    for (size_t i = 0; i + 9 <= soup.size(); i += 9)
    {
        V3 n = tri_normal(&soup[i]);
        f << "  facet normal " << n.x << " " << n.y << " " << n.z << "\n";
        f << "    outer loop\n";
        for (int k = 0; k < 3; ++k)
            f << "      vertex " << soup[i+k*3] << " "
              << soup[i+k*3+1] << " " << soup[i+k*3+2] << "\n";
        f << "    endloop\n  endfacet\n";
    }
    f << "endsolid " << name << "\n";
}

// Brute-force unsigned distance over a loaded mesh (closest point on
// each triangle), used as ground truth for the accelerated query.
float closest_dist2(const float* p, const float* a,
                    const float* b, const float* c)
{
    float ab[3]={b[0]-a[0],b[1]-a[1],b[2]-a[2]};
    float ac[3]={c[0]-a[0],c[1]-a[1],c[2]-a[2]};
    float ap[3]={p[0]-a[0],p[1]-a[1],p[2]-a[2]};
    auto d3 = [](const float* u, const float* v){
        return u[0]*v[0]+u[1]*v[1]+u[2]*v[2]; };
    float d1=d3(ab,ap), d2=d3(ac,ap);
    float q[3];
    auto setq = [&](float px,float py,float pz){ q[0]=px;q[1]=py;q[2]=pz; };
    if (d1<=0 && d2<=0) { setq(a[0],a[1],a[2]); }
    else {
        float bp[3]={p[0]-b[0],p[1]-b[1],p[2]-b[2]};
        float d3_=d3(ab,bp), d4=d3(ac,bp);
        float vc=d1*d4-d3_*d2;
        float cp[3]={p[0]-c[0],p[1]-c[1],p[2]-c[2]};
        float d5=d3(ab,cp), d6=d3(ac,cp);
        float vb=d5*d2-d1*d6;
        float va=d3_*d6-d5*d4;
        if (d3_>=0 && d4<=d3_) { setq(b[0],b[1],b[2]); }
        else if (vc<=0 && d1>=0 && d3_<=0) {
            float v=d1/(d1-d3_);
            setq(a[0]+ab[0]*v,a[1]+ab[1]*v,a[2]+ab[2]*v);
        }
        else if (d6>=0 && d5<=d6) { setq(c[0],c[1],c[2]); }
        else if (vb<=0 && d2>=0 && d6<=0) {
            float w=d2/(d2-d6);
            setq(a[0]+ac[0]*w,a[1]+ac[1]*w,a[2]+ac[2]*w);
        }
        else if (va<=0 && (d4-d3_)>=0 && (d5-d6)>=0) {
            float w=(d4-d3_)/((d4-d3_)+(d5-d6));
            setq(b[0]+(c[0]-b[0])*w,b[1]+(c[1]-b[1])*w,b[2]+(c[2]-b[2])*w);
        }
        else {
            float denom=1.0f/(va+vb+vc);
            float v=vb*denom, w=vc*denom;
            setq(a[0]+ab[0]*v+ac[0]*w, a[1]+ab[1]*v+ac[1]*w,
                 a[2]+ab[2]*v+ac[2]*w);
        }
    }
    float dx=p[0]-q[0], dy=p[1]-q[1], dz=p[2]-q[2];
    return dx*dx+dy*dy+dz*dz;
}

float brute_distance(const TriMesh& m, float x, float y, float z)
{
    float p[3]={x,y,z};
    float best=1e30f;
    for (uint32_t t=0; t<m.tri_count(); ++t)
    {
        const float* a=&m.verts[3*m.tris[3*t]];
        const float* b=&m.verts[3*m.tris[3*t+1]];
        const float* c=&m.verts[3*m.tris[3*t+2]];
        float d=closest_dist2(p,a,b,c);
        if (d<best) best=d;
    }
    return std::sqrt(best);
}

}  // namespace

TEST_CASE("Mesh: binary STL round-trip welds a cube", "[mesh]")
{
    Soup cube = make_cube(1.0f);
    std::string path = tmp_path("cube_bin.stl");
    write_binary_stl(path, cube, "cube");

    TriMesh m;
    std::string err;
    bool ok = load_stl(path, &m, &err);
    CAPTURE(err);
    REQUIRE(ok);
    REQUIRE(m.tri_count() == 12u);
    REQUIRE(m.vert_count() == 8u);
    REQUIRE(m.bbox[0] == Approx(-1.0f));
    REQUIRE(m.bbox[1] == Approx(-1.0f));
    REQUIRE(m.bbox[2] == Approx(-1.0f));
    REQUIRE(m.bbox[3] == Approx(1.0f));
    REQUIRE(m.bbox[4] == Approx(1.0f));
    REQUIRE(m.bbox[5] == Approx(1.0f));
}

TEST_CASE("Mesh: ASCII STL load matches binary counts", "[mesh]")
{
    Soup cube = make_cube(1.0f);
    std::string path = tmp_path("cube_ascii.stl");
    write_ascii_stl(path, cube, "cube");

    TriMesh m;
    std::string err;
    bool ok = load_stl(path, &m, &err);
    CAPTURE(err);
    REQUIRE(ok);
    REQUIRE(m.tri_count() == 12u);
    REQUIRE(m.vert_count() == 8u);
}

TEST_CASE("Mesh: truncated binary STL is rejected", "[mesh]")
{
    Soup cube = make_cube(1.0f);
    std::string path = tmp_path("cube_trunc.stl");
    write_binary_stl(path, cube, "cube");

    // Chop the tail off so the declared triangle count overruns the file.
    std::ifstream in(path, std::ios::binary);
    std::string data((std::istreambuf_iterator<char>(in)),
                      std::istreambuf_iterator<char>());
    in.close();
    data.resize(data.size() - 120);
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    out.write(data.data(), data.size());
    out.close();

    TriMesh m;
    std::string err;
    bool ok = load_stl(path, &m, &err);
    REQUIRE_FALSE(ok);
    bool has_err = !err.empty();
    REQUIRE(has_err);
}

TEST_CASE("Mesh: binary header is preserved", "[mesh]")
{
    Soup cube = make_cube(1.0f);
    std::string path = tmp_path("cube_hdr.stl");
    const std::string hdr = "Stibium binary STL export v1";
    write_binary_stl(path, cube, hdr);

    TriMesh m;
    std::string err;
    bool ok = load_stl(path, &m, &err);
    CAPTURE(err);
    REQUIRE(ok);
    REQUIRE(m.header == hdr);
}

TEST_CASE("Mesh: unsigned distance matches brute force", "[mesh]")
{
    Soup sphere = make_icosphere(2, 1.0f);
    std::string path = tmp_path("sphere.stl");
    write_binary_stl(path, sphere, "icosphere");

    TriMesh m;
    std::string err;
    bool ok = load_stl(path, &m, &err);
    CAPTURE(err);
    REQUIRE(ok);

    MeshQuery q(m);
    std::mt19937 rng(1234);
    std::uniform_real_distribution<float> u(-2.0f, 2.0f);
    float worst = 0;
    for (int i = 0; i < 200; ++i)
    {
        float x = u(rng), y = u(rng), z = u(rng);
        float got = q.unsigned_distance(x, y, z);
        float ref = brute_distance(m, x, y, z);
        worst = std::max(worst, std::fabs(got - ref));
    }
    CAPTURE(worst);
    bool close = worst < 1e-4f;
    REQUIRE(close);
}

TEST_CASE("Mesh: winding number in/out of a closed sphere", "[mesh]")
{
    Soup sphere = make_icosphere(2, 1.0f);
    std::string path = tmp_path("sphere_w.stl");
    write_binary_stl(path, sphere, "icosphere");

    TriMesh m;
    std::string err;
    REQUIRE(load_stl(path, &m, &err));
    MeshQuery q(m);

    std::mt19937 rng(77);
    std::uniform_real_distribution<float> u(-1.0f, 1.0f);
    for (int i = 0; i < 30; ++i)
    {
        float d[3] = {u(rng), u(rng), u(rng)};
        float l = std::sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
        if (l < 1e-3f) continue;
        float win_in  = q.winding_number(d[0]/l*0.5f, d[1]/l*0.5f, d[2]/l*0.5f);
        float win_out = q.winding_number(d[0]/l*1.5f, d[1]/l*1.5f, d[2]/l*1.5f);
        CAPTURE(win_in);
        REQUIRE(win_in > 0.9f);
        CAPTURE(win_out);
        REQUIRE(win_out < 0.1f);
    }
}

TEST_CASE("Mesh: winding number is robust on an open mesh", "[mesh]")
{
    // Cube with one square face (two triangles) removed -> a hole.
    Soup cube = make_cube(1.0f);
    Soup open(cube.begin() + 18, cube.end());   // drop first 2 triangles
    std::string path = tmp_path("open_cube.stl");
    write_binary_stl(path, open, "open");

    TriMesh m;
    std::string err;
    REQUIRE(load_stl(path, &m, &err));
    REQUIRE(m.tri_count() == 10u);
    MeshQuery q(m);

    float center = q.winding_number(0, 0, 0);
    CAPTURE(center);
    REQUIRE(center > 0.5f);

    float outside = q.winding_number(5, 5, 5);
    CAPTURE(outside);
    REQUIRE(outside < 0.3f);
}

TEST_CASE("Mesh: signed distance sign follows inside/outside", "[mesh]")
{
    Soup sphere = make_icosphere(3, 1.0f);
    std::string path = tmp_path("sphere_s.stl");
    write_binary_stl(path, sphere, "icosphere");

    TriMesh m;
    std::string err;
    REQUIRE(load_stl(path, &m, &err));
    MeshQuery q(m);

    float center = q.signed_distance(0, 0, 0);
    CAPTURE(center);
    bool center_ok = center < 0 && std::fabs(center - (-1.0f)) < 0.05f;
    REQUIRE(center_ok);

    float outside = q.signed_distance(1.5f, 0, 0);
    CAPTURE(outside);
    bool outside_ok = outside > 0 && std::fabs(outside - 0.5f) < 0.05f;
    REQUIRE(outside_ok);
}

TEST_CASE("Mesh: signed_distance throughput", "[mesh]")
{
    Soup sphere = make_icosphere(3, 1.0f);
    std::string path = tmp_path("sphere_p.stl");
    write_binary_stl(path, sphere, "icosphere");

    TriMesh m;
    std::string err;
    REQUIRE(load_stl(path, &m, &err));
    MeshQuery q(m);

    std::mt19937 rng(9);
    std::uniform_real_distribution<float> u(-2.0f, 2.0f);
    const int N = 100000;
    std::vector<float> pts(N * 3);
    for (int i = 0; i < N * 3; ++i) pts[i] = u(rng);

    clock_t t0 = clock();
    double acc = 0;
    for (int i = 0; i < N; ++i)
        acc += q.signed_distance(pts[3*i], pts[3*i+1], pts[3*i+2]);
    clock_t t1 = clock();
    double ms = 1000.0 * double(t1 - t0) / CLOCKS_PER_SEC;
    printf("[mesh] 100k signed_distance on icosphere(%u tris): %.1f ms "
           "(%.2f us/query), checksum=%.3f\n",
           m.tri_count(), ms, ms * 1000.0 / N, acc);
    REQUIRE(N > 0);
}

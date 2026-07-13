#include "fab/mesh/mesh_query.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace fab_mesh {

namespace {

// Minimal 3-vector; everything here is stack-local so queries stay
// const and thread-safe.
struct V3 {
    float x, y, z;
};

inline V3 operator-(const V3& a, const V3& b)
{ return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline V3 operator+(const V3& a, const V3& b)
{ return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline V3 operator*(const V3& a, float s)
{ return {a.x * s, a.y * s, a.z * s}; }
inline float dot(const V3& a, const V3& b)
{ return a.x * b.x + a.y * b.y + a.z * b.z; }
inline V3 cross(const V3& a, const V3& b)
{ return {a.y * b.z - a.z * b.y,
          a.z * b.x - a.x * b.z,
          a.x * b.y - a.y * b.x}; }
inline float len(const V3& a) { return std::sqrt(dot(a, a)); }

const float PI = 3.14159265358979323846f;
const float FOUR_PI = 4.0f * PI;

// Closest point on triangle abc to p (Ericson, Real-Time Collision
// Detection 5.1.5); returns the point via voronoi-region tests.
V3 closest_on_tri(const V3& p, const V3& a, const V3& b, const V3& c)
{
    const V3 ab = b - a;
    const V3 ac = c - a;
    const V3 ap = p - a;
    float d1 = dot(ab, ap);
    float d2 = dot(ac, ap);
    if (d1 <= 0.0f && d2 <= 0.0f)
        return a;

    const V3 bp = p - b;
    float d3 = dot(ab, bp);
    float d4 = dot(ac, bp);
    if (d3 >= 0.0f && d4 <= d3)
        return b;

    float vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f)
    {
        float v = d1 / (d1 - d3);
        return a + ab * v;
    }

    const V3 cp = p - c;
    float d5 = dot(ab, cp);
    float d6 = dot(ac, cp);
    if (d6 >= 0.0f && d5 <= d6)
        return c;

    float vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f)
    {
        float w = d2 / (d2 - d6);
        return a + ac * w;
    }

    float va = d3 * d6 - d5 * d4;
    if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f)
    {
        float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return b + (c - b) * w;
    }

    float denom = 1.0f / (va + vb + vc);
    float v = vb * denom;
    float w = vc * denom;
    return a + ab * v + ac * w;
}

// Signed solid angle subtended by triangle abc at p, via the
// van Oosterom-Strackee formula.  Sign follows the vertex winding.
float solid_angle(const V3& p, const V3& a, const V3& b, const V3& c)
{
    const V3 A = a - p, B = b - p, C = c - p;
    float la = len(A), lb = len(B), lc = len(C);
    float num = dot(A, cross(B, C));
    float den = la * lb * lc
              + dot(A, B) * lc
              + dot(B, C) * la
              + dot(C, A) * lb;
    return 2.0f * std::atan2(num, den);
}

struct Node {
    float bmin[3], bmax[3];
    V3 dipole;      // area-weighted normal sum
    V3 center;      // area-weighted centroid (dipole expansion point)
    float radius;   // max distance from center to any node vertex
    int start, count;   // triangle range in the reordered index array
    int left, right;    // children, or -1/-1 for a leaf
};

}  // namespace

struct MeshQuery::Impl {
    std::vector<float> verts;       // xyz triples
    std::vector<uint32_t> tris;     // index triples
    std::vector<int> order;         // triangle ids, reordered by build
    std::vector<Node> nodes;

    // Per-triangle precomputed quantities (indexed by triangle id).
    std::vector<V3> centroid;
    std::vector<V3> area_normal;    // 0.5*cross(b-a,c-a): |.| == area

    V3 vert(uint32_t i) const
    { return {verts[3*i], verts[3*i+1], verts[3*i+2]}; }

    void tri_pts(int t, V3* a, V3* b, V3* c) const
    {
        *a = vert(tris[3*t]);
        *b = vert(tris[3*t+1]);
        *c = vert(tris[3*t+2]);
    }

    int build(int start, int count);
    void aggregate(Node& n);

    float unsigned_distance(const V3& p) const;
    float winding_number(const V3& p) const;
};

// Fill a node's bbox, dipole, expansion center, and radius from its
// triangle range.  Called after children so the range is stable.
void MeshQuery::Impl::aggregate(Node& n)
{
    float bmin[3] = { 1e30f,  1e30f,  1e30f};
    float bmax[3] = {-1e30f, -1e30f, -1e30f};
    V3 dip{0, 0, 0};
    V3 wc{0, 0, 0};
    float area_sum = 0;

    for (int i = 0; i < n.count; ++i)
    {
        int t = order[n.start + i];
        V3 a, b, c;
        tri_pts(t, &a, &b, &c);
        const V3 pts[3] = {a, b, c};
        for (const V3& v : pts)
        {
            bmin[0] = std::min(bmin[0], v.x); bmax[0] = std::max(bmax[0], v.x);
            bmin[1] = std::min(bmin[1], v.y); bmax[1] = std::max(bmax[1], v.y);
            bmin[2] = std::min(bmin[2], v.z); bmax[2] = std::max(bmax[2], v.z);
        }
        float area = len(area_normal[t]);
        dip = dip + area_normal[t];
        wc = wc + centroid[t] * area;
        area_sum += area;
    }
    for (int j = 0; j < 3; ++j) { n.bmin[j] = bmin[j]; n.bmax[j] = bmax[j]; }
    n.dipole = dip;
    if (area_sum > 0)
        n.center = wc * (1.0f / area_sum);
    else
        n.center = {(bmin[0]+bmax[0])*0.5f, (bmin[1]+bmax[1])*0.5f,
                    (bmin[2]+bmax[2])*0.5f};

    // Radius: farthest node vertex from the expansion center.
    float r2 = 0;
    for (int i = 0; i < n.count; ++i)
    {
        int t = order[n.start + i];
        V3 a, b, c;
        tri_pts(t, &a, &b, &c);
        const V3 pts[3] = {a, b, c};
        for (const V3& v : pts)
            r2 = std::max(r2, dot(v - n.center, v - n.center));
    }
    n.radius = std::sqrt(r2);
}

// Recursive median split on the longest centroid-extent axis.
int MeshQuery::Impl::build(int start, int count)
{
    int id = (int)nodes.size();
    nodes.push_back(Node{});
    // Note: nodes may reallocate during recursion, so index by id and
    // re-fetch after children are built rather than holding a Node&.
    {
        Node n{};
        n.start = start;
        n.count = count;
        n.left = n.right = -1;
        nodes[id] = n;
    }

    const int LEAF = 8;
    if (count > LEAF)
    {
        // Split on the longest axis of the centroid bounds.
        float cmin[3] = { 1e30f,  1e30f,  1e30f};
        float cmax[3] = {-1e30f, -1e30f, -1e30f};
        for (int i = 0; i < count; ++i)
        {
            const V3& c = centroid[order[start + i]];
            const float cc[3] = {c.x, c.y, c.z};
            for (int j = 0; j < 3; ++j)
            {
                cmin[j] = std::min(cmin[j], cc[j]);
                cmax[j] = std::max(cmax[j], cc[j]);
            }
        }
        int axis = 0;
        float ext = cmax[0] - cmin[0];
        if (cmax[1] - cmin[1] > ext) { ext = cmax[1] - cmin[1]; axis = 1; }
        if (cmax[2] - cmin[2] > ext) { ext = cmax[2] - cmin[2]; axis = 2; }

        if (ext > 0)
        {
            int mid = count / 2;
            auto first = order.begin() + start;
            std::nth_element(first, first + mid, first + count,
                [&](int lhs, int rhs) {
                    const V3& a = centroid[lhs];
                    const V3& b = centroid[rhs];
                    float av = axis == 0 ? a.x : axis == 1 ? a.y : a.z;
                    float bv = axis == 0 ? b.x : axis == 1 ? b.y : b.z;
                    return av < bv;
                });
            int l = build(start, mid);
            int r = build(start + mid, count - mid);
            nodes[id].left = l;
            nodes[id].right = r;
        }
    }

    aggregate(nodes[id]);
    return id;
}

// Squared distance from p to an axis-aligned box.
static inline float box_dist2(const Node& n, const V3& p)
{
    float d = 0;
    const float pc[3] = {p.x, p.y, p.z};
    for (int j = 0; j < 3; ++j)
    {
        float e = 0;
        if (pc[j] < n.bmin[j]) e = n.bmin[j] - pc[j];
        else if (pc[j] > n.bmax[j]) e = pc[j] - n.bmax[j];
        d += e * e;
    }
    return d;
}

float MeshQuery::Impl::unsigned_distance(const V3& p) const
{
    if (nodes.empty())
        return 0;
    float best = 1e30f;   // squared

    int stack[128];
    int sp = 0;
    stack[sp++] = 0;
    while (sp > 0)
    {
        int ni = stack[--sp];
        const Node& n = nodes[ni];
        if (box_dist2(n, p) >= best)
            continue;
        if (n.left < 0)
        {
            for (int i = 0; i < n.count; ++i)
            {
                int t = order[n.start + i];
                V3 a, b, c;
                tri_pts(t, &a, &b, &c);
                V3 q = closest_on_tri(p, a, b, c);
                V3 d = p - q;
                best = std::min(best, dot(d, d));
            }
        }
        else
        {
            // Push the farther child first so the nearer is popped next.
            float dl = box_dist2(nodes[n.left], p);
            float dr = box_dist2(nodes[n.right], p);
            if (dl < dr)
            {
                stack[sp++] = n.right;
                stack[sp++] = n.left;
            }
            else
            {
                stack[sp++] = n.left;
                stack[sp++] = n.right;
            }
        }
    }
    return std::sqrt(best);
}

float MeshQuery::Impl::winding_number(const V3& p) const
{
    if (nodes.empty())
        return 0;
    const float beta = 2.0f;
    float sum = 0;   // accumulates solid angle; scaled by 1/4pi at end

    int stack[128];
    int sp = 0;
    stack[sp++] = 0;
    while (sp > 0)
    {
        int ni = stack[--sp];
        const Node& n = nodes[ni];
        V3 d = n.center - p;
        float dist = len(d);
        if (n.left >= 0 && dist > beta * n.radius)
        {
            // Far field: first-order dipole approximation.
            sum += dot(n.dipole, d) / (dist * dist * dist);
        }
        else if (n.left < 0)
        {
            // Leaf: exact per-triangle solid angles.
            for (int i = 0; i < n.count; ++i)
            {
                int t = order[n.start + i];
                V3 a, b, c;
                tri_pts(t, &a, &b, &c);
                sum += solid_angle(p, a, b, c);
            }
        }
        else
        {
            stack[sp++] = n.left;
            stack[sp++] = n.right;
        }
    }
    return sum / FOUR_PI;
}

MeshQuery::MeshQuery(const TriMesh& mesh)
    : impl(new Impl)
{
    impl->verts = mesh.verts;
    impl->tris = mesh.tris;

    const int nt = (int)mesh.tri_count();
    impl->centroid.resize(nt);
    impl->area_normal.resize(nt);
    impl->order.resize(nt);
    for (int t = 0; t < nt; ++t)
    {
        impl->order[t] = t;
        V3 a, b, c;
        impl->tri_pts(t, &a, &b, &c);
        impl->centroid[t] = (a + b + c) * (1.0f / 3.0f);
        impl->area_normal[t] = cross(b - a, c - a) * 0.5f;
    }

    if (nt > 0)
    {
        impl->nodes.reserve(2 * nt);
        impl->build(0, nt);
    }
}

MeshQuery::~MeshQuery() = default;

float MeshQuery::unsigned_distance(float x, float y, float z) const
{
    return impl->unsigned_distance({x, y, z});
}

float MeshQuery::winding_number(float x, float y, float z) const
{
    return impl->winding_number({x, y, z});
}

float MeshQuery::signed_distance(float x, float y, float z) const
{
    float d = impl->unsigned_distance({x, y, z});
    return impl->winding_number({x, y, z}) > 0.5f ? -d : d;
}

}  // namespace fab_mesh

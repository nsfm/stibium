#include <algorithm>
#include <cstring>
#include <iostream>
#include <unordered_set>

#include "fab/tree/triangulate/mesher.h"

#include "fab/tree/tape.h"
#include "fab/util/switches.h"

#if MIN_VOLUME < 60
#error "MIN_VOLUME is below minimum for meshing implementation."
#endif

/*
 *  When feature detection is enabled, the Mesher implements the algorithm from
 *
 *  "Feature Sensitive Surface Extraction from Volume Data"
 *
 *  (Kobbelt, Leif P. and Botsch, Mario and
 *   Schwanecke, Ulrich and Seidel, Hans-Peter)
 *
 *  SIGGRAPH 2001
 *
 *  Triangles are stored as vertex indices into an interned vertex table
 *  (12 bytes per triangle) rather than as inline coordinates; vertex
 *  positions are deduplicated on their float bit patterns as they are
 *  created, so the output is naturally an indexed mesh.
 */

static const uint8_t VERTEX_LOOP[] = {6, 4, 5, 1, 3, 2, 6};

// Based on which vertices are filled, this map tells you which
// edges to interpolate between when forming zero, one, or two
// triangles for a tetrahedron.
// (filled vertex is first in the pair, and is given as a tetrahedron vertex
//  so you have to translate into a proper cube vertex).
static const int EDGE_MAP[16][2][3][2] = {
    {{{-1,-1}, {-1,-1}, {-1,-1}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // ----
    {{{ 0, 2}, { 0, 1}, { 0, 3}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // ---0
    {{{ 1, 2}, { 1, 3}, { 1, 0}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // --1-
    {{{ 1, 2}, { 1, 3}, { 0, 3}}, {{ 0, 3}, { 0, 2}, { 1, 2}}}, // --10
    {{{ 2, 0}, { 2, 3}, { 2, 1}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // -2--
    {{{ 0, 3}, { 2, 3}, { 2, 1}}, {{ 2, 1}, { 0, 1}, { 0, 3}}}, // -2-0
    {{{ 2, 0}, { 2, 3}, { 1, 0}}, {{ 2, 3}, { 1, 3}, { 1, 0}}}, // -21-
    {{{ 2, 3}, { 1, 3}, { 0, 3}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // -210

    {{{ 3, 0}, { 3, 1}, { 3, 2}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // 3---
    {{{ 3, 2}, { 0, 2}, { 0, 1}}, {{ 3, 1}, { 3, 2}, { 0, 1}}}, // 3--0
    {{{ 1, 2}, { 3, 2}, { 3, 0}}, {{ 3, 0}, { 1, 0}, { 1, 2}}}, // 3-1-
    {{{ 1, 2}, { 3, 2}, { 0, 2}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // 3-10
    {{{ 3, 0}, { 3, 1}, { 2, 1}}, {{ 2, 1}, { 2, 0}, { 3, 0}}}, // 32--
    {{{ 3, 1}, { 2, 1}, { 0, 1}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // 32-0
    {{{ 3, 0}, { 1, 0}, { 2, 0}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // 321-
    {{{-1,-1}, {-1,-1}, {-1,-1}}, {{-1,-1}, {-1,-1}, {-1,-1}}}, // 3210
};

Mesher::Mesher(const Deck* deck, bool detect_edges, volatile int* halt,
               std::atomic<uint64_t>* progress)
    : deck(deck), ctx(tape_ctx_new(deck)), block_tape(nullptr),
      detect_edges(detect_edges), halt(halt),
      progress(progress),
      data(new float[MIN_VOLUME]), has_data(false),
      X(new float[MIN_VOLUME]),
      Y(new float[MIN_VOLUME]),
      Z(new float[MIN_VOLUME]),
      ex(new float[MIN_VOLUME]),
      ey(new float[MIN_VOLUME]),
      ez(new float[MIN_VOLUME]),
      nx(new float[MIN_VOLUME]),
      ny(new float[MIN_VOLUME]),
      nz(new float[MIN_VOLUME]),
      finalized(false),
      voxel_start(NONE), voxel_end(NONE), fan_start(NONE)
{
    // Nothing to do here
}

Mesher::~Mesher()
{
    tape_ctx_free(ctx);
    for (auto ptr : {data, X, Y, Z, ex, ey, ez, nx, ny, nz})
        delete [] ptr;
}

uint32_t Mesher::intern(float x, float y, float z)
{
    // Normalize -0.0f to +0.0f: keys are compared by bit pattern, but
    // the mesh topology logic needs value equality (the pre-indexed
    // code keyed float values through std::map/std::set).
    if (x == 0)     x = 0;
    if (y == 0)     y = 0;
    if (z == 0)     z = 0;

    std::array<uint32_t, 3> key;
    memcpy(&key[0], &x, sizeof(float));
    memcpy(&key[1], &y, sizeof(float));
    memcpy(&key[2], &z, sizeof(float));

    auto found = vertex_ids.find(key);
    if (found != vertex_ids.end())
        return found->second;

    // UINT32_MAX is the tombstone sentinel, so it can never be handed
    // out as a live vertex id (that would need a ~50 GB vertex table;
    // fail loudly rather than silently dropping triangles).
    if (vertices.size() >= UINT32_MAX)
    {
        std::cerr << "Mesher: vertex count overflow" << std::endl;
        std::abort();
    }

    const uint32_t id = vertices.size();
    vertices.push_back({{x, y, z}});
    vertex_ids.emplace(key, id);
    return id;
}

Vec3f Mesher::pos(uint32_t v) const
{
    return Vec3f(vertices[v][0], vertices[v][1], vertices[v][2]);
}

// Estimate the normals of a set of contour points.
std::vector<Vec3f> Mesher::get_normals(const std::vector<uint32_t>& contour)
{
    // Find epsilon as the single shortest side length divided by 100.
    float epsilon = INFINITY;
    for (size_t i = 0; i + 1 < contour.size(); ++i)
    {
        const Vec3f d = pos(contour[i+1]) - pos(contour[i]);
        epsilon = fmin(epsilon, d.norm() / 100.0f);
    }

    // We'll be evaluating a dummy region to numerically estimate gradients
    Region dummy;
    dummy.voxels = contour.size() * 7;
    if (dummy.voxels >= MIN_VOLUME)
    {
        // The nx/ny/nz scratch buffers hold MIN_VOLUME floats; writing
        // contour.size() * 7 entries past that would overflow them.
        // Returning no normals makes check_feature treat this fan as
        // featureless, which is safe (it just isn't sharpened).
        std::cerr << "Error: too many normals to calculate at once!"
                  << std::endl;
        return {};
    }
    dummy.X = nx;
    dummy.Y = ny;
    dummy.Z = nz;

    // Load position data into the dummy region
    int i=0;
    for (auto c : contour)
    {
        const Vec3f v = pos(c);
        dummy.X[i]   = v[0];
        dummy.X[i+1] = v[0] + epsilon;
        dummy.X[i+2] = v[0];
        dummy.X[i+3] = v[0];
        dummy.X[i+4] = v[0] - epsilon;
        dummy.X[i+5] = v[0];
        dummy.X[i+6] = v[0];

        dummy.Y[i]   = v[1];
        dummy.Y[i+1] = v[1];
        dummy.Y[i+2] = v[1] + epsilon;
        dummy.Y[i+3] = v[1];
        dummy.Y[i+4] = v[1];
        dummy.Y[i+5] = v[1] - epsilon;
        dummy.Y[i+6] = v[1];

        dummy.Z[i]   = v[2];
        dummy.Z[i+1] = v[2];
        dummy.Z[i+2] = v[2];
        dummy.Z[i+3] = v[2] + epsilon;
        dummy.Z[i+4] = v[2];
        dummy.Z[i+5] = v[2];
        dummy.Z[i+6] = v[2] - epsilon;
        i += 7;
    }

    // The +/- epsilon probes can step just outside the packed block,
    // where the block's specialized tape isn't valid (a pruned branch
    // could win out there); walk up the tape stack until the whole
    // sample cloud is covered.
    float lo[3] = { INFINITY, INFINITY, INFINITY };
    float hi[3] = { -INFINITY, -INFINITY, -INFINITY };
    for (unsigned q = 0; q < dummy.voxels; ++q)
    {
        lo[0] = fmin(lo[0], dummy.X[q]);  hi[0] = fmax(hi[0], dummy.X[q]);
        lo[1] = fmin(lo[1], dummy.Y[q]);  hi[1] = fmax(hi[1], dummy.Y[q]);
        lo[2] = fmin(lo[2], dummy.Z[q]);  hi[2] = fmax(hi[2], dummy.Z[q]);
    }
    Tape* tape = tape_base_for_region(
            block_tape ? block_tape : deck_base(deck),
            lo[0], hi[0], lo[1], hi[1], lo[2], hi[2]);

    const float* out = tape_eval_r(tape, ctx, dummy);

    // Extract normals from the evaluated data.
    std::vector<Vec3f> normals;
    normals.reserve(contour.size());
    i = 0;
    for (size_t n = 0; n < contour.size(); ++n)
    {
        const float dx = (out[i+1] - out[i]) - (out[i+4] - out[i]);
        const float dy = (out[i+2] - out[i]) - (out[i+5] - out[i]);
        const float dz = (out[i+3] - out[i]) - (out[i+6] - out[i]);
        normals.push_back(Vec3f(dx, dy, dz).normalized());
        i += 7;
    }

    return normals;
}

// Mark that the first edge of the most recent triangle is swappable
// (as part of feature detection / extraction).
void Mesher::push_swappable_triangle(Tri t)
{
    auto found = swappable.find(edge(t.a, t.b));
    if (found != swappable.end())
    {
        Tri& other = triangles[found->second];
        other.b = t.c;
        t.b = other.c;
        triangles.push_back(t);
        swappable.erase(found);
    }
    else
    {
        triangles.push_back(t);

        // Store the new triangle's position, keyed on its reversed edge.
        swappable[edge(t.b, t.a)] = triangles.size() - 1;
    }

    // Adjust voxel_end so that it points to the first new triangle.
    if (voxel_end == NONE)
        voxel_end = triangles.size() - 1;
}

std::vector<uint32_t> Mesher::get_contour()
{
    const size_t ve = voxel_limit();

    // Find all of the singular edges in this fan
    // (edges that aren't shared between multiple triangles).
    std::unordered_set<uint64_t> valid_edges;
    for (size_t i = voxel_start; i != ve; ++i)
    {
        const Tri& t = triangles[i];

        if (valid_edges.count(edge(t.b, t.a)))
            valid_edges.erase(edge(t.b, t.a));
        else
            valid_edges.insert(edge(t.a, t.b));

        if (valid_edges.count(edge(t.c, t.b)))
            valid_edges.erase(edge(t.c, t.b));
        else
            valid_edges.insert(edge(t.b, t.c));

        if (valid_edges.count(edge(t.a, t.c)))
            valid_edges.erase(edge(t.a, t.c));
        else
            valid_edges.insert(edge(t.c, t.a));
    }

    std::unordered_set<uint32_t> in_fan;

    std::vector<uint32_t> contour = {triangles[voxel_start].a};
    in_fan.insert(triangles[voxel_start].a);
    in_fan.insert(triangles[voxel_start].b);
    in_fan.insert(triangles[voxel_start].c);

    fan_start = voxel_start;
    voxel_start++;

    while (contour.size() == 1 || contour.front() != contour.back())
    {
        size_t i;
        for (i = fan_start; i != ve; ++i)
        {
            const Tri& t = triangles[i];
            if (contour.back() == t.a && valid_edges.count(edge(t.a, t.b)))
            {
                contour.push_back(t.b);
                break;
            }

            if (contour.back() == t.b && valid_edges.count(edge(t.b, t.c)))
            {
                contour.push_back(t.c);
                break;
            }

            if (contour.back() == t.c && valid_edges.count(edge(t.c, t.a)))
            {
                contour.push_back(t.a);
                break;
            }
        }
        // If we broke out of the loop (meaning i points to a relevant
        // triangle which should be moved to the back of the fan region
        // [fan_start, voxel_start)), rotate it into place and update
        // bookmarks appropriately.  Only elements between the triangle's
        // old and new positions move, so no other stored index is
        // disturbed.
        if (i != ve)
        {
            in_fan.insert(triangles[i].a);
            in_fan.insert(triangles[i].b);
            in_fan.insert(triangles[i].c);

            if (i == voxel_start)
            {
                voxel_start++;
            }
            else if (i > voxel_start)
            {
                // Move a triangle from beyond the fan region: it lands at
                // voxel_start's old position and the region grows by one.
                std::rotate(triangles.begin() + voxel_start,
                            triangles.begin() + i,
                            triangles.begin() + i + 1);
                voxel_start++;
            }
            else if (i != fan_start)
            {
                // Reorder within the fan region: the triangle moves to the
                // back of the region, which doesn't change its extent.
                std::rotate(triangles.begin() + i,
                            triangles.begin() + i + 1,
                            triangles.begin() + voxel_start);
            }
        }
    }

    // Special case to catch triangles that are part of a particular fan but
    // don't have any edges in the contour (which can happen!).
    for (size_t i = voxel_start; i != ve; ++i)
    {
        const Tri& t = triangles[i];
        if (in_fan.count(t.a) && in_fan.count(t.b) && in_fan.count(t.c))
        {
            if (i == voxel_start)
            {
                voxel_start++;
            }
            else if (i != fan_start)
            {
                std::rotate(triangles.begin() + voxel_start,
                            triangles.begin() + i,
                            triangles.begin() + i + 1);
                voxel_start++;
            }
        }
    }

    // Remove the last point of the contour, since it's a closed loop.
    contour.pop_back();
    return contour;
}

void Mesher::check_feature()
{
    auto contour = get_contour();

    // A degenerate contour can't contain a feature, and the math below
    // (center /= size, least-squares fit) misbehaves on empty input.
    if (contour.size() < 2)
        return;

    const auto normals = get_normals(contour);

    // Find the largest cone and the normals that enclose
    // the largest angle as n0, n1.
    float theta = 1;
    Vec3f n0, n1;
    for (auto ni : normals)
    {
        for (auto nj : normals)
        {
            float dot = ni.dot(nj);
            if (dot < theta)
            {
                theta = dot;
                n0 = ni;
                n1 = nj;
            }
        }
    }

    // If there isn't a feature in this fan, then return immediately.
    if (theta > 0.9)
        return;

    // Decide whether this is a corner or edge feature.
    const Vec3f nstar = n0.cross(n1);
    float phi = 0;
    for (auto n : normals)
        phi = fmax(phi, fabs(nstar.dot(n)));
    bool edge_feature = phi < 0.7;

    // Find the center of the contour.
    Vec3f center(0, 0, 0);
    for (auto c : contour)
        center += pos(c);
    center /= contour.size();

    // Construct the matrices for use in our least-square fit.
    Eigen::MatrixX3d A(normals.size(), 3);
    {
        int i=0;
        for (auto n : normals)
            A.row(i++) << n.transpose();
    }

    // When building the second matrix, shift position values to be centered
    // about the origin (because that's what the least-squares fit will
    // minimize).
    Eigen::VectorXd B(normals.size(), 1);
    {
        for (size_t i = 0; i < normals.size(); ++i)
            B.row(i) << normals[i].dot(pos(contour[i]) - center);
    }

    // Use singular value decomposition to solve the least-squares fit.
    Eigen::JacobiSVD<Eigen::MatrixX3d> svd(A, Eigen::ComputeFullU |
                                              Eigen::ComputeFullV);

    // Set the smallest singular value to zero to make fitting happier.
    if (edge_feature)
    {
        auto singular = svd.singularValues();
        svd.setThreshold(singular.minCoeff() / singular.maxCoeff() * 1.01);
    }

    // Solve for the new point's position.
    const Vec3f new_pt = svd.solve(B) + center;

    // Erase this triangle fan, as we'll be inserting a vertex in the center.
    // (Tombstoned rather than removed, so every stored index stays valid.)
    for (size_t i = fan_start; i != voxel_start; ++i)
        triangles[i].a = UINT32_MAX;

    // Construct a new triangle fan.
    const uint32_t np = intern(float(new_pt[0]),
                               float(new_pt[1]),
                               float(new_pt[2]));
    contour.push_back(contour.front());
    for (size_t i = 0; i + 1 < contour.size(); ++i)
        push_swappable_triangle(Tri{contour[i], contour[i+1], np});
}

void Mesher::remove_dupes()
{
    // Key each live triangle on its sorted vertex indices, tagged with
    // its position so that the first copy is the one that survives.
    struct Key {
        uint32_t v[3];
        size_t seq;
    };
    std::vector<Key> keys;
    keys.reserve(triangles.size());

    for (size_t i = 0; i < triangles.size(); ++i)
    {
        const Tri& t = triangles[i];
        if (dead(t))
            continue;
        Key k = {{t.a, t.b, t.c}, i};
        std::sort(std::begin(k.v), std::end(k.v));
        keys.push_back(k);
    }

    std::sort(keys.begin(), keys.end(), [](const Key& a, const Key& b) {
        if (a.v[0] != b.v[0])   return a.v[0] < b.v[0];
        if (a.v[1] != b.v[1])   return a.v[1] < b.v[1];
        if (a.v[2] != b.v[2])   return a.v[2] < b.v[2];
        return a.seq < b.seq;
    });

    for (size_t i = 1; i < keys.size(); ++i)
        if (keys[i].v[0] == keys[i-1].v[0] &&
            keys[i].v[1] == keys[i-1].v[1] &&
            keys[i].v[2] == keys[i-1].v[2])
        {
            triangles[keys[i].seq].a = UINT32_MAX;
            // Keep comparing against the run's first (surviving) entry
            keys[i] = keys[i-1];
        }
}

void Mesher::prune_flags()
{
    // Sorted flat list of every directed edge (duplicates are harmless:
    // this is only used for membership tests).
    std::vector<uint64_t> edges;
    edges.reserve(triangles.size() * 3);
    for (const auto& t : triangles)
    {
        if (dead(t))
            continue;
        edges.push_back(edge(t.a, t.b));
        edges.push_back(edge(t.b, t.c));
        edges.push_back(edge(t.c, t.a));
    }
    std::sort(edges.begin(), edges.end());

    const auto has = [&edges](uint64_t e) {
        return std::binary_search(edges.begin(), edges.end(), e);
    };

    for (auto& t : triangles)
    {
        if (dead(t))
            continue;
        if (!has(edge(t.b, t.a)) ||
            !has(edge(t.c, t.b)) ||
            !has(edge(t.a, t.c)))
        {
            t.a = UINT32_MAX;
        }
    }
}

void Mesher::finalize()
{
    if (finalized)
        return;
    finalized = true;

    // Nothing is interned or swapped after meshing; drop the maps.
    vertex_ids = decltype(vertex_ids)();
    swappable = decltype(swappable)();

    if (detect_edges)
    {
        remove_dupes();
        prune_flags();
    }
}

// Loads a vertex into the vertex list.
// If this vertex completes a triangle, check for features.
void Mesher::push_vert(const float x, const float y, const float z)
{
    triangle.push_back(intern(x, y, z));
    if (triangle.size() == 3)
    {
        triangles.push_back(Tri{triangle[0], triangle[1], triangle[2]});
        triangle.clear();

        // If this is the first triangle being constructed (or voxel_start
        // has just been cleared), store its position so that we know where
        // the next triangle fan begins.
        if (voxel_start == NONE)
            voxel_start = triangles.size() - 1;
    }
}


// Evaluates a region voxel-by-voxel, storing the output in the data
// member of the tristate struct.
bool Mesher::load_packed(const Region& r, Tape* tape)
{
    // Only load the packed matrix if we have few enough voxels.
    const unsigned voxels = (r.ni+1) * (r.nj+1) * (r.nk+1);
    if (voxels >= MIN_VOLUME)
        return false;

    // Zero-crossing / normal refinement over this block must run at
    // the same specialization, so remember the tape (borrowed: the
    // recursion frame that loaded the block outlives its unload).
    block_tape = tape;

    // Flatten a 3D region into a 1D list of points that
    // touches every point in the region, one by one.
    int q = 0;
    for (unsigned k=0; k <= r.nk; ++k) {
        for (unsigned j=0; j <= r.nj; ++j) {
            for (unsigned i=0; i <= r.ni; ++i) {
                X[q] = r.X[i];
                Y[q] = r.Y[j];
                Z[q] = r.Z[k];
                q++;
            }
        }
    }

    // Make a dummy region that has the newly-flattened point arrays as the
    // X, Y, Z coordinate data arrays (so that we can run eval_r on it).
    packed.imin = r.imin;
    packed.jmin = r.jmin;
    packed.kmin = r.kmin;
    packed.ni = r.ni;
    packed.nj = r.nj;
    packed.nk = r.nk;
    packed.X = X;
    packed.Y = Y;
    packed.Z = Z;
    packed.voxels = voxels;

    // Run the region evaluator and copy the data out
    memcpy(data, tape_eval_r(tape, ctx, packed), voxels * sizeof(float));
    has_data = true;

    return true;
}

void Mesher::unload_packed()
{
    block_tape = nullptr;
    has_data = false;
}

bool Mesher::get_corner_data(const Region& r, float d[8])
{
    bool has_positive = false;
    bool has_negative = false;
    // Populates an 8-element array with the function evaluation
    // results from the corner of a single-voxel region.
    for (int i=0; i < 8; ++i)
    {
        // Figure out where this bit of data lives in the larger eval_r array.
        const unsigned index =
            (r.imin - packed.imin + ((i & 4) ? r.ni : 0)) +
            (r.jmin - packed.jmin + ((i & 2) ? r.nj : 0))
                * (packed.ni+1) +
            (r.kmin - packed.kmin + ((i & 1) ? r.nk : 0))
                * (packed.ni+1) * (packed.nj+1);

        d[i] = data[index];

        has_negative |= d[i] < 0;
        has_positive |= d[i] >= 0;
    }

    return has_positive && has_negative;
}

void Mesher::eval_zero_crossings(Vec3f* v0, Vec3f* v1, unsigned count)
{
    float p[count];
    for (unsigned i=0; i < count; ++i)
        p[i] = 0.5;

    float step = 0.25;

    Region dummy;
    dummy.X = ex;
    dummy.Y = ey;
    dummy.Z = ez;
    dummy.voxels = count;

    for (int iteration=0; iteration < 8; ++iteration)
    {
        // Load new data into the x, y, z arrays.
        for (unsigned i=0; i < count; i++)
        {
            dummy.X[i] = v0[i][0] * (1 - p[i]) + v1[i][0] * p[i];
            dummy.Y[i] = v0[i][1] * (1 - p[i]) + v1[i][1] * p[i];
            dummy.Z[i] = v0[i][2] * (1 - p[i]) + v1[i][2] * p[i];
        }
        const float* out = tape_eval_r(
                block_tape ? block_tape : deck_base(deck), ctx, dummy);

        for (unsigned i=0; i < count; i++)
            if      (out[i] < 0)    p[i] += step;
            else if (out[i] > 0)    p[i] -= step;

        step /= 2;
    }
}

// Flushes out a queue of interpolation commands
void Mesher::flush_queue()
{
    Vec3f low[MIN_VOLUME];
    Vec3f high[MIN_VOLUME];

    // Go through the list, saving a list of vertex pairs on which
    // interpolation should be run into low and high.
    unsigned count=0;
    for (auto c : queue)
    {
        if (c.cmd == InterpolateCommand::INTERPOLATE)
        {
            low[count] = c.v0;
            high[count] = c.v1;
            count++;
        }
    }

    if (count)
        eval_zero_crossings(low, high, count);

    // Next, go through and actually load vertices
    // (either directly or from the cache)
    count = 0;
    for (auto c : queue)
    {
        if (c.cmd == InterpolateCommand::INTERPOLATE)
        {
            push_vert(ex[count], ey[count], ez[count]);
            count++;
        }
        else if (c.cmd == InterpolateCommand::CACHED)
        {
            unsigned i = c.cached;
            push_vert(ex[i], ey[i], ez[i]);
        }
        else if (c.cmd == InterpolateCommand::END_OF_VOXEL)
        {
            if (detect_edges)
            {
                // Clear voxel_end
                // (it will be reset when the next triangle is pushed)
                voxel_end = NONE;

                // Then, iterate until no more features are found in
                // the current voxel.
                while (voxel_start != NONE &&
                       voxel_start != voxel_limit() &&
                       voxel_start != triangles.size())
                {
                    check_feature();
                }

                // Clear voxel_start
                // (it will be reset when the next triangle is pushed)
                voxel_start = NONE;
            }
        }
    }
    queue.clear();
}

// Schedule an interpolate calculation in the queue.
void Mesher::interpolate_between(const Vec3f& v0, const Vec3f& v1)
{
    InterpolateCommand next = (InterpolateCommand){
        .cmd=InterpolateCommand::INTERPOLATE, .v0=v0, .v1=v1};

    // Walk through the list, looking for duplicates.
    // If we find the same operation, then switch to a CACHED lookup instead.
    unsigned count = 0;
    for (auto c : queue)
    {
        if (c.cmd == InterpolateCommand::INTERPOLATE)
        {
            if ((v0 == c.v0 && v1 == c.v1) || (v0 == c.v1 && v1 == c.v0))
            {
                next.cmd = InterpolateCommand::CACHED;
                next.cached = count;
            }
            count++;
        }
    }

    queue.push_back(next);
    if (next.cmd == InterpolateCommand::INTERPOLATE && count + 1 == MIN_VOLUME)
        flush_queue();
}


void Mesher::triangulate_tet(const Region& r, const float* const d,
                             const int t)
{
    // Find vertex positions for this tetrahedron
    uint8_t vertices[] = {0, 7, VERTEX_LOOP[t], VERTEX_LOOP[t+1]};

    // Figure out which of the sixteen possible combinations
    // we're currently experiencing.
    uint8_t lookup = 0;
    for (int v=3; v>=0; --v) {
        lookup = (lookup << 1) + (d[vertices[v]] < 0);
    }

    // Iterate over (up to) two triangles in this tetrahedron
    for (int i=0; i < 2; ++i)
    {
        if (EDGE_MAP[lookup][i][0][0] == -1)    break;

        // ...and insert vertices into the mesh.
        for (int v=0; v < 3; ++v)
        {
            const uint8_t v0 = vertices[EDGE_MAP[lookup][i][v][0]];
            const uint8_t v1 = vertices[EDGE_MAP[lookup][i][v][1]];

            interpolate_between(
                        (Vec3f){(v0 & 4) ? r.X[1] : r.X[0],
                                (v0 & 2) ? r.Y[1] : r.Y[0],
                                (v0 & 1) ? r.Z[1] : r.Z[0]},
                        (Vec3f){(v1 & 4) ? r.X[1] : r.X[0],
                                (v1 & 2) ? r.Y[1] : r.Y[0],
                                (v1 & 1) ? r.Z[1] : r.Z[0]});
        }
    }
}

void Mesher::triangulate_voxel(const Region& r, const float* const d)
{
    for (int t=0; t < 6; ++t)
        triangulate_tet(r, d, t);
}

void Mesher::triangulate_region(const Region& r)
{
    triangulate_region(r, deck_base(deck));
}

void Mesher::triangulate_region(const Region& r, Tape* tape)
{
    // Early abort if the halt flag is set
    if (*halt)
        return;

    // Do a round of interval evaluation to skip empty regions.
    const Interval X = {r.X[0], r.X[r.ni]},
                   Y = {r.Y[0], r.Y[r.nj]},
                   Z = {r.Z[0], r.Z[r.nk]};
    auto interval = tape_eval_i(tape, ctx, X, Y, Z);
    if (interval.lower > 0 || interval.upper < 0)
    {
        // Count skipped voxels toward progress, unless we're inside a
        // packed block (those are counted wholesale when they finish).
        if (progress && !has_data)
            progress->fetch_add(r.voxels, std::memory_order_relaxed);
        return;
    }

    // Specialize the tape against the interval results: decided
    // min/max branches drop out for the whole subtree below here
    // (values are exact, so the mesh is unchanged - the old code
    // only pruned once per packed block).  Inside a packed block the
    // field is already sampled and the tape has already flattened
    // out (see STIBIUM_TAPE_STATS), so pushing further is pure churn.
    Tape* sub = has_data ? tape_retain(tape)
                         : tape_push(tape, ctx, X, Y, Z,
                                     TAPE_PUSH_STANDARD);

    // If we can calculate all of the points in this region with a single
    // region-eval call, then do so.  This large chunk will be used in
    // future recursive calls to make things more efficient.
    bool loaded_data;
    if (!has_data)
        loaded_data = load_packed(r, sub);
    else
        loaded_data = false;

    // If we have greater than one voxel, subdivide and recurse.
    if (r.voxels > 1)
    {
        Region octants[8];
        const uint8_t split = octsect(r, octants);
        for (int i=0; i < 8; ++i)
            if (split & (1 << i))
                triangulate_region(octants[i], sub);
    }
    else
    {
        // Load corner values from this voxel
        // (from the packed data array)
        float d[8];

        if (get_corner_data(r, d))
        {
            // Triangulate this particular voxel
            triangulate_voxel(r, d);

            // Mark that a voxel has ended
            // (which triggers mesh refinment)
            queue.push_back((InterpolateCommand){
                    .cmd=InterpolateCommand::END_OF_VOXEL});
        }
    }

    // If this stage of the recursion loaded data into the buffer,
    // flush the interpolation queue (which evaluates with the block's
    // tape) and clear the has_data flag, so that future stages will
    // re-run the region evaluator on their portion of the space.
    if (loaded_data)
    {
        flush_queue();
        unload_packed();
        if (progress)
            progress->fetch_add(r.voxels, std::memory_order_relaxed);
    }

    tape_release(sub);
}

float* Mesher::get_verts(unsigned* count)
{
    finalize();

    size_t live = 0;
    for (const auto& t : triangles)
        if (!dead(t))
            live++;

    // There are 9 floats in each triangle
    *count = live * 9;

    float* out = (float*)malloc(sizeof(float) * (*count));

    unsigned i = 0;
    for (const auto& t : triangles)
    {
        if (dead(t))
            continue;
        for (auto v : {t.a, t.b, t.c})
            for (int j=0; j < 3; ++j)
                out[i++] = vertices[v][j];
    }

    return out;
}

void Mesher::get_mesh(std::vector<float>& out_verts,
                      std::vector<uint32_t>& out_indices,
                      bool run_finalize)
{
    if (run_finalize)
    {
        finalize();
    }
    else
    {
        // Still drop the meshing-time maps, but leave dedup / pruning
        // to the caller (see the header comment).
        vertex_ids = decltype(vertex_ids)();
        swappable = decltype(swappable)();
    }

    out_verts.clear();
    out_indices.clear();

    // Compact the vertex table down to referenced vertices only,
    // numbered in first-use order.
    std::vector<uint32_t> remap(vertices.size(), UINT32_MAX);
    for (const auto& t : triangles)
    {
        if (dead(t))
            continue;
        for (auto v : {t.a, t.b, t.c})
        {
            if (remap[v] == UINT32_MAX)
            {
                remap[v] = out_verts.size() / 3;
                out_verts.push_back(vertices[v][0]);
                out_verts.push_back(vertices[v][1]);
                out_verts.push_back(vertices[v][2]);
            }
            out_indices.push_back(remap[v]);
        }
    }
}

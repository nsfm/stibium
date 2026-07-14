#ifndef MESHER_H
#define MESHER_H

#include <array>
#include <atomic>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "fab/tree/triangulate/triangle.h"
#include "fab/tree/tape.h"

#include "fab/util/region.h"

struct InterpolateCommand {
    enum {INTERPOLATE, CACHED, END_OF_VOXEL} cmd;
    Vec3f v0;
    Vec3f v1;
    unsigned cached;
};

class Mesher {
public:
    /*
     *  progress, if given, is incremented by the number of voxels
     *  completed (empty regions count when skipped; occupied regions
     *  when their packed block finishes).  It reaches the region's
     *  total voxel count when triangulate_region returns unhalted.
     *
     *  The deck is borrowed (the caller frees it after meshing); it
     *  is immutable, so any number of Meshers can share one deck -
     *  each brings its own evaluation workspace.
     */
    Mesher(const Deck* deck, bool detect_edges, volatile int* halt,
           std::atomic<uint64_t>* progress=nullptr);
    ~Mesher();

    /*
     *  Recursively triangulates the given voxel region.
     */
    void triangulate_region(const Region& r);

    /*
     *  Allocates memory (using malloc) and returns a flat set of vertices.
     *  count is set to the number of floats allocated
     *  (i.e. number of vertices * 3)
     */
    float* get_verts(unsigned* count);

    /*
     *  Returns the mesh in indexed form: out_verts is a flat list of
     *  unique vertices (3 floats each); out_indices is 3 vertex indices
     *  per triangle.  Vertices are emitted in first-use order.
     *
     *  run_finalize=false skips the dedup / flag-pruning passes; used
     *  by the multithreaded driver, which must run them globally after
     *  merging (pruning per-chunk would drop valid triangles whose
     *  reverse edges live in a neighboring chunk's mesh).
     */
    void get_mesh(std::vector<float>& out_verts,
                  std::vector<uint32_t>& out_indices,
                  bool run_finalize=true);

protected:
    // Sentinel for the triangle-range bookmarks below (the equivalent of
    // a list.end() iterator: "no position yet" / "runs to the tail").
    static constexpr size_t NONE = SIZE_MAX;

    /*
     *  Returns the vertex table index for (x, y, z), creating a new
     *  entry if this position hasn't been seen before.  Positions are
     *  keyed on their exact float bit patterns.
     */
    uint32_t intern(float x, float y, float z);

    /*
     *  Looks up an interned vertex's position.
     */
    Vec3f pos(uint32_t v) const;

    /*
     *  Packs a directed edge into a single map key.
     */
    static uint64_t edge(uint32_t a, uint32_t b)
        { return (uint64_t(a) << 32) | b; }

    /*
     *  An erased (tombstoned) triangle.  Erasing in place keeps every
     *  triangle index stable, so bookmarks and the swappable-edge map
     *  never need fixing up; dead triangles are skipped on output.
     */
    static bool dead(const Tri& t) { return t.a == UINT32_MAX; }

    /*
     *  Upper bound of the current voxel's triangle range
     *  (voxel_end, or the current tail if voxel_end is unset).
     */
    size_t voxel_limit() const
        { return voxel_end == NONE ? triangles.size() : voxel_end; }

    /*
     *  Finds the normals of each vertex on the contour.
     */
    std::vector<Vec3f> get_normals(const std::vector<uint32_t>& contour);

    /*
     *  Records another vertex.
     *  Calls add_triangle if this vertex completes a triangle.
     */
    void push_vert(const float x, const float y, const float z);

    /*
     *  The recursive worker behind triangulate_region: tape is the
     *  (already-specialized) tape for this region's ancestors, and is
     *  pushed further at each ambiguous subdivision.
     */
    void triangulate_region(const Region& r, Tape* tape);

    /*
     *  Attempts to evaluate every voxel in the given region.
     *  Returns false if there are too many voxels; true on success.
     *  On success the tape is remembered as block_tape, so that
     *  zero-crossing / normal refinement over this block evaluates
     *  with the same specialization (the old code disabled nodes for
     *  exactly this scope).
     */
    bool load_packed(const Region& r, Tape* tape);
    void unload_packed();

    /*
     *  Looks up the corner values for the given region, storing them in d.
     *  Returns true if this voxel has anything of interest in it.
     */
    bool get_corner_data(const Region& r, float d[8]);

    /*
     *  Performs binary search on a set of edges.
     *  v0 and v1 must contain 'count' points.
     *  Found x, y, z value are stored in ex, ey, ez.
     */
    void eval_zero_crossings(Vec3f* v0, Vec3f* v1, unsigned count);

    /*
     *  Flushes the command queue.
     *  This will involve calculating a set of interpolated positions
     *  then assembling them into triangles.
     */
    void flush_queue();

    /*
     *  Schedules an interpolation command in the command queue.
     */
    void interpolate_between(const Vec3f& v0, const Vec3f& v1);

    /*
     *  Evaluates the given voxel.
     *      r is the voxel region
     *      d is the corner values
     */
    void triangulate_voxel(const Region& r, const float* const d);

    /*
     *  Evaluates the given tetrahedron in a voxel.
     *      r is the voxel region
     *      d is the corner values
     *      t is the tetrahedron's ID.
     */
    void triangulate_tet(const Region& r, const float* const d, const int t);

    /*
     *  Marks that the first edge of this triangle is swappable,
     *  and performs the swap if a match is found.
     */
    void push_swappable_triangle(Tri t);

    /*
     *  Check the most recent fan (from voxel_start to voxel_limit())
     *  for features and process them if they are found.
     */
    void check_feature();

    /*
     *  Removes duplicate triangles (same three vertices, any order).
     */
    void remove_dupes();

    /*
     *  Removes triangles with edges that aren't connected to the
     *  rest of the mesh (which happens sometimes when refining geometry).
     */
    void prune_flags();

    /*
     *  One-time cleanup before output: drops the interning and
     *  swappable-edge maps, then runs remove_dupes / prune_flags
     *  when feature detection is enabled.
     */
    void finalize();

    /*
     *  Returns a closed contour that traces the most recent fan.
     *
     *  Modifies triangles, voxel_start, and fan_start so that the
     *  most recent fan is stored between fan_start and voxel_start.
     */
    std::vector<uint32_t> get_contour();

    // Compiled evaluation deck (borrowed, shared) and our workspace
    const Deck* deck;
    TapeCtx* ctx;

    // Tape of the packed block currently being refined (borrowed from
    // the recursion frame that loaded it; null outside a block)
    Tape* block_tape;

    bool detect_edges;
    volatile int* halt;
    std::atomic<uint64_t>* progress;

    // Cached region and data from an eval_r call
    Region packed;
    float* data;
    bool has_data;

    // Buffers used for eval_r
    float* X;
    float* Y;
    float* Z;

    // Buffers used in eval_zero_crossings
    float* ex;
    float* ey;
    float* ez;

    // Buffers used in evaluating normals
    float* nx;
    float* ny;
    float* nz;

    // Queue of interpolation commands to be run soon.  A vector on
    // purpose: entries are only appended, scanned, and bulk-cleared,
    // and clear() keeps capacity across packed blocks (the old
    // std::list re-allocated a node per command - visible in export
    // profiles).
    std::vector<InterpolateCommand> queue;

    // Dedup index over the queue's INTERPOLATE entries: canonical
    // vertex pair -> INTERPOLATE ordinal.  Replaces the per-edge
    // linear queue walk (O(edges x queue) with an Eigen compare in
    // the middle - ~30% of export wall time on big models).  Keys
    // are bit-packed coords - Vec3f is Eigen::Vector3d despite the
    // name, so each is a full 8-byte word - with -0 normalized to
    // +0 so bitwise equality matches ==; NaN coords are never
    // entered (they match nothing under ==).
    struct PairHash {
        size_t operator()(const std::array<uint64_t, 6>& k) const {
            // FNV-1a over the six packed doubles
            uint64_t h = 0xcbf29ce484222325ull;
            for (uint64_t w : k) {
                h = (h ^ w) * 0x100000001b3ull;
            }
            return h;
        }
    };
    std::unordered_map<std::array<uint64_t, 6>, unsigned, PairHash>
        interp_cache;
    unsigned interp_count=0;   // INTERPOLATE entries now in queue

    // Triangle that's being constructed (vertex indices, up to 3)
    std::vector<uint32_t> triangle;

    // Interned vertex positions, and a bit-pattern lookup table for them
    std::vector<std::array<float, 3>> vertices;
    struct VertHash {
        size_t operator()(const std::array<uint32_t, 3>& k) const {
            // FNV-1a over the three packed floats
            uint64_t h = 0xcbf29ce484222325ull;
            for (uint32_t w : k) {
                h = (h ^ w) * 0x100000001b3ull;
            }
            return h;
        }
    };
    std::unordered_map<std::array<uint32_t, 3>, uint32_t, VertHash> vertex_ids;

    // Flat list of triangles (12 bytes each); erased ones are tombstoned
    std::vector<Tri> triangles;
    bool finalized;

    // Bookmarks into the triangle list (see get_contour / check_feature)
    size_t voxel_start;
    size_t voxel_end;
    size_t fan_start;
    std::unordered_map<uint64_t, size_t> swappable;
};

#endif

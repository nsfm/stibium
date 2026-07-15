#ifndef DELAUNAY_H
#define DELAUNAY_H

/*  The adaptive-Delaunay meshing campaign (doc/MESH-NEXT.md).
 *
 *  Stage A: point generation from the tape kernel - octree descent
 *  with STANDARD pushes, corner-lattice sampling in leaf blocks,
 *  batched bisection of sign-change lattice edges.  The interval
 *  bound gives the far field for free: a culled box contributes its
 *  corners with a PROVEN sign and no pointwise evaluation at all.
 *
 *  Stage B/C (CGAL side, delaunay_mesh): incremental Delaunay over
 *  all points, opposite-sign tet edges refined to convergence,
 *  triangles emitted between inside/outside tets across surface-
 *  corner faces.  Faces of a tet complex cannot self-intersect, and
 *  at convergence the surface is the boundary of the inside-tet
 *  union - closed by construction.  (Kobbelt-style feature points
 *  and interval-driven hidden-feature drill-down are later stages;
 *  see MESH-NEXT.)
 */

#include <cstdint>
#include <vector>

#include "fab/tree/tape.h"
#include "fab/util/region.h"

/*  A signed sample: pointwise-evaluated lattice corner, or a culled
 *  octree box corner whose sign the interval bound proves.  */
struct DSample
{
    float x, y, z;
    bool inside;      // f < 0, strictly (matches the renderer)
    bool on_surface;  // f == 0 exactly: grid-aligned geometry lands
                      // lattice samples ON the surface; they enter
                      // the triangulation as surface vertices
};

/*  A surface point from bisecting a sign-change lattice edge.  */
struct DSurfPoint
{
    float x, y, z;
};

struct DSoup
{
    std::vector<DSample> samples;
    std::vector<DSurfPoint> surface;

    /*  Diagnostics  */
    uint64_t leaf_blocks = 0;
    uint64_t culled_empty = 0, culled_full = 0;
    /*  Interval says the surface may be here, but every sample in
     *  the block agrees on sign: a candidate hidden thin feature.
     *  Stage-D drill-down hook; counted from day one.  */
    uint64_t hidden_candidates = 0;
    /*  QEF-placed sharp-feature points appended to `surface`.  */
    uint64_t feature_points = 0;
    /*  Crossings replaced by their cell's feature point.  */
    uint64_t suppressed = 0;
    /*  Outlier feature points mediated back onto their chain.  */
    uint64_t mediated = 0;
    /*  Finest lattice pitch (for downstream keep-out radii).  */
    float spacing = 0;

    /*  Crease-tracer output (delaunay_trace): exactly-on-crease
     *  points appended to `surface`, chained in traced order.  When
     *  non-empty the mesher constrains these and skips the QEF
     *  chain extractor.  */
    std::vector<std::vector<uint32_t>> tchains;
    std::vector<uint8_t> tclosed;    // parallel to tchains
    uint64_t traced = 0;
};

/*  Stage A: collect the point soup for a region.  */
DSoup delaunay_sample(const Deck* deck, Region r, volatile int* halt);

/*  The crease tracer: marches every min/max clause's crease
 *  {f_A = 0, f_B = 0} with an SSI predictor-corrector on the tape's
 *  prefix evaluation, trimmed to the boundary by the full oracle,
 *  junction endpoints bisected onto their corners and shared.
 *  Seeds come from the soup's QEF feature points; traced polylines
 *  land in soup->tchains.  Returns false when nothing was traced
 *  (no min/max clauses, no seeds, or nothing converged) - the
 *  caller falls back to delaunay_chains.  */
bool delaunay_trace(const Deck* deck, Region r, DSoup* soup,
                    volatile int* halt);

/*  Crease chains: feature points ordered into polylines by a
 *  radius graph (1.6 cells).  Degree-2 runs are chains; degree>=3
 *  vertices are junctions (cube corners) and terminate chains;
 *  all-degree-2 cycles are closed loops.  Indices refer into
 *  soup.surface (the feature tail).  Feeds the constrained-edge
 *  round (MESH-NEXT).  */
struct DChains
{
    std::vector<std::vector<uint32_t>> chains;
    std::vector<uint8_t> closed;   // parallel to chains
    uint64_t junctions = 0;
    uint64_t stray = 0;            // representatives in no chain
    uint64_t reps = 0;             // features after duplicate merge
};
DChains delaunay_chains(const DSoup& soup);

struct DMesh
{
    std::vector<float> verts;      // xyz triples
    std::vector<uint32_t> tris;    // index triples, CCW seen from outside

    /*  Diagnostics  */
    uint64_t iterations = 0;       // refinement rounds to convergence
    uint64_t inserted = 0;         // surface points added by refinement
    uint64_t open_edges = 0;       // boundary edges not shared by 2 tris
    uint64_t nonmanifold_edges = 0;
    uint64_t repaired = 0;         // wart midpoints projected+inserted
    uint64_t repair_rounds = 0;
    /*  Constrained-crease round (STIBIUM_DMESH_CCDT, default on):
     *  crease-chain segments inserted as constrained edges, and the
     *  Steiner vertices the conforming machinery placed on them
     *  (on-crease surface vertices by construction).  */
    uint64_t constrained = 0;
    uint64_t steiner = 0;
    /*  Manifold pass (STIBIUM_DMESH_MANIFOLD, default on): pinched
     *  fans split into one coincident vertex per surface sheet.  */
    uint64_t split_verts = 0;
};

/*  Stages B+C: point soup -> Delaunay -> refine -> extract.
 *  Returns false when built without CGAL support.  */
bool delaunay_mesh(const Deck* deck, Region r, volatile int* halt,
                   DMesh* out);

/*  Same, but over a caller-supplied soup (testing hook: hand it a
 *  soup with the surface points stripped and the refinement loop
 *  must rebuild the surface from signed samples alone).  */
bool delaunay_mesh_soup(const Deck* deck, const DSoup& soup,
                        volatile int* halt, DMesh* out);

#endif

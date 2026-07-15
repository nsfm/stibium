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
    bool inside;
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
};

/*  Stage A: collect the point soup for a region.  */
DSoup delaunay_sample(const Deck* deck, Region r, volatile int* halt);

struct DMesh
{
    std::vector<float> verts;      // xyz triples
    std::vector<uint32_t> tris;    // index triples, CCW seen from outside

    /*  Diagnostics  */
    uint64_t iterations = 0;       // refinement rounds to convergence
    uint64_t inserted = 0;         // surface points added by refinement
    uint64_t open_edges = 0;       // boundary edges not shared by 2 tris
    uint64_t nonmanifold_edges = 0;
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

#ifndef FAB_MESH_QUERY_H
#define FAB_MESH_QUERY_H

#include <memory>

#include "fab/mesh/trimesh.h"

namespace fab_mesh {

/*
 *  Accelerated geometric queries against a triangle mesh: exact
 *  unsigned distance via an AABB tree, and generalized winding
 *  number via a Barnes-Hut dipole hierarchy (Barill et al. 2018),
 *  so inside/outside classification degrades gracefully on meshes
 *  that aren't watertight.
 *
 *  Queries are const and safe to call concurrently from many
 *  threads once construction has finished.
 */
class MeshQuery {
public:
    /*  Builds the acceleration structures (single-threaded).
     *  The mesh is copied; the TriMesh may be freed afterwards. */
    explicit MeshQuery(const TriMesh& mesh);
    ~MeshQuery();

    MeshQuery(const MeshQuery&) = delete;
    MeshQuery& operator=(const MeshQuery&) = delete;

    /*  Exact distance to the closest point on any triangle. */
    float unsigned_distance(float x, float y, float z) const;

    /*  Generalized winding number: ~1 well inside a closed surface,
     *  ~0 well outside, fractional near holes and open boundaries. */
    float winding_number(float x, float y, float z) const;

    /*  unsigned_distance, negated where winding_number > 0.5 */
    float signed_distance(float x, float y, float z) const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl;
};

}  // namespace fab_mesh

#endif

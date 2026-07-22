/*  Stibium <-> libfive bridge (Stibium's own work, AGPL-3.0-or-later).
 *
 *  Translates a Stibium MathTree (the `Node` graph produced by
 *  deck_from_tree's input tree) into a libfive::Tree, then meshes it with
 *  libfive's flagship Dual Contouring mesher for a head-to-head against the
 *  stibnite adaptive-Delaunay mesher.
 *
 *  This header is deliberately Eigen-free so it can be included from the
 *  app's translation units (which use Stibium's vendored Eigen 3.2.4) while
 *  the .cpp compiles against libfive's Eigen (system, >= 3.3) in an isolated
 *  shared library.  Only plain arrays cross the boundary.
 */
#ifndef STIBIUM_LIBFIVE_BRIDGE_H
#define STIBIUM_LIBFIVE_BRIDGE_H

#include <cstdint>
#include <vector>
#include <string>

#if defined(__GNUC__) || defined(__clang__)
#define STIBIUM_LF_EXPORT __attribute__((visibility("default")))
#else
#define STIBIUM_LF_EXPORT
#endif

struct Node_;

namespace stibium_libfive {

/*  Meshes `root` (a Stibium Node graph) with libfive Dual Contouring.
 *
 *  bounds are the axis-aligned meshing region; min_feature is libfive's
 *  smallest-cell-edge target (use 1/resolution to match a Stibium res);
 *  max_err gates DC cell collapse (libfive default 1e-8; -1 disables).
 *
 *  On success returns 0 and fills verts (xyz float triples) and tris
 *  (uint32 index triples, referencing verts).  On failure returns non-zero
 *  and puts a human-readable reason in `err` (e.g. an unsupported opcode).
 */
STIBIUM_LF_EXPORT
int mesh_shape(const Node_* root,
               double x0, double y0, double z0,
               double x1, double y1, double z1,
               double min_feature, double max_err,
               unsigned workers,
               std::vector<float>& verts,
               std::vector<uint32_t>& tris,
               std::string& err);

}   // namespace stibium_libfive

#endif

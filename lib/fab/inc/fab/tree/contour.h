#ifndef CONTOUR_H
#define CONTOUR_H

#include <array>
#include <atomic>
#include <cstdint>
#include <vector>

struct MathTree_;

/*
 *  A closed 2D contour: xy points in order, implicitly closed
 *  (last connects back to first).
 */
typedef std::vector<std::array<float, 2>> ContourPath;

/*
 *  Traces the zero contours of the tree on the z plane over the box
 *  [xmin, xmax] x [ymin, ymax], sampled on an (nx+1) x (ny+1) grid
 *  (marching squares; saddle cells are resolved by evaluating the
 *  cell center).  The sample grid is padded with empty space, so
 *  shapes clipped by the box produce closed loops along its border.
 *
 *  Loops are oriented with the filled region on the left: outer
 *  boundaries wind counterclockwise, holes clockwise.
 *
 *  detect_features recovers sharp corners: where the field normals at
 *  a cell's two crossing points diverge, the crossing chord is split
 *  at the intersection of the two tangent lines (the 2D analog of the
 *  mesher's feature detection).
 *
 *  progress, if given, is incremented once per sample row; it reaches
 *  ny + 1 on completion.
 */
void contour_field(struct MathTree_* tree,
                   float xmin, float ymin, float xmax, float ymax,
                   uint32_t nx, uint32_t ny, float z,
                   bool detect_features, volatile int* halt,
                   std::vector<ContourPath>& paths,
                   std::atomic<uint64_t>* progress = nullptr);

#endif

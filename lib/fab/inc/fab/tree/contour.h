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
 *  Evaluation is spread across threads (<= 0 means hardware
 *  concurrency), each on its own clone of the tree.
 *
 *  Progress: done counts field evaluations as they complete; total is
 *  maintained by this function and grows as later evaluation passes
 *  (saddles, corner gradients) are sized.  done/total is always an
 *  honest fraction.
 */
void contour_field(struct MathTree_* tree,
                   float xmin, float ymin, float xmax, float ymax,
                   uint32_t nx, uint32_t ny, float z,
                   bool detect_features, volatile int* halt,
                   std::vector<ContourPath>& paths,
                   int threads = -1,
                   std::atomic<uint64_t>* progress_done = nullptr,
                   std::atomic<uint64_t>* progress_total = nullptr);

/*
 *  Simplifies each loop with Douglas-Peucker: points are removed
 *  while the path stays within tolerance of the original.  Recovered
 *  corners survive (they are maximal-deviation points by
 *  construction).  Loops degenerating below 3 points are dropped.
 */
void simplify_contours(std::vector<ContourPath>& paths, float tolerance);

#endif

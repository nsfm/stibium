#ifndef TRIANGULATE_H
#define TRIANGULATE_H

#include <atomic>
#include <cstdint>
#include <vector>

#include "fab/util/region.h"

/*
 *  Meshes the tree over the given region as raw triangle soup:
 *  9 floats (3 xyz corners) per triangle, malloc'd into *verts.
 *  Sets *count to the number of floats.
 */
void triangulate(struct MathTree_* tree, Region r,
                 bool detect_edges, volatile int* halt,
                 float** const verts, unsigned* const count);

/*
 *  Meshes the tree over the given region as an indexed mesh:
 *  verts holds unique vertices (3 floats each) and indices holds
 *  3 vertex indices per triangle.
 *
 *  progress, if given, counts completed voxels up to r.voxels.
 */
void triangulate_indexed(struct MathTree_* tree, Region r,
                         bool detect_edges, volatile int* halt,
                         std::vector<float>& verts,
                         std::vector<uint32_t>& indices,
                         std::atomic<uint64_t>* progress=nullptr);

/*
 *  Multithreaded triangulate_indexed: splits the region into chunks
 *  meshed by a pool of workers (each on its own clone of the tree),
 *  then merges the results.  Output triangles are the same set a
 *  single-threaded run produces when feature detection is off; with
 *  feature detection on, seam-adjacent edge swaps may resolve
 *  differently (both results are valid, watertight meshes).
 *
 *  threads <= 0 means hardware concurrency.
 */
void triangulate_indexed_mt(struct MathTree_* tree, Region r,
                            bool detect_edges, volatile int* halt,
                            std::vector<float>& verts,
                            std::vector<uint32_t>& indices,
                            int threads=-1,
                            std::atomic<uint64_t>* progress=nullptr);

#endif

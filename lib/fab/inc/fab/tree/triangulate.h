#ifndef TRIANGULATE_H
#define TRIANGULATE_H

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
 */
void triangulate_indexed(struct MathTree_* tree, Region r,
                         bool detect_edges, volatile int* halt,
                         std::vector<float>& verts,
                         std::vector<uint32_t>& indices);

#endif

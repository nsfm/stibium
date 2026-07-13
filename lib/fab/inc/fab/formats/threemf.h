#ifndef THREEMF_H
#define THREEMF_H

#include <cstdint>

/** Exports an indexed mesh to a .3mf file.
 *
 *  verts is an array of unique xyz vertices; indices holds three
 *  vertex indices per triangle; tri_count is the number of triangles.
 *  Units are declared as millimeters (models are unitless; slicers
 *  assume mm).
 *
 *  Returns true on success, false if the file couldn't be written.
 */
bool save_3mf_indexed(const float* verts, uint32_t vert_count,
                      const uint32_t* indices, uint32_t tri_count,
                      const char* filename);

#endif

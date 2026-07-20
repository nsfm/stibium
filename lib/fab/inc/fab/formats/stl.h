#ifndef STL_H
#define STL_H

#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Exports to a .stl file.
 *
 *  verts is an array of xyz verts packed into triangles.
 *  count is the number of floats in the array.
 */
void save_stl(float* verts, unsigned count, const char* filename);

/** Exports an indexed mesh to a .stl file.
 *
 *  verts is an array of unique xyz vertices; indices holds three
 *  vertex indices per triangle; tri_count is the number of triangles.
 *  (STL itself stores raw soup, so triangles are expanded on the fly.)
 */
void save_stl_indexed(const float* verts, const uint32_t* indices,
                      uint32_t tri_count, const char* filename);

/** Same, with a caller-supplied 80-byte header stamp (truncated /
 *  zero-padded; keep the "Stibium" prefix so mesh import still
 *  recognizes our exports).  */
void save_stl_indexed_stamped(const float* verts,
                              const uint32_t* indices,
                              uint32_t tri_count,
                              const char* filename,
                              const char* stamp);

#ifdef __cplusplus
}
#endif

#endif

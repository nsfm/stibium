#ifndef FAB_TREE_GRID_H
#define FAB_TREE_GRID_H

#include <stdint.h>

#include "fab/util/interval.h"

/*
 *  A MeshGrid is a dense signed-distance grid sampled from an
 *  imported mesh: the heap payload behind OP_GRID nodes.
 *
 *  Grids live in a process-wide registry keyed by an integer id;
 *  math strings reference them as "g<id>", so the strings stay
 *  plain text and the payload never travels through serialization.
 *  Nodes retain/release their grid, and the registry holds one
 *  reference of its own, so a grid survives re-parses for as long
 *  as it stays registered.
 *
 *  Evaluation is pure interpolation: trilinear for point and
 *  gradient queries, and precomputed 8^3-block sample min/max for
 *  conservative interval queries (so all interval-based pruning
 *  keeps working on imported geometry).  Points outside the sampled
 *  box evaluate to the clamped boundary sample plus the Euclidean
 *  distance to the box, which keeps the field positive and
 *  continuous outside.
 */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct MeshGrid_ MeshGrid;

/*  Borrowed pointer, or NULL if the id is unknown. */
MeshGrid* grid_lookup(uint32_t id);

uint32_t grid_id(const MeshGrid* g);

void grid_retain(MeshGrid* g);
void grid_release(MeshGrid* g);

/*  Point, gradient, and interval evaluation (thread-safe; grids are
 *  immutable after registration). */
float grid_eval_f(const MeshGrid* g, float x, float y, float z);
void  grid_eval_g(const MeshGrid* g, float x, float y, float z,
                  float* v, float* dx, float* dy, float* dz);
Interval grid_eval_i(const MeshGrid* g,
                     Interval X, Interval Y, Interval Z);

void grid_sampled_bounds(const MeshGrid* g,
                         float* xmin, float* ymin, float* zmin,
                         float* xmax, float* ymax, float* zmax);

#ifdef __cplusplus
}  // extern "C"

#include <string>
#include <vector>

/*
 *  Registration (C++ only; called from the import path).
 *
 *  samples is x-major: samples[i + ni*(j + nj*k)], sample (i,j,k) at
 *  (xmin + i*dx, ymin + j*dy, zmin + k*dz) with dx = (xmax-xmin)/(ni-1).
 *  Each axis needs at least 2 samples.  Returns the new grid's id
 *  (ids start at 1; 0 is never valid).
 *
 *  key identifies the content (source hash + sample parameters):
 *  grid_find returns an already-registered grid with the same key,
 *  so re-running an import script reuses the existing grid instead
 *  of re-sampling.
 */
uint32_t grid_register(std::vector<float>&& samples,
                       uint32_t ni, uint32_t nj, uint32_t nk,
                       float xmin, float ymin, float zmin,
                       float xmax, float ymax, float zmax,
                       const std::string& key);

/*  Id of the grid registered under key, or 0 if there is none. */
uint32_t grid_find(const std::string& key);

/*  Frees grids no longer referenced by any tree (refcount has
 *  dropped back to the registry's own reference).  Returns the
 *  number of grids freed.  Called before registering a new grid so
 *  abandoned imports (e.g. superseded resolutions) don't pile up. */
size_t grid_registry_trim();

/*  Bytes held by a grid's sample + block arrays (0 if unknown id). */
size_t grid_memory_bytes(uint32_t id);

#endif  // __cplusplus

#endif

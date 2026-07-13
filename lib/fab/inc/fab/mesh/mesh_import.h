#ifndef FAB_MESH_IMPORT_H
#define FAB_MESH_IMPORT_H

#include <cstdint>
#include <string>

namespace fab_mesh {

/*
 *  The import pipeline: STL file -> signed distance grid -> grid
 *  registry (see fab/tree/grid.h).  Sampling is the expensive step
 *  (BVH distance + winding number per voxel, parallel over planes),
 *  so results are reused at two levels:
 *
 *  - the registry, keyed by content-hash + parameters, catches
 *    re-runs of an import script within a session;
 *  - a cache file (optional, regenerable, safe to gitignore)
 *    catches re-opens of a project across sessions.
 */
struct ImportResult {
    uint32_t grid_id = 0;       // 0 means failure; see error
    std::string error;

    bool stibium_stamp = false; // source is one of our own exports
    bool from_cache = false;    // grid loaded, not sampled

    float bounds[6] = {0, 0, 0, 0, 0, 0};   // sampled (padded) box
    uint32_t dims[3] = {0, 0, 0};
    uint64_t tri_count = 0;
    std::string sha256;         // content hash of the source file
};

/*
 *  Imports path as a signed distance grid sampled at voxels_per_unit.
 *  cache_dir may be empty (no file caching).  threads <= 0 means
 *  half the hardware threads (physical cores, by this codebase's
 *  convention).
 */
ImportResult import_mesh_grid(const std::string& path,
                              float voxels_per_unit,
                              const std::string& cache_dir,
                              int threads = -1);

/*  Streaming SHA-256 of a file's contents (empty string on I/O
 *  error).  Exposed for hash-pinning checks in the import node. */
std::string sha256_file(const std::string& path);

}  // namespace fab_mesh

#endif

#ifndef FAB_MESH_TRIMESH_H
#define FAB_MESH_TRIMESH_H

#include <cstdint>
#include <string>
#include <vector>

namespace fab_mesh {

/*
 *  An indexed triangle mesh as loaded from a file, plus enough
 *  provenance (the raw header) for the importer to recognize
 *  Stibium's own exports and advise opening the source .sb instead.
 */
struct TriMesh {
    std::vector<float> verts;       // xyz triples
    std::vector<uint32_t> tris;     // vertex-index triples

    // xmin ymin zmin xmax ymax zmax; valid whenever tris is non-empty
    float bbox[6];

    // Binary STL: the raw 80-byte header (trailing NULs stripped).
    // ASCII STL: the text after "solid" on the first line.
    std::string header;

    uint32_t tri_count() const { return tris.size() / 3; }
    uint32_t vert_count() const { return verts.size() / 3; }
};

/*
 *  Loads a binary or ASCII STL file (sniffed by content, not name:
 *  files starting with "solid" that parse as ASCII are ASCII, all
 *  else is binary).  Vertices are deduplicated by exact coordinate
 *  match so shared corners weld into real indexed topology.
 *  Degenerate triangles (repeated vertex index) are dropped.
 *
 *  Returns true on success; on failure returns false and describes
 *  the problem in *err.
 */
bool load_stl(const std::string& path, TriMesh* out, std::string* err);

}  // namespace fab_mesh

#endif

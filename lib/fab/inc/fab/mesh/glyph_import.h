#ifndef FAB_MESH_GLYPH_IMPORT_H
#define FAB_MESH_GLYPH_IMPORT_H

#include <cstdint>
#include <string>

struct MathTree_;   // defined in fab/tree; used by redistance_grid below

namespace fab_mesh {

/*
 *  Bake one glyph of a TrueType font into a 3D signed-distance grid: the
 *  glyph outline is flattened to line segments, an exact 2D signed distance
 *  (min segment distance + winding-number sign) is combined with a z-slab to
 *  form an extruded solid, and the field is sampled onto a dense grid and
 *  registered in the OP_GRID registry (see fab/tree/grid.h) -- the same
 *  machinery mesh import uses. cap = capital-letter height in mm; the solid is
 *  extruded 0..thickness in z; voxels_per_mm sets sample density.
 *
 *  Result grids are content-keyed, so re-baking the same (font, glyph, cap,
 *  thickness, res) reuses the registered grid instead of re-sampling.
 */
struct GlyphResult {
    uint32_t grid_id = 0;   // 0 means failure; see error
    std::string error;
    float bounds[6] = {0, 0, 0, 0, 0, 0};
    uint32_t dims[3] = {0, 0, 0};
    bool from_cache = false;
};

GlyphResult bake_glyph_grid(const std::string& font_path, int codepoint,
                            float cap, float thickness, float voxels_per_mm);

/*  Advance width of a glyph in mm at the given cap height (for layout). */
float glyph_advance(const std::string& font_path, int codepoint, float cap);

/*
 *  Redistance: sample any Shape's field sign on a grid over the given (finite)
 *  bounds and rebuild it as a TRUE euclidean distance field (Felzenszwalb EDT),
 *  registered as OP_GRID. The general grid-bake keystone -- the Antimony font
 *  rides it (baking Matt Keeter's extruded glyph fields), it fixes offset/shell
 *  on non-exact fields, and it is the freeze-a-subtree performance primitive.
 */
GlyphResult redistance_grid(::MathTree_* tree,
                            float xmin, float ymin, float zmin,
                            float xmax, float ymax, float zmax,
                            float voxels_per_mm, const std::string& key);

/*
 *  Bake a 2D Shape (a laid-out glyph/word, e.g. Matt's Antimony font) into an
 *  extruded SDF grid: the field's zero level is marching-squares contoured for
 *  a smooth boundary, exact distance is taken to that contour, and the sign
 *  comes from the field itself. Clean walls (no EDT striping), same quality as
 *  the TTF path. (x0,y0,x1,y1) are the shape's finite XY bounds.
 */
GlyphResult bake_shape_glyph(::MathTree_* tree, float x0, float y0,
                             float x1, float y1, float thickness,
                             float voxels_per_mm, const std::string& key);

}  // namespace fab_mesh

#endif

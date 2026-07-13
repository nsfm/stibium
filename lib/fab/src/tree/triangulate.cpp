#include "fab/tree/triangulate.h"
#include "fab/tree/triangulate/mesher.h"

#include "fab/tree/tree.h"

// Finds an array of vertices (as x,y,z float triplets).
// Sets *count to the number of vertices returned.
void triangulate(MathTree* tree, const Region r,
                 bool detect_edges, volatile int* halt,
                 float** const verts, unsigned* const count)
{
    Mesher t(tree, detect_edges, halt);

    // Top-level call to the recursive triangulation function.
    t.triangulate_region(r);

    // Copy data from the mesher to the output pointers.
    *verts = t.get_verts(count);
}

void triangulate_indexed(MathTree* tree, const Region r,
                         bool detect_edges, volatile int* halt,
                         std::vector<float>& verts,
                         std::vector<uint32_t>& indices)
{
    Mesher t(tree, detect_edges, halt);
    t.triangulate_region(r);
    t.get_mesh(verts, indices);
}

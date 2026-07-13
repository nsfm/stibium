#include <stdlib.h>
#include <stdio.h>

#include "fab/formats/stl.h"

static void write_header(FILE* stl, uint32_t tris)
{
    // 80-character header
    fprintf(stl, "This is a binary STL file made in Antimony      \n(github.com/mkeeter/antimony)\n\n");

    // Little-endian uint32 triangle count
    fwrite(&tris, sizeof(tris), 1, stl);
}

void save_stl(float* verts, unsigned count, const char* filename)
{
    FILE* stl = fopen(filename, "wb");

    const uint32_t tris = count / 9;
    write_header(stl, tris);

    for (uint32_t t=0; t < tris; ++t)
    {
        // Write the face normal (which we'll keep empty)
        for (unsigned j=0; j < 3*sizeof(float); ++j)    fputc(0, stl);

        // Write out all of the vertices.
        fwrite(&verts[t*9], sizeof(float), 9, stl);

        // Attribute byte count
        fputc(0, stl);
        fputc(0, stl);
    }

    fclose(stl);
}

void save_stl_indexed(const float* verts, const uint32_t* indices,
                      uint32_t tri_count, const char* filename)
{
    FILE* stl = fopen(filename, "wb");

    write_header(stl, tri_count);

    for (uint32_t t=0; t < tri_count; ++t)
    {
        // Write the face normal (which we'll keep empty)
        for (unsigned j=0; j < 3*sizeof(float); ++j)    fputc(0, stl);

        // Expand the triangle's three vertices from the vertex table.
        for (int v=0; v < 3; ++v)
            fwrite(&verts[indices[t*3 + v] * 3], sizeof(float), 3, stl);

        // Attribute byte count
        fputc(0, stl);
        fputc(0, stl);
    }

    fclose(stl);
}

#include <cstdio>

#include "fab/formats/dxf.h"

/*
 *  DXF is a sequence of (group code, value) line pairs.  This writes
 *  the minimal R12 skeleton: a HEADER with version + extents, then an
 *  ENTITIES section of closed POLYLINEs (66/1 = "vertices follow",
 *  70/1 = closed).  No TABLES section is needed for layer 0.
 */

bool save_dxf(const std::vector<ContourPath>& paths,
              float xmin, float ymin, float xmax, float ymax,
              const char* filename)
{
    FILE* f = fopen(filename, "wb");
    if (!f)
        return false;

    fprintf(f,
        "0\nSECTION\n"
        "2\nHEADER\n"
        "9\n$ACADVER\n1\nAC1009\n"
        "9\n$EXTMIN\n10\n%g\n20\n%g\n30\n0\n"
        "9\n$EXTMAX\n10\n%g\n20\n%g\n30\n0\n"
        "0\nENDSEC\n",
        xmin, ymin, xmax, ymax);

    fprintf(f,
        "0\nSECTION\n"
        "2\nENTITIES\n");

    for (const auto& path : paths)
    {
        fprintf(f,
            "0\nPOLYLINE\n"
            "8\n0\n"        // layer 0
            "66\n1\n"       // vertices follow
            "70\n1\n");     // closed polyline

        for (const auto& pt : path)
            fprintf(f,
                "0\nVERTEX\n"
                "8\n0\n"
                "10\n%.6g\n"
                "20\n%.6g\n",
                pt[0], pt[1]);

        fprintf(f, "0\nSEQEND\n");
    }

    fprintf(f,
        "0\nENDSEC\n"
        "0\nEOF\n");

    return fclose(f) == 0;
}

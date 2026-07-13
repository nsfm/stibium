#include <cstdio>

#include "fab/formats/svg.h"

bool save_svg(const std::vector<ContourPath>& paths,
              float xmin, float ymin, float xmax, float ymax,
              const char* filename)
{
    FILE* f = fopen(filename, "wb");
    if (!f)
        return false;

    const float w = xmax - xmin;
    const float h = ymax - ymin;

    fprintf(f,
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        "<svg xmlns=\"http://www.w3.org/2000/svg\"\n"
        "     width=\"%gmm\" height=\"%gmm\"\n"
        "     viewBox=\"%g %g %g %g\">\n",
        w, h, xmin, ymin, w, h);

    fprintf(f, " <path fill=\"black\" fill-rule=\"evenodd\" d=\"\n");
    for (const auto& path : paths)
    {
        for (size_t i = 0; i < path.size(); ++i)
        {
            // Flip y: SVG's y axis points down
            fprintf(f, "%c%.6g %.6g\n", i == 0 ? 'M' : 'L',
                    path[i][0], (ymax + ymin) - path[i][1]);
        }
        fprintf(f, "Z\n");
    }
    fprintf(f, " \"/>\n</svg>\n");

    return fclose(f) == 0;
}

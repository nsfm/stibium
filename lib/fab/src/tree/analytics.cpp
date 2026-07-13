#include <cmath>
#include <cstring>
#include <vector>

#include <algorithm>

#include "fab/tree/analytics.h"
#include "fab/tree/parallel_eval.h"
#include "fab/tree/tree.h"

bool analyze_field(MathTree* tree,
                   float xmin, float ymin, float zmin,
                   float xmax, float ymax, float zmax,
                   float resolution, bool flat,
                   int threads, volatile int* halt,
                   FieldStats* out)
{
    memset(out, 0, sizeof(*out));

    const double w = xmax - xmin, h = ymax - ymin;
    const double d = flat ? 0 : zmax - zmin;

    double cell = 0;
    if (resolution > 0)
    {
        cell = 1.0 / resolution;
    }
    else
    {
        // Aim for ~4M samples over the box
        const double content = flat ? w * h : w * h * d;
        cell = pow(content / 4e6, flat ? 1.0 / 2 : 1.0 / 3);
        if (cell <= 0)
            return false;
    }

    const uint64_t nx = uint64_t(fmax(1, ceil(w / cell)));
    const uint64_t ny = uint64_t(fmax(1, ceil(h / cell)));
    const uint64_t nz = flat ? 1 : uint64_t(fmax(1, ceil(d / cell)));
    out->cell = cell;

    // Sample plane-by-plane to keep the coordinate buffers modest
    const uint64_t plane = nx * ny;
    std::vector<float> xs(plane), ys(plane), zs(plane), vals(plane);
    for (uint64_t j = 0; j < ny; ++j)
        for (uint64_t i = 0; i < nx; ++i)
        {
            xs[j * nx + i] = xmin + (i + 0.5) * cell;
            ys[j * nx + i] = ymin + (j + 0.5) * cell;
        }

    double sx = 0, sy = 0, sz = 0;
    float tight[6] = {INFINITY, INFINITY, INFINITY,
                      -INFINITY, -INFINITY, -INFINITY};

    for (uint64_t k = 0; k < nz && !*halt; ++k)
    {
        const float z = flat ? 0 : zmin + (k + 0.5) * cell;
        std::fill(zs.begin(), zs.end(), z);
        parallel_eval(tree, xs.data(), ys.data(), zs.data(), vals.data(),
                      plane, threads, halt);

        for (uint64_t p = 0; p < plane; ++p)
        {
            if (vals[p] < 0)
            {
                out->inside++;
                sx += xs[p];
                sy += ys[p];
                sz += z;
                tight[0] = fmin(tight[0], xs[p]);
                tight[1] = fmin(tight[1], ys[p]);
                tight[2] = fmin(tight[2], z);
                tight[3] = fmax(tight[3], xs[p]);
                tight[4] = fmax(tight[4], ys[p]);
                tight[5] = fmax(tight[5], z);
            }
        }
        out->samples += plane;
    }
    if (*halt)
        return false;

    const double cell_measure = flat ? cell * cell : cell * cell * cell;
    out->volume = double(out->inside) * cell_measure;
    if (out->inside)
    {
        out->com[0] = sx / out->inside;
        out->com[1] = sy / out->inside;
        out->com[2] = sz / out->inside;
        const float half = cell / 2;
        for (int i = 0; i < 3; ++i)
        {
            out->tight[i] = tight[i] - half;
            out->tight[i + 3] = tight[i + 3] + half;
        }
        if (flat)
            out->tight[2] = out->tight[5] = 0;
    }
    return true;
}

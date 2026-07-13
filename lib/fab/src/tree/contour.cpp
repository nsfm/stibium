#include <cmath>
#include <cstring>
#include <unordered_map>

#include "fab/tree/contour.h"

#include "fab/tree/tree.h"
#include "fab/tree/eval.h"
#include "fab/util/region.h"
#include "fab/util/switches.h"

/*
 *  Marching squares over one plane of the field, with exact
 *  edge-indexed chaining (segments meet at shared grid-edge crossings
 *  identified by edge id, not by comparing floats) and optional sharp
 *  corner recovery.
 *
 *  Filled space is f < 0, matching the 3D mesher.  Segments are
 *  emitted with the filled region on the left of the direction of
 *  travel, so outer boundaries wind CCW and holes wind CW.
 */

namespace {

/*
 *  Batched field evaluation: xs/ys/zs are parallel coordinate arrays.
 *  Results land in out.  Points are fed through eval_r in chunks of
 *  MIN_VOLUME (the per-node result buffers' capacity).
 */
void eval_points(MathTree* tree,
                 std::vector<float>& xs, std::vector<float>& ys,
                 std::vector<float>& zs, std::vector<float>& out,
                 volatile int* halt)
{
    out.resize(xs.size());
    size_t done = 0;
    while (done < xs.size() && !*halt)
    {
        const size_t count = std::min(size_t(MIN_VOLUME) - 1,
                                      xs.size() - done);
        Region dummy = {};
        dummy.X = xs.data() + done;
        dummy.Y = ys.data() + done;
        dummy.Z = zs.data() + done;
        dummy.voxels = count;

        const float* result = eval_r(tree, dummy);
        memcpy(out.data() + done, result, count * sizeof(float));
        done += count;
    }
}

struct Segment {
    uint32_t p0, p1;        // crossing-point indices
    uint64_t e0, e1;        // grid-edge ids of the crossings
    uint32_t ci, cj;        // owning cell (for the corner sanity box)
    int32_t corner;         // recovered corner point index, or -1
};

}  // namespace

void contour_field(MathTree* tree,
                   float xmin, float ymin, float xmax, float ymax,
                   uint32_t nx, uint32_t ny, float z,
                   bool detect_features, volatile int* halt,
                   std::vector<ContourPath>& paths,
                   std::atomic<uint64_t>* progress)
{
    paths.clear();
    if (nx == 0 || ny == 0)
        return;

    const float dx = (xmax - xmin) / nx;
    const float dy = (ymax - ymin) / ny;

    // Sample grid, padded with one ring of empty space so shapes
    // clipped by the bounds still produce closed loops.
    const uint32_t NX = nx + 3;
    const uint32_t NY = ny + 3;
    const float EMPTY = 1e9f;

    std::vector<float> gx(NX), gy(NY);
    for (uint32_t i = 0; i < NX; ++i)
        gx[i] = xmin + (int(i) - 1) * dx;
    for (uint32_t j = 0; j < NY; ++j)
        gy[j] = ymin + (int(j) - 1) * dy;

    std::vector<float> values(size_t(NX) * NY, EMPTY);

    {   // Evaluate the interior samples, row by row
        std::vector<float> xs(nx + 1), ys(nx + 1), zs(nx + 1), row;
        for (uint32_t i = 0; i <= nx; ++i)
            xs[i] = gx[i + 1];
        std::fill(zs.begin(), zs.end(), z);

        for (uint32_t j = 1; j + 1 < NY && !*halt; ++j)
        {
            std::fill(ys.begin(), ys.end(), gy[j]);
            eval_points(tree, xs, ys, zs, row, halt);
            memcpy(&values[size_t(j) * NX + 1], row.data(),
                   (nx + 1) * sizeof(float));
            if (progress)
                progress->fetch_add(1, std::memory_order_relaxed);
        }
    }
    if (*halt)
        return;

    const auto val = [&](uint32_t i, uint32_t j) -> float {
        return values[size_t(j) * NX + i];
    };
    const auto inside = [&](uint32_t i, uint32_t j) -> bool {
        return val(i, j) < 0;
    };

    // Grid-edge ids: horizontal edges join (i,j)-(i+1,j), vertical
    // edges join (i,j)-(i,j+1).
    const auto edge_h = [&](uint32_t i, uint32_t j) -> uint64_t {
        return (uint64_t(j) * NX + i) * 2;
    };
    const auto edge_v = [&](uint32_t i, uint32_t j) -> uint64_t {
        return (uint64_t(j) * NX + i) * 2 + 1;
    };

    // Crossing points, deduplicated by grid edge (neighboring cells
    // compute bit-identical positions, but the map makes sharing
    // explicit and gives each crossing a stable index).
    std::vector<std::array<float, 2>> points;
    std::unordered_map<uint64_t, uint32_t> point_by_edge;

    const auto crossing = [&](uint64_t eid,
                              uint32_t ai, uint32_t aj,
                              uint32_t bi, uint32_t bj) -> uint32_t {
        auto found = point_by_edge.find(eid);
        if (found != point_by_edge.end())
            return found->second;

        const float va = val(ai, aj);
        const float vb = val(bi, bj);
        const float t = va / (va - vb);
        const float x = gx[ai] + t * (gx[bi] - gx[ai]);
        const float y = gy[aj] + t * (gy[bj] - gy[aj]);

        const uint32_t id = points.size();
        points.push_back({{x, y}});
        point_by_edge.emplace(eid, id);
        return id;
    };

    // Saddle cells (cases 5 and 10) are resolved by the sign of the
    // field at the cell center; collect them first and batch-evaluate.
    std::vector<std::pair<uint32_t, uint32_t>> saddles;
    for (uint32_t j = 0; j + 1 < NY; ++j)
        for (uint32_t i = 0; i + 1 < NX; ++i)
        {
            const int c = (inside(i,   j)   ? 1 : 0) |
                          (inside(i+1, j)   ? 2 : 0) |
                          (inside(i+1, j+1) ? 4 : 0) |
                          (inside(i,   j+1) ? 8 : 0);
            if (c == 5 || c == 10)
                saddles.push_back({i, j});
        }

    std::unordered_map<uint64_t, bool> saddle_inside;
    if (!saddles.empty())
    {
        std::vector<float> xs, ys, zs, out;
        xs.reserve(saddles.size());
        for (auto& s : saddles)
        {
            xs.push_back((gx[s.first] + gx[s.first + 1]) / 2);
            ys.push_back((gy[s.second] + gy[s.second + 1]) / 2);
            zs.push_back(z);
        }
        eval_points(tree, xs, ys, zs, out, halt);
        if (*halt)
            return;
        for (size_t k = 0; k < saddles.size(); ++k)
            saddle_inside[uint64_t(saddles[k].second) * NX +
                          saddles[k].first] = out[k] < 0;
    }

    // March the cells, emitting directed segments (filled on the left)
    std::vector<Segment> segments;
    std::unordered_map<uint64_t, uint32_t> seg_by_start;

    for (uint32_t j = 0; j + 1 < NY; ++j)
    {
        for (uint32_t i = 0; i + 1 < NX; ++i)
        {
            const int c = (inside(i,   j)   ? 1 : 0) |
                          (inside(i+1, j)   ? 2 : 0) |
                          (inside(i+1, j+1) ? 4 : 0) |
                          (inside(i,   j+1) ? 8 : 0);
            if (c == 0 || c == 15)
                continue;

            const uint64_t eB = edge_h(i, j);
            const uint64_t eT = edge_h(i, j + 1);
            const uint64_t eL = edge_v(i, j);
            const uint64_t eR = edge_v(i + 1, j);

            // Each entry is a (from-edge, to-edge) directed pair
            uint64_t segs[2][2];
            int n = 0;
            switch (c)
            {
                case 1:  segs[n][0] = eB; segs[n][1] = eL; n++; break;
                case 2:  segs[n][0] = eR; segs[n][1] = eB; n++; break;
                case 3:  segs[n][0] = eR; segs[n][1] = eL; n++; break;
                case 4:  segs[n][0] = eT; segs[n][1] = eR; n++; break;
                case 6:  segs[n][0] = eT; segs[n][1] = eB; n++; break;
                case 7:  segs[n][0] = eT; segs[n][1] = eL; n++; break;
                case 8:  segs[n][0] = eL; segs[n][1] = eT; n++; break;
                case 9:  segs[n][0] = eB; segs[n][1] = eT; n++; break;
                case 11: segs[n][0] = eR; segs[n][1] = eT; n++; break;
                case 12: segs[n][0] = eL; segs[n][1] = eR; n++; break;
                case 13: segs[n][0] = eB; segs[n][1] = eR; n++; break;
                case 14: segs[n][0] = eL; segs[n][1] = eB; n++; break;

                case 5:
                    if (saddle_inside[uint64_t(j) * NX + i])
                    {
                        segs[n][0] = eB; segs[n][1] = eR; n++;
                        segs[n][0] = eT; segs[n][1] = eL; n++;
                    }
                    else
                    {
                        segs[n][0] = eB; segs[n][1] = eL; n++;
                        segs[n][0] = eT; segs[n][1] = eR; n++;
                    }
                    break;

                case 10:
                    if (saddle_inside[uint64_t(j) * NX + i])
                    {
                        segs[n][0] = eL; segs[n][1] = eB; n++;
                        segs[n][0] = eR; segs[n][1] = eT; n++;
                    }
                    else
                    {
                        segs[n][0] = eR; segs[n][1] = eB; n++;
                        segs[n][0] = eL; segs[n][1] = eT; n++;
                    }
                    break;
            }

            const auto point_of = [&](uint64_t e) -> uint32_t {
                if (e == eB) return crossing(e, i, j,     i+1, j);
                if (e == eT) return crossing(e, i, j+1,   i+1, j+1);
                if (e == eL) return crossing(e, i, j,     i,   j+1);
                else         return crossing(e, i+1, j,   i+1, j+1);
            };

            for (int s = 0; s < n; ++s)
            {
                Segment seg;
                seg.e0 = segs[s][0];
                seg.e1 = segs[s][1];
                seg.p0 = point_of(seg.e0);
                seg.p1 = point_of(seg.e1);
                seg.ci = i;
                seg.cj = j;
                seg.corner = -1;
                seg_by_start.emplace(seg.e0, segments.size());
                segments.push_back(seg);
            }
        }
    }

    // Sharp corner recovery: where the field normals at a segment's
    // two crossings diverge, split the chord at the intersection of
    // the two tangent lines (2D analog of the mesher's Kobbelt pass).
    if (detect_features && !segments.empty())
    {
        // Batch-evaluate central-difference gradients at every
        // crossing point (four evaluations per point).
        const float eps = 0.01f * fmin(dx, dy);
        std::vector<float> xs, ys, zs, out;
        xs.reserve(points.size() * 4);
        for (const auto& p : points)
        {
            xs.push_back(p[0] + eps); ys.push_back(p[1]);
            xs.push_back(p[0] - eps); ys.push_back(p[1]);
            xs.push_back(p[0]);       ys.push_back(p[1] + eps);
            xs.push_back(p[0]);       ys.push_back(p[1] - eps);
        }
        zs.assign(xs.size(), z);
        eval_points(tree, xs, ys, zs, out, halt);
        if (*halt)
            return;

        std::vector<std::array<float, 2>> normals(points.size());
        for (size_t k = 0; k < points.size(); ++k)
        {
            const float gxv = out[k*4] - out[k*4 + 1];
            const float gyv = out[k*4 + 2] - out[k*4 + 3];
            const float len = std::sqrt(gxv*gxv + gyv*gyv);
            normals[k] = len > 0
                ? std::array<float, 2>{{gxv / len, gyv / len}}
                : std::array<float, 2>{{0, 0}};
        }

        for (auto& seg : segments)
        {
            const auto& n0 = normals[seg.p0];
            const auto& n1 = normals[seg.p1];
            if ((n0[0] == 0 && n0[1] == 0) || (n1[0] == 0 && n1[1] == 0))
                continue;

            // Same threshold as the 3D mesher: cone angle > ~25 deg
            if (n0[0]*n1[0] + n0[1]*n1[1] > 0.9f)
                continue;

            // Intersect the tangent lines through each crossing
            const auto& p0 = points[seg.p0];
            const auto& p1 = points[seg.p1];
            const float t0x = -n0[1], t0y = n0[0];
            const float t1x = -n1[1], t1y = n1[0];
            const float det = t1x*t0y - t0x*t1y;
            if (std::fabs(det) < 1e-12f)
                continue;
            const float s = ((p1[0] - p0[0]) * (-t1y) +
                             (p1[1] - p0[1]) * t1x) / det;
            const float qx = p0[0] + s * t0x;
            const float qy = p0[1] + s * t0y;

            // Only accept corners near the owning cell; a far-away
            // intersection means nearly-parallel tangents or noise.
            if (qx < gx[seg.ci] - dx || qx > gx[seg.ci + 1] + dx ||
                qy < gy[seg.cj] - dy || qy > gy[seg.cj + 1] + dy)
                continue;

            seg.corner = int32_t(points.size());
            points.push_back({{qx, qy}});
        }
    }

    // Chain segments into closed loops by following edge ids
    std::vector<char> visited(segments.size(), 0);
    for (size_t start = 0; start < segments.size(); ++start)
    {
        if (visited[start])
            continue;

        ContourPath path;
        size_t s = start;
        while (!visited[s])
        {
            visited[s] = 1;
            const Segment& seg = segments[s];

            if (path.empty() || points[seg.p0] != path.back())
                path.push_back(points[seg.p0]);
            if (seg.corner >= 0)
                path.push_back(points[seg.corner]);

            auto next = seg_by_start.find(seg.e1);
            if (next == seg_by_start.end())
                break;  // shouldn't happen on a padded grid
            s = next->second;
        }

        // Drop the closing duplicate (paths are implicitly closed)
        if (path.size() > 1 && path.front() == path.back())
            path.pop_back();
        if (path.size() >= 3)
            paths.push_back(std::move(path));
    }
}

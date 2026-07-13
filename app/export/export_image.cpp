#include <Python.h>

#include <boost/python.hpp>

#include <cmath>
#include <cstring>
#include <memory>
#include <vector>

#include <QImage>
#include <QVector3D>

#include "app/settings.h"

#include "export/export_image.h"

#include "graph/graph.h"
#include "graph/node.h"
#include "graph/datum.h"

#include "viewport/render/task.h"

#include "fab/fab.h"
#include "fab/types/shape.h"
#include "fab/tree/render.h"
#include "fab/tree/render_mt.h"
#include "fab/util/region.h"

namespace ImageExport
{

QString render(Graph* graph, const Options& opt)
{

    // Union the file's Shape outputs, 3D and 2D separately: when 3D
    // geometry exists, 2D shapes are construction profiles and are
    // left out of the thumbnail rather than flattening it.
    const std::string want = opt.node_name.toStdString();
    bool node_found = false;
    std::unique_ptr<Shape> u3d, u2d;
    for (auto n : graph->childNodes())
    {
        if (!want.empty() && n->getName() != want)
            continue;
        node_found = true;
        for (auto d : n->childDatums())
        {
            // Match the viewport: only terminal shape outputs are
            // drawn (an output feeding another node is construction
            // geometry, e.g. a cutter feeding a difference).  In
            // single-node mode, that node's shape outputs render
            // even when they feed other nodes - that's the point.
            if (!d->isValid() || d->getType() != fab::ShapeType ||
                !d->isOutput() || !d->currentValue())
                continue;
            if (want.empty() && !d->outgoingLinks().empty())
                continue;
            boost::python::extract<Shape> es(d->currentValue());
            if (!es.check())
                continue;
            const Shape s = es();
            if (std::isinf(s.bounds.xmin) || std::isinf(s.bounds.xmax) ||
                std::isinf(s.bounds.ymin) || std::isinf(s.bounds.ymax))
                continue;
            auto& u = (std::isinf(s.bounds.zmin) ||
                       std::isinf(s.bounds.zmax)) ? u2d : u3d;
            u.reset(u ? new Shape(*u | s) : new Shape(s));
        }
    }

    if (!want.empty() && !node_found)
        return "no node named '" + opt.node_name + "'";

    const bool flat = !u3d;
    const std::unique_ptr<Shape>& total = flat ? u2d : u3d;
    if (!total)
        return "no renderable shapes";
    const Bounds b = total->bounds;

    // 2D content always renders top-down
    QMatrix4x4 M;
    if (!flat)
        M = opt.M;

    Shape src = flat
        ? Shape(total->math, Bounds(b.xmin, b.ymin, 0,
                                    b.xmax, b.ymax, 0.0001f))
        : *total;
    Shape transformed = src.map(RenderTask::getTransform(M));
    const Bounds tb = transformed.bounds;

    float res = opt.resolution;
    if (res <= 0)
    {
        const float extent = fmax(tb.xmax - tb.xmin, tb.ymax - tb.ymin);
        res = extent > 0 ? opt.fit_px / extent : 1;
    }

    // Section view: cut away the near side of the depth range
    const float tzmax = tb.zmax - (1 - opt.section) * (tb.zmax - tb.zmin);

    Region r = {};
    r.ni = uint32_t(fmax(1, (tb.xmax - tb.xmin) * res));
    r.nj = uint32_t(fmax(1, (tb.ymax - tb.ymin) * res));
    r.nk = uint32_t(fmax(1, (tzmax - tb.zmin) * res));
    build_arrays(&r, tb.xmin, tb.ymin, tb.zmin, tb.xmax, tb.ymax, tzmax);

    const int W = r.ni, H = r.nj;
    std::vector<uint16_t> d16(size_t(W) * H, 0);
    std::vector<uint16_t*> d16_rows(H);
    auto s8 = new uint8_t[size_t(W) * H][3];
    memset(s8, 0, size_t(W) * H * 3);
    auto s8_rows = new decltype(s8)[H];
    for (int i = 0; i < H; ++i)
    {
        d16_rows[i] = &d16[size_t(W) * i];
        s8_rows[i] = &s8[size_t(W) * i];
    }

    volatile int halt = 0;
    render16_mt(transformed.tree.get(), r, d16_rows.data(), &halt);
    if (flat)
    {
        // 2D fields have no meaningful z gradient; fill with a
        // straight-on normal (the viewport does the same)
        for (int j = 0; j < H; ++j)
            for (int i = 0; i < W; ++i)
                if (d16_rows[j][i])
                {
                    s8_rows[j][i][0] = 128;
                    s8_rows[j][i][1] = 128;
                    s8_rows[j][i][2] = 255;
                }
    }
    else
    {
        shaded8_mt(transformed.tree.get(), r, d16_rows.data(), s8_rows,
                   &halt);
    }
    free_arrays(&r);

    // Compose an ARGB image: decode the packed normals (n*127+128)
    // and light them (key + ambient, warm material) where depth
    // exists; transparent elsewhere.  Flip vertically (rows are y-up).
    // Key light: honor the user's configured direction (the same
    // one the viewport gizmo drags), falling back to the default
    float Lx = 0.33f, Ly = 0.32f, Lz = 0.89f;
    {
        const auto parts = Settings::get(
                "render/key_light", "0.57,-0.57,0.57").toString().split(',');
        if (parts.size() == 3)
        {
            QVector3D l(parts[0].toFloat(), parts[1].toFloat(),
                        parts[2].toFloat());
            l.normalize();
            // The gizmo's y points down in view coords; image space
            // is y-up here, so flip to match the viewport's lighting
            Lx = l.x();
            Ly = -l.y();
            Lz = l.z();
        }
    }
    const float base[3] = {228, 219, 205};            // warm gray
    QImage img(W, H, QImage::Format_ARGB32);
    img.fill(opt.transparent ? QColor(0, 0, 0, 0) : QColor(28, 26, 24));
    for (int j = 0; j < H; ++j)
        for (int i = 0; i < W; ++i)
            if (d16_rows[j][i])
            {
                const auto& px = s8_rows[j][i];
                const float nx = (px[0] - 128) / 127.f;
                const float ny = (px[1] - 128) / 127.f;
                const float nz = (px[2] - 128) / 127.f;
                const float diff = fmax(0.f, nx*Lx + ny*Ly + nz*Lz);
                const float lum = 0.32f + 0.68f * diff;
                img.setPixel(i, H - 1 - j, qRgba(
                    int(fmin(255.f, base[0] * lum)),
                    int(fmin(255.f, base[1] * lum)),
                    int(fmin(255.f, base[2] * lum)), 255));
            }

    delete [] s8;
    delete [] s8_rows;

    if (!img.save(opt.filename))
        return "could not write " + opt.filename;

    return QString();
}

}  // namespace ImageExport

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

QString collectShapeList(Graph* graph, const QString& node_name,
                         std::vector<Shape>& shapes3d,
                         std::vector<Shape>& shapes2d)
{
    // Gather the file's Shape outputs, 3D and 2D separately: when 3D
    // geometry exists, 2D shapes are construction profiles and are
    // left out of the render rather than flattening it.
    const std::string want = node_name.toStdString();
    bool node_found = false;
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
            auto& list = (std::isinf(s.bounds.zmin) ||
                          std::isinf(s.bounds.zmax)) ? shapes2d
                                                     : shapes3d;
            list.push_back(s);
        }
    }

    if (!want.empty() && !node_found)
        return "no node named '" + node_name + "'";

    return QString();
}

QString collectShapes(Graph* graph, const QString& node_name,
                      std::unique_ptr<Shape>& u3d,
                      std::unique_ptr<Shape>& u2d)
{
    std::vector<Shape> shapes3d, shapes2d;
    const QString err = collectShapeList(graph, node_name,
                                         shapes3d, shapes2d);
    if (!err.isEmpty())
        return err;
    for (const auto& s : shapes3d)
        u3d.reset(u3d ? new Shape(*u3d | s) : new Shape(s));
    for (const auto& s : shapes2d)
        u2d.reset(u2d ? new Shape(*u2d | s) : new Shape(s));
    return QString();
}

QImage renderShapesImage(const std::vector<Shape>& shapes, bool flat,
                         const Options& opt, QString* err)
{
    // 2D content always renders top-down
    QMatrix4x4 M;
    if (!flat)
        M = opt.M;
    const Transform view = RenderTask::getTransform(M);

    // World-space circumsphere for rotation-stable framing
    QVector3D sphere_center;
    float sphere_r = 0;
    if (opt.fit_sphere)
    {
        Bounds wb(INFINITY, INFINITY, INFINITY,
                  -INFINITY, -INFINITY, -INFINITY);
        for (const auto& s : shapes)
        {
            if (s.bounds.xmin > s.bounds.xmax)
                continue;
            wb = Bounds(fmin(wb.xmin, s.bounds.xmin),
                        fmin(wb.ymin, s.bounds.ymin),
                        fmin(wb.zmin, flat ? 0 : s.bounds.zmin),
                        fmax(wb.xmax, s.bounds.xmax),
                        fmax(wb.ymax, s.bounds.ymax),
                        fmax(wb.zmax, flat ? 0 : s.bounds.zmax));
        }
        sphere_center = QVector3D((wb.xmin + wb.xmax) / 2,
                                  (wb.ymin + wb.ymax) / 2,
                                  (wb.zmin + wb.zmax) / 2);
        sphere_r = 0.5f * QVector3D(wb.xmax - wb.xmin,
                                    wb.ymax - wb.ymin,
                                    wb.zmax - wb.zmin).length();
    }

    // Transform each shape into view space, carrying its color, and
    // accumulate the view-space bounding box
    std::vector<Shape> ts;
    Bounds tb(INFINITY, INFINITY, INFINITY,
              -INFINITY, -INFINITY, -INFINITY);
    for (const auto& s : shapes)
    {
        if (s.bounds.xmin > s.bounds.xmax ||
            s.bounds.ymin > s.bounds.ymax)
            continue;       // empty (e.g. a diff layer with no volume)

        const Shape src = flat
            ? Shape(s.math, Bounds(s.bounds.xmin, s.bounds.ymin, 0,
                                   s.bounds.xmax, s.bounds.ymax,
                                   0.0001f))
            : s;
        Shape t = src.map(view);
        t.r = s.r;  t.g = s.g;  t.b = s.b;    // map() drops color

        tb = Bounds(fmin(tb.xmin, t.bounds.xmin),
                    fmin(tb.ymin, t.bounds.ymin),
                    fmin(tb.zmin, t.bounds.zmin),
                    fmax(tb.xmax, t.bounds.xmax),
                    fmax(tb.ymax, t.bounds.ymax),
                    fmax(tb.zmax, t.bounds.zmax));
        ts.push_back(t);
    }
    if (ts.empty())
    {
        *err = "no renderable shapes";
        return QImage();
    }

    // Rotation-stable framing: override the fitted xy box with the
    // circumsphere's square around the view-space center (depth
    // keeps the per-view fit so the z budget stays tight)
    if (opt.fit_sphere && sphere_r > 0)
    {
        const QVector3D c = M.map(sphere_center);
        tb.xmin = c.x() - sphere_r;
        tb.xmax = c.x() + sphere_r;
        tb.ymin = c.y() - sphere_r;
        tb.ymax = c.y() + sphere_r;
    }

    const int ss = opt.supersample > 1 ? opt.supersample : 1;
    float res = opt.resolution;
    if (res <= 0)
    {
        const float extent = fmax(tb.xmax - tb.xmin, tb.ymax - tb.ymin);
        res = extent > 0 ? opt.fit_px / extent : 1;
    }
    res *= ss;

    // Section view: cut away the near side of the depth range
    const float tzmax = tb.zmax - (1 - opt.section) * (tb.zmax - tb.zmin);

    Region r = {};
    r.ni = uint32_t(fmax(1, (tb.xmax - tb.xmin) * res));
    r.nj = uint32_t(fmax(1, (tb.ymax - tb.ymin) * res));
    r.nk = uint32_t(fmax(1, (tzmax - tb.zmin) * res));
    build_arrays(&r, tb.xmin, tb.ymin, tb.zmin, tb.xmax, tb.ymax, tzmax);

    const int W = r.ni, H = r.nj;
    const size_t px_count = size_t(W) * H;

    // Master buffers: nearest depth wins, carrying its normal + color
    std::vector<uint16_t> d16(px_count, 0);
    std::vector<uint8_t> n8(px_count * 3, 0);
    std::vector<uint8_t> c8(px_count * 3, 0);

    // Scratch buffers, reused per shape
    std::vector<uint16_t> sd16(px_count);
    std::vector<uint16_t*> sd16_rows(H);
    auto ss8 = new uint8_t[px_count][3];
    auto ss8_rows = new decltype(ss8)[H];
    for (int i = 0; i < H; ++i)
    {
        sd16_rows[i] = &sd16[size_t(W) * i];
        ss8_rows[i] = &ss8[size_t(W) * i];
    }

    const uint8_t warm[3] = {228, 219, 205};        // default material
    volatile int halt = 0;
    for (const auto& t : ts)
    {
        memset(sd16.data(), 0, px_count * sizeof(uint16_t));
        memset(ss8, 0, px_count * 3);

        render16_mt(t.tree.get(), r, sd16_rows.data(), &halt);
        if (flat)
        {
            // 2D fields have no meaningful z gradient; fill with a
            // straight-on normal (the viewport does the same)
            for (size_t p = 0; p < px_count; ++p)
                if (sd16[p])
                {
                    ss8[p][0] = 128;
                    ss8[p][1] = 128;
                    ss8[p][2] = 255;
                }
        }
        else
        {
            shaded8_mt(t.tree.get(), r, sd16_rows.data(), ss8_rows,
                       &halt);
        }

        const uint8_t rgb[3] = {
            uint8_t(t.r >= 0 ? t.r : warm[0]),
            uint8_t(t.g >= 0 ? t.g : warm[1]),
            uint8_t(t.b >= 0 ? t.b : warm[2])};
        for (size_t p = 0; p < px_count; ++p)
        {
            if (sd16[p] > d16[p])
            {
                d16[p] = sd16[p];
                memcpy(&n8[p * 3], ss8[p], 3);
                memcpy(&c8[p * 3], rgb, 3);
            }
        }
    }
    free_arrays(&r);
    delete [] ss8;
    delete [] ss8_rows;

    // Compose an ARGB image: decode the packed normals (n*127+128)
    // and light them (key + ambient) with each pixel's material
    // color; transparent elsewhere.  Flip vertically (rows are y-up).
    // Key light: honor the user's configured direction (the same
    // one the viewport gizmo drags), falling back to the default
    float Lx = 0.33f, Ly = 0.32f, Lz = 0.89f;
    if (opt.light_override)
    {
        Lx = opt.light[0];
        Ly = opt.light[1];
        Lz = opt.light[2];
    }
    else
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
    QImage img(W, H, QImage::Format_ARGB32);
    img.fill(opt.transparent ? QColor(0, 0, 0, 0) : QColor(28, 26, 24));
    for (int j = 0; j < H; ++j)
        for (int i = 0; i < W; ++i)
        {
            const size_t p = size_t(W) * j + i;
            if (!d16[p])
                continue;
            const float nx = (n8[p*3 + 0] - 128) / 127.f;
            const float ny = (n8[p*3 + 1] - 128) / 127.f;
            const float nz = (n8[p*3 + 2] - 128) / 127.f;
            const float diff = fmax(0.f, nx*Lx + ny*Ly + nz*Lz);
            const float lum = 0.32f + 0.68f * diff;
            img.setPixel(i, H - 1 - j, qRgba(
                int(fmin(255.f, c8[p*3 + 0] * lum)),
                int(fmin(255.f, c8[p*3 + 1] * lum)),
                int(fmin(255.f, c8[p*3 + 2] * lum)), 255));
        }

    // Antialiasing: the render above happened at ss x the target
    // resolution; a smooth downscale averages the subpixels
    if (ss > 1)
        img = img.scaled(W / ss, H / ss, Qt::IgnoreAspectRatio,
                         Qt::SmoothTransformation);

    return img;
}

QString renderShapes(const std::vector<Shape>& shapes, bool flat,
                     const Options& opt)
{
    QString err;
    const QImage img = renderShapesImage(shapes, flat, opt, &err);
    if (!err.isEmpty())
        return err;

    // Fall back to PNG when the filename has no recognized suffix
    // (freedesktop thumbnailers pass bare temp paths)
    if (!img.save(opt.filename) && !img.save(opt.filename, "PNG"))
        return "could not write " + opt.filename;

    return QString();
}

QString render(Graph* graph, const Options& opt)
{
    std::vector<Shape> shapes3d, shapes2d;
    const QString err = collectShapeList(graph, opt.node_name,
                                         shapes3d, shapes2d);
    if (!err.isEmpty())
        return err;

    const bool flat = shapes3d.empty();
    return renderShapes(flat ? shapes2d : shapes3d, flat, opt);
}

}  // namespace ImageExport

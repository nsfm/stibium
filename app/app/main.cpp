#include <Python.h>

#include <QDebug>

#include <QCommandLineParser>
#include <QStandardPaths>
#include <QMainWindow>
#include <QCoreApplication>
#include <QSurfaceFormat>
#include <QStringList>
#include <QMessageBox>

#include "app/app.h"
#include "graph/hooks/hooks.h"

#include "fab/fab.h"
#include "graph/graph.h"

#include "graph/node.h"
#include "graph/script_node.h"
#include "graph/datum.h"
#include "graph/proxy/graph.h"
#include "graph/proxy/node.h"
#include "export/export_worker.h"
#include "graph/serialize/serializer.h"

#include <QFile>
#include <QJsonDocument>
#include "viewport/render/task.h"
#include "fab/fab.h"
#include "fab/types/shape.h"
#include "fab/tree/render.h"
#include "fab/util/region.h"

#include <QImage>
#include <QMatrix4x4>

#include <cstring>
#include <memory>

/*
 *  Implements --render: composites every shape in the file (what the
 *  viewport would show, unioned) into a shaded PNG.  resolution <= 0
 *  fits the longest side to `size` pixels.  Pure-C render path; no GL.
 */
static int renderHeadless(App& app, const QString& out, float resolution,
                          const QString& view, int size, float section)
{
    // Union the file's Shape outputs, 3D and 2D separately: when 3D
    // geometry exists, 2D shapes are construction profiles and are
    // left out of the thumbnail rather than flattening it.
    std::unique_ptr<Shape> u3d, u2d;
    for (auto n : app.getGraph()->childNodes())
        for (auto d : n->childDatums())
        {
            // Match the viewport: only terminal shape outputs are
            // drawn (an output feeding another node is construction
            // geometry, e.g. a cutter feeding a difference)
            if (!d->isValid() || d->getType() != fab::ShapeType ||
                !d->isOutput() || !d->outgoingLinks().empty() ||
                !d->currentValue())
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

    const bool flat = !u3d;
    const std::unique_ptr<Shape>& total = flat ? u2d : u3d;
    if (!total)
    {
        fprintf(stderr, "render: no renderable shapes in this file\n");
        return 1;
    }
    const Bounds b = total->bounds;

    // View transform (2D content always renders top-down)
    QMatrix4x4 M;
    if (!flat && view != "top")
    {
        if (view == "front")
        {
            M.rotate(-90, 1, 0, 0);
        }
        else  // iso
        {
            M.rotate(-55, 1, 0, 0);
            M.rotate(25, 0, 0, 1);
        }
    }

    Shape src = flat
        ? Shape(total->math, Bounds(b.xmin, b.ymin, 0,
                                    b.xmax, b.ymax, 0.0001f))
        : *total;
    Shape transformed = src.map(RenderTask::getTransform(M));
    const Bounds tb = transformed.bounds;

    float res = resolution;
    if (res <= 0)
    {
        const float extent = fmax(tb.xmax - tb.xmin, tb.ymax - tb.ymin);
        res = extent > 0 ? size / extent : 1;
    }

    // Section view: cut away the near side of the depth range
    const float tzmax = tb.zmax - (1 - section) * (tb.zmax - tb.zmin);

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
    render16(transformed.tree.get(), r, d16_rows.data(), &halt, nullptr);
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
        shaded8(transformed.tree.get(), r, d16_rows.data(), s8_rows,
                &halt, nullptr);
    }
    free_arrays(&r);

    // Compose an ARGB image: decode the packed normals (n*127+128)
    // and light them (key + ambient, warm material) where depth
    // exists; transparent elsewhere.  Flip vertically (rows are y-up).
    const float Lx = 0.33f, Ly = 0.32f, Lz = 0.89f;   // key light
    const float base[3] = {228, 219, 205};            // warm gray
    QImage img(W, H, QImage::Format_ARGB32);
    img.fill(Qt::transparent);
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

    if (!img.save(out))
    {
        fprintf(stderr, "render: could not write %s\n",
                out.toLocal8Bit().constData());
        return 1;
    }
    return 0;
}

/*
 *  Implements --validate and --export: reports script/datum errors,
 *  then (for export) runs the file's export worker synchronously.
 *  Returns the process exit code.
 */
static int runHeadless(App& app, bool validate_only,
                       const QString& export_file, float resolution,
                       int detect_features, bool skip_export)
{
    int errors = 0;
    for (auto n : app.getGraph()->childNodes())
    {
        if (auto sn = dynamic_cast<ScriptNode*>(n))
        {
            const auto err = sn->getError();
            if (!err.empty())
            {
                fprintf(stderr, "[%s] script error (line %i):\n%s\n",
                        n->getName().c_str(), sn->getErrorLine(),
                        err.c_str());
                errors++;
            }
        }
        for (auto d : n->childDatums())
        {
            if (!d->isValid())
            {
                fprintf(stderr, "[%s.%s] error: %s\n",
                        n->getName().c_str(), d->getName().c_str(),
                        d->getError().c_str());
                errors++;
            }
        }
    }

    if (validate_only)
    {
        printf("%s\n", errors ? "invalid" : "valid");
        return errors ? 1 : 0;
    }

    if (errors)
    {
        fprintf(stderr, "aborted: %i error(s) in the file\n", errors);
        return 1;
    }

    if (skip_export)
        return 0;

    // Find the export worker registered by the file's node scripts
    QList<ExportWorker*> workers;
    for (auto np : app.getProxy()->nodeProxies())
        if (auto w = np->getExportWorker())
            workers << w;

    if (workers.isEmpty())
    {
        fprintf(stderr, "export: no export node in this file\n");
        return 1;
    }
    if (workers.size() > 1)
        fprintf(stderr, "warning: %lld export nodes; using the first\n",
                (long long)workers.size());

    return workers.first()->runHeadless(export_file, resolution,
                                        detect_features) ? 0 : 1;
}

/*
 *  Implements --resave: writes the loaded graph back out in the
 *  current protocol.  Load + resave is also the canonicalizing
 *  round-trip: a second resave of the output is byte-identical.
 */
static int resaveHeadless(App& app, const QString& out)
{
    QFile file(out);
    if (!file.open(QIODevice::WriteOnly))
    {
        fprintf(stderr, "resave: could not write %s\n",
                out.toLocal8Bit().constData());
        return 1;
    }
    auto i = app.getProxy()->canvasInfo();
    file.write(QJsonDocument(
                SceneSerializer::run(app.getGraph(), &i)).toJson());
    return 0;
}

int main(int argc, char *argv[])
{
    {   // Set the default OpenGL version to be 2.1 with sample buffers
        QSurfaceFormat format;
        format.setVersion(2, 1);
        QSurfaceFormat::setDefaultFormat(format);
    }

    // Headless verbs run without a display: pick the offscreen
    // platform before Qt initializes (unless the user chose one)
    for (int i = 1; i < argc; ++i)
    {
        const QByteArray arg(argv[i]);
        if (arg == "--export" || arg == "--validate" ||
            arg == "--render" || arg == "--resave")
        {
            if (!qEnvironmentVariableIsSet("QT_QPA_PLATFORM"))
                qputenv("QT_QPA_PLATFORM", "offscreen");
            break;
        }
    }

    // Create the Application object
    App app(argc, argv);

    // Initialize various Python modules and the interpreter itself
    fab::preInit();
    Graph::preInit();
    AppHooks::preInit();
    Py_Initialize();

    // Set locale to C to make atof correctly parse floats
    setlocale(LC_NUMERIC, "C");

    {   // Modify Python's default search path to include the application's
        // directory (as this doesn't happen on Linux by default)
#if defined Q_OS_MAC
        QStringList path = QCoreApplication::applicationDirPath().split("/");
        path.removeLast();
        path << "Resources";
        fab::postInit({path.join("/").toStdString()});
#elif defined Q_OS_LINUX
        auto dir = QCoreApplication::applicationDirPath();
        std::vector<std::string> fab_paths =
            {(dir + "/sb").toStdString(),
             (dir + "/../share/stibium/").toStdString(),
             (dir + "/../share/antimony/").toStdString()};
        for (auto p : QStandardPaths::standardLocations(
                QStandardPaths::AppDataLocation))
        {
            fab_paths.push_back(p.toStdString());
        }
        fab::postInit(fab_paths);
#elif defined Q_OS_OPENBSD
        auto dir = QCoreApplication::applicationDirPath();
        std::vector<std::string> fab_paths =
            {(dir + "/sb").toStdString(),
             (dir + "/../share/stibium/").toStdString(),
             (dir + "/../share/antimony/").toStdString()};
        for (auto p : QStandardPaths::standardLocations(
                QStandardPaths::AppDataLocation))
        {
            fab_paths.push_back(p.toStdString());
        }
        fab::postInit(fab_paths);	
#elif defined Q_OS_WIN32
        auto dir = QCoreApplication::applicationDirPath();
        fab::postInit({(dir + "/sb").toStdString()});
#else
#error "Unknown OS!"
#endif
    }

    {   // Install operator.or_ as a reducer for shapes
        auto op = PyImport_ImportModule("operator");
        Datum::installReducer(fab::ShapeType, PyObject_GetAttrString(op, "or_"));
        Py_DECREF(op);
    }

    {   // Check to make sure that the fab module exists
        PyObject* fab = PyImport_ImportModule("fab");
        if (!fab)
        {
            PyErr_Print();
            QMessageBox::critical(NULL, "Import error",
                    "Import Error:<br><br>"
                    "Could not find <tt>fab</tt> Python module.<br>"
                    "Stibium will now exit.");
            exit(1);
        }
        Py_DECREF(fab);
    }

    {   // Parse command-line arguments
        QCommandLineParser parser;
        parser.setApplicationDescription("CAD from a parallel universe");
        parser.addHelpOption();
        QCommandLineOption forceHeightmap("heightmap",
                "Open 3D windows in heightmap mode");
        parser.addOption(forceHeightmap);
        QCommandLineOption exportOpt("export",
                "Run the file's export node headlessly, writing FILE "
                "(format follows the extension), then exit", "FILE");
        parser.addOption(exportOpt);
        QCommandLineOption resolutionOpt("resolution",
                "Resolution for --export (voxels per unit); overrides "
                "the script's value", "R");
        parser.addOption(resolutionOpt);
        QCommandLineOption detectOpt("detect-features",
                "Force feature detection on for --export");
        parser.addOption(detectOpt);
        QCommandLineOption renderOpt("render",
                "Render every shape in the file to an image (PNG "
                "recommended) and exit; iso view for 3D, top-down "
                "for 2D", "FILE");
        parser.addOption(renderOpt);
        QCommandLineOption viewOpt("view",
                "View for --render: iso (default), top, or front",
                "VIEW", "iso");
        parser.addOption(viewOpt);
        QCommandLineOption sectionOpt("section",
                "Cross-section for --render: fraction of the model's "
                "depth to keep, cutting from the viewer's side "
                "(0-1; default 1 = no cut)", "F", "1");
        parser.addOption(sectionOpt);
        QCommandLineOption sizeOpt("size",
                "Longest image side in pixels for --render when "
                "--resolution isn't given (default 512)", "N", "512");
        parser.addOption(sizeOpt);
        QCommandLineOption resaveOpt("resave",
                "Load the file and save it back out as FILE in the "
                "current protocol (batch migration / canonicalizing), "
                "then exit", "FILE");
        parser.addOption(resaveOpt);
        QCommandLineOption validateOpt("validate",
                "Load the file, report script and datum errors to "
                "stderr, and exit (0 = valid)");
        parser.addOption(validateOpt);
        parser.addPositionalArgument("file", "File to open", "[file]");

        parser.process(app);

        const bool headless = parser.isSet(exportOpt) ||
                              parser.isSet(renderOpt) ||
                              parser.isSet(resaveOpt) ||
                              parser.isSet(validateOpt);

        auto args = parser.positionalArguments();
        if (args.length() > 1)
        {
            qCritical("Too many command-line arguments");
            exit(1);
        }
        else if (headless && args.length() != 1)
        {
            fprintf(stderr,
                    "--export/--render/--validate require a file\n");
            exit(1);
        }
        else if (args.length() == 1)
        {
            app.setHeadless(headless);
            app.loadFile(args[0]);
        }

        if (headless)
        {
            const float res = parser.isSet(resolutionOpt)
                ? parser.value(resolutionOpt).toFloat() : -1;

            int code = runHeadless(app, parser.isSet(validateOpt),
                                   parser.value(exportOpt), res,
                                   parser.isSet(detectOpt) ? 1 : -1,
                                   !parser.isSet(exportOpt));
            if (code == 0 && parser.isSet(renderOpt))
                code = renderHeadless(app, parser.value(renderOpt), res,
                                      parser.value(viewOpt),
                                      parser.value(sizeOpt).toInt(),
                                      parser.value(sectionOpt).toFloat());
            if (code == 0 && parser.isSet(resaveOpt))
                code = resaveHeadless(app, parser.value(resaveOpt));
            // Return (rather than exit()) so Qt tears down cleanly
            return code;
        }
    }

    app.makeDefaultWindows();
    return app.exec();
}

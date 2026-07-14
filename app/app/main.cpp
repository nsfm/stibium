#include <Python.h>

#include <QDebug>

#include <QCommandLineParser>
#include <QStandardPaths>
#include <QMainWindow>
#include <QCoreApplication>
#include <QElapsedTimer>
#include <QProgressDialog>
#include <QSocketNotifier>
#include <QSurfaceFormat>
#include <QStringList>
#include <QMessageBox>

#include "app/app.h"
#include "app/vm_import.h"
#include "graph/hooks/hooks.h"

#include "fab/fab.h"
#include "graph/graph.h"

#include "graph/node.h"
#include "graph/script_node.h"
#include "graph/datum.h"
#include "graph/proxy/graph.h"
#include "graph/proxy/node.h"
#include "export/export_worker.h"
#include "export/export_image.h"
#include "export/gif_writer.h"
#include "graph/serialize/serializer.h"

#include <QFile>
#include <QJsonDocument>
#include "viewport/render/task.h"
#include "fab/fab.h"
#include "fab/types/shape.h"
#include "fab/tree/analytics.h"
#include "fab/tree/render.h"
#include "fab/tree/render_mt.h"
#include "fab/util/region.h"

#include <QDirIterator>
#include <QImage>
#include <QJsonArray>
#include <QJsonObject>
#include <QMatrix4x4>

#include <cmath>
#include <cstring>
#include <memory>

/*
 *  Implements --describe-nodes: dumps the node library as JSON
 *  (name, title, category, inputs with types and defaults, outputs).
 *  Static text parse of the .node scripts - nothing is evaluated.
 */
static int describeNodes(App& app)
{
    static const QList<QRegularExpression> title_res = {
        QRegularExpression("title\\('+([^']+)'+\\)"),
        QRegularExpression("title\\(\"+([^\"]+)\"+\\)")};
    static const QRegularExpression input_re(
        "^\\s*input\\('([^']+)'\\s*,\\s*([A-Za-z_][A-Za-z0-9_.]*)"
        "(?:\\s*,\\s*(.*))?\\)\\s*$");
    static const QRegularExpression output_re(
        "^\\s*output\\('([^']+)'");

    QJsonArray out;
    for (auto path : app.nodePaths())
    {
        QDirIterator itr(path, QDirIterator::Subdirectories);
        while (itr.hasNext())
        {
            const auto fname = itr.next();
            if (!fname.endsWith(".node"))
                continue;
            QFile file(fname);
            if (!file.open(QIODevice::ReadOnly))
                continue;
            const QString txt = QTextStream(&file).readAll();

            auto split = fname.split('/');
            while (!split.isEmpty() && split.first() != "nodes")
                split.removeFirst();
            if (!split.isEmpty())
                split.removeFirst();
            const QString base = split.isEmpty()
                ? fname : QString(split.takeLast()).replace(".node", "");
            const QString category = split.join("/");

            QString title = base;
            for (const auto& re : title_res)
            {
                const auto match = re.match(txt);
                if (match.hasMatch())
                    title = match.captured(1);
            }

            QJsonArray inputs, outputs;
            for (const auto& line : txt.split('\n'))
            {
                if (auto match = input_re.match(line); match.hasMatch())
                {
                    QJsonObject in;
                    in["name"] = match.captured(1);
                    in["type"] = match.captured(2);
                    auto def = match.captured(3).trimmed();
                    if (!def.isEmpty())
                        in["default"] = def;
                    inputs.append(in);
                }
                else if (auto o = output_re.match(line); o.hasMatch())
                {
                    outputs.append(o.captured(1));
                }
            }

            QJsonObject node;
            node["name"] = base;
            node["title"] = title;
            node["category"] = category;
            node["inputs"] = inputs;
            node["outputs"] = outputs;
            out.append(node);
        }
    }

    printf("%s", QJsonDocument(out).toJson().constData());
    return 0;
}

/*
 *  Implements --analyze: grid-integrates the model's shapes and
 *  prints volume, center of mass, and tight bounds as JSON.
 */
static int analyzeHeadless(App& app, float resolution,
                           const QString& node)
{
    std::unique_ptr<Shape> u3d, u2d;
    const QString err = ImageExport::collectShapes(
            app.getGraph(), node, u3d, u2d);
    if (!err.isEmpty())
    {
        fprintf(stderr, "analyze: %s\n", err.toLocal8Bit().constData());
        return 1;
    }

    const bool flat = !u3d;
    const Shape& s = flat ? *u2d : *u3d;
    const Bounds b = s.bounds;

    volatile int halt = 0;
    FieldStats stats;
    if (!analyze_field(s.tree.get(),
                       b.xmin, b.ymin, flat ? 0 : b.zmin,
                       b.xmax, b.ymax, flat ? 0 : b.zmax,
                       resolution, flat, -1, &halt, &stats))
    {
        fprintf(stderr, "analyze: field integration failed\n");
        return 1;
    }

    QJsonObject out;
    out[flat ? "area" : "volume"] = stats.volume;
    out["center_of_mass"] = QJsonArray(
            {stats.com[0], stats.com[1], stats.com[2]});
    out["bounds_declared"] = QJsonArray(
            {b.xmin, b.ymin, b.zmin, b.xmax, b.ymax, b.zmax});
    out["bounds_tight"] = QJsonArray(
            {stats.tight[0], stats.tight[1], stats.tight[2],
             stats.tight[3], stats.tight[4], stats.tight[5]});
    out["flat"] = flat;
    out["samples"] = double(stats.samples);
    out["cell"] = stats.cell;
    printf("%s", QJsonDocument(out).toJson().constData());
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

/*
 *  Implements --render: canned views over the shared image exporter.
 */
/*
 *  Implements --diff: loads a second file, splits the two models
 *  into unchanged / removed / added fields (set algebra on the
 *  distance fields - meshes can't do this, math can), prints their
 *  measures as JSON, and optionally renders a composite image
 *  (unchanged in gray, removed in red, added in green).
 */
static int diffHeadless(App& app, const QString& other_file,
                        const QString& render_file, float resolution,
                        const QString& view, int size, float section,
                        int aa)
{
    std::unique_ptr<Shape> a3, a2;
    QString err = ImageExport::collectShapes(app.getGraph(), "", a3, a2);
    if (err.isEmpty() && !a3 && !a2)
        err = "no renderable shapes";
    if (!err.isEmpty())
    {
        fprintf(stderr, "diff: base file: %s\n",
                err.toLocal8Bit().constData());
        return 1;
    }

    app.loadFile(other_file);
    std::unique_ptr<Shape> b3, b2;
    err = ImageExport::collectShapes(app.getGraph(), "", b3, b2);
    if (err.isEmpty() && !b3 && !b2)
        err = "no renderable shapes";
    if (!err.isEmpty())
    {
        fprintf(stderr, "diff: '%s': %s\n",
                other_file.toLocal8Bit().constData(),
                err.toLocal8Bit().constData());
        return 1;
    }

    const bool a_flat = !a3, b_flat = !b3;
    if (a_flat != b_flat)
    {
        fprintf(stderr, "diff: one model is 2D and the other 3D\n");
        return 1;
    }
    const bool flat = a_flat;
    const Shape& A = flat ? *a2 : *a3;
    const Shape& B = flat ? *b2 : *b3;

    Shape removed = A & ~B;
    Shape added = B & ~A;
    Shape unchanged = A & B;
    removed.r = 220;  removed.g = 60;   removed.b = 50;
    added.r = 70;     added.g = 190;    added.b = 90;

    // Integrate the changed regions.  The exact-identity shortcut
    // skips integration when the fields are literally the same math.
    const bool identical = A.math == B.math;
    double vol_removed = 0, vol_added = 0;
    if (!identical)
    {
        volatile int halt = 0;
        for (auto* layer : {&removed, &added})
        {
            const Bounds& lb = layer->bounds;
            if (lb.xmin > lb.xmax || lb.ymin > lb.ymax)
                continue;               // provably empty
            FieldStats stats;
            if (analyze_field(layer->tree.get(),
                              lb.xmin, lb.ymin, flat ? 0 : lb.zmin,
                              lb.xmax, lb.ymax, flat ? 0 : lb.zmax,
                              resolution, flat, -1, &halt, &stats))
            {
                (layer == &removed ? vol_removed : vol_added)
                        = stats.volume;
            }
        }
    }

    QJsonObject out;
    out["identical"] = identical ||
            (vol_removed == 0 && vol_added == 0);
    out[flat ? "removed_area" : "removed_volume"] = vol_removed;
    out[flat ? "added_area" : "added_volume"] = vol_added;
    out["flat"] = flat;
    printf("%s", QJsonDocument(out).toJson().constData());

    if (!render_file.isEmpty())
    {
        ImageExport::Options opt;
        if (view == "front")
        {
            opt.M.rotate(-90, 1, 0, 0);
        }
        else if (view != "top")  // iso
        {
            opt.M.rotate(-55, 1, 0, 0);
            opt.M.rotate(25, 0, 0, 1);
        }
        opt.section = section;
        opt.resolution = resolution;
        opt.fit_px = size;
        opt.supersample = aa;
        opt.filename = render_file;

        err = ImageExport::renderShapes({unchanged, removed, added},
                                        flat, opt);
        if (!err.isEmpty())
        {
            fprintf(stderr, "diff render: %s\n",
                    err.toLocal8Bit().constData());
            return 1;
        }
    }
    return 0;
}

/*
 *  Implements --turntable / --wiggle / --stereo: the model rendered
 *  from a sweep of yaw angles with rotation-stable framing (the
 *  circumsphere), assembled into an animated GIF or a side-by-side
 *  stereo pair.
 */
static int animationHeadless(App& app, const QString& turntable,
                             const QString& wiggle,
                             const QString& stereo,
                             const QString& lightsweep, int frames,
                             float resolution, int size, int aa,
                             const QString& view)
{
    std::vector<Shape> s3, s2;
    QString err = ImageExport::collectShapeList(app.getGraph(), "",
                                                s3, s2);
    if (err.isEmpty() && s3.empty())
        err = s2.empty() ? "no renderable shapes"
                         : "2D models have no depth to animate";
    if (!err.isEmpty())
    {
        fprintf(stderr, "animate: %s\n", err.toLocal8Bit().constData());
        return 1;
    }

    const auto frame = [&](float yaw, float light_az, QString* e) {
        ImageExport::Options opt;
        if (view == "front")
            opt.M.rotate(-90, 1, 0, 0);
        else if (view != "top")     // iso pitch
            opt.M.rotate(-55, 1, 0, 0);
        opt.M.rotate(yaw, 0, 0, 1);
        opt.resolution = resolution;
        opt.fit_px = size;
        opt.supersample = aa;
        opt.fit_sphere = true;      // stable framing across the sweep
        opt.transparent = false;    // GIF alpha is binary; use theme
        if (light_az == light_az)   // NaN = keep configured light
        {
            // Key light circles the model at ~35 degrees elevation
            const float el = 0.6f;
            opt.light_override = true;
            opt.light[0] = cosf(light_az) * cosf(el);
            opt.light[1] = sinf(light_az) * cosf(el);
            opt.light[2] = sinf(el);
        }
        return ImageExport::renderShapesImage(s3, false, opt, e);
    };
    const float KEEP_LIGHT = NAN;

    if (!stereo.isEmpty())
    {
        // Parallel-view pair (cross-eyed viewers: swap panels)
        QImage L = frame(2.5f, KEEP_LIGHT, &err);
        QImage R = err.isEmpty() ? frame(-2.5f, KEEP_LIGHT, &err)
                                 : QImage();
        if (!err.isEmpty())
        {
            fprintf(stderr, "stereo: %s\n",
                    err.toLocal8Bit().constData());
            return 1;
        }
        const int gap = 8;
        QImage pair(L.width() + R.width() + gap, L.height(),
                    QImage::Format_ARGB32);
        pair.fill(QColor(28, 26, 24));
        for (int j = 0; j < L.height(); ++j)
        {
            memcpy(pair.scanLine(j), L.constScanLine(j),
                   size_t(L.width()) * 4);
            memcpy(pair.scanLine(j) + (L.width() + gap) * 4,
                   R.constScanLine(j), size_t(R.width()) * 4);
        }
        if (!pair.save(stereo) && !pair.save(stereo, "PNG"))
        {
            fprintf(stderr, "stereo: could not write '%s'\n",
                    stereo.toLocal8Bit().constData());
            return 1;
        }
    }

    if (!wiggle.isEmpty())
    {
        QImage a = frame(-3, KEEP_LIGHT, &err);
        QImage b = err.isEmpty() ? frame(3, KEEP_LIGHT, &err) : QImage();
        if (err.isEmpty() &&
            !save_gif({a, b}, {12, 12}, wiggle))
            err = "could not write '" + wiggle + "'";
        if (!err.isEmpty())
        {
            fprintf(stderr, "wiggle: %s\n",
                    err.toLocal8Bit().constData());
            return 1;
        }
    }

    if (!lightsweep.isEmpty())
    {
        const int n = frames < 4 ? 36 : frames;
        const int delay = std::max(2, 400 / n);
        QList<QImage> seq;
        for (int i = 0; i < n; ++i)
        {
            seq.append(frame(0, 6.2831853f * i / n, &err));
            if (!err.isEmpty())
            {
                fprintf(stderr, "lightsweep: %s\n",
                        err.toLocal8Bit().constData());
                return 1;
            }
            fprintf(stderr, "\rlightsweep: frame %d/%d", i + 1, n);
        }
        fprintf(stderr, "\n");
        if (!save_gif(seq, {delay}, lightsweep))
        {
            fprintf(stderr, "lightsweep: could not write '%s'\n",
                    lightsweep.toLocal8Bit().constData());
            return 1;
        }
    }

    if (!turntable.isEmpty())
    {
        const int n = frames < 4 ? 36 : frames;
        const int delay = std::max(2, 400 / n);   // ~4s per revolution
        QList<QImage> seq;
        for (int i = 0; i < n; ++i)
        {
            seq.append(frame(360.f * i / n, KEEP_LIGHT, &err));
            if (!err.isEmpty())
            {
                fprintf(stderr, "turntable: %s\n",
                        err.toLocal8Bit().constData());
                return 1;
            }
            fprintf(stderr, "\rturntable: frame %d/%d", i + 1, n);
        }
        fprintf(stderr, "\n");
        if (!save_gif(seq, {delay}, turntable))
        {
            fprintf(stderr, "turntable: could not write '%s'\n",
                    turntable.toLocal8Bit().constData());
            return 1;
        }
    }
    return 0;
}

static int renderHeadless(App& app, const QString& out, float resolution,
                          const QString& view, int size, float section,
                          const QString& node, int aa)
{
    ImageExport::Options opt;
    opt.node_name = node;
    if (view == "front")
    {
        opt.M.rotate(-90, 1, 0, 0);
    }
    else if (view != "top")  // iso
    {
        opt.M.rotate(-55, 1, 0, 0);
        opt.M.rotate(25, 0, 0, 1);
    }
    opt.section = section;
    opt.resolution = resolution;
    opt.fit_px = size;
    opt.supersample = aa;
    opt.filename = out;

    const QString err = ImageExport::render(app.getGraph(), opt);
    if (!err.isEmpty())
    {
        fprintf(stderr, "render: %s\n", err.toLocal8Bit().constData());
        return 1;
    }
    return 0;
}

/*  Progress pump for long operations inside script evaluation (mesh
 *  import sampling).  Evaluation blocks the GUI thread, so this
 *  repaints a modal progress dialog by hand, with user input
 *  excluded — the graph can't safely change mid-eval. */
static void guiLongOpHook(const char* label, uint64_t done, uint64_t total)
{
    static QProgressDialog* dialog = NULL;
    static QElapsedTimer elapsed;

    if (done >= total)
    {
        if (dialog)
        {
            dialog->close();
            dialog->deleteLater();
            dialog = NULL;
            QCoreApplication::processEvents(
                    QEventLoop::ExcludeUserInputEvents);
        }
        elapsed.invalidate();
        return;
    }

    if (!elapsed.isValid())
    {
        elapsed.start();
        return;             // quick operations never show a dialog
    }
    if (!dialog)
    {
        if (elapsed.elapsed() < 350)
            return;

        // Deterministic show (QProgressDialog's minimumDuration
        // estimation is unreliable when we're also the event pump).
        // No cancel button: sampling has no safe abort point today.
        dialog = new QProgressDialog(
                QString::fromUtf8(label), QString(), 0, 1000);
        dialog->setWindowModality(Qt::ApplicationModal);
        dialog->setMinimumDuration(0);
        dialog->setWindowTitle("Stibium");
        dialog->show();
    }
    dialog->setValue(int(1000.0 * done / (total ? total : 1)));
    QCoreApplication::processEvents(QEventLoop::ExcludeUserInputEvents);
}

#ifndef Q_OS_WIN
#include <csignal>
#include <sys/socket.h>
#include <unistd.h>

/*
 *  Ctrl-C from the launching terminal should run the normal quit
 *  path (save prompt included), not kill the process.  Signal
 *  handlers can't touch Qt, so the handler writes one byte to a
 *  socketpair and a QSocketNotifier picks it up on the event loop.
 */
static int sigint_fds[2] = {-1, -1};

static void sigintHandler(int)
{
    const char byte = 1;
    const auto r = ::write(sigint_fds[0], &byte, 1);
    (void)r;
}
#endif

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
            arg == "--render" || arg == "--resave" ||
            arg == "--analyze" || arg == "--describe-nodes" ||
            arg == "--import-vm" || arg == "--diff" ||
            arg == "--turntable" || arg == "--wiggle" ||
            arg == "--stereo" || arg == "--lightsweep")
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
    // Pass 0 so CPython does NOT install its own signal handlers:
    // we want SIGINT to run the app's quit path (save prompt), and
    // reactive script evaluation would otherwise let Python reclaim
    // the handler after ours is set.
    Py_InitializeEx(0);

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
        QCommandLineOption describeOpt("describe-nodes",
                "Dump the node library (names, categories, inputs, "
                "outputs) as JSON and exit");
        parser.addOption(describeOpt);
        QCommandLineOption nodeOpt("node",
                "With --render: render only the named node's shape "
                "outputs (visual bisection)", "NAME");
        parser.addOption(nodeOpt);
        QCommandLineOption analyzeOpt("analyze",
                "Integrate the model's shapes and print volume, "
                "center of mass, and tight bounds as JSON, then exit");
        parser.addOption(analyzeOpt);
        QCommandLineOption resaveOpt("resave",
                "Load the file and save it back out as FILE in the "
                "current protocol (batch migration / canonicalizing), "
                "then exit", "FILE");
        parser.addOption(resaveOpt);
        QCommandLineOption aaOpt("aa",
                "Antialiasing for rendered images: supersample factor "
                "(1 disables; default 2)", "N", "2");
        parser.addOption(aaOpt);
        QCommandLineOption turntableOpt("turntable",
                "Render a looping turntable GIF of the model to FILE "
                "(--frames, --size, --view and --aa apply)", "FILE");
        parser.addOption(turntableOpt);
        QCommandLineOption wiggleOpt("wiggle",
                "Render a two-frame wigglegram GIF to FILE (depth "
                "through motion)", "FILE");
        parser.addOption(wiggleOpt);
        QCommandLineOption stereoOpt("stereo",
                "Render a side-by-side parallel-view stereo pair to "
                "FILE", "FILE");
        parser.addOption(stereoOpt);
        QCommandLineOption lightsweepOpt("lightsweep",
                "Render a GIF where the key light circles the still "
                "model (shadow study) to FILE", "FILE");
        parser.addOption(lightsweepOpt);
        QCommandLineOption framesOpt("frames",
                "Frame count for --turntable (default 36)", "N", "36");
        parser.addOption(framesOpt);
        QCommandLineOption diffOpt("diff",
                "Compare the loaded model against FILE: print the "
                "removed/added material measures as JSON, and with "
                "--render, write a composite image (gray unchanged, "
                "red removed, green added)", "FILE");
        parser.addOption(diffOpt);
        QCommandLineOption importVmOpt("import-vm",
                "Translate a Fidget .vm math tape (given as the "
                "positional argument) into FILE as a Stibium project, "
                "then exit", "FILE");
        parser.addOption(importVmOpt);
        QCommandLineOption validateOpt("validate",
                "Load the file, report script and datum errors to "
                "stderr, and exit (0 = valid)");
        parser.addOption(validateOpt);
        parser.addPositionalArgument("file", "File to open", "[file]");

        parser.process(app);

        if (parser.isSet(describeOpt))
            return describeNodes(app);

        const bool headless = parser.isSet(exportOpt) ||
                              parser.isSet(renderOpt) ||
                              parser.isSet(resaveOpt) ||
                              parser.isSet(analyzeOpt) ||
                              parser.isSet(validateOpt) ||
                              parser.isSet(importVmOpt) ||
                              parser.isSet(diffOpt) ||
                              parser.isSet(turntableOpt) ||
                              parser.isSet(wiggleOpt) ||
                              parser.isSet(stereoOpt) ||
                              parser.isSet(lightsweepOpt);

        auto args = parser.positionalArguments();
        // GUI sessions get a progress dialog for long script-eval
        // operations; register before any file loads so imports in
        // command-line-opened files report too
        if (!headless)
            fab::longOpHook = guiLongOpHook;

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
        // .vm import never loads the positional as a project: it IS
        // the tape to translate
        if (parser.isSet(importVmOpt))
        {
            if (args.length() != 1)
            {
                fprintf(stderr, "--import-vm needs the input .vm as "
                                "the file argument\n");
                return 1;
            }
            return importVmHeadless(args[0],
                                    parser.value(importVmOpt));
        }

        if (args.length() == 1)
        {
            app.setHeadless(headless);
            app.loadFile(args[0]);
        }

        if (headless)
        {
            const float res = parser.isSet(resolutionOpt)
                ? parser.value(resolutionOpt).toFloat() : -1;

            const int aa = std::max(1, parser.value(aaOpt).toInt());

            if (parser.isSet(diffOpt))
                return diffHeadless(app, parser.value(diffOpt),
                                    parser.value(renderOpt), res,
                                    parser.value(viewOpt),
                                    parser.value(sizeOpt).toInt(),
                                    parser.value(sectionOpt).toFloat(),
                                    aa);

            if (parser.isSet(turntableOpt) || parser.isSet(wiggleOpt) ||
                parser.isSet(stereoOpt) || parser.isSet(lightsweepOpt))
                return animationHeadless(app, parser.value(turntableOpt),
                                         parser.value(wiggleOpt),
                                         parser.value(stereoOpt),
                                         parser.value(lightsweepOpt),
                                         parser.value(framesOpt).toInt(),
                                         res,
                                         parser.value(sizeOpt).toInt(),
                                         aa, parser.value(viewOpt));

            int code = runHeadless(app, parser.isSet(validateOpt),
                                   parser.value(exportOpt), res,
                                   parser.isSet(detectOpt) ? 1 : -1,
                                   !parser.isSet(exportOpt));
            if (code == 0 && parser.isSet(renderOpt))
                code = renderHeadless(app, parser.value(renderOpt), res,
                                      parser.value(viewOpt),
                                      parser.value(sizeOpt).toInt(),
                                      parser.value(sectionOpt).toFloat(),
                                      parser.value(nodeOpt), aa);
            if (code == 0 && parser.isSet(analyzeOpt))
                code = analyzeHeadless(app, res, parser.value(nodeOpt));
            if (code == 0 && parser.isSet(resaveOpt))
                code = resaveHeadless(app, parser.value(resaveOpt));
            // Return (rather than exit()) so Qt tears down cleanly
            return code;
        }
    }

#ifndef Q_OS_WIN
    if (::socketpair(AF_UNIX, SOCK_STREAM, 0, sigint_fds) == 0)
    {
        auto notifier = new QSocketNotifier(sigint_fds[1],
                                            QSocketNotifier::Read, &app);
        QObject::connect(notifier, &QSocketNotifier::activated, &app,
                [&app](QSocketDescriptor, QSocketNotifier::Type){
                    char byte;
                    const auto r = ::read(sigint_fds[1], &byte, 1);
                    (void)r;
                    app.onQuit();
                });
        struct sigaction sa = {};
        sa.sa_handler = sigintHandler;
        ::sigaction(SIGINT, &sa, NULL);
    }
#endif

    app.makeDefaultWindows();
    return app.exec();
}

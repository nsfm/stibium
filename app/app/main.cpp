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

/*
 *  Implements --validate and --export: reports script/datum errors,
 *  then (for export) runs the file's export worker synchronously.
 *  Returns the process exit code.
 */
static int runHeadless(App& app, bool validate_only,
                       const QString& export_file, float resolution,
                       int detect_features)
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
        fprintf(stderr, "export aborted: %i error(s) in the file\n",
                errors);
        return 1;
    }

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
        if (arg == "--export" || arg == "--validate")
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
        QCommandLineOption validateOpt("validate",
                "Load the file, report script and datum errors to "
                "stderr, and exit (0 = valid)");
        parser.addOption(validateOpt);
        parser.addPositionalArgument("file", "File to open", "[file]");

        parser.process(app);

        const bool headless = parser.isSet(exportOpt) ||
                              parser.isSet(validateOpt);

        auto args = parser.positionalArguments();
        if (args.length() > 1)
        {
            qCritical("Too many command-line arguments");
            exit(1);
        }
        else if (headless && args.length() != 1)
        {
            fprintf(stderr, "--export/--validate require a file\n");
            exit(1);
        }
        else if (args.length() == 1)
        {
            app.setHeadless(headless);
            app.loadFile(args[0]);
        }

        if (headless)
            exit(runHeadless(app, parser.isSet(validateOpt),
                             parser.value(exportOpt),
                             parser.isSet(resolutionOpt)
                                 ? parser.value(resolutionOpt).toFloat()
                                 : -1,
                             parser.isSet(detectOpt) ? 1 : -1));
    }

    app.makeDefaultWindows();
    return app.exec();
}

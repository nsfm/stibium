#include <Python.h>

#include <QMessageBox>
#include <QStandardPaths>
#include <QDir>
#include <QFileDialog>
#include <QFileOpenEvent>
#include <QJsonDocument>

#include "app/app.h"

#include "graph/proxy/graph.h"
#include "graph/serialize/serializer.h"
#include "graph/serialize/deserializer.h"

#include "undo/undo_stack.h"

#include "graph/graph.h"
#include "app/theme.h"
#include "app/settings.h"
#include <QGraphicsView>
#include <QFileInfo>
#include "canvas/scene.h"
#include "canvas/canvas_view.h"

App::App(int& argc, char** argv)
    : QApplication(argc, argv),
      graph(new Graph()), proxy(new GraphProxy(graph)),
      undo_stack(new UndoStack(this))
{
    connect(undo_stack, &QUndoStack::cleanChanged,
            this, &App::cleanChanged);

    autosave_timer = new QTimer(this);
    connect(autosave_timer, &QTimer::timeout,
        this, &App::onAutosave);

    setApplicationName("Stibium");
    setOrganizationDomain("github.com/nsfm/stibium");

    Settings::init();
    autosave_interval = Settings::get(
            "autosave/interval_ms", autosave_interval).toInt();

    Theme::apply(this);
}

App::~App()
{
    delete graph;
    delete proxy;
    delete undo_stack;
}

App* App::instance()
{
    Q_ASSERT(dynamic_cast<App*>(QApplication::instance()));
    return static_cast<App*>(QApplication::instance());
}

void App::makeDefaultWindows()
{
    newCanvasWindow();
    newViewportWindow();
}

QStringList App::nodePaths() const
{
    QStringList paths;
#if defined Q_OS_MAC
    // On Mac, the 'nodes' folder must be at
    // Antimony.app/Contents/Resources/nodes
    auto path = applicationDirPath().split("/");
    path.removeLast(); // Trim the MacOS folder from the path
    paths << path.join("/") + "/Resources/nodes";
#elif defined Q_OS_LINUX
    // If we're running Antimony from the build folder, use sb/nodes
    paths << applicationDirPath() + "/sb/nodes";
    paths << applicationDirPath() + "/../share/antimony/nodes";
#elif defined Q_OS_OPENBSD
    // If we're running Antimony from the build folder, use sb/nodes
    paths << applicationDirPath() + "/sb/nodes";
    paths << applicationDirPath() + "/../share/antimony/nodes";
#elif defined Q_OS_WIN32
    // Windows only supports running from the build directory
    paths << applicationDirPath() + "/sb/nodes";
#else
#error "Unknown OS!"
#endif

    for (auto p : QStandardPaths::standardLocations(
            QStandardPaths::AppDataLocation))
    {
        paths << p + "/nodes";
    }

    // Filter paths to canonical forms, keeping only the paths that actually
    // exist as real folders
    QSet<QString> existing_paths;
    for (auto p : paths)
    {
        auto d = QDir(p);
        if (d.exists())
        {
            existing_paths.insert(d.canonicalPath());
        }
    }

    return existing_paths.values();
}

////////////////////////////////////////////////////////////////////////////////

void App::onNew()
{
    if (undo_stack->isClean() || QMessageBox::question(
                NULL, "Discard unsaved changes?",
                "Discard unsaved changes?") == QMessageBox::Yes)
    {
        graph->clear();
        filename.clear();
        undo_stack->clear();

        emit(filenameChanged(""));
    }
}

void App::onSave()
{
    if (filename.isEmpty())
        return onSaveAs();

    QFile file(filename);
    file.open(QIODevice::WriteOnly);

    {
        auto i = proxy->canvasInfo();
        file.write(QJsonDocument(SceneSerializer::run(graph, &i)).toJson());
    }

    undo_stack->setClean();
}

void App::onAutosave()
{
    if (filename.isEmpty())
        return;
    return onSave();
}

void App::onToggleAutosave(bool enabled)
{
    if (enabled) {
      autosave_timer->start(autosave_interval);
    } else {
      autosave_timer->stop();
    }
}

void App::onSaveAs()
{
    QString f = QFileDialog::getSaveFileName(NULL, "Save as",
            Settings::get("files/last_dir", "").toString(), "*.sb");
    if (!f.isEmpty())
    {
#ifdef Q_OS_LINUX
        if (!f.endsWith(".sb"))
            f += ".sb";
#elif defined Q_OS_OPENBSD
	if (!f.endsWith(".sb"))
            f += ".sb";
#endif
        if (!QFileInfo(QFileInfo(f).path()).isWritable())
        {
            QMessageBox::critical(NULL, "Save As error",
                    "<b>Save As error:</b><br>"
                    "Target file is not writable.");
            return;
        }
        filename = f;
        Settings::set("files/last_dir", QFileInfo(f).absolutePath());
        emit(filenameChanged(filename));
        return onSave();
    }
}

void App::onOpen()
{
    if (undo_stack->isClean() || QMessageBox::question(
                NULL, "Discard unsaved changes?",
                "Discard unsaved changes?") == QMessageBox::Yes)
    {
        QString f = QFileDialog::getOpenFileName(NULL, "Open",
            Settings::get("files/last_dir", "").toString(), "*.sb");
        if (!f.isEmpty())
            loadFile(f);
    }
}

void App::onQuit()
{
    if (undo_stack->isClean() || QMessageBox::question(
                NULL, "Discard unsaved changes?",
                "Discard unsaved changes?") == QMessageBox::Yes)
    {
        quit();
    }
}

////////////////////////////////////////////////////////////////////////////////

void App::newCanvasWindow()
{
    proxy->newCanvasWindow();
}

void App::newViewportWindow()
{
    proxy->newViewportWindow();
}

void App::newQuadWindow()
{
    proxy->newQuadWindow();
}

////////////////////////////////////////////////////////////////////////////////

void App::onAbout()
{
    QString txt(
            "<i>Stibium</i><br><br>"
            "CAD from a parallel universe.<br>"
            "<a href=\"https://github.com/nsfm/stibium\">https://github.com/nsfm/stibium</a><br><br>"
            "A continuation of <a href=\"https://github.com/mkeeter/antimony\">Antimony</a><br>"
            "© 2013-2015 Matthew Keeter<br>"
            "_________________________________________________<br><br>"
            "Includes code from <a href=\"https://github.com/mkeeter/kokopelli\">kokopelli</a>, which is <br>"
            "© 2012-2013 MIT<br>"
            "© 2013-2014 Matthew Keeter<br><br>"
            "Inspired by the <a href=\"http://kokompe.cba.mit.edu\">fab modules</a><br>"
            "_________________________________________________<br><br>"
    );
    QString tag(GITTAG);
    QString branch(GITBRANCH);
    QString rev(GITREV);

    if (!tag.isEmpty())
        txt += "Release: <tt>" + tag + "</tt>";
    else
        txt += "Branch: <tt>" + branch + "</tt>";
    txt += "<br>Git revision: <tt>" + rev + "</tt>";

    QMessageBox::about(NULL, "Stibium", txt);
}

////////////////////////////////////////////////////////////////////////////////

bool App::event(QEvent *event)
{
    switch (event->type()) {
        case QEvent::FileOpen:
            loadFile(static_cast<QFileOpenEvent*>(event)->file());
            return true;
        default:
            return QApplication::event(event);
    }
}

////////////////////////////////////////////////////////////////////////////////

void App::loadFile(QString f)
{
    Settings::set("files/last_dir", QFileInfo(f).absolutePath());
    filename = f;
    graph->clear();

    // XXX disable rendering

    QFile file(f);
    if (!file.open(QIODevice::ReadOnly))
    {
        QMessageBox::critical(NULL, "Loading error",
                "<b>Loading error:</b><br>"
                "File does not exist.");
        onNew();
        return;
    }

    SceneDeserializer::Info ds;
    const bool success = SceneDeserializer::run(
            QJsonDocument::fromJson(file.readAll()).object(),
            graph, &ds);

    if (!success)
    {
        QMessageBox::critical(NULL, "Loading error",
                "<b>Loading error:</b><br>" +
                ds.error_message);
        onNew();
    } else {
        // If there's a warning message, show it in a box.
        if (!ds.warning_message.isNull())
            QMessageBox::information(NULL, "Loading information",
                    "<b>Loading information:</b><br>" +
                    ds.warning_message);

        proxy->setPositions(ds.frames);
        emit(filenameChanged(filename));

        // Center canvas views on the loaded graph, which may live far
        // from the origin. Deferred one event-loop turn so windows are
        // laid out and inspector positions applied; uses the view's own
        // zoom/center properties (the same machinery as zoom-to-node).
        QTimer::singleShot(0, this, [this]{
            auto scene = proxy->canvasScene();
            if (!scene)
                return;
            const auto r = scene->itemsBoundingRect();
            if (r.isNull())
                return;
            for (auto v : scene->views())
                if (auto c = dynamic_cast<CanvasView*>(v))
                    c->zoomToFit(r);
        });
    }
}

////////////////////////////////////////////////////////////////////////////////

void App::onNewWindow()
{
    emit(filenameChanged(filename));
    emit(cleanChanged(undo_stack->isClean()));
}

////////////////////////////////////////////////////////////////////////////////

void App::undo()
{
    undo_stack->undo();
}

void App::redo()
{
    undo_stack->redo();
}

void App::pushUndoStack(UndoCommand* c)
{
    undo_stack->push(c);
}

void App::beginUndoMacro(QString text)
{
    undo_stack->beginMacro(text);
}

void App::endUndoMacro()
{
    undo_stack->endMacro();
}

QAction* App::getUndoAction()
{
    auto a = undo_stack->createUndoAction(this);
    a->setShortcuts(QKeySequence::Undo);
    return a;
}

QAction* App::getRedoAction()
{
    auto a = undo_stack->createRedoAction(this);
    a->setShortcuts(QKeySequence::Redo);
    return a;
}

bool App::isUndoStackClean() const
{
    return undo_stack->isClean();
}

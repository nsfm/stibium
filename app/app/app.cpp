#include <Python.h>

#include <QFileSystemWatcher>
#include <QFutureWatcher>
#include <QtConcurrent>
#include <QMessageBox>
#include <QPainter>
#include <QSplashScreen>
#include <QStandardPaths>
#include <QDir>
#include <QFileDialog>
#include <QFileOpenEvent>
#include <QJsonDocument>

#include <memory>

#include "app/app.h"
#include "export/export_image.h"
#include "dialog/exporting.h"

#include "graph/proxy/graph.h"
#include "graph/serialize/serializer.h"
#include "graph/serialize/deserializer.h"

#include "undo/undo_stack.h"

#include "graph/graph.h"
#include "app/theme.h"
#include "app/settings.h"
#include "fab/fab.h"
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

////////////////////////////////////////////////////////////////////////////////

/*  The shared splash artwork: wordmark, tagline, and a status line. */
static QPixmap splashPixmap(const QString& subtitle)
{
    const qreal dpr = qApp->devicePixelRatio();
    QPixmap pm(int(460 * dpr), int(150 * dpr));
    pm.setDevicePixelRatio(dpr);
    pm.fill(QColor(28, 26, 24));

    QPainter painter(&pm);
    painter.setRenderHint(QPainter::Antialiasing);
    painter.setRenderHint(QPainter::TextAntialiasing);

    auto font = painter.font();
    font.setPointSize(26);
    font.setBold(true);
    painter.setFont(font);
    painter.setPen(QColor(228, 219, 205));
    painter.drawText(QRect(0, 26, 460, 44), Qt::AlignCenter, "Stibium");

    font.setPointSize(11);
    font.setBold(false);
    font.setItalic(true);
    painter.setFont(font);
    painter.setPen(QColor(150, 142, 130));
    painter.drawText(QRect(0, 74, 460, 24), Qt::AlignCenter,
                     "CAD from a parallel universe");

    font.setItalic(false);
    font.setPointSize(9);
    painter.setFont(font);
    painter.setPen(QColor(120, 113, 104));
    painter.drawText(QRect(0, 104, 460, 22), Qt::AlignCenter, subtitle);

    painter.setPen(QColor(70, 64, 58));
    painter.drawRect(0, 0, 459, 149);
    return pm;
}

void App::showStartupSplash(const QString& subtitle)
{
    if (headless || startup_splash)
        return;
    startup_splash = new QSplashScreen(splashPixmap(subtitle));
    startup_splash->show();
    startup_splash->raise();

    // Pump events long enough for the compositor to actually map AND
    // paint the splash NOW.  A single processEvents flushes the map
    // request but not the expose that paints it; the blocking startup
    // that follows (Python init, graph evaluation) runs no event loop,
    // so the splash would otherwise paint - as a one-frame flash -
    // only when that work finishes.  The maxtime overload waits on the
    // round-trip instead of spinning past it.
    processEvents(QEventLoop::ExcludeUserInputEvents, 150);
    startup_splash->repaint();
}

void App::finishStartupSplash(QWidget* w)
{
    if (!startup_splash)
        return;
    startup_splash->finish(w);
    delete startup_splash;
    startup_splash = nullptr;
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
        fab::setProjectDir("");
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
    ignore_next_file_change = true;
    watchCurrentFile();
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
        fab::setProjectDir(QFileInfo(f).absolutePath().toStdString());
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
            "<a href=\"https://github.com/nsfm/stibium\">https://github.com/nsfm/stibium</a><br>"
            "© 2022-2026 Nate Dube<br>"
            "Licensed under the GNU Affero General Public License v3.0 or later.<br>"
            "_________________________________________________<br><br>"
            "A continuation of <a href=\"https://github.com/mkeeter/antimony\">Antimony</a>, which is<br>"
            "© 2013-2022 Matthew Keeter and other contributors (MIT)<br>"
            "_________________________________________________<br><br>"
            "Includes code from <a href=\"https://github.com/mkeeter/kokopelli\">kokopelli</a>, which is <br>"
            "© 2012-2013 MIT<br>"
            "© 2013 Matthew Keeter<br><br>"
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

void App::watchCurrentFile()
{
    if (filename.isEmpty())
        return;

    if (!file_watcher)
    {
        file_watcher = new QFileSystemWatcher(this);
        connect(file_watcher, &QFileSystemWatcher::fileChanged,
                this, &App::onFileChangedOnDisk);
    }

    const auto watched = file_watcher->files();
    if (!watched.isEmpty())
        file_watcher->removePaths(watched);
    file_watcher->addPath(filename);
}

void App::onFileChangedOnDisk(const QString& path)
{
    if (ignore_next_file_change)
    {
        ignore_next_file_change = false;
        watchCurrentFile();
        return;
    }

    // Editors and tools often replace files (write + rename), which
    // drops the watch and can fire before the new content lands, so
    // wait a beat, re-arm, and only then decide whether to reload.
    QTimer::singleShot(150, this, [this, path]{
        if (path != filename || !QFile::exists(path))
            return;
        watchCurrentFile();

        if (!undo_stack->isClean())
        {
            qWarning("%s changed on disk; not reloading "
                     "(session has unsaved changes)",
                     path.toLocal8Bit().constData());
            return;
        }
        loadFile(path);
    });
}

void App::openFile(QString f)
{
    if (undo_stack->isClean() || QMessageBox::question(
                NULL, "Discard unsaved changes?",
                "Discard unsaved changes?") == QMessageBox::Yes)
    {
        loadFile(f);
    }
}

void App::touchRecentFile(const QString& f)
{
    auto recent = Settings::get("files/recent", QStringList())
            .toStringList();
    recent.removeAll(f);
    recent.prepend(f);
    while (recent.size() > 10)
        recent.removeLast();
    Settings::set("files/recent", recent);
}

bool App::runAnalytics()
{
    std::unique_ptr<Shape> u3d, u2d;
    const QString err = ImageExport::collectShapes(
            graph, QString(), u3d, u2d);
    if (!err.isEmpty())
    {
        QMessageBox::information(NULL, "Analytics",
                "<b>Analytics:</b><br>" + err);
        return false;
    }

    const bool flat = !u3d;
    const Shape s = flat ? *u2d : *u3d;
    const Bounds b = s.bounds;

    auto dialog = new ExportingDialog();
    dialog->setStatus("Analyzing...");

    volatile int halt = 0;
    FieldStats stats;
    auto future = QtConcurrent::run([&]{
        return analyze_field(s.tree.get(),
                b.xmin, b.ymin, flat ? 0 : b.zmin,
                b.xmax, b.ymax, flat ? 0 : b.zmax,
                -1, flat, -1, &halt, &stats);
    });
    QFutureWatcher<bool> watcher;
    watcher.setFuture(future);
    connect(&watcher, &QFutureWatcher<bool>::finished,
            dialog, &QDialog::accept);
    if (dialog->exec() == QDialog::Rejected)
    {
        halt = 1;
        future.waitForFinished();
    }
    delete dialog;

    if (!future.result())
        return false;

    analytics_stats = stats;
    analytics_flat = flat;
    analytics_valid = true;
    return true;
}

void App::loadFile(QString f)
{
    Settings::set("files/last_dir", QFileInfo(f).absolutePath());
    filename = f;

    // Relative mesh-import paths resolve against the project dir,
    // so it must be current before deserialization runs any scripts
    fab::setProjectDir(QFileInfo(f).absolutePath().toStdString());

    // Big graphs take seconds to evaluate on load with no feedback;
    // show a splash card for the duration.  On cold start main()
    // already owns one (startup_splash), so only make a mid-session
    // splash here for File > Open on an already-running window.
    std::unique_ptr<QSplashScreen> splash;
    if (!headless && !startup_splash)
    {
        splash.reset(new QSplashScreen(
                splashPixmap("Loading " + QFileInfo(f).fileName() +
                             "...")));
        splash->show();
        QCoreApplication::processEvents(
                QEventLoop::ExcludeUserInputEvents);
    }

    // Drop undo history before demolishing the graph: stale commands
    // hold pointers into the old nodes, and undoing across a load
    // would replay them into freed memory. (File > New always did
    // this; File > Open and the live-reload path never did.)
    undo_stack->clear();

    graph->clear();

    // XXX disable rendering

    QFile file(f);
    if (!file.open(QIODevice::ReadOnly))
    {
        if (headless)
        {
            fprintf(stderr, "loading error: cannot read %s\n",
                    f.toLocal8Bit().constData());
            exit(1);
        }
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
        if (headless)
        {
            fprintf(stderr, "loading error: %s\n",
                    ds.error_message.toLocal8Bit().constData());
            exit(1);
        }
        QMessageBox::critical(NULL, "Loading error",
                "<b>Loading error:</b><br>" +
                ds.error_message);
        onNew();
    } else {
        // If there's a warning message, show it in a box.
        if (!ds.warning_message.isNull())
        {
            if (headless)
                fprintf(stderr, "loading warning: %s\n",
                        ds.warning_message.toLocal8Bit().constData());
            else
                QMessageBox::information(NULL, "Loading information",
                        "<b>Loading information:</b><br>" +
                        ds.warning_message);
        }

        proxy->setPositions(ds.frames);
        // A freshly loaded file matches disk exactly: mark this the
        // clean baseline so quitting doesn't prompt to save it.
        undo_stack->setClean();
        emit(filenameChanged(filename));
        touchRecentFile(filename);

        if (!headless)
            watchCurrentFile();

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

#pragma once

#include <Python.h>

#include <QApplication>

#include "fab/tree/analytics.h"
#include <QAction>
#include <QTimer>

class Graph;
class GraphProxy;
class QFileSystemWatcher;
class UndoCommand;
class UndoStack;

class App : public QApplication
{
    Q_OBJECT
public:
    explicit App(int& argc, char **argv);
    ~App();

    /*
     *  Constructs a canvas and viewport window for the root graph.
     */
    void makeDefaultWindows();

    /*
     *  Helper function to get running instance.
     */
    static App* instance();

    /*
     *  Returns the path to the node directories
     *  (which varies depending on OS).
     */
    QStringList nodePaths() const;

    /*
     *  Returns the top-level graph proxy object
     */
    GraphProxy* getProxy() const { return proxy; }

    /*
     *  Returns the root graph
     */
    Graph* getGraph() const { return graph; }

    /*
     *  In headless mode, load errors go to stderr (and exit) instead
     *  of popping modal dialogs that nothing can dismiss.
     */
    void setHeadless(bool h) { headless = h; }

    /*
     *  Global Undo and Redo operations
     */
    void undo();
    void redo();

    /*
     *  Pushes an undo command to the stack
     */
    void pushUndoStack(UndoCommand* c);

    /*
     *  Begins a multi-command undo macro with the given description
     */
    void beginUndoMacro(QString text);

    /*
     *  Ends a multi-command undo macro
     */
    void endUndoMacro();

    /*
     *  Get undo and redo actions to populate in menus
     */
    QAction* getUndoAction();
    QAction* getRedoAction();

    /*
     *  Loads a file specified by name
     */
    void loadFile(QString f);

    /*
     *  loadFile with the same discard-unsaved-changes confirmation
     *  File > Open uses (for the recent-files menu)
     */
    void openFile(QString f);

    /*
     *  Moves f to the front of the recent-files list
     */
    static void touchRecentFile(const QString& f);

    /*
     *  Integrates the model's shapes (busy dialog + worker thread),
     *  storing the result for viewport overlays.  Returns false if
     *  there was nothing to analyze.
     */
    bool runAnalytics();

    bool analyticsValid() const { return analytics_valid; }
    bool analyticsFlat() const { return analytics_flat; }
    const FieldStats& analyticsStats() const { return analytics_stats; }

    /*
     *  (Re)arms the file watcher on the current filename
     */
    void watchCurrentFile();

    /*
     *  Live reload: handles the open file changing on disk
     */
    void onFileChangedOnDisk(const QString& path);

    /*
     *  Emits the signals used by windows to set their titles
     *  (to give a new window a correct title bar)
     */
    void onNewWindow();

    /*
     *  Checks if undo stack is clean
     */
    bool isUndoStackClean() const;

public slots:
    /*
     *  Commands from File menu
     */
    void onNew();
    void onSave();
    void onAutosave();
    void onToggleAutosave(bool);
    void onSaveAs();
    void onOpen();
    void onQuit();

    /*
     *  Open new windows for the root graph.
     */
    void newCanvasWindow();
    void newViewportWindow();
    void newQuadWindow();

    /*
     *  Help menu
     */
    void onAbout();

signals:
    /*
     *  Emitted when the file's name changes
     */
    void filenameChanged(QString f);

    /*
     *  Emitted when the unsaved state changes
     */
    void cleanChanged(bool u);

protected:
    /*
     *  Overload event handling to detect drag-and-drop events
     */
    bool event(QEvent* event);

    Graph* graph;
    GraphProxy* proxy;
    UndoStack* undo_stack;

    /*  Latest field analytics (viewport overlays read these)  */
    FieldStats analytics_stats;
    bool analytics_valid=false;
    bool analytics_flat=false;

    /*
     *  Watches the open file and reloads when it changes on disk
     *  (only while the session has no unsaved edits), so external
     *  external tools can edit the model live.
     */
    QFileSystemWatcher* file_watcher=nullptr;
    bool ignore_next_file_change=false;
    bool headless=false;

    QString filename;
    QTimer* autosave_timer;
    int autosave_interval = 60000;
};

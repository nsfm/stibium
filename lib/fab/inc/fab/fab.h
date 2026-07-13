#ifndef FAB_H
#define FAB_H

#include <Python.h>

#include <cstdint>
#include <vector>
#include <string>

namespace fab
{
    struct ParseError {};
    struct ShapeError {};

    /** Callback that raises a Python runtime exception. */
    void onParseError(ParseError const& e);

    /** Callback that raises a Python runtime exception. */
    void onShapeError(ShapeError const& e);

    /** Loads the fab module in Python's namespace.
     *
     *  Must be called before Py_Initialize().
     */
    void preInit();

    /** Loads the fab module in Python's namespace.
     *
     *  fab_paths are added to Python's search path.
     */
    void postInit(std::vector<std::string> fab_paths);

    /** Directory against which relative mesh-import paths resolve
     *  (the directory holding the current .sb file; empty when the
     *  project has never been saved).  The app updates this on
     *  open/save/new, before any scripts run. */
    void setProjectDir(std::string dir);
    std::string projectDir();

    /** Optional hook for long blocking operations inside script
     *  evaluation (today: mesh-import sampling).  Called on the
     *  thread that is blocked, every few tens of ms, with
     *  done < total while running and done == total exactly once at
     *  the end (close your dialog there).  The GUI registers a
     *  progress-dialog pump; headless leaves it NULL. */
    extern void (*longOpHook)(const char* label,
                              uint64_t done, uint64_t total);

    extern PyTypeObject* ShapeType;
}

#endif // FAB_H

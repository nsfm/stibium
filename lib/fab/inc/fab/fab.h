#ifndef FAB_H
#define FAB_H

#include <Python.h>

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

    extern PyTypeObject* ShapeType;
}

#endif // FAB_H

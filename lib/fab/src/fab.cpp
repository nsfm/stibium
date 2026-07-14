#include <boost/python.hpp>
#include <boost/python/raw_function.hpp>

#include "fab/fab.h"
#include "fab/types/shape.h"
#include "fab/types/transform.h"
#include "fab/mesh/mesh_import.h"
#include "fab/tree/analytics.h"
#include "fab/tree/grid.h"

using namespace boost::python;

////////////////////////////////////////////////////////////////////////////////

void (*fab::longOpHook)(const char*, uint64_t, uint64_t) = NULL;

static std::string _project_dir;

void fab::setProjectDir(std::string dir)
{
    _project_dir = dir;
}

std::string fab::projectDir()
{
    return _project_dir;
}

static std::string project_dir_py()
{
    return _project_dir;
}

/*  Grid-integrates a shape's field over its bounds: the measurement
 *  primitive behind the Checks node family and scripted analysis.
 *  Returns a dict; raises on unbounded shapes.  resolution is
 *  samples per unit (<= 0 targets ~4M samples); 2D shapes report
 *  area in 'volume' with z fields zeroed. */
static dict measure_py(const Shape& s, float resolution)
{
    const Bounds& b = s.bounds;
    const bool flat = std::isinf(b.zmin) || std::isinf(b.zmax);
    if (std::isinf(b.xmin) || std::isinf(b.xmax) ||
        std::isinf(b.ymin) || std::isinf(b.ymax))
    {
        PyErr_SetString(PyExc_RuntimeError,
                "measure: shape has unbounded x/y extents");
        throw_error_already_set();
    }

    volatile int halt = 0;
    FieldStats stats;
    if (!analyze_field(s.tree.get(),
                       b.xmin, b.ymin, flat ? 0 : b.zmin,
                       b.xmax, b.ymax, flat ? 0 : b.zmax,
                       resolution, flat, -1, &halt, &stats))
    {
        PyErr_SetString(PyExc_RuntimeError,
                "measure: field integration failed");
        throw_error_already_set();
    }

    dict out;
    out[flat ? "area" : "volume"] = stats.volume;
    out["com"] = boost::python::make_tuple(
            stats.com[0], stats.com[1], stats.com[2]);
    out["tight"] = boost::python::make_tuple(
            stats.tight[0], stats.tight[1], stats.tight[2],
            stats.tight[3], stats.tight[4], stats.tight[5]);
    out["inside"] = (unsigned long long)stats.inside;
    out["samples"] = (unsigned long long)stats.samples;
    out["cell"] = stats.cell;
    out["flat"] = flat;
    return out;
}

/*  Backend for fab.shapes.import_mesh (see py/fab/shapes.py, which
 *  owns path resolution and user-facing checks).
 *  Returns (Shape, stibium_stamp, sha256, from_cache). */
static tuple import_mesh_py(std::string path, float voxels_per_unit,
                            std::string cache_dir)
{
    /*  Grids abandoned by superseded imports (e.g. an earlier
     *  resolution of the same file) are only referenced by the
     *  registry now; drop them before allocating another. */
    grid_registry_trim();

    const auto res = fab_mesh::import_mesh_grid(
            path, voxels_per_unit, cache_dir);
    if (!res.grid_id)
    {
        PyErr_SetString(PyExc_RuntimeError,
                ("import_mesh: " + res.error).c_str());
        throw_error_already_set();
    }

    Shape s("g" + std::to_string(res.grid_id),
            Bounds(res.bounds[0], res.bounds[1], res.bounds[2],
                   res.bounds[3], res.bounds[4], res.bounds[5]));
    return boost::python::make_tuple(s, res.stibium_stamp,
                                     res.sha256, res.from_cache);
}

void fab::onParseError(const fab::ParseError &e)
{
    (void)e;
    PyErr_SetString(PyExc_RuntimeError, "Failed to parse math expression");
}

void fab::onShapeError(const fab::ShapeError& e)
{
    (void)e;
    PyErr_SetString(PyExc_RuntimeError, "Could not construct Shape object.");
}

BOOST_PYTHON_MODULE(_fabtypes)
{
    class_<Bounds>("Bounds", init<>())
            .def(init<>())
            .def(init<float, float, float, float>())
            .def(init<float, float, float, float, float, float>())
            .def_readonly("xmin", &Bounds::xmin)
            .def_readonly("ymin", &Bounds::ymin)
            .def_readonly("zmin", &Bounds::zmin)
            .def_readonly("xmax", &Bounds::xmax)
            .def_readonly("ymax", &Bounds::ymax)
            .def_readonly("zmax", &Bounds::zmax)
            .def("is_bounded_xy", &Bounds::is_bounded_xy,
                 "Returns True if both x and y bounds are non-infinite.")
            .def("is_bounded_xyz", &Bounds::is_bounded_xyz,
                 "Returns True if both all bounds are non-infinite.");

    class_<Shape>("Shape", no_init)
            .def("__init__", raw_function(&Shape::init), "raw constructor")
            .def(init<std::string, Bounds, int, int, int>())
            .def_readonly("math", &Shape::math)
            .def_readonly("bounds", &Shape::bounds)
            .def_readwrite("_r", &Shape::r)
            .def_readwrite("_g", &Shape::g)
            .def_readwrite("_b", &Shape::b)
            .def("map", &Shape::map)
            .def("__repr__", &Shape::repr)
            .def(self & self)
            .def(self | self)
            .def(~self);


    class_<Transform>("Transform",
        init<std::string, std::string, std::string, std::string>())
            .def(init<std::string, std::string, std::string,
                      std::string, std::string, std::string>())
            .def_readonly("x_forward", &Transform::x_forward)
            .def_readonly("y_forward", &Transform::y_forward)
            .def_readonly("z_forward", &Transform::z_forward)
            .def_readonly("z_reverse", &Transform::z_reverse)
            .def_readonly("y_reverse", &Transform::y_reverse)
            .def_readonly("z_reverse", &Transform::z_reverse);

    def("project_dir", &project_dir_py,
        "Directory of the current project file ('' if unsaved)");
    def("measure", &measure_py,
        "Grid-integrate a Shape: volume/area, center of mass, tight "
        "bounds. measure(shape, resolution=-1)");
    def("_import_mesh", &import_mesh_py,
        "Backend for fab.shapes.import_mesh; "
        "returns (Shape, stibium_stamp, sha256, from_cache)");

    register_exception_translator<fab::ParseError>(fab::onParseError);
    register_exception_translator<fab::ShapeError>(fab::onShapeError);
}


PyTypeObject* fab::ShapeType = NULL;

void fab::preInit()
{
    PyImport_AppendInittab("_fabtypes", PyInit__fabtypes);
}

void fab::postInit(std::vector<std::string> fab_paths)
{
    PyObject* fabtypes = PyImport_ImportModule("_fabtypes");
    ShapeType = (PyTypeObject*)PyObject_GetAttrString(fabtypes, "Shape");

    // Modify the default search path to include the application's directory
    // (as this doesn't happen on Linux by default)
    for (auto p : fab_paths)
    {
        PyList_Insert(PySys_GetObject("path"), 0, PyUnicode_FromString(p.c_str()));
    }
}

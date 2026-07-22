#include <boost/python.hpp>
#include <boost/python/raw_function.hpp>

#include "fab/fab.h"
#include "fab/types/shape.h"
#include "fab/types/transform.h"
#include "fab/mesh/mesh_import.h"
#include "fab/tree/analytics.h"
#include "fab/tree/grid.h"
#include "fab/mesh/glyph_import.h"

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

/*  Backend for fab.shapes.glyph (see py/fab/shapes.py). Bakes one TTF glyph
 *  into an OP_GRID-backed Shape. Returns (Shape, from_cache). */
static tuple bake_glyph_py(std::string font_path, int codepoint,
                           float cap, float thickness, float voxels_per_mm)
{
    grid_registry_trim();
    const auto res = fab_mesh::bake_glyph_grid(
            font_path, codepoint, cap, thickness, voxels_per_mm);
    if (!res.grid_id)
    {
        PyErr_SetString(PyExc_RuntimeError,
                ("glyph: " + res.error).c_str());
        throw_error_already_set();
    }
    Shape s("g" + std::to_string(res.grid_id),
            Bounds(res.bounds[0], res.bounds[1], res.bounds[2],
                   res.bounds[3], res.bounds[4], res.bounds[5]));
    return boost::python::make_tuple(s, res.from_cache);
}

/*  Advance width (mm) of a glyph at cap height, for fab.shapes.text_font. */
static float glyph_advance_py(std::string font_path, int codepoint, float cap)
{
    return fab_mesh::glyph_advance(font_path, codepoint, cap);
}

/*  Backend for fab.shapes.redistance: rebuild a Shape as a true-distance
 *  OP_GRID. Returns (Shape, from_cache). The shape must be 3D-bounded. */
static tuple redistance_py(const Shape& s, float res, std::string key)
{
    grid_registry_trim();
    const Bounds& b = s.bounds;
    if (std::isinf(b.xmin) || std::isinf(b.xmax) || std::isinf(b.ymin) ||
        std::isinf(b.ymax) || std::isinf(b.zmin) || std::isinf(b.zmax))
    {
        PyErr_SetString(PyExc_RuntimeError,
                "redistance: shape has unbounded extents (extrude 2D first)");
        throw_error_already_set();
    }
    const auto r = fab_mesh::redistance_grid(s.tree.get(),
            b.xmin, b.ymin, b.zmin, b.xmax, b.ymax, b.zmax, res, key);
    if (!r.grid_id)
    {
        PyErr_SetString(PyExc_RuntimeError,
                ("redistance: " + r.error).c_str());
        throw_error_already_set();
    }
    Shape out("g" + std::to_string(r.grid_id),
              Bounds(r.bounds[0], r.bounds[1], r.bounds[2],
                     r.bounds[3], r.bounds[4], r.bounds[5]));
    return boost::python::make_tuple(out, r.from_cache);
}

/*  Backend for the Antimony font: bake a 2D Shape into an extruded grid via
 *  marching-squares contour + exact distance. Returns (Shape, from_cache). */
static tuple bake_shape_glyph_py(const Shape& s, float thickness,
                                 float res, std::string key)
{
    grid_registry_trim();
    const Bounds& b = s.bounds;
    if (std::isinf(b.xmin) || std::isinf(b.xmax) ||
        std::isinf(b.ymin) || std::isinf(b.ymax))
    {
        PyErr_SetString(PyExc_RuntimeError,
                "text_font: 2D shape has unbounded x/y extents");
        throw_error_already_set();
    }
    const auto r = fab_mesh::bake_shape_glyph(s.tree.get(),
            b.xmin, b.ymin, b.xmax, b.ymax, thickness, res, key);
    if (!r.grid_id)
    {
        PyErr_SetString(PyExc_RuntimeError,
                ("text_font (antimony): " + r.error).c_str());
        throw_error_already_set();
    }
    Shape out("g" + std::to_string(r.grid_id),
              Bounds(r.bounds[0], r.bounds[1], r.bounds[2],
                     r.bounds[3], r.bounds[4], r.bounds[5]));
    return boost::python::make_tuple(out, r.from_cache);
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
    def("_bake_glyph", &bake_glyph_py,
        "Backend for fab.shapes.glyph; bakes a TTF glyph into an "
        "OP_GRID Shape. returns (Shape, from_cache)");
    def("_glyph_advance", &glyph_advance_py,
        "Advance width (mm) of a TTF glyph at the given cap height");
    def("_redistance", &redistance_py,
        "Rebuild a Shape as a true-distance OP_GRID (Felzenszwalb EDT). "
        "returns (Shape, from_cache)");
    def("_bake_shape_glyph", &bake_shape_glyph_py,
        "Bake a 2D Shape into an extruded OP_GRID via marching-squares "
        "contour + exact distance. returns (Shape, from_cache)");

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

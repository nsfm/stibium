#include <Python.h>

#include <QVector>
#include <QPair>

#include "app/colors.h"

#include "graph/datum.h"
#include "fab/fab.h"

namespace Colors
{
QColor red("#e05252");
QColor orange("#e8883f");
QColor yellow("#f2c14e");
QColor green("#97bf70");
QColor teal("#4fc1b0");
QColor blue("#5a9fe0");
QColor violet("#b183d6");
QColor brown("#c08552");
QColor amber("#f0a63c");    // UI accent: selection, focus, hover

// Warm charcoal ramp (dark -> light)
QColor base00("#16130e");
QColor base01("#1e1a14");
QColor base02("#2a241b");
QColor base03("#4d4536");
QColor base04("#b5aa96");
QColor base05("#d8cfbe");
QColor base06("#eae3d5");
QColor base07("#f8f4eb");

QColor adjust(QColor c, float scale)
{
    // Lightness-based shift: multiplicative RGB scaling can't brighten
    // dark colors (0 * anything is 0), so hover highlights vanished on
    // dark fills.
    return scale >= 1 ? c.lighter(int(scale * 100))
                      : c.darker(int(100 / scale));
}

QColor highlight(QColor c)
{
    return adjust(c, 1.4);
}

QColor dim(QColor c)
{
    return adjust(c, 1/1.4);
}

QColor getColor(const Datum *d)
{
    auto t = d->getType();
    if (t == &PyUnicode_Type)
        return brown;
    else if (t == &PyFloat_Type)
        return yellow;
    else if (t == &PyLong_Type)
        return orange;
    else if (t == fab::ShapeType)
        return green;
    else
        return red;
}

PyObject* PyColors()
{
    // Here are all of our standard colors and their names.
    static QVector<QPair<QString, QColor>> colors = {
        {"red", red},
        {"orange", orange},
        {"yellow", yellow},
        {"green", green},
        {"teal", teal},
        {"blue", blue},
        {"violet", violet},
        {"brown", brown},
        {"amber", amber},
        {"base00", base00},
        {"base01", base01},
        {"base02", base02},
        {"base03", base03},
        {"base04", base04},
        {"base05", base05},
        {"base06", base06},
        {"base07", base07}};

    // Lazy initialization of NamedTuple constructor
    static PyObject* colors_tuple = NULL;
    if (colors_tuple == NULL)
    {
        PyObject* tuple_constructor;

        // Build a namedtuple constructor that has all of the colors
        // as arguments.
        {
            auto list = PyList_New(colors.size());
            size_t i=0;
            for (auto c : colors)
                PyList_SetItem(list, i++, Py_BuildValue(
                            "s", c.first.toStdString().c_str()));

            auto collections = PyImport_ImportModule("collections");
            auto nt = PyObject_GetAttrString(collections, "namedtuple");
            auto args = Py_BuildValue("(sO)", "SbColors", list);
            tuple_constructor = PyObject_Call(nt, args, NULL);

            for (auto o : {collections, nt, args, list})
                Py_DECREF(o);
        }

        // Then, call this constructor on a list of color tuples.
        {
            auto list = PyList_New(colors.size());
            size_t i=0;
            for (auto c : colors)
                PyList_SetItem(list, i++, Py_BuildValue(
                        "(iii)", c.second.red(),
                        c.second.green(), c.second.blue()));
            auto args = PyList_AsTuple(list);
            colors_tuple = PyObject_Call(tuple_constructor, args, NULL);

            for (auto o : {list, args})
                Py_DECREF(o);
        }
        Q_ASSERT(!PyErr_Occurred());
    }

    return colors_tuple;
}

} // end of Colors namespace

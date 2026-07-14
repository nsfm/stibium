#pragma once

#include <Python.h>

#include <QImage>
#include <QMatrix4x4>
#include <QString>

#include <memory>

#include "fab/types/shape.h"

class Graph;

namespace ImageExport
{

struct Options
{
    QMatrix4x4 M;           /*  view transform (rotation is what counts;
                                the model is auto-fit to the image)  */
    float section = 1;      /*  screen-parallel cut (1 = no cut)  */
    float resolution = -1;  /*  voxels per unit; <= 0 fits fit_px  */
    int fit_px = 512;       /*  longest image side when fitting  */
    int supersample = 2;    /*  antialiasing: render at N x then
                                smooth-downscale (1 = off)  */
    bool fit_sphere = false;/*  frame the model's circumsphere
                                instead of its per-view bounds, so
                                rotating views keep a fixed framing
                                (turntables)  */
    bool transparent = true;/*  transparent vs. theme background  */
    QString filename;       /*  output image (format by extension)  */
    QString node_name;      /*  non-empty: render only this node's
                                shape outputs (even non-terminal) -
                                the visual-bisection mode  */
};

/*
 *  Renders the graph's terminal shapes (what the viewport shows) to
 *  an image, shaded, with no GL and no display needed.  Each shape
 *  renders separately and composites by depth, so per-shape colors
 *  (set_color) survive into the image.
 *  Returns an empty string on success, else an error message.
 */
QString render(Graph* graph, const Options& opt);

/*
 *  Renders an explicit list of shapes (all 3D or all 2D, flagged by
 *  flat) with the same compositing pipeline; used by --diff.
 */
QString renderShapes(const std::vector<Shape>& shapes, bool flat,
                     const Options& opt);

/*
 *  Same, but returns the image instead of saving it (animation
 *  frames).  On failure returns a null image and sets *err.
 */
QImage renderShapesImage(const std::vector<Shape>& shapes, bool flat,
                         const Options& opt, QString* err);

/*
 *  Collects the graph's renderable shapes (terminal outputs only,
 *  unless node_name selects one node), 3D and 2D separately, without
 *  unioning them.  Returns an error message, or empty on success.
 */
QString collectShapeList(Graph* graph, const QString& node_name,
                         std::vector<Shape>& shapes3d,
                         std::vector<Shape>& shapes2d);

/*
 *  Unions the graph's renderable shapes, 3D and 2D separately
 *  (terminal outputs only, unless node_name selects one node).
 *  Returns an error message, or empty on success.
 */
QString collectShapes(Graph* graph, const QString& node_name,
                      std::unique_ptr<Shape>& u3d,
                      std::unique_ptr<Shape>& u2d);

}  // namespace ImageExport

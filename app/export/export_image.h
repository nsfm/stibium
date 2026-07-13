#pragma once

#include <Python.h>

#include <QMatrix4x4>
#include <QString>

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
    bool transparent = true;/*  transparent vs. theme background  */
    QString filename;       /*  output image (format by extension)  */
    QString node_name;      /*  non-empty: render only this node's
                                shape outputs (even non-terminal) -
                                the visual-bisection mode  */
};

/*
 *  Renders the graph's terminal shapes (what the viewport shows) to
 *  an image, shaded, with no GL and no display needed.
 *  Returns an empty string on success, else an error message.
 */
QString render(Graph* graph, const Options& opt);

}  // namespace ImageExport

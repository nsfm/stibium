#ifndef SVG_H
#define SVG_H

#include "fab/tree/contour.h"

/*
 *  Exports closed contour loops to an SVG file.
 *
 *  Loops become subpaths of one filled path (fill-rule evenodd, so
 *  holes render as holes regardless of orientation).  Coordinates are
 *  interpreted as millimeters; the viewBox spans the given bounds and
 *  the y axis is flipped (SVG is y-down, models are y-up).
 *
 *  Returns true on success, false if the file couldn't be written.
 */
bool save_svg(const std::vector<ContourPath>& paths,
              float xmin, float ymin, float xmax, float ymax,
              const char* filename);

#endif

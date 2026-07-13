#ifndef DXF_H
#define DXF_H

#include "fab/tree/contour.h"

/*
 *  Exports closed contour loops to a DXF file (R12 dialect: one
 *  closed POLYLINE entity per loop, layer 0).  R12 is the most widely
 *  read DXF flavor across laser-cutter and CAM toolchains.
 *
 *  Coordinates are written as-is (y-up, model units; by convention
 *  laser tools read them as millimeters).
 *
 *  Returns true on success, false if the file couldn't be written.
 */
bool save_dxf(const std::vector<ContourPath>& paths,
              float xmin, float ymin, float xmax, float ymax,
              const char* filename);

#endif

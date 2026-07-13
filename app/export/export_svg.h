#pragma once

#include <Python.h>

#include "export/export_worker.h"
#include "fab/types/shape.h"

////////////////////////////////////////////////////////////////////////////////

class ExportSvgWorker : public ExportWorker
{
public:
    explicit ExportSvgWorker(Shape s, Bounds b, QString f, float r,
                             bool detect=true, float sim=-1)
        : ExportWorker(s, b, f, r), detect_features(detect),
          simplify(sim) {}

    /*
     *  Top-level function that accepts user input and start the export
     */
    void run() override;

    /*
     *  Asynchronous function that traces contours and writes the file
     */
    void async() override;

protected:
    /*
     *  Recover sharp corners on the traced contours
     *  (on by default: vector output exists to be exact)
     */
    const bool detect_features;

    /*
     *  Maximum deviation (model units) for path simplification.
     *  Negative means automatic (a quarter cell); 0 disables.
     */
    const float simplify;
};

#pragma once

#include <Python.h>

#include "export/export_worker.h"
#include "fab/types/shape.h"

////////////////////////////////////////////////////////////////////////////////

class ExportSvgWorker : public ExportWorker
{
public:
    explicit ExportSvgWorker(Shape s, Bounds b, QString f, float r,
                             bool detect=true)
        : ExportWorker(s, b, f, r), detect_features(detect) {}

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
};

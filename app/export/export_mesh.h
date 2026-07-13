#pragma once

#include <Python.h>

#include <cstdint>
#include <vector>

#include "export/export_worker.h"
#include "fab/types/shape.h"

////////////////////////////////////////////////////////////////////////////////

class ExportMeshWorker : public ExportWorker
{
public:
    explicit ExportMeshWorker(Shape s, Bounds b, QString f, float r, bool d,
                              float sim=0)
        : ExportWorker(s, b, f, r), detect_features(d), simplify(sim) {}

    /*
     *  Top-level function that accepts user input and start the export
     */
    void run() override;

    /*
     *  Asynchronous function that actually performs rendering
     */
    void async() override;

    bool runHeadless(const QString& fname, float res,
                     int detect) override;

protected:
    /*
     *  Simplifies the indexed mesh to within _simplify model units of
     *  deviation, replacing verts / indices with the (much smaller)
     *  result.
     */
    void simplifyMesh(std::vector<float>& verts,
                      std::vector<uint32_t>& indices) const;

    /*
     *  Call-time settings
     */
    const bool detect_features;
    const float simplify;

    /*
     *  Run-time, set by dialogs
     */
    bool _detect_features;

    /*
     *  Maximum geometric deviation (in model units) allowed when
     *  simplifying the mesh after triangulation; <= 0 disables it.
     */
    float _simplify;
};

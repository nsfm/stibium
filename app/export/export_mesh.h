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

    int facedevHeadless(const QString& mesh, float res) override;

protected:
    /*
     *  Simplifies the indexed mesh to within _simplify model units of
     *  deviation, replacing verts / indices with the (much smaller)
     *  result.
     */
    void simplifyMesh(float deviation, std::vector<float>& verts,
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

    /*  Stibnite integration (2026-07-18): mesher routing + the
     *  curated advanced knobs.  The dialog path applies these to
     *  the process environment before meshing (one export at a
     *  time is law, so setenv is safe); the headless path leaves
     *  the environment untouched and env-driven as ever.  */
    int _mesher = 0;                // 0 Stibnite, 1 Classic
    bool _from_dialog = false;
    bool _adv_autodense = true;
    int _adv_density_cap = 2;
    bool _adv_decimate = true;
    bool _adv_snap = true;
    int _adv_stall = 1;

    /*  Post-export report, filled by async(), shown by run() on
     *  the GUI thread (empty = nothing to show).  */
    QString _stats;
};

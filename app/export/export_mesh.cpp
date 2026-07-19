#include <Python.h>

#include <atomic>
#include <chrono>
#include <cstring>
#include <thread>
#include <vector>

#include <QFileDialog>
#include <QFileInfo>

#include "app/settings.h"
#include <QMessageBox>

#include "export/export_mesh.h"
#include "vendor/meshoptimizer/meshoptimizer.h"

#include "dialog/resolution.h"

#include "fab/util/region.h"
#include "fab/tree/triangulate.h"
#include "fab/tree/triangulate/delaunay.h"
#include "fab/formats/stl.h"
#include "fab/formats/threemf.h"

////////////////////////////////////////////////////////////////////////////////

void ExportMeshWorker::run()
{
    // Sanity-check bounds
    if (std::isinf(bounds.xmin) || std::isinf(bounds.xmax) ||
        std::isinf(bounds.ymin) || std::isinf(bounds.ymax) ||
        std::isinf(bounds.zmin) || std::isinf(bounds.zmax))
    {
        QMessageBox::critical(NULL, "Export error",
                "<b>Export error:</b><br>"
                "Target shape has invalid (infinite) bounds");
        return;
    }

    // Get resolution, either hardcoded or from the user
    if (resolution == -1)
    {
        auto resolution_dialog = new ResolutionDialog(
                bounds, RESOLUTION_DIALOG_3D, UNITLESS, 1 << 22, NULL, 1);
        if (!resolution_dialog->exec())
            return;
        _resolution = resolution_dialog->getResolution();
        _detect_features = resolution_dialog->getDetectFeatures();
        _simplify = resolution_dialog->getSimplifyDeviation();
        _mesher = resolution_dialog->getMesher();
        _adv_autodense = resolution_dialog->getAutodense();
        _adv_density_cap = resolution_dialog->getDensityCap();
        _adv_decimate = resolution_dialog->getDecimate();
        _adv_snap = resolution_dialog->getSnap();
        _adv_stall = resolution_dialog->getStallPatience();
        _from_dialog = true;
        delete resolution_dialog;
    }
    else
    {
        _resolution = resolution;
        _detect_features = detect_features;
        _simplify = simplify;
        _mesher = getenv("STIBIUM_EXPORT_CLASSIC") ? 1 : 0;
    }

    if (_resolution == 0)
    {
        QMessageBox::critical(NULL, "Export error",
                "<b>Export error:</b><br>"
                "Resolution cannot be set to 0");
        return;
    }

    //  Get a target filename, either hardcoded or from the user
    if (filename.isEmpty())
    {
        QString filter = "3MF (*.3mf)";
        _filename = QFileDialog::getSaveFileName(
                NULL, "Export mesh",
                Settings::get("files/last_export_dir", "").toString(),
                "3MF (*.3mf);;STL (*.stl)", &filter);

        // If no recognized extension was typed, take it from the filter
        if (!_filename.isEmpty() &&
            !_filename.endsWith(".3mf", Qt::CaseInsensitive) &&
            !_filename.endsWith(".stl", Qt::CaseInsensitive))
        {
            _filename += filter.contains("*.stl") ? ".stl" : ".3mf";
        }
    }
    else
    {
        _filename = filename;
    }
    if (_filename.isEmpty())
        return;
    Settings::set("files/last_export_dir",
                  QFileInfo(_filename).absolutePath());

    if (checkWritable())
    {
        runAsync();
        /*  Post-export report (Nate's stats-dialog ask,
         *  2026-07-18): async() fills _stats; we are back on the
         *  GUI thread here.  Silent when cancelled or classic.  */
        if (!_stats.isEmpty() && !halt)
            QMessageBox::information(NULL, "Export complete",
                                     _stats);
    }
}

////////////////////////////////////////////////////////////////////////////////

bool ExportMeshWorker::runHeadless(const QString& fname, float res,
                                   int detect)
{
    if (std::isinf(bounds.xmin) || std::isinf(bounds.xmax) ||
        std::isinf(bounds.ymin) || std::isinf(bounds.ymax) ||
        std::isinf(bounds.zmin) || std::isinf(bounds.zmax))
    {
        fprintf(stderr, "export: shape has infinite bounds\n");
        return false;
    }

    _filename = fname.isEmpty() ? filename : fname;
    _resolution = res > 0 ? res : resolution;
    _detect_features = detect < 0 ? detect_features : bool(detect);
    _simplify = simplify;
    /*  Headless default is Stibnite (integration 2026-07-18);
     *  STIBIUM_EXPORT_CLASSIC=1 routes to the classic sampler.
     *  STIBIUM_EXPORT_DMESH=1 remains accepted (now a no-op) so
     *  every documented harness command keeps working.  Knobs
     *  stay env-driven on this path.  */
    _mesher = getenv("STIBIUM_EXPORT_CLASSIC") ? 1 : 0;
    _from_dialog = false;

    if (_filename.isEmpty())
    {
        fprintf(stderr, "export: no filename (pass --export FILE or "
                        "set filename= in the script)\n");
        return false;
    }
    if (_resolution <= 0)
    {
        fprintf(stderr, "export: no resolution (pass --resolution R or "
                        "set resolution= in the script)\n");
        return false;
    }

    async();
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void ExportMeshWorker::async()
{
    Region r = {};
    r.ni = uint32_t((bounds.xmax - bounds.xmin) * _resolution);
    r.nj = uint32_t((bounds.ymax - bounds.ymin) * _resolution);
    r.nk = uint32_t((bounds.zmax - bounds.zmin) * _resolution);
    r.voxels = uint64_t(r.ni) * r.nj * r.nk;

    build_arrays(
            &r, bounds.xmin, bounds.ymin, bounds.zmin,
                bounds.xmax, bounds.ymax, bounds.zmax);

    // The mesher produces an indexed mesh directly (unique vertices +
    // three indices per triangle), so no welding pass is needed.
    progress_phase = PHASE_MESHING;
    progress_total = r.voxels;
    std::vector<float> verts;
    std::vector<uint32_t> indices;
    _stats.clear();
    /*  Stibnite (the adaptive-Delaunay mesher, doc/MESH-WAR.md) is
     *  the default; Classic is the dialog's other choice or
     *  STIBIUM_EXPORT_CLASSIC=1 headless.  */
    if (_mesher == 0)
    {
        /*  Dialog knobs -> process environment.  Safe: one export
         *  at a time is law, and the headless path never enters
         *  here with _from_dialog set.  */
        if (_from_dialog)
        {
            qputenv("STIBIUM_DMESH_AUTODENSE",
                    _adv_autodense ? "1" : "0");
            qputenv("STIBIUM_DMESH_AUTODENSE_MAX",
                    QByteArray::number(_adv_density_cap));
            qputenv("STIBIUM_DMESH_DECIMATE",
                    _adv_decimate ? "1" : "0");
            qputenv("STIBIUM_DMESH_SNAP", _adv_snap ? "1" : "0");
            qputenv("STIBIUM_DMESH_STALL",
                    QByteArray::number(_adv_stall));
        }

        /*  Mirror the mesher's progress sink into the dialog's
         *  bar until meshing returns: overall fraction, stage
         *  index, and a remaining-time estimate.
         *
         *  The raw estimate (elapsed x (1-f)/f) whipsaws while the
         *  early stages burn through their phase weights, so the
         *  display runs on a LOCKED DEADLINE: past 15% overall the
         *  smoothed estimate sets a finish time and the label just
         *  counts down toward it (smooth by construction).  The
         *  deadline only moves for sustained disagreement - pushed
         *  out when the smoothed estimate exceeds the countdown by
         *  25%, pulled in when it drops below half.  Overshooting
         *  the deadline shows "wrapping up" (eta sentinel -2)
         *  instead of a lie.  */
        progress_total = 1000;
        std::atomic<bool> mesh_done{ false };
        std::thread mirror([&]() {
            const auto m0 = std::chrono::steady_clock::now();
            double ema = -1;
            double deadline = -1;   // seconds since m0
            while (!mesh_done.load())
            {
                const float f = dmesh_progress()->overall.load();
                progress_done = uint64_t(f * 1000);
                progress_stage = dmesh_progress()->stage.load();
                const double el = std::chrono::duration<double>(
                        std::chrono::steady_clock::now() - m0)
                        .count();
                if (f > 0.15f && f < 1.f)
                {
                    const double raw = el * (1.f - f) / f;
                    ema = ema < 0 ? raw : 0.9 * ema + 0.1 * raw;
                    const double left0 = deadline - el;
                    if (deadline < 0 || ema > left0 * 1.25 ||
                        ema < left0 * 0.5)
                        deadline = el + ema;
                    const double left = deadline - el;
                    progress_eta_s = left >= 1 ? int(left) : -2;
                }
                else
                    progress_eta_s = -1;
                std::this_thread::sleep_for(
                        std::chrono::milliseconds(100));
            }
            progress_stage = -1;
            progress_eta_s = -1;
        });

        const auto t0 = std::chrono::steady_clock::now();
        Deck* deck = deck_from_tree(shape.tree.get());
        DMesh m;
        const bool ok = delaunay_mesh(deck, r, &halt, &m);
        deck_free(deck);
        mesh_done = true;
        mirror.join();
        const double wall = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();

        fprintf(stderr, "stibnite mesh report: %s, %zu tris, "
                "%llu open, "
                "%llu non-manifold, %llu constrained, %llu steiner, "
                "%llu repaired, %llu split, %llu snapped\n",
                ok ? "ok" : "FAILED (built without CGAL?)",
                m.tris.size() / 3,
                (unsigned long long)m.open_edges,
                (unsigned long long)m.nonmanifold_edges,
                (unsigned long long)m.constrained,
                (unsigned long long)m.steiner,
                (unsigned long long)m.repaired,
                (unsigned long long)m.split_verts,
                (unsigned long long)m.snapped);
        if (ok && !halt)
            _stats = QString(
                    "<b>Stibnite mesh report</b><br><br>"
                    "%1 triangles, %2 vertices<br>"
                    "open edges: %3<br>"
                    "feature edges constrained: %4<br>"
                    "repairs: %5 &nbsp; crease snaps: %6<br>"
                    "meshing time: %7 s")
                    .arg(m.tris.size() / 3)
                    .arg(m.verts.size() / 3)
                    .arg(m.open_edges)
                    .arg(m.constrained)
                    .arg(m.repaired)
                    .arg(m.snapped)
                    .arg(wall, 0, 'f', 1);
        verts = std::move(m.verts);
        indices = std::move(m.tris);
    }
    else
        triangulate_indexed_mt(shape.tree.get(), r, _detect_features,
                               &halt, verts, indices, -1,
                               &progress_done);

    // Simplification and file writing can't report granular progress;
    // drop the bar back to its busy animation instead of a full bar.
    progress_total = 0;

    /*  QEM on the dmesh path (quadric assessment 2026-07-17):
     *  STIBIUM_DMESH_SIMPLIFY=<mm> forces the vendored meshopt
     *  pass on dmesh exports - deviation-bounded, in model
     *  units.  Nate's standalone experiment (50% off, no visual
     *  defects, twice) green-lit QEM for this geometry; this
     *  delivers it from code already in the tree.  Phase 2
     *  (crease vertex_lock) is designed in doc/reviews/.  */
    float simplify = _simplify;
    if (_mesher == 0)
        if (const char* se = getenv("STIBIUM_DMESH_SIMPLIFY"))
            simplify = float(atof(se));
    if (simplify > 0 && indices.size() >= 3 && !halt)
    {
        progress_phase = PHASE_SIMPLIFYING;
        simplifyMesh(simplify, verts, indices);
        if (!_stats.isEmpty())
            _stats += QString("<br>after simplify: %1 triangles, "
                              "%2 vertices")
                    .arg(indices.size() / 3)
                    .arg(verts.size() / 3);
    }

    /*  Pinch split LAST (strictly after QEM): meshopt tears
     *  coincident sheet copies into real boundary holes if the
     *  split runs first, and its collapses mint fresh pinches of
     *  their own (bino 274 -> 325 nm edges, measured) - the tail
     *  position cures both.  Zero vertex motion; index-level
     *  only, which 3MF preserves.  */
    if (_mesher == 0 && indices.size() >= 3 && !halt)
    {
        const uint64_t copies = dmesh_split_pinches(verts, indices);
        if (copies && !_stats.isEmpty())
            _stats += QString("<br>pinch seams split: %1 sheet "
                              "copies").arg(copies);
    }

    progress_phase = PHASE_WRITING;

    // Format follows the file extension; STL is the fallback for
    // scripted exports with unrecognized names (the legacy behavior).
    if (_filename.endsWith(".3mf", Qt::CaseInsensitive))
        save_3mf_indexed(verts.data(), verts.size() / 3,
                         indices.data(), indices.size() / 3,
                         _filename.toStdString().c_str());
    else
        save_stl_indexed(verts.data(), indices.data(), indices.size() / 3,
                         _filename.toStdString().c_str());
    free_arrays(&r);
}

////////////////////////////////////////////////////////////////////////////////

void ExportMeshWorker::simplifyMesh(float deviation,
                                    std::vector<float>& verts,
                                    std::vector<uint32_t>& indices) const
{
    const size_t vertex_count = verts.size() / 3;

    // target_error is relative to mesh extent; deviation is in model units
    const float scale = meshopt_simplifyScale(
            verts.data(), vertex_count, 3 * sizeof(float));
    const float target_error = scale > 0 ? deviation / scale : deviation;

    std::vector<uint32_t> simplified(indices.size());
    const size_t out_index_count = meshopt_simplify(
            simplified.data(), indices.data(), indices.size(),
            verts.data(), vertex_count, 3 * sizeof(float),
            0, target_error, meshopt_SimplifyPrune, NULL);

    simplified.resize(out_index_count);
    indices = std::move(simplified);

    // Drop the vertices orphaned by simplification (renumbering the
    // rest in first-use order), so indexed formats don't carry dead
    // entries.
    std::vector<uint32_t> remap(vertex_count, UINT32_MAX);
    std::vector<float> packed;
    packed.reserve(verts.size());
    for (auto& i : indices)
    {
        if (remap[i] == UINT32_MAX)
        {
            remap[i] = packed.size() / 3;
            packed.insert(packed.end(),
                          verts.begin() + i*3, verts.begin() + i*3 + 3);
        }
        i = remap[i];
    }
    verts = std::move(packed);
}

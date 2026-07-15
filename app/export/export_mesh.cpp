#include <Python.h>

#include <cstring>
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
                bounds, RESOLUTION_DIALOG_3D, UNITLESS, 1 << 22, NULL, 7);
        if (!resolution_dialog->exec())
            return;
        _resolution = resolution_dialog->getResolution();
        _detect_features = resolution_dialog->getDetectFeatures();
        _simplify = resolution_dialog->getSimplifyDeviation();
        delete resolution_dialog;
    }
    else
    {
        _resolution = resolution;
        _detect_features = detect_features;
        _simplify = simplify;
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
        runAsync();
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
    // Meshing is chunked across all cores; the progress counters feed
    // the export dialog's bar.
    progress_phase = PHASE_MESHING;
    progress_total = r.voxels;
    std::vector<float> verts;
    std::vector<uint32_t> indices;
    /*  STIBIUM_EXPORT_DMESH=1 routes the export through the
     *  adaptive-Delaunay mesher (doc/MESH-NEXT.md) instead of the
     *  production Kobbelt path - the final-quality pipeline, single
     *  threaded, no progress reporting.  Experimental.  */
    if (getenv("STIBIUM_EXPORT_DMESH"))
    {
        Deck* deck = deck_from_tree(shape.tree.get());
        DMesh m;
        const bool ok = delaunay_mesh(deck, r, &halt, &m);
        deck_free(deck);
        fprintf(stderr, "dmesh export: %s, %zu tris, %llu open, "
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

    if (_simplify > 0 && indices.size() >= 3 && !halt)
    {
        progress_phase = PHASE_SIMPLIFYING;
        simplifyMesh(verts, indices);
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

void ExportMeshWorker::simplifyMesh(std::vector<float>& verts,
                                    std::vector<uint32_t>& indices) const
{
    const size_t vertex_count = verts.size() / 3;

    // target_error is relative to mesh extent; _simplify is in model units
    const float scale = meshopt_simplifyScale(
            verts.data(), vertex_count, 3 * sizeof(float));
    const float target_error = scale > 0 ? _simplify / scale : _simplify;

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

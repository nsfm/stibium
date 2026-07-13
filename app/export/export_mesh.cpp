#include <Python.h>

#include <cstring>
#include <vector>

#include <QFileDialog>
#include <QMessageBox>

#include "export/export_mesh.h"
#include "vendor/meshoptimizer/meshoptimizer.h"

#include "dialog/resolution.h"

#include "fab/util/region.h"
#include "fab/tree/triangulate.h"
#include "fab/formats/stl.h"

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
                bounds, RESOLUTION_DIALOG_3D, UNITLESS);
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
        _filename = QFileDialog::getSaveFileName(
                NULL, "Export STL", "", "*.stl");
    else
        _filename = filename;
    if (_filename.isEmpty())
        return;

    if (checkWritable())
        runAsync();
}

////////////////////////////////////////////////////////////////////////////////

void ExportMeshWorker::async()
{
    Region r = (Region){
        .imin=0, .jmin=0, .kmin=0,
        .ni=uint32_t((bounds.xmax - bounds.xmin) * _resolution),
        .nj=uint32_t((bounds.ymax - bounds.ymin) * _resolution),
        .nk=uint32_t((bounds.zmax - bounds.zmin) * _resolution),
    };
    r.voxels = r.ni * r.nj * r.nk;

    build_arrays(
            &r, bounds.xmin, bounds.ymin, bounds.zmin,
                bounds.xmax, bounds.ymax, bounds.zmax);

    // The mesher produces an indexed mesh directly (unique vertices +
    // three indices per triangle), so no welding pass is needed.
    std::vector<float> verts;
    std::vector<uint32_t> indices;
    triangulate_indexed(shape.tree.get(), r, _detect_features, &halt,
                        verts, indices);

    if (_simplify > 0 && indices.size() >= 3 && !halt)
        simplifyMesh(verts, indices);

    save_stl_indexed(verts.data(), indices.data(), indices.size() / 3,
                     _filename.toStdString().c_str());
    free_arrays(&r);
}

////////////////////////////////////////////////////////////////////////////////

void ExportMeshWorker::simplifyMesh(const std::vector<float>& verts,
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
}

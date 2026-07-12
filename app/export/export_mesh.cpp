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
    float* verts;
    unsigned count;

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

    triangulate(shape.tree.get(), r, _detect_features, &halt, &verts, &count);

    if (_simplify > 0 && count >= 9 && !halt)
        simplifyMesh(&verts, &count);

    save_stl(verts, count, _filename.toStdString().c_str());
    free_arrays(&r);
    free(verts);
}

////////////////////////////////////////////////////////////////////////////////

void ExportMeshWorker::simplifyMesh(float** verts, unsigned* count) const
{
    // The mesher emits raw triangle soup (three floats per vertex, three
    // vertices per triangle, with every vertex duplicated per face), so
    // weld it into an indexed mesh before simplification.
    const size_t index_count = *count / 3;
    std::vector<unsigned int> remap(index_count);
    const size_t vertex_count = meshopt_generateVertexRemap(
            remap.data(), NULL, index_count,
            *verts, index_count, 3 * sizeof(float));

    std::vector<float> vertices(vertex_count * 3);
    meshopt_remapVertexBuffer(vertices.data(), *verts, index_count,
                              3 * sizeof(float), remap.data());
    // Soup has implicit sequential indices, so the remap table is
    // already the welded index buffer.

    // target_error is relative to mesh extent; _simplify is in model units
    const float scale = meshopt_simplifyScale(
            vertices.data(), vertex_count, 3 * sizeof(float));
    const float target_error = scale > 0 ? _simplify / scale : _simplify;

    std::vector<unsigned int> simplified(index_count);
    const size_t out_index_count = meshopt_simplify(
            simplified.data(), remap.data(), index_count,
            vertices.data(), vertex_count, 3 * sizeof(float),
            0, target_error, meshopt_SimplifyPrune, NULL);

    // Expand back into flat triangle soup for save_stl
    const unsigned out_count = out_index_count * 3;
    float* out = static_cast<float*>(malloc(out_count * sizeof(float)));
    for (size_t i=0; i < out_index_count; ++i)
        memcpy(&out[i*3], &vertices[simplified[i]*3], 3 * sizeof(float));

    free(*verts);
    *verts = out;
    *count = out_count;
}

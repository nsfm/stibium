#include <Python.h>

#include <cmath>

#include <QFileDialog>
#include <QMessageBox>

#include "export/export_svg.h"

#include "dialog/resolution.h"

#include "fab/formats/svg.h"
#include "fab/tree/contour.h"

////////////////////////////////////////////////////////////////////////////////

void ExportSvgWorker::run()
{
    if (std::isinf(bounds.xmin) || std::isinf(bounds.xmax) ||
        std::isinf(bounds.ymin) || std::isinf(bounds.ymax))
    {
        QMessageBox::critical(NULL, "Export error",
                "<b>Export error:</b><br>"
                "Target shape has invalid (infinite) X or Y bounds");
        return;
    }

    // Get resolution, either hardcoded or from the user
    if (resolution == -1)
    {
        auto resolution_dialog = new ResolutionDialog(
                bounds, RESOLUTION_DIALOG_2D, UNITLESS);
        if (!resolution_dialog->exec())
            return;
        _resolution = resolution_dialog->getResolution();
        delete resolution_dialog;
    }
    else
    {
        _resolution = resolution;
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
        _filename = QFileDialog::getSaveFileName(
                NULL, "Export SVG", "", "SVG (*.svg)");
        if (!_filename.isEmpty() &&
            !_filename.endsWith(".svg", Qt::CaseInsensitive))
        {
            _filename += ".svg";
        }
    }
    else
    {
        _filename = filename;
    }
    if (_filename.isEmpty())
        return;

    if (checkWritable())
        runAsync();
}

////////////////////////////////////////////////////////////////////////////////

void ExportSvgWorker::async()
{
    const uint32_t nx = uint32_t((bounds.xmax - bounds.xmin) * _resolution);
    const uint32_t ny = uint32_t((bounds.ymax - bounds.ymin) * _resolution);

    // 2D shapes have unbounded Z and evaluate ignoring it; slicing a
    // 3D shape takes the cross-section at its Z midpoint.
    const float z =
        (std::isinf(bounds.zmin) || std::isinf(bounds.zmax))
            ? 0
            : (bounds.zmin + bounds.zmax) / 2;

    progress_phase = PHASE_CONTOURING;
    progress_total = uint64_t(ny) + 1;

    std::vector<ContourPath> paths;
    contour_field(shape.tree.get(),
                  bounds.xmin, bounds.ymin, bounds.xmax, bounds.ymax,
                  nx, ny, z, detect_features, &halt, paths,
                  &progress_done);

    if (halt)
        return;

    progress_phase = PHASE_WRITING;
    progress_total = 0;

    save_svg(paths, bounds.xmin, bounds.ymin, bounds.xmax, bounds.ymax,
             _filename.toStdString().c_str());
}

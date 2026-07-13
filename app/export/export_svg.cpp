#include <Python.h>

#include <cmath>

#include <QFileDialog>
#include <QMessageBox>

#include "export/export_svg.h"

#include "dialog/resolution.h"

#include "fab/formats/dxf.h"
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
        QString filter = "SVG (*.svg)";
        _filename = QFileDialog::getSaveFileName(
                NULL, "Export vector", "",
                "SVG (*.svg);;DXF (*.dxf)", &filter);

        // If no recognized extension was typed, take it from the filter
        if (!_filename.isEmpty() &&
            !_filename.endsWith(".svg", Qt::CaseInsensitive) &&
            !_filename.endsWith(".dxf", Qt::CaseInsensitive))
        {
            _filename += filter.contains("*.dxf") ? ".dxf" : ".svg";
        }
    }
    else
    {
        _filename = filename;
    }
    if (_filename.isEmpty())
        return;

    _detect_features = detect_features;

    if (checkWritable())
        runAsync();
}

////////////////////////////////////////////////////////////////////////////////

bool ExportSvgWorker::runHeadless(const QString& fname, float res,
                                  int detect)
{
    if (std::isinf(bounds.xmin) || std::isinf(bounds.xmax) ||
        std::isinf(bounds.ymin) || std::isinf(bounds.ymax))
    {
        fprintf(stderr, "export: shape has infinite XY bounds\n");
        return false;
    }

    _filename = fname.isEmpty() ? filename : fname;
    _resolution = res > 0 ? res : resolution;

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

    // detect < 0 keeps the script's setting (default on for vectors)
    _detect_features = detect < 0 ? detect_features : bool(detect);

    async();
    return true;
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

    std::vector<ContourPath> paths;
    contour_field(shape.tree.get(),
                  bounds.xmin, bounds.ymin, bounds.xmax, bounds.ymax,
                  nx, ny, z, _detect_features, &halt, paths,
                  -1, &progress_done, &progress_total);

    if (halt)
        return;

    progress_phase = PHASE_WRITING;
    progress_total = 0;

    // Drop redundant collinear points (marching squares emits one
    // chord per cell even along straight runs); corners are maximal-
    // deviation points, so they always survive.
    const float tolerance = simplify < 0 ? 0.25f / _resolution : simplify;
    simplify_contours(paths, tolerance);

    // Format follows the file extension; SVG is the fallback for
    // scripted exports with unrecognized names.
    if (_filename.endsWith(".dxf", Qt::CaseInsensitive))
        save_dxf(paths, bounds.xmin, bounds.ymin, bounds.xmax,
                 bounds.ymax, _filename.toStdString().c_str());
    else
        save_svg(paths, bounds.xmin, bounds.ymin, bounds.xmax,
                 bounds.ymax, _filename.toStdString().c_str());
}

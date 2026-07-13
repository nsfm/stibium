#include <Python.h>

#include <cmath>

#include "ui_resolution_dialog.h"
#include "dialog/resolution.h"

ResolutionDialog::ResolutionDialog(Bounds bounds, bool dimensions, bool has_units,
                                   long max_voxels, QWidget* parent,
                                   float min_resolution)
    : QDialog(parent), bounds(bounds), ui(new Ui::ResolutionDialog),
      z_bounded(!std::isinf(bounds.zmax) && !std::isinf(bounds.zmin))
{
    ui->setupUi(this);

    if (!has_units)
    {
        ui->units->hide();
        ui->unit_label->hide();
    }

    if (dimensions == RESOLUTION_DIALOG_2D)
    {
        ui->detect_features->hide();
        ui->simplify_mesh->hide();
        ui->max_deviation_label->hide();
        ui->max_deviation->hide();
    }

    connect(ui->simplify_mesh, &QCheckBox::toggled,
            ui->max_deviation, &QDoubleSpinBox::setEnabled);
    connect(ui->simplify_mesh, &QCheckBox::toggled,
            ui->max_deviation_label, &QLabel::setEnabled);

    // Re-do the layout, since things may have just been hidden
    layout()->invalidate();
    adjustSize();

    {   // Real model dimensions (voxel counts were never useful)
        const float w = bounds.xmax - bounds.xmin;
        const float h = bounds.ymax - bounds.ymin;
        if (z_bounded)
            ui->export_size->setText(QString("model  %1 \u00d7 %2 \u00d7 %3 mm")
                    .arg(w, 0, 'g', 4).arg(h, 0, 'g', 4)
                    .arg(bounds.zmax - bounds.zmin, 0, 'g', 4));
        else
            ui->export_size->setText(QString("model  %1 \u00d7 %2 mm")
                    .arg(w, 0, 'g', 4).arg(h, 0, 'g', 4));
    }

    connect(ui->max_deviation,
            static_cast<void (QDoubleSpinBox::*)(double)>(
                &QDoubleSpinBox::valueChanged),
            [this](double){
                if (!updating_deviation)
                    deviation_touched = true;
            });

    // This connection is awkward because of function overloading.
    connect(ui->export_res,
            static_cast<void (QSpinBox::*)(int)>(&QSpinBox::valueChanged),
            this, &ResolutionDialog::onValueChanged);
    onValueChanged(ui->export_res->value());

    if (dimensions == RESOLUTION_DIALOG_2D || !z_bounded)
    {
        float area = (bounds.xmax - bounds.xmin) *
                     (bounds.ymax - bounds.ymin);
        ui->export_res->setValue(std::max(double(min_resolution),
                pow(max_voxels / area, 1/2.) / 4.));
    }
    else
    {
        float volume = (bounds.xmax - bounds.xmin) *
                       (bounds.ymax - bounds.ymin) *
                       (bounds.zmax - bounds.zmin);
        ui->export_res->setValue(std::max(double(min_resolution),
                pow(max_voxels / volume, 1/3.) / 2.52));
    }
}

void ResolutionDialog::onValueChanged(int i)
{
    ui->export_size->setText(QString("%1 x %2 x %3")
            .arg(int((bounds.xmax - bounds.xmin) * i))
            .arg(int((bounds.ymax - bounds.ymin) * i))
            .arg(z_bounded
                    ? int((bounds.zmax - bounds.zmin) * i)
                    : 1));

    // Sampling can miss details thinner than ~2 voxels; surface the
    // implied minimum feature size so the tradeoff is visible
    const double feature = 2.0 / i;
    ui->feature_size->setText(feature < 0.1
            ? QString("min feature \u2248 %1 \u00b5m  \u24d8")
                    .arg(feature * 1000, 0, 'g', 3)
            : QString("min feature \u2248 %1 mm  \u24d8")
                    .arg(feature, 0, 'g', 3));

    // Simplification deviation follows the sampling error (half the
    // feature size = one voxel) until the user edits it themselves,
    // making the resolution the gold standard down the whole chain
    if (!deviation_touched)
    {
        updating_deviation = true;
        ui->max_deviation->setValue(1.0 / i);
        updating_deviation = false;
    }
}

float ResolutionDialog::getResolution() const
{
    return ui->export_res->value();
}

float ResolutionDialog::getMMperUnit() const
{
    QString u = ui->units->currentText();

    if (u == "mm")              return 1;
    else if (u == "cm")         return 10;
    else if (u == "inches")     return 25.4;

    return 1;
}

bool ResolutionDialog::getDetectFeatures() const
{
    return ui->detect_features->isChecked();
}

float ResolutionDialog::getSimplifyDeviation() const
{
    return ui->simplify_mesh->isChecked() ? ui->max_deviation->value() : 0;
}

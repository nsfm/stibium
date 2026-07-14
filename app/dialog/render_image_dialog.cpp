#include <QCheckBox>
#include <QDialogButtonBox>
#include <QFormLayout>
#include <QSpinBox>

#include "dialog/render_image_dialog.h"

RenderImageDialog::RenderImageDialog(bool section_active, QWidget* parent)
    : QDialog(parent)
{
    setWindowTitle("Export image");

    auto layout = new QFormLayout(this);

    size_box = new QSpinBox(this);
    size_box->setRange(64, 16384);
    size_box->setValue(2048);
    size_box->setSuffix(" px");
    size_box->setToolTip("Longest side of the exported image");
    layout->addRow("Image size", size_box);

    aa_box = new QSpinBox(this);
    aa_box->setRange(1, 4);
    aa_box->setValue(2);
    aa_box->setPrefix("x");
    aa_box->setToolTip("Supersampling factor: renders at Nx the size "
                       "and smooth-downscales (1 disables)");
    layout->addRow("Antialiasing", aa_box);

    transparent_box = new QCheckBox(this);
    transparent_box->setChecked(true);
    layout->addRow("Transparent background", transparent_box);

    section_box = new QCheckBox(this);
    section_box->setChecked(section_active);
    section_box->setEnabled(section_active);
    section_box->setToolTip(section_active
            ? "Keep the viewport's section cut in the export"
            : "No section cut is active in the viewport");
    layout->addRow("Apply section cut", section_box);

    auto buttons = new QDialogButtonBox(
            QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
    connect(buttons, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
    layout->addRow(buttons);
}

int RenderImageDialog::imageSize() const
{
    return size_box->value();
}

int RenderImageDialog::antialias() const
{
    return aa_box->value();
}

bool RenderImageDialog::transparentBackground() const
{
    return transparent_box->isChecked();
}

bool RenderImageDialog::applySection() const
{
    return section_box->isChecked();
}

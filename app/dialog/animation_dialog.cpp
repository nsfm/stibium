#include <QDialogButtonBox>
#include <QFormLayout>
#include <QSpinBox>

#include "dialog/animation_dialog.h"

AnimationDialog::AnimationDialog(const QString& title, bool show_frames,
                                 QWidget* parent)
    : QDialog(parent), frames_box(NULL)
{
    setWindowTitle(title);

    auto layout = new QFormLayout(this);

    size_box = new QSpinBox(this);
    size_box->setRange(64, 8192);
    size_box->setValue(512);
    size_box->setSuffix(" px");
    size_box->setToolTip("Longest side of each frame");
    layout->addRow("Image size", size_box);

    if (show_frames)
    {
        frames_box = new QSpinBox(this);
        frames_box->setRange(4, 360);
        frames_box->setValue(36);
        frames_box->setToolTip("Frames per revolution");
        layout->addRow("Frames", frames_box);
    }

    aa_box = new QSpinBox(this);
    aa_box->setRange(1, 4);
    aa_box->setValue(2);
    aa_box->setPrefix("x");
    aa_box->setToolTip("Supersampling factor: renders at Nx the size "
                       "and smooth-downscales (1 disables)");
    layout->addRow("Antialiasing", aa_box);

    auto buttons = new QDialogButtonBox(
            QDialogButtonBox::Ok | QDialogButtonBox::Cancel, this);
    connect(buttons, &QDialogButtonBox::accepted, this, &QDialog::accept);
    connect(buttons, &QDialogButtonBox::rejected, this, &QDialog::reject);
    layout->addRow(buttons);
}

int AnimationDialog::imageSize() const
{
    return size_box->value();
}

int AnimationDialog::frameCount() const
{
    return frames_box ? frames_box->value() : 2;
}

int AnimationDialog::antialias() const
{
    return aa_box->value();
}

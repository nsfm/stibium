#pragma once

#include <QDialog>

class QSpinBox;
class QCheckBox;

/*
 *  Options dialog for View > Render > Export image
 */
class RenderImageDialog : public QDialog
{
public:
    explicit RenderImageDialog(bool section_active, QWidget* parent=0);

    int imageSize() const;
    int antialias() const;
    bool transparentBackground() const;
    bool applySection() const;

protected:
    QSpinBox* size_box;
    QSpinBox* aa_box;
    QCheckBox* transparent_box;
    QCheckBox* section_box;
};

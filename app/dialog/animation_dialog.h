#pragma once

#include <QDialog>

class QSpinBox;

/*
 *  Options dialog for the Render menu's animation exports
 *  (turntable / light sweep / wigglegram / stereo pair).
 */
class AnimationDialog : public QDialog
{
public:
    /*  show_frames: turntable and light sweep take a frame count;
     *  wigglegrams and stereo pairs are fixed two-frame renders. */
    AnimationDialog(const QString& title, bool show_frames,
                    QWidget* parent=0);

    int imageSize() const;
    int frameCount() const;
    int antialias() const;

protected:
    QSpinBox* size_box;
    QSpinBox* frames_box;
    QSpinBox* aa_box;
};

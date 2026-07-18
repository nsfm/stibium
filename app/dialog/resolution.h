#pragma once

#include <QDialog>

#include "fab/types/bounds.h"

namespace Ui {
class ResolutionDialog;
}

#define RESOLUTION_DIALOG_2D 0
#define RESOLUTION_DIALOG_3D 1

#define HAS_UNITS 1
#define UNITLESS 0

/*
 *  General-purpose dialog for getting 2D and 3D resolution
 */
class ResolutionDialog : public QDialog
{
    Q_OBJECT
public:
    /*
     *  Construct a new resolution dialog.
     *
     *  b is the export bounds (used for printing pixel / voxel counts)
     *  dimensions is RESOLUTION_DIALOG_2D or RESOLUTION_DIALOG_3D
     *  has_units is HAS_UNITS or UNITLESS
     *  max_voxels is the maximum allowable voxel count
     */
    explicit ResolutionDialog(Bounds b, bool dimensions, bool has_units,
                              long max_voxels=(1<<22), QWidget* parent=0,
                     float min_resolution=1);

    /*
     *  Returns the current resolution (from the UI)
     */
    float getResolution() const;

    /*
     *  Returns the current mm/unit setting
     */
    float getMMperUnit() const;

    /*
     *  Returns the state of the "Detect Features (experimental)" checkbox
     */
    bool getDetectFeatures() const;

    /*
     *  Returns the maximum deviation (in model units) for mesh
     *  simplification, or 0 if simplification is disabled
     */
    float getSimplifyDeviation() const;

    /*  Stibnite integration (2026-07-18): mesher choice and the
     *  curated advanced knobs.  Only meaningful for 3D dialogs;
     *  2D callers never see the widgets.  */
    enum Mesher { MESHER_STIBNITE = 0, MESHER_CLASSIC = 1 };
    int getMesher() const;
    bool getAutodense() const;
    int getDensityCap() const;      // 2 (standard) or 3 (high)
    bool getDecimate() const;
    bool getSnap() const;
    int getStallPatience() const;   // 1 (fast) or 2 (patient)

protected slots:
    /*
     *  When a value changes, update pixel / voxel count
     */
    void onValueChanged(int i);

    /*
     *  Mesher / quality selection changed: swap the dialog between
     *  Stibnite mode (quality presets, advanced group) and Classic
     *  mode (voxels-per-unit spin, feature detection).
     */
    void onModeChanged();

protected:
    Bounds bounds;
    Ui::ResolutionDialog* ui;
    bool deviation_touched=false;
    bool updating_deviation=false;
    const bool z_bounded;
    bool is_3d=false;
};

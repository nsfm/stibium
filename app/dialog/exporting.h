#pragma once

#include <cstdint>

#include <QDialog>

namespace Ui {
class ExportingDialog;
}

/*
 *  Dialog used by the export system
 */
class ExportingDialog : public QDialog
{
public:
    explicit ExportingDialog(QWidget* parent=0);

    /*
     *  Updates the progress bar.  A total of 0 means "unknown" and
     *  leaves the bar in its indeterminate (busy) state.
     */
    void setProgress(uint64_t done, uint64_t total);

    /*
     *  Updates the status line ("Meshing...", "Simplifying...", ...).
     */
    void setStatus(const QString& status);

    /*
     *  Updates the right-aligned, fixed-width remaining-time label
     *  (empty string clears it).  Separate from the status text so
     *  the countdown the user is staring at never shifts position.
     */
    void setEta(const QString& eta);

protected:
    Ui::ExportingDialog* ui;
};

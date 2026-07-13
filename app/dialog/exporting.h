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

protected:
    Ui::ExportingDialog* ui;
};

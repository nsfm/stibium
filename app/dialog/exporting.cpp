// cmake doesn't find ui_exporting_dialog without this dummy line
#include "ui_exporting_dialog.h"
#include "dialog/exporting.h"

ExportingDialog::ExportingDialog(QWidget* parent)
    : QDialog(parent, Qt::CustomizeWindowHint|Qt::WindowTitleHint),
      ui(new Ui::ExportingDialog)
{
    ui->setupUi(this);
}

void ExportingDialog::setProgress(uint64_t done, uint64_t total)
{
    // A total of 0 means "working, duration unknown": show Qt's busy
    // animation rather than a bar pinned at its last value (phases
    // like mesh decimation can't report granular progress).
    if (total == 0)
    {
        if (ui->progressBar->maximum() != 0)
        {
            ui->progressBar->setMaximum(0);
            ui->progressBar->setValue(-1);
        }
        return;
    }

    if (ui->progressBar->maximum() == 0)
        ui->progressBar->setMaximum(1000);

    ui->progressBar->setValue(int((done * 1000) / total));
}

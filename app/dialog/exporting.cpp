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
    if (total == 0)
        return;

    // Leave the .ui default (maximum == 0, Qt's busy indicator) until
    // the first real progress report arrives.
    if (ui->progressBar->maximum() == 0)
        ui->progressBar->setMaximum(1000);

    ui->progressBar->setValue(int((done * 1000) / total));
}

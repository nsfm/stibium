#include <Python.h>

#include <QCoreApplication>
#include <QtConcurrent>
#include <QFuture>
#include <QFutureWatcher>
#include <QMessageBox>
#include <QDialog>
#include <QTimer>

#include "export/export_worker.h"
#include "dialog/exporting.h"

void ExportWorker::runAsync()
{
    auto exporting_dialog = new ExportingDialog();
    progress_done = 0;
    progress_total = 0;
    progress_phase = PHASE_DEFAULT;

    auto future = QtConcurrent::run(&ExportWorker::async, this);
    QFutureWatcher<void> watcher;
    watcher.setFuture(future);

    connect(&watcher, &decltype(watcher)::finished,
            exporting_dialog, &ExportingDialog::accept);

    // Poll the worker's progress counters into the dialog's bar
    QTimer timer;
    connect(&timer, &QTimer::timeout, exporting_dialog, [=, this]{
        exporting_dialog->setStatus(phaseLabel(progress_phase.load()));
        exporting_dialog->setProgress(progress_done.load(),
                                      progress_total.load());
    });
    timer.start(50);

    // Run until the dialog closes, either because it was accepeted
    // (which indicates that the future has finished) or cancelled
    // (which indicates that someone his escape).
    //
    // If the dialog was cancelled, set the halt flag and wait for the
    // thread to finish (processing events all the while).
    if (exporting_dialog->exec() == QDialog::Rejected)
    {
        halt = 1;
        while (future.isRunning())
            QCoreApplication::processEvents();
    }

    // Reset halt so that we don't halt immediately on future runs
    halt = 0;
    delete exporting_dialog;
}

QString ExportWorker::phaseLabel(int phase)
{
    switch (phase)
    {
        case PHASE_MESHING:     return "Meshing...";
        case PHASE_SIMPLIFYING: return "Simplifying...";
        case PHASE_WRITING:     return "Writing...";
        case PHASE_CONTOURING:  return "Contouring...";
        default:                return "Exporting...";
    }
}

bool ExportWorker::checkWritable() const
{
    if (!QFileInfo(QFileInfo(_filename).path()).isWritable())
    {
        QMessageBox::critical(NULL, "Export error",
                "<b>Export error:</b><br>"
                "Target file is not writable.");
        return false;
    }
    return true;
}

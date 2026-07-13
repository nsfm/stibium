#pragma once

#include <Python.h>
#include <QObject>

#include <atomic>
#include <cstdint>

#include "fab/types/shape.h"
#include "fab/types/bounds.h"

/*
 *  Abstract base class used to export files from a Shape
 */
class ExportWorker : public QObject
{
    Q_OBJECT
public:
    explicit ExportWorker(Shape s, Bounds b, QString f, float r)
        : shape(s), bounds(b), filename(f), resolution(r), halt(0) {}

    /*
     *  Synchronous part of export (e.g. getting resolution from dialog)
     */
    virtual void run()=0;

    /*
     *  Asynchronous part of export (e.g. actually meshing the model)
     */
    virtual void async()=0;

    /*
     *  Handles async running and stopping
     */
    void runAsync();

    /*
     *  Checks if _filename is writable.
     *  If so, returns true; otherwise, shows a warning and returns false.
     */
    bool checkWritable() const;

protected:
    /*
     *  These are constants for the export worker and set at call-time
     */
    const Shape shape;
    const Bounds bounds;
    const QString filename;
    const float resolution;

    /*
     *  These are variables that are set by dialogs (and don't persist
     *  between calls to the worker).
     */
    float _resolution;
    QString _filename;

    /*
     *  Flag used to abort rendering
     */
    volatile int halt;

    /*
     *  Progress reporting: async() sets progress_total to the number
     *  of work units (0 leaves the dialog's bar indeterminate) and
     *  advances progress_done; the GUI thread polls both.
     */
    std::atomic<uint64_t> progress_done{0};
    std::atomic<uint64_t> progress_total{0};

    /*
     *  Which step async() is currently in, shown as the dialog's
     *  status line.
     */
    enum Phase { PHASE_DEFAULT=0, PHASE_MESHING, PHASE_SIMPLIFYING,
                 PHASE_WRITING, PHASE_CONTOURING };
    std::atomic<int> progress_phase{PHASE_DEFAULT};

    static QString phaseLabel(int phase);
};

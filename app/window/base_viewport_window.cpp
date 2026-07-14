#include "app/app.h"
#include "app/settings.h"

#include "viewport/scene.h"
#include "viewport/view.h"

#include <QActionGroup>

#include <QToolButton>
#include <QFileDialog>
#include <QFileInfo>
#include <QFutureWatcher>
#include <QInputDialog>
#include <QMessageBox>
#include <QTimer>
#include <QtConcurrent>

#include <atomic>
#include <memory>

#include "export/export_image.h"
#include "export/gif_writer.h"
#include "dialog/exporting.h"
#include "dialog/render_image_dialog.h"

#include "window/base_viewport_window.h"

BaseViewportWindow::BaseViewportWindow(QList<ViewportView*> views)
    : BaseWindow("Viewport")
{
    ui->menuAdd->deleteLater();
    ui->menuReference->deleteLater();

    // Make heightmap and shaded options mutually exclusive
    QActionGroup* view_actions = new QActionGroup(this);
    view_actions->addAction(ui->actionEnhanced);
    view_actions->addAction(ui->actionShaded);
    view_actions->addAction(ui->actionHeightmap);
    view_actions->setExclusive(true);

    // Accept the global command-line argument '--heightmap'
    // to always open scenes in height-map view.
    if (App::instance()->arguments().contains("--heightmap"))
        ui->actionHeightmap->setChecked(true);

    for (auto view : views)
    {
        connect(ui->actionShaded, &QAction::triggered,
                [=]{ view->scene()->invalidate(); });
        connect(ui->actionEnhanced, &QAction::triggered,
                [=]{ view->scene()->invalidate(); });
        connect(ui->actionHeightmap, &QAction::triggered,
                [=]{ view->scene()->invalidate(); });
        // Checked means visible (the action reads "Bounding boxes")
        connect(ui->actionHideUI, &QAction::toggled,
                [=](bool b){ view->hideUI(!b); });
        connect(ui->actionSectionSlider, &QAction::toggled,
                [=](bool b){ view->setSectionVisible(b); });
        connect(ui->actionLightControl, &QAction::toggled,
                [=](bool b){ view->setLightGizmoVisible(b); });

        // The viewport's eye button and the menu action stay in sync
        // (setChecked only emits toggled on change, so no feedback loop)
        connect(view->hideUIButton(), &QToolButton::toggled,
                ui->actionHideUI, &QAction::setChecked);
        connect(ui->actionHideUI, &QAction::toggled,
                view->hideUIButton(), &QToolButton::setChecked);
    }

    // Render > Analytics overlay: integrate the model, then show the
    // stats card + center-of-mass marker in every pane
    connect(ui->actionAnalytics, &QAction::toggled, this, [=](bool on){
        if (on && !App::instance()->analyticsValid())
        {
            if (!App::instance()->runAnalytics())
            {
                ui->actionAnalytics->setChecked(false);
                return;
            }
        }
        else if (on)
        {
            // Re-toggling refreshes: a check while already computed
            // means the user wants fresh numbers
            App::instance()->runAnalytics();
        }
        for (auto view : views)
            view->setAnalyticsVisible(on);
    });
    for (auto view : views)
        connect(view, &ViewportView::analyticsDismissed,
                [=]{ ui->actionAnalytics->setChecked(false); });

    // Render > Export image: view-matched shaded render to a file
    connect(ui->actionExportImage, &QAction::triggered, this, [=]{
        auto view = views.first();

        RenderImageDialog dialog(view->getSection() < 1, this);
        if (!dialog.exec())
            return;

        const QString filename = QFileDialog::getSaveFileName(
                this, "Export image",
                Settings::get("files/last_export_dir", "").toString(),
                "Images (*.png)");
        if (filename.isEmpty())
            return;
        Settings::set("files/last_export_dir",
                      QFileInfo(filename).absolutePath());

        ImageExport::Options opt;
        opt.M = view->getMatrix(ViewportView::ROT);
        opt.section = dialog.applySection() ? view->getSection() : 1;
        opt.fit_px = dialog.imageSize();
        opt.transparent = dialog.transparentBackground();
        opt.filename = filename.endsWith(".png", Qt::CaseInsensitive)
                ? filename : filename + ".png";

        // Render in the background behind the busy dialog
        auto exporting = new ExportingDialog(this);
        exporting->setStatus("Rendering...");
        auto graph = App::instance()->getGraph();
        auto future = QtConcurrent::run([=]{
            return ImageExport::render(graph, opt);
        });
        QFutureWatcher<QString> watcher;
        watcher.setFuture(future);
        connect(&watcher, &QFutureWatcher<QString>::finished,
                exporting, &QDialog::accept);
        exporting->exec();
        delete exporting;

        const QString err = future.result();
        if (!err.isEmpty())
            QMessageBox::critical(this, "Export error",
                    "<b>Export error:</b><br>" + err);
    });

    /*
     *  Render > Turntable / Wigglegram / Stereo pair.
     *
     *  Shapes are collected on the GUI thread (touching the graph
     *  and Python is main-thread business); the frame loop runs on a
     *  worker with pure Shape copies, reporting progress through an
     *  atomic the busy dialog polls.
     */
    const auto animate = [=](const QString& kind) {
        auto view = views.first();

        // 3D only: a spinning top-down plane is just a static image
        std::vector<Shape> s3, s2;
        const QString collect_err = ImageExport::collectShapeList(
                App::instance()->getGraph(), "", s3, s2);
        if (!collect_err.isEmpty() || s3.empty())
        {
            QMessageBox::critical(this, "Render error",
                    "<b>Render error:</b><br>" +
                    (collect_err.isEmpty()
                        ? "no 3D shapes to animate"
                        : collect_err));
            return;
        }

        int frames = 36;
        if (kind == "turntable")
        {
            bool ok = true;
            frames = QInputDialog::getInt(this, "Turntable",
                    "Frames per revolution:", 36, 4, 360, 1, &ok);
            if (!ok)
                return;
        }

        const bool gif = kind != "stereo";
        const QString filename = QFileDialog::getSaveFileName(
                this, "Render " + kind,
                Settings::get("files/last_export_dir", "").toString(),
                gif ? "GIF (*.gif)" : "Images (*.png)");
        if (filename.isEmpty())
            return;
        Settings::set("files/last_export_dir",
                      QFileInfo(filename).absolutePath());
        const QString suffix = gif ? ".gif" : ".png";
        const QString out = filename.endsWith(suffix, Qt::CaseInsensitive)
                ? filename : filename + suffix;

        // Spin the model about world z, seen from the current view
        const QMatrix4x4 base = view->getMatrix(ViewportView::ROT);
        const auto frame_opt = [base](float yaw) {
            ImageExport::Options opt;
            opt.M = base;
            opt.M.rotate(yaw, 0, 0, 1);
            opt.fit_sphere = true;
            opt.transparent = false;
            return opt;
        };

        auto exporting = new ExportingDialog(this);
        exporting->setStatus(kind == "stereo" ? "Rendering stereo pair..."
                                              : "Rendering frames...");
        auto done = std::make_shared<std::atomic<int>>(0);
        const int total = kind == "turntable" ? frames : 2;

        auto future = QtConcurrent::run([=]() -> QString {
            QString err;
            if (kind == "stereo")
            {
                const QImage L = ImageExport::renderShapesImage(
                        s3, false, frame_opt(2.5f), &err);
                done->store(1);
                const QImage R = err.isEmpty()
                    ? ImageExport::renderShapesImage(
                            s3, false, frame_opt(-2.5f), &err)
                    : QImage();
                if (!err.isEmpty())
                    return err;
                done->store(2);
                const int gap = 8;
                QImage pair(L.width() + R.width() + gap, L.height(),
                            QImage::Format_ARGB32);
                pair.fill(QColor(28, 26, 24));
                for (int j = 0; j < L.height(); ++j)
                {
                    memcpy(pair.scanLine(j), L.constScanLine(j),
                           size_t(L.width()) * 4);
                    memcpy(pair.scanLine(j) + (L.width() + gap) * 4,
                           R.constScanLine(j), size_t(R.width()) * 4);
                }
                return pair.save(out) ? QString()
                                      : "could not write " + out;
            }

            QList<QImage> seq;
            QList<int> delays;
            if (kind == "wiggle")
            {
                for (const float yaw : {-3.f, 3.f})
                {
                    seq.append(ImageExport::renderShapesImage(
                            s3, false, frame_opt(yaw), &err));
                    if (!err.isEmpty())
                        return err;
                    done->fetch_add(1);
                }
                delays = {12, 12};
            }
            else
            {
                for (int i = 0; i < frames; ++i)
                {
                    seq.append(ImageExport::renderShapesImage(
                            s3, false, frame_opt(360.f * i / frames),
                            &err));
                    if (!err.isEmpty())
                        return err;
                    done->fetch_add(1);
                }
                delays = {std::max(2, 400 / frames)};
            }
            return save_gif(seq, delays, out)
                    ? QString() : "could not write " + out;
        });

        auto poll = new QTimer(exporting);
        connect(poll, &QTimer::timeout, exporting, [=]{
            exporting->setProgress(done->load(), total);
        });
        poll->start(50);

        QFutureWatcher<QString> watcher;
        watcher.setFuture(future);
        connect(&watcher, &QFutureWatcher<QString>::finished,
                exporting, &QDialog::accept);
        exporting->exec();
        delete exporting;

        const QString err = future.result();
        if (!err.isEmpty())
            QMessageBox::critical(this, "Render error",
                    "<b>Render error:</b><br>" + err);
    };
    connect(ui->actionTurntable, &QAction::triggered, this,
            [=]{ animate("turntable"); });
    connect(ui->actionWiggle, &QAction::triggered, this,
            [=]{ animate("wiggle"); });
    connect(ui->actionStereo, &QAction::triggered, this,
            [=]{ animate("stereo"); });

    show();
}

BaseViewportWindow::ShadingMode BaseViewportWindow::shadingMode() const
{
    if (ui->actionHeightmap->isChecked())
        return SHADE_HEIGHTMAP;
    else if (ui->actionShaded->isChecked())
        return SHADE_SHADED;
    return SHADE_ENHANCED;
}

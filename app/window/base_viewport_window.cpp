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
#include <cmath>
#include <memory>

#include "export/export_image.h"
#include "export/gif_writer.h"
#include "dialog/exporting.h"
#include "dialog/animation_dialog.h"
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
        opt.supersample = dialog.antialias();
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

        const bool sweep = kind == "turntable" || kind == "lightsweep";
        AnimationDialog dialog(
                kind == "turntable" ? "Turntable"
                : kind == "lightsweep" ? "Light sweep"
                : kind == "wiggle" ? "Wigglegram" : "Stereo pair",
                sweep, this);
        if (!dialog.exec())
            return;
        const int frames = dialog.frameCount();
        const int img_size = dialog.imageSize();
        const int aa = dialog.antialias();

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
        const auto frame_opt = [base, img_size, aa](float yaw,
                                                    float light_az) {
            ImageExport::Options opt;
            opt.M = base;
            opt.M.rotate(yaw, 0, 0, 1);
            opt.fit_px = img_size;
            opt.supersample = aa;
            opt.fit_sphere = true;
            opt.transparent = false;
            if (light_az == light_az)   // NaN = configured light
            {
                const float el = 0.6f;  // ~35 degrees elevation
                opt.light_override = true;
                opt.light[0] = cosf(light_az) * cosf(el);
                opt.light[1] = sinf(light_az) * cosf(el);
                opt.light[2] = sinf(el);
            }
            return opt;
        };
        const float KEEP = NAN;

        auto exporting = new ExportingDialog(this);
        exporting->setStatus(kind == "stereo" ? "Rendering stereo pair..."
                                              : "Rendering frames...");
        auto done = std::make_shared<std::atomic<int>>(0);
        const int total = sweep ? frames : 2;

        auto future = QtConcurrent::run([=]() -> QString {
            QString err;
            if (kind == "stereo")
            {
                const QImage L = ImageExport::renderShapesImage(
                        s3, false, frame_opt(2.5f, KEEP), &err);
                done->store(1);
                const QImage R = err.isEmpty()
                    ? ImageExport::renderShapesImage(
                            s3, false, frame_opt(-2.5f, KEEP), &err)
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
                            s3, false, frame_opt(yaw, KEEP), &err));
                    if (!err.isEmpty())
                        return err;
                    done->fetch_add(1);
                }
                delays = {12, 12};
            }
            else
            {
                const auto sweep_opt = [&](int i) {
                    return kind == "lightsweep"
                        ? frame_opt(0, 6.2831853f * i / frames)
                        : frame_opt(360.f * i / frames, KEEP);
                };
                for (int i = 0; i < frames; ++i)
                {
                    seq.append(ImageExport::renderShapesImage(
                            s3, false, sweep_opt(i),
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
    connect(ui->actionLightSweep, &QAction::triggered, this,
            [=]{ animate("lightsweep"); });
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

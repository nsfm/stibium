#include "app/app.h"

#include "viewport/scene.h"
#include "viewport/view.h"

#include <QActionGroup>

#include <QToolButton>
#include <QFileDialog>
#include <QFutureWatcher>
#include <QMessageBox>
#include <QtConcurrent>

#include "export/export_image.h"
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

    // Render > Export image: view-matched shaded render to a file
    connect(ui->actionExportImage, &QAction::triggered, this, [=]{
        auto view = views.first();

        RenderImageDialog dialog(view->getSection() < 1, this);
        if (!dialog.exec())
            return;

        const QString filename = QFileDialog::getSaveFileName(
                this, "Export image", "", "Images (*.png)");
        if (filename.isEmpty())
            return;

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

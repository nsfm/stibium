#include "app/app.h"

#include "viewport/scene.h"
#include "viewport/view.h"

#include <QActionGroup>

#include <QToolButton>

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
        connect(ui->actionHideUI, &QAction::toggled,
                [=](bool b){ view->hideUI(b); });
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

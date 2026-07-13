#include <Python.h>

#include <QFontDatabase>

#include <QOpenGLWidget>
#include <QMatrix4x4>
#include <QMouseEvent>
#include <QMenu>
#include <QPropertyAnimation>

#include <QSlider>
#include <QToolButton>
#include <QToolTip>

#include "viewport/view.h"
#include "app/app.h"
#include "viewport/scene.h"
#include "viewport/image.h"
#include "viewport/control/control.h"
#include "viewport/control/control_instance.h"

#include "app/colors.h"
#include "app/settings.h"

#include "graph/proxy/datum.h"
#include "graph/constructor/populate.h"

ViewportView::ViewportView(QWidget* parent, ViewportScene* scene)
    : QGraphicsView(new QGraphicsScene(), parent), gl(new QOpenGLWidget(this)),
      scale(100), pitch(0), yaw(0), view_scene(scene)
{
    setStyleSheet("QGraphicsView { border-style: none; }");
    setRenderHints(QPainter::Antialiasing);

    setViewport(gl.context);

    setSceneRect(-width()/2, -height()/2, width(), height());
    setMouseTracking(true);

    QAbstractScrollArea::setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    QAbstractScrollArea::setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

    {   // Section slider (left edge, on by default)
        section_slider = new QSlider(Qt::Vertical, this);
        section_slider->setRange(1, 100);
        section_slider->setValue(100);
        section_slider->setToolTip(
                "Section view: cuts the model on a screen-parallel plane.\n"
                "Rotate the view to aim the cut; slide to set its depth.");
        connect(section_slider, &QSlider::valueChanged,
                [this](int v){ this->setSection(v / 100.0f); });
    }

    {   // Eye button (top right): show / hide bounding boxes
        hide_ui_button = new QToolButton(this);
        hide_ui_button->setCheckable(true);
        hide_ui_button->setChecked(true);   // lit = boxes visible
        hide_ui_button->setText(QString::fromUtf8("\xf0\x9f\x91\x81"));
        hide_ui_button->setToolTip("Show / hide bounding boxes");
        hide_ui_button->setCursor(Qt::PointingHandCursor);
        hide_ui_button->setStyleSheet(
                "QToolButton {"
                "  border: 1px solid rgba(120, 120, 120, 120);"
                "  border-radius: 12px;"
                "  background: rgba(40, 40, 40, 140);"
                "  padding: 1px 8px;"
                "}"
                "QToolButton:checked {"
                "  background: rgba(200, 150, 60, 160);"
                "}");
    }
}

////////////////////////////////////////////////////////////////////////////////

QMatrix4x4 ViewportView::getMatrix(int params) const
{
    QMatrix4x4 M;

    // Remember that these operations are applied in reverse order.
    if (params & SCALE)
    {
        M.scale(scale, -scale, scale);
    }

    if (params & ROT)
    {
        M.rotate(pitch * 180 / M_PI, QVector3D(1, 0, 0));
        M.rotate(yaw  *  180 / M_PI, QVector3D(0, 0, 1));
    }

    if (params & MOVE)
    {
        M.translate(-center.x(), -center.y(), -center.z());
    }

    return M;
}

QVector3D ViewportView::sceneToWorld(QPointF p) const
{
    return getMatrix().inverted() * QVector3D(p.x(), p.y(), 0);
}

void ViewportView::installImage(DepthImage* d)
{
    connect(this, &ViewportView::getDepth,
            d, &DepthImage::getDepth);
    connect(this, &ViewportView::paintImage,
            d, &DepthImage::paint);
}

////////////////////////////////////////////////////////////////////////////////

void ViewportView::lockAngle(float y, float p)
{
    yaw = y;
    pitch = p;
    angle_locked = true;
}

////////////////////////////////////////////////////////////////////////////////

void ViewportView::drawBackground(QPainter* painter, const QRectF& rect)
{
    QGraphicsView::drawBackground(painter, rect);

    painter->beginNativePainting();
    if (!gl_initialized)
    {
        initializeOpenGLFunctions();
        gl_initialized = true;
    }
    glClear(GL_DEPTH_BUFFER_BIT);
    painter->endNativePainting();

    {   // Vertical gradient backdrop in the warm palette
        QLinearGradient grad(rect.topLeft(), rect.bottomLeft());
        grad.setColorAt(0, QColor("#12100b"));
        grad.setColorAt(1, QColor("#080706"));
        painter->fillRect(rect, grad);
    }

    painter->beginNativePainting();

    // Get bounds from all child images
    float zmin = INFINITY;
    float zmax = -INFINITY;
    auto m = getMatrix();
    emit(getDepth(m, &zmin, &zmax));

    // Paint all images
    emit(paintImage(m, zmin, zmax));

    painter->endNativePainting();
}

void ViewportView::drawAxes(QPainter* painter) const
{
    // First, draw the axes.
    auto m = getMatrix();
    QVector3D o = m * QVector3D(0, 0, 0);
    QVector3D x = m * QVector3D(1, 0, 0);
    QVector3D y = m * QVector3D(0, 1, 0);
    QVector3D z = m * QVector3D(0, 0, 1);

    QList<QPair<QVector3D, QColor>> pts = {
        {x, Colors::red},
        {y, Colors::green},
        {z, Colors::blue}};

    // Sort the axes to fake proper z clipping
    std::sort(pts.begin(), pts.end(),
            [](QPair<QVector3D, QColor> a, QPair<QVector3D, QColor> b)
            { return a.first.z() < b.first.z(); });

    QList<QPair<QVector3D, QString>> labels = {
        {x, "X"}, {y, "Y"}, {z, "Z"}};

    int i = 0;
    for (auto p : pts)
    {
        painter->setPen(QPen(p.second, 2));
        painter->drawLine(o.toPointF(), p.first.toPointF());
    }
    for (auto l : labels)
    {
        const auto c = QList<QColor>({Colors::red, Colors::green,
                                      Colors::blue})[i++];
        painter->setPen(c);
        const auto tip = o.toPointF() +
            (l.first.toPointF() - o.toPointF()) * 1.12;
        painter->drawText(QRectF(tip.x() - 6, tip.y() - 6, 12, 12),
                          Qt::AlignCenter, l.second);
    }
}

void ViewportView::drawCoords(QPainter* painter) const
{
    QPointF mouse_pos = mapToScene(mapFromGlobal(QCursor::pos()));
    if (!sceneRect().contains(mouse_pos))
    {
        return;
    }

    // Get rotate-only transform matrix
    QMatrix4x4 M = getMatrix(ROT);
    const float threshold = 0.98;

    const auto a = M.inverted() * QVector3D(0, 0, 1);

    QList<QPair<char, QVector3D>> axes = {
        {'x', QVector3D(1, 0, 0)},
        {'y', QVector3D(0, 1, 0)},
        {'z', QVector3D(0, 0, 1)}};

    char axis = 0;
    float opacity = 0;
    for (const auto v : axes)
    {
        float dot = fabs(QVector3D::dotProduct(a, v.second));
        if (dot > threshold)
        {
            axis = v.first;
            opacity = (dot - threshold) / (1 - threshold);
        }
    }

    auto p = sceneToWorld(mouse_pos);

    int value = opacity * 200;

    auto readout = Colors::base05;
    readout.setAlpha(value);
    painter->setPen(QPen(readout));
    painter->setFont(QFontDatabase::systemFont(QFontDatabase::FixedFont));
    QString txt;
    if (axis == 'z')
    {
        txt = QString("X: %1\nY: %2").arg(p.x()).arg(p.y());
    }
    else if (axis == 'y')
    {
        txt = QString("X: %1\nZ: %2").arg(p.x()).arg(p.z());
    }
    else if (axis == 'x')
    {
        txt = QString("Y: %1\nZ: %2").arg(p.y()).arg(p.z());
    }

    painter->drawText({-width()/2.0 + 10, -height()/2.0 + 10, 300, 200}, txt);
}

void ViewportView::drawForeground(QPainter* painter, const QRectF& rect)
{
    QGraphicsView::drawForeground(painter, rect);
    drawAxes(painter);
    drawCoords(painter);
    drawLightGizmo(painter);
    drawAnalytics(painter);
}

QPointF ViewportView::lightGizmoCenter() const
{
    return QPointF(width()/2.0 - 58, height()/2.0 - 58);
}

void ViewportView::drawLightGizmo(QPainter* painter) const
{
    if (!light_gizmo_visible)
        return;

    const float R = 34;
    const auto c = lightGizmoCenter();

    painter->setRenderHint(QPainter::Antialiasing);

    auto fill = Colors::base01;
    fill.setAlpha(160);
    painter->setBrush(fill);
    painter->setPen(QPen(Colors::base03, 1));
    painter->drawEllipse(c, R, R);

    // Read the light direction and place the sun dot
    const auto parts = Settings::get(
            "render/key_light", "0.57,-0.57,0.57").toString().split(',');
    QVector3D l(0.57, -0.57, 0.57);
    if (parts.size() == 3)
        l = QVector3D(parts[0].toFloat(), parts[1].toFloat(),
                      parts[2].toFloat());
    l.normalize();
    const QPointF dot = c + QPointF(l.x(), l.y()) * R;

    painter->setPen(QPen(Colors::base03, 1));
    painter->drawLine(c, dot);
    painter->setPen(QPen(Colors::base00, 1));
    painter->setBrush(Colors::amber);
    painter->drawEllipse(dot, 5, 5);
}

void ViewportView::updateLightFromGizmo(QPointF scene_pos)
{
    const float R = 34;
    auto d = (scene_pos - lightGizmoCenter()) / R;
    float x = d.x(), y = d.y();
    const float len = sqrt(x*x + y*y);
    if (len > 0.95)
    {
        x *= 0.95 / len;
        y *= 0.95 / len;
    }
    const float z = sqrt(fmax(0.0475, 1 - x*x - y*y));

    Settings::set("render/key_light", QString("%1,%2,%3")
            .arg(x, 0, 'f', 3).arg(y, 0, 'f', 3).arg(z, 0, 'f', 3));
    scene()->invalidate();
}

////////////////////////////////////////////////////////////////////////////////

void ViewportView::update()
{
    scene()->invalidate();
    emit(changed(getMatrix(), geometry(), section));

    emit(centerChanged(center));
    emit(scaleChanged(scale));
}

////////////////////////////////////////////////////////////////////////////////

void ViewportView::installControl(Control* c)
{
    c->makeInstanceFor(this);
}

void ViewportView::installDatum(BaseDatumProxy* d)
{
    d->addViewport(this);
}

////////////////////////////////////////////////////////////////////////////////

void ViewportView::setScale(float s)
{
    scale = s;

    // Update the view without calling ViewportView::update
    // (because that would lead to infinite recursion)
    scene()->invalidate();
    emit(changed(getMatrix(), geometry(), section));
}

void ViewportView::setCenter(QVector3D c)
{
    center = c;

    // Update the view without calling ViewportView::update
    // (because that would lead to infinite recursion)
    scene()->invalidate();
    emit(changed(getMatrix(), geometry(), section));
}

////////////////////////////////////////////////////////////////////////////////

void ViewportView::mousePressEvent(QMouseEvent* event)
{
    // The analytics card's close box swallows its click
    if (show_analytics && event->button() == Qt::LeftButton &&
        App::instance()->analyticsValid() &&
        analytics_close.contains(event->position()))
    {
        emit(analyticsDismissed());
        event->accept();
        return;
    }

    // Clicks on the light gizmo aim the key light instead of the view
    if (event->button() == Qt::LeftButton)
    {
        const auto p = mapToScene(event->position().toPoint());
        if (light_gizmo_visible &&
            QLineF(p, lightGizmoCenter()).length() < 42)
        {
            gizmo_drag = true;
            updateLightFromGizmo(p);
            event->accept();
            return;
        }
    }

    QGraphicsView::mousePressEvent(event);

    // If the event hasn't been accepted, record click position for
    // panning / rotation on mouse drag.
    if (!event->isAccepted())
    {
        if (event->button() == Qt::LeftButton)
        {
            click_pos = mapToScene(event->pos());
            click_pos_world = sceneToWorld(click_pos);
        }
        else
        {
            click_pos = event->pos();
        }
        dragged = false;
    }
}

void ViewportView::mouseMoveEvent(QMouseEvent* event)
{
    if (gizmo_drag)
    {
        updateLightFromGizmo(mapToScene(event->position().toPoint()));
        return;
    }

    QGraphicsView::mouseMoveEvent(event);

    current_pos = event->pos();
    if (scene()->mouseGrabberItem() == NULL)
    {
        if (event->buttons() == Qt::LeftButton)
        {
            center += click_pos_world - sceneToWorld(mapToScene(event->pos()));
            update();
        }
        else if (event->buttons() == Qt::RightButton && !angle_locked)
        {
            QPointF d = click_pos - event->pos();
            pitch = fmin(0, fmax(-M_PI, pitch - 0.01 * d.y()));
            yaw = fmod(yaw + M_PI - 0.01 * d.x(), M_PI*2) - M_PI;

            click_pos = event->pos();
            update();
        }
        dragged = true;
    }

    // Redraw to update cursor position
    scene()->invalidate(QRect(), QGraphicsScene::ForegroundLayer);
}

void ViewportView::mouseReleaseEvent(QMouseEvent* event)
{
    if (gizmo_drag)
    {
        gizmo_drag = false;
        event->accept();
        return;
    }

    if (event->button() == Qt::RightButton && !dragged)
    {
        auto is = items(event->pos());
        if (is.size() > 1)
        {
            openRaiseMenu(is);
        }
        else if (is.size())
        {
            if (auto c = dynamic_cast<ControlInstance*>(is.front()))
            {
                c->openContextMenu();
            }
        }
        else
        {
            openAddMenu(true);
        }
    }
    else
    {
        QGraphicsView::mouseReleaseEvent(event);
    }
}

void ViewportView::wheelEvent(QWheelEvent* event)
{
    QVector3D a = sceneToWorld(mapToScene(event->position().toPoint()));
    scale *= pow(1.001, -event->angleDelta().y());
    QVector3D b = sceneToWorld(mapToScene(event->position().toPoint()));
    center += a - b;
    update();
}

void ViewportView::resizeEvent(QResizeEvent* e)
{
    if (section_slider)
        section_slider->setGeometry(10, 44, 18, height() - 54);
    if (hide_ui_button)
        hide_ui_button->setGeometry(width() - 50, 10, 40, 24);

    Q_UNUSED(e);
    setSceneRect(-width()/2, -height()/2, width(), height());
}

void ViewportView::keyPressEvent(QKeyEvent* event)
{
    QGraphicsView::keyPressEvent(event);
    if (event->isAccepted())
    {
        return;
    }
    if (event->key() == Qt::Key_A &&
        (event->modifiers() & Qt::ShiftModifier))
    {
        openAddMenu();
    }
}

void ViewportView::setSection(float s)
{
    section = s;
    scene()->invalidate();
    emit(changed(getMatrix(), geometry(), section));
}

void ViewportView::setSectionVisible(bool visible)
{
    section_slider->setVisible(visible);
    if (!visible)
        setSection(1);
    else
        setSection(section_slider->value() / 100.0f);
}

void ViewportView::setLightGizmoVisible(bool visible)
{
    light_gizmo_visible = visible;
    scene()->invalidate();
}

void ViewportView::setAnalyticsVisible(bool visible)
{
    show_analytics = visible;
    scene()->invalidate();
}

void ViewportView::drawAnalytics(QPainter* painter) const
{
    if (!show_analytics || !App::instance()->analyticsValid())
        return;

    const auto& s = App::instance()->analyticsStats();
    const bool flat = App::instance()->analyticsFlat();

    // Center-of-mass marker, projected into the scene
    {
        const auto m = getMatrix();
        const QVector3D com = m * QVector3D(s.com[0], s.com[1], s.com[2]);
        const QPointF p(com.x(), com.y());
        painter->setPen(QPen(Colors::amber, 1.5));
        painter->setBrush(Qt::NoBrush);
        painter->drawEllipse(p, 7, 7);
        painter->drawLine(p - QPointF(12, 0), p + QPointF(12, 0));
        painter->drawLine(p - QPointF(0, 12), p + QPointF(0, 12));
    }

    // Stats card, fixed to the viewport's top-left
    painter->save();
    painter->resetTransform();

    QFont font = painter->font();
    font.setPixelSize(11);
    painter->setFont(font);

    QStringList lines;
    lines << QString(flat ? "area    %1" : "volume  %1")
                .arg(s.volume, 0, 'g', 6);
    lines << QString("center  %1, %2, %3")
                .arg(s.com[0], 0, 'g', 4)
                .arg(s.com[1], 0, 'g', 4)
                .arg(s.com[2], 0, 'g', 4);
    lines << QString("size    %1 x %2 x %3")
                .arg(s.tight[3] - s.tight[0], 0, 'g', 4)
                .arg(s.tight[4] - s.tight[1], 0, 'g', 4)
                .arg(s.tight[5] - s.tight[2], 0, 'g', 4);

    const QFontMetricsF fm(font);
    float wmax = 0;
    for (const auto& l : lines)
        wmax = fmax(wmax, fm.horizontalAdvance(l));

    // Top-center, clear of the section slider on the left
    const float card_w = wmax + 38;
    const QRectF card((width() - card_w) / 2, 10, card_w,
                      lines.size() * fm.height() + 14);

    auto bg = Colors::base01;
    bg.setAlphaF(0.85);
    painter->setBrush(bg);
    painter->setPen(QPen(Colors::amber, 1));
    painter->drawRoundedRect(card, 5, 5);

    painter->setPen(Colors::base06);
    float y = card.top() + 7 + fm.ascent();
    for (const auto& l : lines)
    {
        painter->drawText(QPointF(card.left() + 10, y), l);
        y += fm.height();
    }

    // Close box, top-right of the card
    analytics_close = QRectF(card.right() - 20, card.top() + 4, 16, 16);
    painter->setPen(QPen(Colors::base04, 1.4));
    const auto x = analytics_close.adjusted(4.5, 4.5, -4.5, -4.5);
    painter->drawLine(x.topLeft(), x.bottomRight());
    painter->drawLine(x.topRight(), x.bottomLeft());

    painter->restore();
}

bool ViewportView::viewportEvent(QEvent* e)
{
    if (e->type() == QEvent::ToolTip && light_gizmo_visible)
    {
        auto he = static_cast<QHelpEvent*>(e);
        const auto p = mapToScene(he->pos());
        if (QLineF(p, lightGizmoCenter()).length() < 42)
        {
            QToolTip::showText(he->globalPos(),
                    "Key light: drag the sun to relight the\n"
                    "Enhanced render mode", this);
            return true;
        }
    }
    return QGraphicsView::viewportEvent(e);
}

////////////////////////////////////////////////////////////////////////////////

void ViewportView::openAddMenu(bool view_commands)
{
    QMenu* m = new QMenu(this);

    if (view_commands && !angle_locked)
    {
        auto sub = new QMenu("View");

        connect(sub->addAction("Top"), &QAction::triggered,
                [=]{ this->spinTo(0, 0); });
        connect(sub->addAction("Bottom"), &QAction::triggered,
                [=]{ this->spinTo(0, -M_PI); });
        connect(sub->addAction("Left"), &QAction::triggered,
                [=]{ this->spinTo(M_PI/2, -M_PI/2); });
        connect(sub->addAction("Right"), &QAction::triggered,
                [=]{ this->spinTo(-M_PI/2, -M_PI/2); });
        connect(sub->addAction("Front"), &QAction::triggered,
                [=]{ this->spinTo(0, -M_PI/2); });
        connect(sub->addAction("Back"), &QAction::triggered,
                [=]{ this->spinTo(-M_PI, -M_PI/2); });

        m->addMenu(sub);
        m->addSeparator();
    }

    populateNodeMenu(m, view_scene->getGraph());

    m->exec(QCursor::pos());
    m->deleteLater();
}

void ViewportView::openRaiseMenu(QList<QGraphicsItem*> items)
{
    QScopedPointer<QMenu> m(new QMenu(this));

    int found = 0;
    for (auto i : items)
    {
        if (auto c = dynamic_cast<ControlInstance*>(i))
        {
            auto a = new QAction(c->getName(), m.data());
            a->setData(QVariant::fromValue(c));
            m->addAction(a);
            found++;
        }
    }

    QAction* chosen = (found > 1) ? m->exec(QCursor::pos())
                                  : NULL;

    if (chosen)
    {
        if (raised)
        {
            raised->setZValue(0);
        }
        raised = chosen->data().value<ControlInstance*>();
        raised->setZValue(0.1);
    }
}

////////////////////////////////////////////////////////////////////////////////

void ViewportView::setYaw(float y)
{
    yaw = y;
    update();
}

void ViewportView::setPitch(float p)
{
    pitch = p;
    update();
}

void ViewportView::spinTo(float new_yaw, float new_pitch)
{
    QPropertyAnimation* a = new QPropertyAnimation(this, "_yaw", this);
    a->setDuration(100);
    a->setStartValue(yaw);
    a->setEndValue(new_yaw);

    QPropertyAnimation* b = new QPropertyAnimation(this, "_pitch", this);
    b->setDuration(100);
    b->setStartValue(pitch);
    b->setEndValue(new_pitch);

    a->start(QPropertyAnimation::DeleteWhenStopped);
    b->start(QPropertyAnimation::DeleteWhenStopped);
}

////////////////////////////////////////////////////////////////////////////////

void ViewportView::zoomTo(Node* n)
{
    // Find all ControlInstances that are declared by this node
    QList<ControlInstance*> instances;
    for (auto i : items())
    {
        if (auto c = dynamic_cast<ControlInstance*>(i))
        {
            if (c->getNode() == n)
            {
                instances.push_back(c);
            }
        }
    }

    // Find a weighted sum of central points
    QVector3D pos;
    float area_sum = 0;
    for (auto i : instances)
    {
        const float area = i->boundingRect().width() *
                           i->boundingRect().height();
        pos += i->getControl()->pos() * area;
        area_sum += area;
    }
    pos /= area_sum;

    auto a = new QPropertyAnimation(this, "_center");
    a->setDuration(100);
    a->setStartValue(center);
    a->setEndValue(pos);

    a->start(QPropertyAnimation::DeleteWhenStopped);
}

void ViewportView::hideUI(bool b)
{
    ui_hidden = b;
    for (auto i : scene()->items())
    {
        if (auto c = dynamic_cast<ControlInstance*>(i))
        {
            if (b)
            {
                c->hide();
            }
            else
            {
                c->show();
            }
        }
    }
}

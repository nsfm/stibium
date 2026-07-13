#include <Python.h>
#include <QPainter>

#include "canvas/scene.h"
#include "canvas/connection/base.h"
#include <QStyleOptionGraphicsItem>
#include <QVariantAnimation>

#include "app/colors.h"

BaseConnection::BaseConnection(QColor color)
    : base_color(color), hover(false)
{
    setAcceptHoverEvents(true);
    setZValue(1);
}

QRectF BaseConnection::boundingRect() const
{
    QPainterPathStroker s;
    s.setWidth(20);
    return s.createStroke(path()).boundingRect();
}

QColor BaseConnection::color() const
{
    return isSelected() ? Colors::highlight(base_color) : base_color;
}

void BaseConnection::paint(QPainter *painter,
                           const QStyleOptionGraphicsItem *option,
                           QWidget *widget)
{
    Q_UNUSED(widget);

    painter->setRenderHint(QPainter::Antialiasing);

    const float lod = option->levelOfDetailFromTransform(
            painter->worldTransform());

    // Zoomed way out: constant screen-width stroke (never sub-pixel)
    if (lod < CANVAS_LOD_THRESHOLD)
    {
        QPen pen(color(), 1.5, Qt::SolidLine, Qt::RoundCap);
        pen.setCosmetic(true);
        painter->setPen(pen);
        painter->drawPath(path());
        return;
    }

    // Fade the expensive dressing in over a band above the threshold,
    // so the transition reads as gradual rather than a snap.
    const float fade = fmin(1.0f,
            (lod - CANVAS_LOD_THRESHOLD) / 0.13f);

    // Soft accent glow when hovered or selected, eased in and out
    const qreal glow_t = isSelected() ? 1.0 : hover_glow;
    if (glow_t > 0.01)
    {
        QColor glow = Colors::amber;
        glow.setAlpha(40 * fade * glow_t);
        painter->setPen(QPen(glow, 11, Qt::SolidLine, Qt::RoundCap));
        painter->drawPath(path(true));
        glow.setAlpha(70 * fade * glow_t);
        painter->setPen(QPen(glow, 7, Qt::SolidLine, Qt::RoundCap));
        painter->drawPath(path(true));
    }

    // Drop shadow pass
    painter->save();
    painter->translate(0, 1.5);
    painter->setPen(QPen(QColor(0, 0, 0, int(90 * fade)), 4,
                         Qt::SolidLine, Qt::RoundCap));
    painter->drawPath(path());
    painter->restore();

    // Main stroke: subtle gradient along the wire to imply flow
    QLinearGradient grad(startPos(), endPos());
    grad.setColorAt(0, Colors::dim(color()));
    grad.setColorAt(1, Colors::highlight(color()));
    painter->setPen(QPen(QBrush(grad), 1.5 + 1.5 * fade,
                         Qt::SolidLine, Qt::RoundCap));
    painter->drawPath(path());
}

QPainterPath BaseConnection::path(bool only_bezier) const
{
    QPointF start = startPos();
    QPointF end = endPos();

    float length = 50;
    if (end.x() <= start.x())
    {
        length += (start.x() - end.x()) / 2;
    }

    QPainterPath p;
    p.moveTo(start);
    if (only_bezier)
        p.moveTo(start + QPointF(15, 0));
    else
        p.lineTo(start + QPointF(15, 0));

    p.cubicTo(QPointF(start.x() + length, start.y()),
              QPointF(end.x() - length, end.y()),
              QPointF(end.x() - 15, end.y()));

    if (!only_bezier)
        p.lineTo(end);

   return p;
}

QPainterPath BaseConnection::shape() const
{
    QPainterPathStroker s;
    s.setWidth(20);
    return s.createStroke(path(true));
}

void BaseConnection::animateGlow(qreal target)
{
    auto anim = new QVariantAnimation(this);
    anim->setDuration(110);
    anim->setStartValue(hover_glow);
    anim->setEndValue(target);
    connect(anim, &QVariantAnimation::valueChanged,
            [this](const QVariant& v){
                hover_glow = v.toReal();
                update(); });
    anim->start(QAbstractAnimation::DeleteWhenStopped);
}

void BaseConnection::hoverEnterEvent(QGraphicsSceneHoverEvent *event)
{
    Q_UNUSED(event);
    if (!hover)
    {
        hover = true;
        animateGlow(1.0);
    }
}

void BaseConnection::hoverLeaveEvent(QGraphicsSceneHoverEvent *event)
{
    Q_UNUSED(event);
    if (hover)
    {
        hover = false;
        animateGlow(0.0);
    }
}

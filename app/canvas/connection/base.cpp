#include <Python.h>
#include <QPainter>

#include "canvas/connection/base.h"
#include <QStyleOptionGraphicsItem>

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

    // Zoomed way out: single cheap stroke, no shadow or glow
    if (option->levelOfDetailFromTransform(painter->worldTransform()) < 0.4)
    {
        painter->setPen(QPen(color(), 2, Qt::SolidLine, Qt::RoundCap));
        painter->drawPath(path());
        return;
    }

    // Soft accent glow when hovered or selected
    if (hover || isSelected())
    {
        QColor glow = Colors::amber;
        glow.setAlpha(40);
        painter->setPen(QPen(glow, 11, Qt::SolidLine, Qt::RoundCap));
        painter->drawPath(path(true));
        glow.setAlpha(70);
        painter->setPen(QPen(glow, 7, Qt::SolidLine, Qt::RoundCap));
        painter->drawPath(path(true));
    }

    // Drop shadow pass
    painter->save();
    painter->translate(0, 1.5);
    painter->setPen(QPen(QColor(0, 0, 0, 90), 4, Qt::SolidLine, Qt::RoundCap));
    painter->drawPath(path());
    painter->restore();

    // Main stroke: subtle gradient along the wire to imply flow
    QLinearGradient grad(startPos(), endPos());
    grad.setColorAt(0, Colors::dim(color()));
    grad.setColorAt(1, Colors::highlight(color()));
    painter->setPen(QPen(QBrush(grad), 3, Qt::SolidLine, Qt::RoundCap));
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

void BaseConnection::hoverEnterEvent(QGraphicsSceneHoverEvent *event)
{
    Q_UNUSED(event);
    if (!hover)
    {
        hover = true;
        update();
    }
}

void BaseConnection::hoverLeaveEvent(QGraphicsSceneHoverEvent *event)
{
    Q_UNUSED(event);
    if (hover)
    {
        hover = false;
        update();
    }
}

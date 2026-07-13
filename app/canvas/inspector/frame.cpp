#include <Python.h>

#include <algorithm>

#include <QPropertyAnimation>
#include <QVariantAnimation>

#include <QGraphicsScene>
#include <QPainter>
#include <QGraphicsSceneMouseEvent>
#include <QGraphicsSceneContextMenuEvent>
#include <QMenu>

#include "canvas/inspector/frame.h"
#include "canvas/datum_row.h"
#include "canvas/datum_editor.h"
#include "canvas/inspector/title.h"
#include "canvas/inspector/export.h"
#include "canvas/scene.h"

#include "graph/node.h"
#include "graph/datum.h"
#include "fab/fab.h"

#include "app/colors.h"

const float InspectorFrame::PADDING_ROWS = 3;

InspectorFrame::InspectorFrame(Node* node, QGraphicsScene* scene)
    : QGraphicsObject(), node(node), title_row(new InspectorTitle(node, this)),
      export_button(new InspectorExportButton(this)),
      show_hidden(false), dragging(false)
{
    setFlags(QGraphicsItem::ItemIsMovable |
             QGraphicsItem::ItemIsSelectable);
    setAcceptHoverEvents(true);

    scene->addItem(this);
    if (auto cs = dynamic_cast<CanvasScene*>(scene))
        cs->registerInspector(this);
    redoLayout();
}

InspectorFrame::~InspectorFrame()
{
    if (auto cs = dynamic_cast<CanvasScene*>(scene()))
        cs->unregisterInspector(this);
}

QRectF InspectorFrame::tightBoundingRect() const
{
    // In low-detail mode the children are hidden; use the rect
    // captured when the mode was entered.
    if (low_detail)
        return low_detail_rect;

    QRectF b;
    for (auto c : childItems())
    {
        if (c->isVisible())
        {
            b = b.united(c->boundingRect().translated(c->pos()));
        }
    }
    b.setBottom(b.bottom() + PADDING_ROWS);
    return b;
}

QList<DatumRow*> InspectorFrame::visibleRows() const
{
    QList<DatumRow*> rows;
    for (auto c : childItems())
        if (auto row = dynamic_cast<DatumRow*>(c))
        {
            if (show_hidden || !row->shouldBeHidden())
            {
                rows.append(row);
            }
        }

    // Sort datums by row order
    std::sort(rows.begin(), rows.end(),
          [](const DatumRow* a, const DatumRow* b)
          { return a->getIndex() < b->getIndex(); });

    return rows;
}

QRectF InspectorFrame::boundingRect() const
{
    auto r = tightBoundingRect();
    r += {10, 10, 10, 10};
    return r;
}

void InspectorFrame::paint(QPainter *painter,
                           const QStyleOptionGraphicsItem *option,
                           QWidget *widget)
{
    Q_UNUSED(option);
    Q_UNUSED(widget);

    const auto r = tightBoundingRect();

    painter->setRenderHint(QPainter::Antialiasing);

    // Zoomed way out: draw a skeleton card echoing the full node's
    // anatomy - type-tinted title band (slightly enlarged), body in a
    // much darker shade of the same tint, light stub lines where the
    // datum rows live. Shape without detail, react-skeleton style;
    // naming is the floating labels' job out here.
    if (low_detail)
    {
        const float px = painter->worldTransform().m11();  // screen px per unit
        const bool exportish = export_button->getWorker() != nullptr;
        const QColor tint = exportish ? Colors::amber : typeTint();

        // Sub-8px cards: a single flat rect is all anyone can see
        if (r.height() * px < 8)
        {
            painter->setRenderHint(QPainter::Antialiasing, false);
            painter->setPen(Qt::NoPen);
            painter->setBrush(tint);
            painter->drawRect(r);
            if (isSelected())
            {
                painter->setBrush(Qt::NoBrush);
                painter->setPen(QPen(Colors::amber, 2 / px));
                painter->drawRect(r);
            }
            return;
        }

        // (No shadows out here: at 150+ cards they're a per-frame
        // fill pass nobody can see)
        painter->setPen(Qt::NoPen);

        // Body: the tint, pulled way down
        painter->setBrush(Colors::adjust(tint, 1 / 2.6f));
        painter->drawRoundedRect(r, 10, 10);

        // Title band, rounded only on top, a touch taller than life
        const float band_h = fmin(r.height() * 0.5f,
                                  fmax(24.f, r.height() * 0.22f));
        {
            const QRectF band_r(r.left(), r.top(), r.width(), band_h);
            QPainterPath band;
            band.setFillRule(Qt::WindingFill);
            band.addRoundedRect(band_r, 10, 10);
            band.addRect(band_r.adjusted(0, band_h / 2, 0, 0));
            painter->setBrush(tint);
            painter->drawPath(band.simplified());
        }

        // Datum-row stubs: light tint bars on the dark body,
        // skipped when they'd land under ~3 screen pixels
        int rows = 0;
        for (auto d : node->childDatums())
            if (!d->getName().empty() && d->getName().front() != '_')
                rows++;
        const float row_h = 8;
        const float row_gap = 6;
        const float body_top = r.top() + band_h + 10;
        rows = std::min({rows, 6,
                         int((r.bottom() - 8 - body_top) / (row_h + row_gap))});

        // Stubs fade in between 3 and 6 screen pixels of row pitch
        // instead of snapping at a visibility threshold
        const float pitch_px = (row_h + row_gap) * px;
        const float stub_fade = fmin(1.f, fmax(0.f, (pitch_px - 3) / 3));
        if (rows > 0 && stub_fade > 0)
        {
            auto stub = Colors::adjust(tint, 1.6f);
            stub.setAlphaF(0.55 * stub_fade);
            painter->setBrush(stub);
            const float x0 = r.left() + 10;
            float y = body_top;
            for (int i = 0; i < rows; ++i)
            {
                // Vary widths a touch so it reads as content, not bars
                const float w = (r.width() - 20) * (i % 2 ? 0.55f : 0.8f);
                painter->drawRoundedRect(QRectF(x0, y, w, row_h), 3, 3);
                y += row_h + row_gap;
            }
        }

        painter->setBrush(Qt::NoBrush);
        if (exportish && !isSelected())
        {
            // Findable glow for render/mesh/export nodes
            auto halo = Colors::amber;
            halo.setAlphaF(0.5);
            painter->setPen(QPen(halo, 5));
            painter->drawRoundedRect(r.adjusted(-4, -4, 4, 4), 12, 12);
        }
        painter->setPen(isSelected() ? QPen(Colors::amber, 2)
                                     : QPen(Colors::base02, 1));
        painter->drawRoundedRect(r, 10, 10);
        return;
    }

    // Drop shadow: stacked translucent rects (deliberately not a
    // QGraphicsDropShadowEffect, which rasterizes items on every repaint)
    painter->setPen(Qt::NoPen);
    for (int i=0; i < 3; ++i)
    {
        painter->setBrush(QColor(0, 0, 0, 34 - i*10));
        painter->drawRoundedRect(r.translated(0, 1.5 + i), 8 + i, 8 + i);
    }

    // Body with a subtle vertical gradient
    QLinearGradient grad(r.topLeft(), r.bottomLeft());
    grad.setColorAt(0, Colors::base01);
    grad.setColorAt(1, Colors::adjust(Colors::base01, 1/1.2f));
    painter->setBrush(grad);
    painter->drawRoundedRect(r, 8, 8);

    {   // Title band, rounded only at the top corners
        const auto tr = title_row->boundingRect();
        QPainterPath band;
        band.setFillRule(Qt::WindingFill);
        band.addRoundedRect(tr, 8, 8);
        band.addRect(tr.adjusted(0, tr.height()/2, 0, 0));
        painter->setBrush(typeTint());
        painter->drawPath(band.simplified());
    }

    // Soft accent glow when hovered (focus without selection)
    painter->setBrush(Qt::NoBrush);
    if (focus_glow > 0.01 && !isSelected())
    {
        QColor glow = Colors::amber;
        glow.setAlpha(70 * focus_glow);
        painter->setPen(QPen(glow, 4));
        painter->drawRoundedRect(r, 8, 8);
    }

    // Outline: amber when selected, quiet otherwise
    painter->setPen(isSelected() ? QPen(Colors::amber, 1.5)
                                 : QPen(Colors::base02, 1));
    painter->drawRoundedRect(r, 8, 8);
}

////////////////////////////////////////////////////////////////////////////////

QColor InspectorFrame::typeTint() const
{
    QColor type = Colors::base03;
    bool shape_input = false;
    for (auto d : node->childDatums())
    {
        if (d->isOutput())
            type = Colors::getColor(d);
        else if (d->getType() == fab::ShapeType)
            shape_input = true;
    }

    // Shape-consuming, shape-producing nodes are operators (CSG,
    // deforms, transforms); tint them violet so they read differently
    // from shape sources.
    if (shape_input && type == Colors::green)
        type = Colors::violet;

    const auto base = Colors::base02;
    return QColor(base.red()   * 0.72 + type.red()   * 0.28,
                  base.green() * 0.72 + type.green() * 0.28,
                  base.blue()  * 0.72 + type.blue()  * 0.28);
}

void InspectorFrame::setLowDetail(bool low)
{
    if (low == low_detail)
        return;

    if (low)
        low_detail_rect = tightBoundingRect();  // capture while visible

    prepareGeometryChange();
    low_detail = low;

    for (auto c : childItems())
        c->setVisible(!low);
    if (!low)
        redoLayout();   // restores per-datum hidden rules

    // Soften the swap with a quick opacity pop
    auto anim = new QPropertyAnimation(this, "opacity");
    anim->setDuration(130);
    anim->setStartValue(0.5);
    anim->setEndValue(1.0);
    anim->start(QAbstractAnimation::DeleteWhenStopped);

    update();
}

////////////////////////////////////////////////////////////////////////////////

void InspectorFrame::setNameValid(bool valid)
{
    title_row->setNameValid(valid);
}

////////////////////////////////////////////////////////////////////////////////

void InspectorFrame::setTitle(QString title)
{
    title_row->setTitle(title);
}

////////////////////////////////////////////////////////////////////////////////

void InspectorFrame::setShowHidden(bool h)
{
    if (h != show_hidden)
    {
        show_hidden = h;
        redoLayout();
    }
}

////////////////////////////////////////////////////////////////////////////////

void InspectorFrame::redoLayout()
{
    QList<DatumRow*> rows;
    for (auto c : childItems())
        if (auto row = dynamic_cast<DatumRow*>(c))
        {
            if (show_hidden || !row->shouldBeHidden())
            {
                rows.append(row);
            }
            else
            {
                row->hide();
            }
        }

    // Show all rows that made it into our list
    for (auto r : rows)
    {
        r->show();
    }

    // Sort datums by row order
    std::sort(rows.begin(), rows.end(),
          [](const DatumRow* a, const DatumRow* b)
          { return a->getIndex() < b->getIndex(); });

    {   // Pad the row labels for alignment
        float max_label = 0;
        for (auto row : rows)
            max_label = std::max(max_label, row->labelWidth());
        for (auto row : rows)
            row->padLabel(max_label);
    }

    {   // Pad all of the rows (including the title) to the same width
        float max_width = title_row->minWidth();
        for (auto row : rows)
            max_width = std::max(max_width, row->minWidth());
        title_row->setWidth(max_width);
        for (auto row : rows)
            row->setWidth(max_width);
    }

    {   // Spread out the rows along the Y axis
        float y = title_row->boundingRect().height() + PADDING_ROWS;
        for (auto row : rows)
        {
            row->setPos(0, y);
            y += row->boundingRect().height() + PADDING_ROWS;
        }

        // Set the export button's size and position
        if (export_button->isVisible()) {
            auto w = boundingRect().width() / 2;
            export_button->setPos(w/2, y);
            export_button->setWidth(w);
        }
    }

    prepareGeometryChange();
}

////////////////////////////////////////////////////////////////////////////////

void InspectorFrame::setFocus(bool focus)
{
    if (focus != has_focus)
    {
        has_focus = focus;

        // Ease the hover glow in and out
        auto anim = new QVariantAnimation(this);
        anim->setDuration(120);
        anim->setStartValue(focus_glow);
        anim->setEndValue(focus ? 1.0 : 0.0);
        connect(anim, &QVariantAnimation::valueChanged,
                [this](const QVariant& v){
                    focus_glow = v.toReal();
                    update(); });
        anim->start(QAbstractAnimation::DeleteWhenStopped);

        prepareGeometryChange();
    }
}

void InspectorFrame::focusNext(DatumEditor* prev)
{
    bool next = false;

    for (auto row : visibleRows())
    {
        if (prev == row->editor)
        {
            next = true;
        }
        else if (next && row->editor->isEnabled())
        {
            prev->clearFocus();
            row->editor->setFocus();
            return;
        }
    }
}

void InspectorFrame::focusPrev(DatumEditor* next)
{
    DatumRow* prev = NULL;

    for (auto row : visibleRows())
    {
        if (next == row->editor)
        {
            if (prev)
            {
                prev->editor->setFocus();
                next->clearFocus();
            }
            return;
        }
        prev = row;
    }
}

////////////////////////////////////////////////////////////////////////////////

void InspectorFrame::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
    if (dragging)
    {
        setPos(event->scenePos());
        event->accept();
    }
    else
    {
        QGraphicsItem::mouseMoveEvent(event);
    }
}

void InspectorFrame::contextMenuEvent(QGraphicsSceneContextMenuEvent* e)
{
    Q_UNUSED(e);
    QString desc = QString::fromStdString(node->getName());

    QScopedPointer<QMenu> menu(new QMenu());
    auto jump_to = new QAction("Zoom to " + desc, menu.data());

    menu->addAction(jump_to);
    connect(jump_to, &QAction::triggered, this, &InspectorFrame::onZoomTo);

    menu->exec(QCursor::pos());
}

void InspectorFrame::mouseReleaseEvent(QGraphicsSceneMouseEvent* event)
{
    QGraphicsItem::mouseReleaseEvent(event);

    if (dragging)
    {
        ungrabMouse();
    }
    else if (event->button() == Qt::LeftButton)
    {
        // Store an Undo command for this drag
        const auto delta = event->scenePos() -
                     event->buttonDownScenePos(Qt::LeftButton);
        if (delta != QPointF(0,0))
        {
            static_cast<CanvasScene*>(scene())->endDrag(delta);
        }
    }
    dragging = false;
}

////////////////////////////////////////////////////////////////////////////////

void InspectorFrame::hoverEnterEvent(QGraphicsSceneHoverEvent *event)
{
    QGraphicsItem::hoverEnterEvent(event);
    emit(onFocus(true));
}

void InspectorFrame::hoverLeaveEvent(QGraphicsSceneHoverEvent *event)
{
    QGraphicsItem::hoverLeaveEvent(event);
    emit(onFocus(false));
}

////////////////////////////////////////////////////////////////////////////////

QString InspectorFrame::getTitle() const
{
    return title_row->getTitle();
}

ExportWorker* InspectorFrame::getExportWorker() const
{
    return export_button->getWorker();
}

void InspectorFrame::setExportWorker(ExportWorker* worker)
{
    if (export_button->setWorker(worker))
    {
        redoLayout();
    }
}

void InspectorFrame::clearExportWorker()
{
    if (export_button->clearWorker())
    {
        redoLayout();
    }
}

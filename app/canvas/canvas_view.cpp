#include <Python.h>

#include <QMenu>
#include <QMouseEvent>
#include <QClipboard>
#include <QMimeData>
#include <QTextCursor>
#include <QJsonDocument>
#include <QPropertyAnimation>

#include "app/app.h"
#include "app/colors.h"
#include "canvas/canvas_view.h"
#include "canvas/scene.h"
#include "canvas/add_node_popup.h"

#include "canvas/connection/connection.h"
#include "canvas/inspector/frame.h"
#include "canvas/subdatum/subdatum_frame.h"

#include "undo/undo_delete_multi.h"

#include "graph/serialize/serializer.h"
#include "graph/serialize/deserializer.h"
#include "graph/proxy/graph.h"

#include "graph/constructor/populate.h"

////////////////////////////////////////////////////////////////////////////////

CanvasView::CanvasView(CanvasScene* scene, QWidget* parent)
    : QGraphicsView(scene, parent), selecting(false)
{
    setStyleSheet("QGraphicsView { border-style: none; }");
    setRenderHints(QPainter::Antialiasing);
    setSceneRect(-width()/2, -height()/2, width(), height());

    QAbstractScrollArea::setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    QAbstractScrollArea::setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
}

////////////////////////////////////////////////////////////////////////////////

void CanvasView::mousePressEvent(QMouseEvent* event)
{
    QGraphicsView::mousePressEvent(event);
    if (!event->isAccepted()) {
            if (event->button() == Qt::LeftButton) {
                click_pos = mapToScene(event->pos());
            }
            else if (event->button() == Qt::RightButton)
            {
                openAddMenu();
            }
    }
}

void CanvasView::mouseReleaseEvent(QMouseEvent* event)
{
    QGraphicsView::mouseReleaseEvent(event);

    if (selecting)
    {
        QPainterPath p;
        p.addRect(QRectF(click_pos, drag_pos));
        scene()->setSelectionArea(p);
        selecting = false;
        scene()->invalidate(QRectF(), QGraphicsScene::ForegroundLayer);
    }
}

void CanvasView::mouseMoveEvent(QMouseEvent* event)
{
    QGraphicsView::mouseMoveEvent(event);
    if (scene()->mouseGrabberItem() == NULL && event->buttons() == Qt::LeftButton)
    {
        if (event->modifiers() & Qt::ShiftModifier)
        {
            drag_pos = mapToScene(event->pos());
            selecting = true;
            scene()->invalidate(QRectF(), QGraphicsScene::ForegroundLayer);
        }
        else
        {
            auto d = click_pos - mapToScene(event->pos());
            setSceneRect(sceneRect().translated(d.x(), d.y()));
        }
    }
}

void CanvasView::wheelEvent(QWheelEvent* event)
{
    QPointF a = mapToScene(event->position().toPoint());
    auto s = pow(1.001, -event->angleDelta().y());
    // Clamp cumulative zoom to a sane range
    const float zoom = transform().m11();
    s = fmax(0.02 / zoom, fmin(4 / zoom, s));
    scale(s, s);
    updateLOD();
    auto d = a - mapToScene(event->position().toPoint());
    setSceneRect(sceneRect().translated(d.x(), d.y()));
}

////////////////////////////////////////////////////////////////////////////////

void CanvasView::keyPressEvent(QKeyEvent* event)
{
    QGraphicsView::keyPressEvent(event);
    if (event->isAccepted())
    {
        return;
    }
    else if (event->key() == Qt::Key_Delete || event->key() == Qt::Key_Backspace)
    {
        deleteSelected();
    }
    else if (event->key() == Qt::Key_A &&
                (event->modifiers() & Qt::ShiftModifier))
    {
        openAddMenu();
    }
    else if (event->key() == Qt::Key_Tab)
    {
        openAddPalette();
    }
}

void CanvasView::mouseDoubleClickEvent(QMouseEvent* event)
{
    QGraphicsView::mouseDoubleClickEvent(event);
    if (!event->isAccepted() && event->button() == Qt::LeftButton)
        openAddPalette();
}

void CanvasView::openAddPalette()
{
    auto p = new AddNodePopup(
            static_cast<CanvasScene*>(scene())->getGraph(),
            [this](Node* n){ this->grabNode(n); },
            this);
    p->popUp(QCursor::pos());
}

void CanvasView::openAddMenu()
{
    QMenu* m = new QMenu(this);
    populateNodeMenu(
            m, static_cast<CanvasScene*>(scene())->getGraph(),
            [&](Node* n){ this->grabNode(n); },
            [&](Datum* d){ this->grabDatum(d); });

    m->exec(QCursor::pos());
    m->deleteLater();
}

////////////////////////////////////////////////////////////////////////////////

void CanvasView::deleteSelected()
{
    QSet<Node*> nodes;
    QSet<QPair<const Datum*, Datum*>> links;
    QSet<Datum*> datums;

     // Find all selected links
    for (auto i : scene()->selectedItems())
    {
        if (auto c = dynamic_cast<Connection*>(i))
        {
            links.insert(QPair<const Datum*, Datum*>(
                        c->sourceDatum(),
                        c->targetDatum()));
        }
        else if (auto p = dynamic_cast<InspectorFrame*>(i))
        {
            nodes.insert(p->getNode());
        }
        else if (auto s = dynamic_cast<SubdatumFrame*>(i))
        {
            datums.insert(s->getDatum());
        }
    }

    App::instance()->pushUndoStack(new UndoDeleteMulti(nodes, datums, links));
}

////////////////////////////////////////////////////////////////////////////////

void CanvasView::drawBackground(QPainter* painter, const QRectF& rect)
{
    painter->setBrush(Colors::base00);
    painter->setPen(Qt::NoPen);
    painter->drawRect(rect);

    // Two-level dot grid whose opacity fades with zoom instead of
    // vanishing at a hard threshold.
    const float zoom = painter->transform().m11();
    const int minor = 20, major = 100;

    auto drawGrid = [&](int d, float alpha){
        if (alpha <= 0.02f)
            return;
        QVector<QPointF> pts;
        pts.reserve(int((rect.width() / d + 2) * (rect.height() / d + 2)));
        for (int i = int(rect.left() / d) * d; i < rect.right(); i += d)
            for (int j = int(rect.top() / d) * d; j < rect.bottom(); j += d)
                pts.append(QPointF(i, j));
        auto c = Colors::base02;
        c.setAlphaF(alpha);
        painter->setPen(QPen(c, 2 / zoom));
        painter->drawPoints(pts.data(), pts.size());
    };

    // Minor dots fade out as you zoom away; major dots persist longer.
    drawGrid(minor, fmin(1.f, (zoom - 0.25f) * 2));
    drawGrid(major, fmin(1.f, zoom * 2));
}

void CanvasView::drawForeground(QPainter* painter, const QRectF& rect)
{
    Q_UNUSED(rect);

    if (selecting)
    {
        auto fill = Colors::amber;
        fill.setAlpha(25);
        painter->setBrush(fill);
        painter->setPen(QPen(Colors::amber, 1));
        painter->drawRect(QRectF(click_pos, drag_pos));
    }
}

////////////////////////////////////////////////////////////////////////////////

void CanvasView::onCopy()
{
    if (auto i = dynamic_cast<QGraphicsTextItem*>(scene()->focusItem()))
    {
        QApplication::clipboard()->setText(i->textCursor().selectedText());
    }
    else
    {
        // Find all selected nodes
        QList<Node*> selected;
        for (auto i : scene()->selectedItems())
            if (auto r = dynamic_cast<InspectorFrame*>(i))
                selected << r->getNode();

        if (!selected.isEmpty())
        {
            QJsonArray out;
            auto i = App::instance()->getProxy()->canvasInfo();
            for (auto n : selected)
                out << SceneSerializer::serializeNode(n, &i);

            auto data = new QMimeData();
            data->setData("sb::canvas", QJsonDocument(out).toJson());
            QApplication::clipboard()->setMimeData(data);
        }
    }
}

void CanvasView::onCut()
{
    if (auto i = dynamic_cast<QGraphicsTextItem*>(scene()->focusItem()))
    {
        QApplication::clipboard()->setText(i->textCursor().selectedText());
        i->textCursor().insertText("");
    }
}

void CanvasView::onPaste()
{
    if (auto i = dynamic_cast<QGraphicsTextItem*>(scene()->focusItem()))
    {
        i->textCursor().insertText(QApplication::clipboard()->text());
    }
    else
    {
        auto data = QApplication::clipboard()->mimeData();
        if (data->hasFormat("sb::canvas"))
            pasteNodes(QJsonDocument::fromJson(
                        data->data("sb::canvas")).array());
    }
}

void CanvasView::pasteNodes(QJsonArray array)
{
    QMap<uint64_t, uint64_t> uid_map;
    auto g = static_cast<CanvasScene*>(scene())->getGraph();

    {   // Get the next n UIDs and update the nodes in the array.
        std::list<uint64_t> uids = g->getUIDs(array.size());
        auto itr = uids.begin();
        for (int i=0; i < array.size(); ++i)
        {
            // Update this node's UID and store the change in uid_map
            auto node = array[i].toObject();
            uid_map[node["uid"].toInt()] = *itr;
            node["uid"] = int((*itr)++);
            array[i] = node;
        }

        // For every connection datum, remap UIDs in the expression
        for (int i=0; i < array.size(); ++i)
        {
            auto node = array[i].toObject();
            auto datums = node["datums"].toArray();
            for (int j=0; j < datums.size(); ++j)
            {
                auto d = datums[j].toObject();
                auto expr = d["expr"].toString();

                if (expr.startsWith(Datum::SIGIL_CONNECTION))
                {
                    // Process the list of connections, remapping any
                    // connection that is within the pasted group.
                    QStringList ds = expr.mid(2, expr.size() - 3).split(",");
                    QString out = Datum::SIGIL_CONNECTION + QString("[");
                    for (auto d : ds)
                    {
                        bool ok = true;
                        auto split = d.replace("__", "").split(".");
                        int n = split[0].toInt(&ok);
                        if (ok)
                        {
                            if (out.size() > 2)
                            {
                                out += ",";
                            }

                            if (uid_map.contains(n))
                            {
                                out += "__" + QString::number(uid_map[n]);
                            }
                            else
                            {
                                out += "__" + QString::number(n);
                            }
                            out += ".__" + split[1];
                        }
                    }
                    expr = out + "]";
                }
                d["expr"] = expr;
                datums[j] = d;
            }
            node["datums"] = datums;
            array[i] = node;
        }
    }

    SceneDeserializer::Info ds;
    for (auto n_ : array)
    {
        auto n = n_.toObject();

        // Update this node's name to make it unique
        auto name = n["name"].toString();
        if (!g->isNameUnique(name.toStdString()))
        {
            // Trim trailing numbers from the node's name
            while (name.at(name.size() - 1).isNumber())
                name = name.left(name.size() - 1);
            if (name.isEmpty())
                name = "n";
            // Then use the remaining string as a prefix
            n["name"] = QString::fromStdString(g->nextName(name.toStdString()));
        }
        SceneDeserializer::deserializeNode(n, g, &ds);
    }

    // Update the inspector positions by shifting a bit down and over
    for (auto& i : ds.frames.inspector)
        i += QPointF(10, 10);

    // Pull out the nodes that were just appended to the graph
    QSet<Node*> nodes;
    {
        auto all_nodes = g->childNodes();
        auto itr = all_nodes.rbegin();
        for (int i=0; i < array.size(); ++i)
            nodes.insert(*(itr++));
    }

    // Clear all of the currently selected nodes
    scene()->clearSelection();

    // Load new inspector positions and select them
    App::instance()->getProxy()->setPositions(ds.frames, true);
}

////////////////////////////////////////////////////////////////////////////////

template <class T>
void CanvasView::grab(T* t)
{
    t->setSelected(true);
    t->setPos(mapToScene(mapFromGlobal(QCursor::pos())));
    t->setDragging(true);
    t->grabMouse();
}

InspectorFrame* CanvasView::inspectorFor(Node* n)
{
    for (auto c : scene()->items())
    {
        if (auto i = dynamic_cast<InspectorFrame*>(c))
        {
            if (i->getNode() == n)
            {
                return i;
            }
        }
    }
    return nullptr;
}

void CanvasView::grabNode(Node* n)
{
    if (auto i = inspectorFor(n))
    {
        grab(i);
    }
}

void CanvasView::grabDatum(Datum* d)
{
    for (auto c : scene()->items())
    {
        if (auto s = dynamic_cast<SubdatumFrame*>(c))
        {
            if (s->getDatum() == d)
            {
                grab(s);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void CanvasView::zoomTo(Node* n)
{
    auto inspector = inspectorFor(n);
    Q_ASSERT(inspector);

    auto a = new QPropertyAnimation(this, "_center");
    a->setDuration(100);
    a->setStartValue(getCenter());
    a->setEndValue(inspector->sceneBoundingRect().center());

    auto b = new QPropertyAnimation(this, "_zoom");
    b->setDuration(100);
    b->setStartValue(getZoom());
    b->setEndValue(1);
    b->setEasingCurve(QEasingCurve::InQuart);

    a->start(QPropertyAnimation::DeleteWhenStopped);
    b->start(QPropertyAnimation::DeleteWhenStopped);
}

void CanvasView::updateLOD()
{
    // Below the readability threshold, nodes collapse to name cards
    const bool low = transform().m11() < CANVAS_LOD_THRESHOLD;
    if (low != low_detail_mode)
    {
        low_detail_mode = low;

        // Flag the scene first: connections consult it when their
        // ports' visibility changes, so they survive the collapse.
        if (auto s = dynamic_cast<CanvasScene*>(scene()))
            s->setLowDetail(low);

        for (auto i : scene()->items())
            if (auto f = dynamic_cast<InspectorFrame*>(i))
                f->setLowDetail(low);
    }
}

void CanvasView::zoomToFit(QRectF r)
{
    const auto vp = viewport()->rect();
    float fit = 1;
    if (vp.width() > 0 && r.width() > 0 && r.height() > 0)
        fit = fmin(1.0, 0.85 * fmin(vp.width() / r.width(),
                                    vp.height() / r.height()));
    setZoom(fit);
    setCenter(r.center());
}

void CanvasView::setCenter(QPointF p)
{
    auto t = sceneRect();
    setSceneRect(sceneRect().translated(p.x() - t.center().x(),
                                        p.y() - t.center().y()));
}

QPointF CanvasView::getCenter() const
{
    auto t = sceneRect();
    return QPointF(t.center().x(), t.center().y());
}

void CanvasView::setZoom(float z)
{
    auto t = transform();
    scale(z / t.m11(), z / t.m22());
    updateLOD();
}

float CanvasView::getZoom() const
{
    return transform().m11();
}

#pragma once

#include <QGraphicsScene>
#include <QSet>

class Graph;
class CanvasView;
class DummyConnection;
class InputPort;
class Datum;
class InspectorFrame;

// Zoom level below which the canvas switches to low-detail rendering
static constexpr float CANVAS_LOD_THRESHOLD = 0.32f;

class CanvasScene : public QGraphicsScene
{
public:
    CanvasScene(Graph* g, QObject* parent);

    /*
     *  Canvas-wide low-detail (zoomed out) mode; set by views, consulted
     *  by connections deciding whether hidden ports should hide them.
     */
    void setLowDetail(bool low) { low_detail = low; }
    bool lowDetail() const { return low_detail; }

    /*
     *  The in-progress sticky wire, if any (see DummyConnection)
     */
    void setActiveDummy(DummyConnection* d) { active_dummy = d; }
    DummyConnection* activeDummy() const { return active_dummy; }

    /*
     *  Returns a new CanvasView object looking at this scene
     */
    CanvasView* getView(QWidget* parent=NULL);

    /*
     *  Inspector registry, maintained by InspectorFrame's ctor/dtor
     *  so per-frame passes (floating labels) don't have to walk and
     *  dynamic_cast every item in the scene.
     */
    void registerInspector(InspectorFrame* f) { inspectors.insert(f); }
    void unregisterInspector(InspectorFrame* f) { inspectors.remove(f); }
    const QSet<InspectorFrame*>& allInspectors() const { return inspectors; }

private:
    QSet<InspectorFrame*> inspectors;
    bool low_detail=false;
    DummyConnection* active_dummy=nullptr;

public:

    /*
     *  Returns this scene's graph object
     */
    Graph* getGraph() const { return g; }

    /*
     *  Creates an undo command that undoes a drag operation
     */
    void endDrag(QPointF delta);

    /*
     *  Returns the top item of the given class at the given position
     */
    template <class T> T* itemAt(QPointF pos) const;

    /*
     *  Returns the input port nearest to the given point
     *  that accepts the given datum as an incoming link.
     */
    InputPort* inputPortNear(QPointF pos, Datum* source=NULL) const;

    /*
     *  Returns the top input port at the given position
     */
    InputPort* inputPortAt(QPointF pos) const;

protected:
    Graph* const g;
};

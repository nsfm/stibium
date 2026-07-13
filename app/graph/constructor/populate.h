#pragma once

#include <QMenu>

#include <functional>

#include "graph/graph.h"
#include "graph/constructor/constructor.h"

/*
 *  One addable node: menu category path, display title, constructor.
 */
struct NodeEntry
{
    QStringList category;
    QString title;
    NodeConstructorFunction constructor;
};

/*
 *  Returns the list of all .node entries (scanned once, then cached).
 */
const QList<NodeEntry>& nodeEntries();

void emptyNodeCallback(Node* n);
void emptyDatumCallback(Datum* d);

void populateNodeMenu(
        QMenu* menu, Graph* g,
        std::function<void(Node*)> callback=emptyNodeCallback,
        std::function<void(Datum*)> datum_callback=emptyDatumCallback);

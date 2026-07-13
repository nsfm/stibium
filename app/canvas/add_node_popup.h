#pragma once

#include <Python.h>

#include <QFrame>
#include <QLineEdit>
#include <QListWidget>

#include <functional>

class Graph;
class Node;

/*
 *  Fuzzy-search node palette: type a few characters, hit enter, get a
 *  node. Opened by double-clicking the canvas or pressing Tab.
 */
class AddNodePopup : public QFrame
{
    Q_OBJECT
public:
    AddNodePopup(Graph* g, std::function<void(Node*)> callback,
                 QWidget* parent=nullptr);

    /*
     *  Shows the popup near the given global position and focuses
     *  the search field.
     */
    void popUp(QPoint global_pos);

protected:
    bool eventFilter(QObject* obj, QEvent* event) override;

    void refilter(QString q);
    void accept();

    Graph* const graph;
    const std::function<void(Node*)> callback;

    QLineEdit* edit;
    QListWidget* list;
};

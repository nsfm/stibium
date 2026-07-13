#include <Python.h>

#include <QVBoxLayout>
#include <QKeyEvent>

#include "canvas/add_node_popup.h"

#include "app/app.h"
#include "app/colors.h"
#include "graph/constructor/populate.h"
#include "undo/undo_add_node.h"

namespace
{
// True if every character of q appears in order within t
bool subsequence(const QString& q, const QString& t)
{
    int i = 0;
    for (auto c : t)
        if (i < q.size() && c.toLower() == q[i].toLower())
            i++;
    return i == q.size();
}

// Lower is better; -1 excludes the entry
int score(const QString& q, const NodeEntry& e)
{
    const auto title = e.title;
    if (q.isEmpty())
        return 5;
    if (title.startsWith(q, Qt::CaseInsensitive))
        return 0;
    for (auto word : title.split(' '))
        if (word.startsWith(q, Qt::CaseInsensitive))
            return 1;
    if (title.contains(q, Qt::CaseInsensitive))
        return 2;
    if (subsequence(q, title))
        return 3;
    if (e.category.join(' ').contains(q, Qt::CaseInsensitive))
        return 4;
    return -1;
}
}

AddNodePopup::AddNodePopup(Graph* g, std::function<void(Node*)> cb,
                           QWidget* parent)
    : QFrame(parent, Qt::Popup | Qt::FramelessWindowHint),
      graph(g), callback(cb),
      edit(new QLineEdit(this)), list(new QListWidget(this))
{
    setAttribute(Qt::WA_DeleteOnClose);
    setFixedWidth(340);

    auto layout = new QVBoxLayout(this);
    layout->setContentsMargins(6, 6, 6, 6);
    layout->setSpacing(4);
    layout->addWidget(edit);
    layout->addWidget(list);

    edit->setPlaceholderText("Add node...");
    edit->installEventFilter(this);
    list->setFocusPolicy(Qt::NoFocus);

    setStyleSheet(QString(
        "AddNodePopup { background: %1; border: 1px solid %2;"
        "               border-radius: 6px; }"
        "QLineEdit { background: %3; border: 1px solid %2;"
        "            border-radius: 4px; padding: 4px; }"
        "QListWidget { background: transparent; border: none; }")
        .arg(Colors::base01.name())
        .arg(Colors::base02.name())
        .arg(Colors::base00.name()));

    connect(edit, &QLineEdit::textChanged,
            this, &AddNodePopup::refilter);
    connect(edit, &QLineEdit::returnPressed,
            this, &AddNodePopup::accept);
    connect(list, &QListWidget::itemActivated,
            this, &AddNodePopup::accept);
    connect(list, &QListWidget::itemClicked,
            this, &AddNodePopup::accept);

    refilter("");
}

void AddNodePopup::popUp(QPoint global_pos)
{
    move(global_pos);
    show();
    edit->setFocus();
}

void AddNodePopup::refilter(QString q)
{
    list->clear();

    QList<QPair<int, const NodeEntry*>> matches;
    for (const auto& e : nodeEntries())
    {
        const int s = score(q, e);
        if (s >= 0)
            matches.append({s, &e});
    }
    std::sort(matches.begin(), matches.end(),
            [](const QPair<int, const NodeEntry*>& a,
               const QPair<int, const NodeEntry*>& b){
        return a.first != b.first ? a.first < b.first
                                  : a.second->title < b.second->title; });

    for (const auto& m : matches)
    {
        auto item = new QListWidgetItem(
                QString("%1      %2").arg(m.second->title)
                                     .arg(m.second->category.join(" / ")));
        item->setData(Qt::UserRole,
                      QVariant::fromValue((void*)m.second));
        list->addItem(item);
    }

    if (list->count())
        list->setCurrentRow(0);

    const int rows = std::min(list->count(), 12);
    list->setFixedHeight(rows ? rows * list->sizeHintForRow(0) + 8 : 24);
    adjustSize();
}

void AddNodePopup::accept()
{
    auto item = list->currentItem();
    if (!item)
        return;
    auto entry = static_cast<const NodeEntry*>(
            item->data(Qt::UserRole).value<void*>());

    Node* n = entry->constructor(graph);
    App::instance()->pushUndoStack(new UndoAddNode(n));
    close();
    callback(n);
}

bool AddNodePopup::eventFilter(QObject* obj, QEvent* event)
{
    if (obj == edit && event->type() == QEvent::KeyPress)
    {
        auto k = static_cast<QKeyEvent*>(event);
        if (k->key() == Qt::Key_Down || k->key() == Qt::Key_Up)
        {
            const int d = (k->key() == Qt::Key_Down) ? 1 : -1;
            const int row = list->currentRow() + d;
            if (row >= 0 && row < list->count())
                list->setCurrentRow(row);
            return true;
        }
    }
    return QFrame::eventFilter(obj, event);
}

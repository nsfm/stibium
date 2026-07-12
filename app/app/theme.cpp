#include <Python.h>

#include <QApplication>
#include <QPalette>
#include <QStyleFactory>

#include "app/theme.h"
#include "app/colors.h"

void Theme::apply(QApplication* app)
{
    app->setStyle(QStyleFactory::create("fusion"));

    QPalette p;
    p.setColor(QPalette::Window, Colors::base01);
    p.setColor(QPalette::WindowText, Colors::base05);
    p.setColor(QPalette::Base, Colors::base00);
    p.setColor(QPalette::AlternateBase, Colors::base01);
    p.setColor(QPalette::Text, Colors::base05);
    p.setColor(QPalette::Button, Colors::base01);
    p.setColor(QPalette::ButtonText, Colors::base05);
    p.setColor(QPalette::BrightText, Colors::base07);
    p.setColor(QPalette::Highlight, Colors::amber);
    p.setColor(QPalette::HighlightedText, Colors::base00);
    p.setColor(QPalette::ToolTipBase, Colors::base02);
    p.setColor(QPalette::ToolTipText, Colors::base05);
    p.setColor(QPalette::PlaceholderText, Colors::base03);
    p.setColor(QPalette::Link, Colors::amber);
    p.setColor(QPalette::LinkVisited, Colors::brown);

    for (auto role : {QPalette::WindowText, QPalette::Text,
                      QPalette::ButtonText})
        p.setColor(QPalette::Disabled, role, Colors::base03);
    p.setColor(QPalette::Disabled, QPalette::Highlight, Colors::base02);

    app->setPalette(p);

    app->setStyleSheet(QString(
        "QMenu {"
        "    background-color: %1;"
        "    border: 1px solid %2;"
        "    padding: 4px;"
        "}"
        "QMenu::item {"
        "    padding: 4px 24px 4px 12px;"
        "    border-radius: 4px;"
        "}"
        "QMenu::item:selected {"
        "    background-color: %3;"
        "    color: %4;"
        "}"
        "QMenu::separator {"
        "    height: 1px;"
        "    background: %2;"
        "    margin: 4px 8px;"
        "}"
        "QToolTip {"
        "    background-color: %1;"
        "    color: %5;"
        "    border: 1px solid %2;"
        "    padding: 3px;"
        "}"
        "QScrollBar:vertical {"
        "    background: %6;"
        "    width: 10px;"
        "    margin: 0;"
        "}"
        "QScrollBar::handle:vertical {"
        "    background: %2;"
        "    border-radius: 4px;"
        "    min-height: 24px;"
        "}"
        "QScrollBar::handle:vertical:hover { background: %7; }"
        "QScrollBar:horizontal {"
        "    background: %6;"
        "    height: 10px;"
        "    margin: 0;"
        "}"
        "QScrollBar::handle:horizontal {"
        "    background: %2;"
        "    border-radius: 4px;"
        "    min-width: 24px;"
        "}"
        "QScrollBar::handle:horizontal:hover { background: %7; }"
        "QScrollBar::add-line, QScrollBar::sub-line { height: 0; width: 0; }"
        "QScrollBar::add-page, QScrollBar::sub-page { background: none; }")
        .arg(Colors::base01.name())
        .arg(Colors::base02.name())
        .arg(Colors::amber.name())
        .arg(Colors::base00.name())
        .arg(Colors::base05.name())
        .arg(Colors::base00.name())
        .arg(Colors::base03.name()));
}

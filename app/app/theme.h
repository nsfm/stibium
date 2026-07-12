#pragma once

class QApplication;

namespace Theme
{
    /*
     *  Applies the Stibium dark theme: Fusion style, a QPalette built
     *  from the Colors namespace, and a small stylesheet for widgets
     *  the palette can't reach (menus, scrollbars, tooltips).
     */
    void apply(QApplication* app);
}

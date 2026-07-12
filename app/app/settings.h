#pragma once

#include <QVariant>
#include <QString>

/*
 *  Thin wrapper over QSettings, INI format: ~/.config/Stibium/Stibium.ini
 *  on Linux, platform-appropriate locations elsewhere. Keys are grouped
 *  with slashes, e.g. "autosave/interval_ms", "editor/font_family".
 */
namespace Settings
{
    /* Selects INI format; call once before the first get/set. */
    void init();

    QVariant get(const QString& key, const QVariant& fallback);
    void set(const QString& key, const QVariant& value);
}

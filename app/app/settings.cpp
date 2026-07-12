#include <QSettings>

#include "app/settings.h"

void Settings::init()
{
    QSettings::setDefaultFormat(QSettings::IniFormat);
}

QVariant Settings::get(const QString& key, const QVariant& fallback)
{
    return QSettings().value(key, fallback);
}

void Settings::set(const QString& key, const QVariant& value)
{
    QSettings().setValue(key, value);
}

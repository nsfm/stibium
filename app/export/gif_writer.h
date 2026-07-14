#pragma once

#include <QImage>
#include <QList>
#include <QString>

/*
 *  Minimal animated GIF89a writer: LZW compression, per-frame local
 *  color tables (Qt's Indexed8 quantization + dithering), infinite
 *  loop.  delays_cs holds one per-frame delay in centiseconds (a
 *  single entry repeats for all frames).
 *
 *  Returns false on I/O failure or bad input.
 */
bool save_gif(const QList<QImage>& frames,
              const QList<int>& delays_cs,
              const QString& filename);

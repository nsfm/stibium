#include <QFile>

#include <cstdint>
#include <cstring>
#include <unordered_map>
#include <vector>

#include "export/gif_writer.h"

namespace {

/*
 *  LZW bit-packer: codes are written LSB-first and flushed through
 *  255-byte sub-blocks (the GIF data-block framing).
 */
struct BlockWriter {
    QFile& f;
    uint32_t bits = 0;
    int nbits = 0;
    std::vector<uint8_t> block;

    explicit BlockWriter(QFile& file) : f(file) {}

    void put(uint16_t code, int width)
    {
        bits |= uint32_t(code) << nbits;
        nbits += width;
        while (nbits >= 8)
        {
            block.push_back(uint8_t(bits & 0xff));
            bits >>= 8;
            nbits -= 8;
            if (block.size() == 255)
                flush_block();
        }
    }

    void finish()
    {
        if (nbits > 0)
        {
            block.push_back(uint8_t(bits & 0xff));
            bits = 0;
            nbits = 0;
        }
        if (!block.empty())
            flush_block();
        const uint8_t terminator = 0;
        f.write(reinterpret_cast<const char*>(&terminator), 1);
    }

private:
    void flush_block()
    {
        const uint8_t len = uint8_t(block.size());
        f.write(reinterpret_cast<const char*>(&len), 1);
        f.write(reinterpret_cast<const char*>(block.data()),
                block.size());
        block.clear();
    }
};

/*  Standard GIF LZW over 8-bit pixels: clear=256, end=257, codes
 *  grow 9 -> 12 bits, dictionary reset at 4096. */
void lzw_encode(QFile& f, const uint8_t* px, size_t count)
{
    const uint16_t CLEAR = 256, END = 257;
    const uint8_t min_code_size = 8;
    f.write(reinterpret_cast<const char*>(&min_code_size), 1);

    BlockWriter out(f);
    std::unordered_map<uint32_t, uint16_t> dict;
    uint16_t next = 258;
    int width = 9;

    out.put(CLEAR, width);
    uint16_t cur = px[0];
    for (size_t i = 1; i < count; ++i)
    {
        const uint8_t c = px[i];
        const uint32_t key = (uint32_t(cur) << 8) | c;
        const auto it = dict.find(key);
        if (it != dict.end())
        {
            cur = it->second;
            continue;
        }

        out.put(cur, width);
        dict[key] = next++;
        // The decoder runs one dictionary entry behind the encoder
        // (it cannot add on its first code), so it widens one code
        // later than a naive encoder would: match that, not the
        // local counter
        if (next == (1u << width) + 1 && width < 12)
            ++width;
        if (next == 4096)
        {
            out.put(CLEAR, width);
            dict.clear();
            next = 258;
            width = 9;
        }
        cur = c;
    }
    out.put(cur, width);
    out.put(END, width);
    out.finish();
}

void put16(QFile& f, uint16_t v)
{
    const uint8_t b[2] = {uint8_t(v & 0xff), uint8_t(v >> 8)};
    f.write(reinterpret_cast<const char*>(b), 2);
}

}  // namespace

bool save_gif(const QList<QImage>& frames,
              const QList<int>& delays_cs,
              const QString& filename)
{
    if (frames.isEmpty() || delays_cs.isEmpty())
        return false;
    const int W = frames[0].width(), H = frames[0].height();
    if (W <= 0 || H <= 0 || W > 0xffff || H > 0xffff)
        return false;

    QFile f(filename);
    if (!f.open(QIODevice::WriteOnly))
        return false;

    f.write("GIF89a", 6);
    put16(f, uint16_t(W));
    put16(f, uint16_t(H));
    // No global color table; each frame carries its own palette
    const uint8_t screen[3] = {0x70, 0, 0};
    f.write(reinterpret_cast<const char*>(screen), 3);

    // NETSCAPE looping extension (0 = forever)
    const uint8_t loop[] = {0x21, 0xff, 0x0b,
                            'N','E','T','S','C','A','P','E','2','.','0',
                            0x03, 0x01, 0x00, 0x00, 0x00};
    f.write(reinterpret_cast<const char*>(loop), sizeof(loop));

    for (int fi = 0; fi < frames.size(); ++fi)
    {
        // Quantize: flatten alpha, then let Qt pick a palette with
        // diffusion dithering (per-frame local color table)
        QImage q = frames[fi];
        if (q.width() != W || q.height() != H)
            q = q.scaled(W, H);
        q = q.convertToFormat(QImage::Format_RGB32)
             .convertToFormat(QImage::Format_Indexed8,
                              Qt::DiffuseDither);
        const auto palette = q.colorTable();
        if (palette.isEmpty())
            return false;

        // Graphic control: per-frame delay, no transparency
        const int delay = delays_cs[fi < delays_cs.size() ? fi : 0];
        const uint8_t gce[] = {0x21, 0xf9, 0x04, 0x04,
                               uint8_t(delay & 0xff),
                               uint8_t((delay >> 8) & 0xff),
                               0x00, 0x00};
        f.write(reinterpret_cast<const char*>(gce), sizeof(gce));

        // Image descriptor with a 256-entry local color table
        const uint8_t sep = 0x2c;
        f.write(reinterpret_cast<const char*>(&sep), 1);
        put16(f, 0);
        put16(f, 0);
        put16(f, uint16_t(W));
        put16(f, uint16_t(H));
        const uint8_t flags = 0x87;     // local table, 2^8 entries
        f.write(reinterpret_cast<const char*>(&flags), 1);

        uint8_t table[256][3];
        for (int i = 0; i < 256; ++i)
        {
            const QRgb c = i < palette.size() ? palette[i] : 0;
            table[i][0] = uint8_t(qRed(c));
            table[i][1] = uint8_t(qGreen(c));
            table[i][2] = uint8_t(qBlue(c));
        }
        f.write(reinterpret_cast<const char*>(table), sizeof(table));

        // Pixel rows may be padded; pack them tight for the encoder
        std::vector<uint8_t> px(size_t(W) * H);
        for (int j = 0; j < H; ++j)
            memcpy(&px[size_t(W) * j], q.constScanLine(j), W);
        lzw_encode(f, px.data(), px.size());
    }

    const uint8_t trailer = 0x3b;
    f.write(reinterpret_cast<const char*>(&trailer), 1);
    return f.flush();
}

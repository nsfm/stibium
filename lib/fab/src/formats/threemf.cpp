#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>

#include <zlib.h>

#include "fab/formats/threemf.h"

/*
 *  A 3MF file is an OPC package: a ZIP archive holding a content-types
 *  manifest, a relationships file, and the model XML.  The writer below
 *  is a minimal streaming ZIP implementation on top of zlib's raw
 *  deflate - entries are compressed as they are generated, then the
 *  local header is patched with the final CRC and sizes (the file is
 *  seekable, so no data descriptors are needed).
 *
 *  Limits: no ZIP64, so any single entry (or the archive) must stay
 *  under 4 GiB compressed and uncompressed.  A mesh big enough to hit
 *  that (roughly 60M+ triangles of XML) fails cleanly instead.
 */

namespace {

// Fixed DOS timestamp for deterministic output (2026-01-01 00:00:00)
const uint16_t DOS_DATE = uint16_t(((2026 - 1980) << 9) | (1 << 5) | 1);
const uint16_t DOS_TIME = 0;

class ZipWriter {
public:
    explicit ZipWriter(const char* path)
        : f(fopen(path, "wb")), in_entry(false), failed(f == NULL) {}

    ~ZipWriter() {
        if (f)
            fclose(f);
    }

    bool ok() const { return !failed; }

    /*  Starts a new entry; write() streams data into it. */
    void begin(const std::string& name)
    {
        if (failed)
            return;

        Entry e;
        e.name = name;
        e.offset = ftell(f);
        e.crc = crc32(0, NULL, 0);
        e.csize = 0;
        e.usize = 0;
        entries.push_back(e);

        // Local file header, with zeroed CRC/sizes (patched in end())
        put32(0x04034b50);
        put16(20);              // version needed
        put16(0);               // flags
        put16(8);               // method: deflate
        put16(DOS_TIME);
        put16(DOS_DATE);
        put32(0);               // crc32 (patched)
        put32(0);               // compressed size (patched)
        put32(0);               // uncompressed size (patched)
        put16(name.size());
        put16(0);               // extra length
        raw(name.data(), name.size());

        memset(&strm, 0, sizeof(strm));
        if (deflateInit2(&strm, Z_DEFAULT_COMPRESSION, Z_DEFLATED,
                         -15, 8, Z_DEFAULT_STRATEGY) != Z_OK)
            failed = true;
        in_entry = true;
    }

    void write(const void* p, size_t n)
    {
        if (failed || !in_entry)
            return;

        Entry& e = entries.back();
        e.crc = crc32(e.crc, static_cast<const Bytef*>(p), n);
        if (e.usize + n < e.usize || e.usize + n > UINT32_MAX)
        {
            failed = true;  // would need ZIP64
            return;
        }
        e.usize += n;

        strm.next_in = const_cast<Bytef*>(static_cast<const Bytef*>(p));
        strm.avail_in = n;
        pump(Z_NO_FLUSH);
    }

    void write(const std::string& s) { write(s.data(), s.size()); }

    /*  printf-style write */
    void writef(const char* fmt, ...)
    {
        char buf[512];
        va_list args;
        va_start(args, fmt);
        const int n = vsnprintf(buf, sizeof(buf), fmt, args);
        va_end(args);
        if (n > 0)
            write(buf, size_t(n) < sizeof(buf) ? size_t(n)
                                               : sizeof(buf) - 1);
    }

    /*  Finishes the current entry and patches its local header. */
    void end()
    {
        if (!in_entry)
            return;
        in_entry = false;

        strm.next_in = NULL;
        strm.avail_in = 0;
        pump(Z_FINISH);
        deflateEnd(&strm);
        if (failed)
            return;

        Entry& e = entries.back();
        const long here = ftell(f);
        e.csize = uint32_t(here - long(e.offset)) - 30 - e.name.size();

        // Patch crc32 + sizes at offset 14 of the local header
        if (fseek(f, long(e.offset) + 14, SEEK_SET) != 0)
        {
            failed = true;
            return;
        }
        put32(e.crc);
        put32(e.csize);
        put32(e.usize);
        if (fseek(f, here, SEEK_SET) != 0)
            failed = true;
    }

    /*  Adds a complete entry in one call. */
    void add(const std::string& name, const std::string& contents)
    {
        begin(name);
        write(contents);
        end();
    }

    /*  Writes the central directory; returns overall success. */
    bool finish()
    {
        if (failed)
            return false;

        const long cd_start = ftell(f);
        for (const auto& e : entries)
        {
            put32(0x02014b50);
            put16(20);          // version made by
            put16(20);          // version needed
            put16(0);           // flags
            put16(8);           // method
            put16(DOS_TIME);
            put16(DOS_DATE);
            put32(e.crc);
            put32(e.csize);
            put32(e.usize);
            put16(e.name.size());
            put16(0);           // extra length
            put16(0);           // comment length
            put16(0);           // disk number
            put16(0);           // internal attributes
            put32(0);           // external attributes
            put32(e.offset);
            raw(e.name.data(), e.name.size());
        }
        const long cd_end = ftell(f);

        put32(0x06054b50);      // end of central directory
        put16(0);               // this disk
        put16(0);               // central directory disk
        put16(entries.size()); // entries on this disk
        put16(entries.size()); // entries total
        put32(cd_end - cd_start);
        put32(cd_start);
        put16(0);               // comment length

        if (fflush(f) != 0)
            failed = true;
        return !failed;
    }

private:
    struct Entry {
        std::string name;
        uint32_t offset;
        uint32_t crc;
        uint32_t csize;
        uint32_t usize;
    };

    void raw(const void* p, size_t n)
    {
        if (!failed && fwrite(p, 1, n, f) != n)
            failed = true;
    }
    void put16(uint16_t v) { raw(&v, 2); }   // little-endian host assumed
    void put32(uint32_t v) { raw(&v, 4); }

    /*  Runs deflate until the input is consumed, writing output. */
    void pump(int flush)
    {
        unsigned char out[1 << 16];
        do {
            strm.next_out = out;
            strm.avail_out = sizeof(out);
            const int ret = deflate(&strm, flush);
            if (ret == Z_STREAM_ERROR)
            {
                failed = true;
                return;
            }
            raw(out, sizeof(out) - strm.avail_out);
            if (failed)
                return;
            if (flush == Z_FINISH && ret == Z_STREAM_END)
                return;
        } while (strm.avail_out == 0 || strm.avail_in > 0);
    }

    FILE* f;
    z_stream strm;
    std::vector<Entry> entries;
    bool in_entry;
    bool failed;
};

const char* CONTENT_TYPES =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/"
        "content-types\">\n"
    " <Default Extension=\"rels\" ContentType=\"application/vnd."
        "openxmlformats-package.relationships+xml\"/>\n"
    " <Default Extension=\"model\" ContentType=\"application/vnd."
        "ms-package.3dmanufacturing-3dmodel+xml\"/>\n"
    "</Types>\n";

const char* RELS =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
    "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/"
        "2006/relationships\">\n"
    " <Relationship Target=\"/3D/3dmodel.model\" Id=\"rel0\" "
        "Type=\"http://schemas.microsoft.com/3dmanufacturing/2013/01/"
        "3dmodel\"/>\n"
    "</Relationships>\n";

}  // namespace

bool save_3mf_indexed(const float* verts, uint32_t vert_count,
                      const uint32_t* indices, uint32_t tri_count,
                      const char* filename)
{
    ZipWriter zip(filename);
    if (!zip.ok())
        return false;

    zip.add("[Content_Types].xml", CONTENT_TYPES);
    zip.add("_rels/.rels", RELS);

    zip.begin("3D/3dmodel.model");
    zip.write(
        "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        "<model unit=\"millimeter\" xml:lang=\"en-US\" "
            "xmlns=\"http://schemas.microsoft.com/3dmanufacturing/core/"
            "2015/02\">\n"
        " <metadata name=\"Application\">Stibium</metadata>\n"
        " <resources>\n"
        "  <object id=\"1\" type=\"model\">\n"
        "   <mesh>\n"
        "    <vertices>\n");

    // %.9g round-trips single-precision floats exactly
    for (uint32_t v = 0; v < vert_count; ++v)
        zip.writef("     <vertex x=\"%.9g\" y=\"%.9g\" z=\"%.9g\"/>\n",
                   verts[v*3], verts[v*3 + 1], verts[v*3 + 2]);

    zip.write(
        "    </vertices>\n"
        "    <triangles>\n");

    for (uint32_t t = 0; t < tri_count; ++t)
        zip.writef("     <triangle v1=\"%u\" v2=\"%u\" v3=\"%u\"/>\n",
                   indices[t*3], indices[t*3 + 1], indices[t*3 + 2]);

    zip.write(
        "    </triangles>\n"
        "   </mesh>\n"
        "  </object>\n"
        " </resources>\n"
        " <build>\n"
        "  <item objectid=\"1\"/>\n"
        " </build>\n"
        "</model>\n");
    zip.end();

    return zip.finish();
}

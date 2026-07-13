#include "fab/mesh/trimesh.h"

#include <cctype>
#include <cstring>
#include <fstream>
#include <sstream>
#include <unordered_map>

namespace fab_mesh {

namespace {

// Read an entire file as raw bytes.  STL is binary-or-text, so we
// always open in binary mode and sniff the content ourselves.
bool slurp(const std::string& path, std::string* data)
{
    std::ifstream f(path, std::ios::binary);
    if (!f)
        return false;
    std::ostringstream ss;
    ss << f.rdbuf();
    *data = ss.str();
    return true;
}

// -0.0f and +0.0f compare equal but have different bit patterns, so
// normalize before hashing to weld corners that differ only in sign
// of zero.
inline uint32_t coord_bits(float f)
{
    if (f == 0.0f)      // catches -0.0f too
        f = 0.0f;
    uint32_t bits;
    std::memcpy(&bits, &f, sizeof(bits));
    return bits;
}

struct VertKey {
    uint32_t x, y, z;
    bool operator==(const VertKey& o) const
    { return x == o.x && y == o.y && z == o.z; }
};

struct VertKeyHash {
    size_t operator()(const VertKey& k) const
    {
        // A cheap 96-bit mix; collisions only cost a bucket compare.
        size_t h = k.x;
        h = h * 0x9e3779b1u + k.y;
        h = h * 0x9e3779b1u + k.z;
        return h;
    }
};

// Weld a flat triangle soup (9 floats per triangle) into indexed
// topology, dropping triangles that collapse to an edge or point.
// Fails (returns false) if nothing survives.
bool finalize(const std::vector<float>& soup, std::string header,
              TriMesh* out, std::string* err)
{
    out->verts.clear();
    out->tris.clear();
    out->header = std::move(header);

    std::unordered_map<VertKey, uint32_t, VertKeyHash> index;
    index.reserve(soup.size() / 9 + 1);

    auto intern = [&](const float* p) -> uint32_t {
        VertKey k{coord_bits(p[0]), coord_bits(p[1]), coord_bits(p[2])};
        auto it = index.find(k);
        if (it != index.end())
            return it->second;
        uint32_t id = out->vert_count();
        out->verts.push_back(p[0]);
        out->verts.push_back(p[1]);
        out->verts.push_back(p[2]);
        index.emplace(k, id);
        return id;
    };

    for (size_t i = 0; i + 9 <= soup.size(); i += 9)
    {
        uint32_t a = intern(&soup[i]);
        uint32_t b = intern(&soup[i + 3]);
        uint32_t c = intern(&soup[i + 6]);
        if (a == b || b == c || a == c)
            continue;   // degenerate after welding
        out->tris.push_back(a);
        out->tris.push_back(b);
        out->tris.push_back(c);
    }

    if (out->tris.empty())
    {
        *err = "mesh has no non-degenerate triangles";
        return false;
    }

    // Bounding box over the surviving vertices.
    out->bbox[0] = out->bbox[3] = out->verts[0];
    out->bbox[1] = out->bbox[4] = out->verts[1];
    out->bbox[2] = out->bbox[5] = out->verts[2];
    for (size_t i = 0; i < out->verts.size(); i += 3)
        for (int j = 0; j < 3; ++j)
        {
            float v = out->verts[i + j];
            if (v < out->bbox[j])     out->bbox[j] = v;
            if (v > out->bbox[3 + j]) out->bbox[3 + j] = v;
        }
    return true;
}

// Strip trailing NULs and spaces from a fixed-width binary header.
std::string trim_header(const char* p, size_t n)
{
    while (n > 0 && (p[n - 1] == '\0' || p[n - 1] == ' '))
        --n;
    return std::string(p, n);
}

// True if the first non-whitespace token is "solid".
bool sniff_ascii(const std::string& d)
{
    size_t i = 0;
    while (i < d.size() && std::isspace((unsigned char)d[i]))
        ++i;
    return d.compare(i, 5, "solid") == 0;
}

// Parse ASCII STL.  Returns false (leaving soup untouched) if the text
// does not actually contain the facet/vertex grammar -- this is what
// lets a binary file whose 80-byte header happens to start with
// "solid" fall through to the binary path.
bool parse_ascii(const std::string& d, std::vector<float>* soup,
                 std::string* header)
{
    // Header: the remainder of the first line after "solid".
    {
        size_t nl = d.find('\n');
        size_t start = 5;   // past "solid"
        std::string line = d.substr(start,
                                    (nl == std::string::npos ? d.size() : nl) - start);
        size_t a = line.find_first_not_of(" \t\r");
        size_t b = line.find_last_not_of(" \t\r");
        *header = (a == std::string::npos) ? "" : line.substr(a, b - a + 1);
    }

    std::istringstream in(d);
    std::string tok;
    std::vector<float> local;
    bool saw_facet = false;
    while (in >> tok)
    {
        if (tok == "facet")
        {
            saw_facet = true;
        }
        else if (tok == "vertex")
        {
            float x, y, z;
            if (!(in >> x >> y >> z))
                return false;   // malformed -> not really ASCII STL
            local.push_back(x);
            local.push_back(y);
            local.push_back(z);
        }
    }

    // Need real facet grammar and a whole number of triangles.
    if (!saw_facet || local.empty() || local.size() % 9 != 0)
        return false;

    *soup = std::move(local);
    return true;
}

bool parse_binary(const std::string& d, std::vector<float>* soup,
                  std::string* header, std::string* err)
{
    if (d.size() < 84)
    {
        *err = "binary STL truncated: shorter than 84-byte preamble";
        return false;
    }
    *header = trim_header(d.data(), 80);

    uint32_t count;
    std::memcpy(&count, d.data() + 80, 4);

    const size_t expected = 84 + size_t(count) * 50;
    if (d.size() < expected)
    {
        std::ostringstream ss;
        ss << "binary STL truncated: header claims " << count
           << " triangles (" << expected << " bytes) but file is "
           << d.size() << " bytes";
        *err = ss.str();
        return false;
    }

    soup->clear();
    soup->reserve(size_t(count) * 9);
    for (uint32_t t = 0; t < count; ++t)
    {
        const char* rec = d.data() + 84 + size_t(t) * 50;
        // Skip the 3-float normal (rec[0..11]); keep the 9 vertex floats.
        for (int f = 3; f < 12; ++f)
        {
            float v;
            std::memcpy(&v, rec + f * 4, 4);
            soup->push_back(v);
        }
    }
    return true;
}

}  // namespace

bool load_stl(const std::string& path, TriMesh* out, std::string* err)
{
    std::string data;
    if (!slurp(path, &data))
    {
        *err = "could not open '" + path + "'";
        return false;
    }
    if (data.empty())
    {
        *err = "file is empty";
        return false;
    }

    std::vector<float> soup;
    std::string header;

    if (sniff_ascii(data) && parse_ascii(data, &soup, &header))
        return finalize(soup, std::move(header), out, err);

    if (!parse_binary(data, &soup, &header, err))
        return false;
    return finalize(soup, std::move(header), out, err);
}

}  // namespace fab_mesh

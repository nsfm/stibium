#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <thread>
#include <vector>

#include "fab/mesh/mesh_import.h"
#include "fab/mesh/trimesh.h"
#include "fab/mesh/mesh_query.h"
#include "fab/tree/grid.h"

namespace fab_mesh {

////////////////////////////////////////////////////////////////////////////////
// SHA-256 (compact implementation of FIPS 180-4; no external deps)

namespace {

struct Sha256 {
    uint32_t h[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a,
                     0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
    uint8_t buf[64];
    uint64_t total = 0;
    size_t fill = 0;

    static uint32_t rotr(uint32_t x, int n)
    {
        return (x >> n) | (x << (32 - n));
    }

    void block(const uint8_t* p)
    {
        static const uint32_t K[64] = {
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
            0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
            0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
            0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
            0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
            0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
            0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
            0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
            0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
            0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
            0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
            0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
            0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
            0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
            0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
            0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};

        uint32_t w[64];
        for (int i = 0; i < 16; ++i)
            w[i] = uint32_t(p[i*4]) << 24 | uint32_t(p[i*4 + 1]) << 16 |
                   uint32_t(p[i*4 + 2]) << 8 | p[i*4 + 3];
        for (int i = 16; i < 64; ++i) {
            const uint32_t s0 = rotr(w[i-15], 7) ^ rotr(w[i-15], 18) ^
                                (w[i-15] >> 3);
            const uint32_t s1 = rotr(w[i-2], 17) ^ rotr(w[i-2], 19) ^
                                (w[i-2] >> 10);
            w[i] = w[i-16] + s0 + w[i-7] + s1;
        }

        uint32_t a = h[0], b = h[1], c = h[2], d = h[3],
                 e = h[4], f = h[5], g = h[6], hh = h[7];
        for (int i = 0; i < 64; ++i) {
            const uint32_t S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
            const uint32_t ch = (e & f) ^ (~e & g);
            const uint32_t t1 = hh + S1 + ch + K[i] + w[i];
            const uint32_t S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
            const uint32_t mj = (a & b) ^ (a & c) ^ (b & c);
            const uint32_t t2 = S0 + mj;
            hh = g;  g = f;  f = e;  e = d + t1;
            d = c;  c = b;  b = a;  a = t1 + t2;
        }
        h[0] += a;  h[1] += b;  h[2] += c;  h[3] += d;
        h[4] += e;  h[5] += f;  h[6] += g;  h[7] += hh;
    }

    void update(const void* data, size_t n)
    {
        const uint8_t* p = static_cast<const uint8_t*>(data);
        total += n;
        while (n) {
            const size_t take = n < 64 - fill ? n : 64 - fill;
            memcpy(buf + fill, p, take);
            fill += take;  p += take;  n -= take;
            if (fill == 64) {
                block(buf);
                fill = 0;
            }
        }
    }

    std::string hex()
    {
        const uint64_t bits = total * 8;
        const uint8_t one = 0x80;
        update(&one, 1);
        const uint8_t zero = 0;
        while (fill != 56)
            update(&zero, 1);
        uint8_t len[8];
        for (int i = 0; i < 8; ++i)
            len[i] = uint8_t(bits >> (56 - 8*i));
        update(len, 8);

        char out[65];
        for (int i = 0; i < 8; ++i)
            snprintf(out + i*8, 9, "%08x", h[i]);
        return std::string(out, 64);
    }
};

}  // namespace

std::string sha256_file(const std::string& path)
{
    FILE* f = fopen(path.c_str(), "rb");
    if (!f)
        return "";
    Sha256 s;
    char buf[1 << 16];
    size_t n;
    while ((n = fread(buf, 1, sizeof(buf), f)) > 0)
        s.update(buf, n);
    const bool ok = !ferror(f);
    fclose(f);
    return ok ? s.hex() : "";
}

////////////////////////////////////////////////////////////////////////////////
// Cache files

namespace {

const char CACHE_MAGIC[8] = {'S', 'b', 'G', 'r', 'i', 'd', '1', 0};

/*  <cache_dir>/<first 16 hash chars>-<res>.sbgrid  */
std::string cache_path(const std::string& dir, const std::string& sha,
                       float vpu)
{
    char res[32];
    snprintf(res, sizeof(res), "%g", vpu);
    return dir + "/" + sha.substr(0, 16) + "-" + res + ".sbgrid";
}

bool load_cache(const std::string& path, std::vector<float>* data,
                uint32_t dims[3], float bounds[6])
{
    FILE* f = fopen(path.c_str(), "rb");
    if (!f)
        return false;

    char magic[8];
    bool ok = fread(magic, 1, 8, f) == 8 &&
              memcmp(magic, CACHE_MAGIC, 8) == 0 &&
              fread(dims, sizeof(uint32_t), 3, f) == 3 &&
              fread(bounds, sizeof(float), 6, f) == 6;
    if (ok) {
        const size_t count = size_t(dims[0]) * dims[1] * dims[2];
        // Sanity bound: reject absurd headers before allocating
        ok = count > 0 && count <= (size_t(1) << 31);
        if (ok) {
            data->resize(count);
            ok = fread(data->data(), sizeof(float), count, f) == count;
            // Must be exactly at EOF
            ok = ok && fgetc(f) == EOF;
        }
    }
    fclose(f);
    if (!ok)
        data->clear();
    return ok;
}

void save_cache(const std::string& path, const std::vector<float>& data,
                const uint32_t dims[3], const float bounds[6])
{
    // Write-then-rename so a crash never leaves a torn cache file
    const std::string tmp = path + ".tmp";
    FILE* f = fopen(tmp.c_str(), "wb");
    if (!f)
        return;
    bool ok = fwrite(CACHE_MAGIC, 1, 8, f) == 8 &&
              fwrite(dims, sizeof(uint32_t), 3, f) == 3 &&
              fwrite(bounds, sizeof(float), 6, f) == 6 &&
              fwrite(data.data(), sizeof(float), data.size(), f)
                  == data.size();
    ok = (fclose(f) == 0) && ok;
    if (ok)
        rename(tmp.c_str(), path.c_str());
    else
        remove(tmp.c_str());
}

}  // namespace

////////////////////////////////////////////////////////////////////////////////
// Import

ImportResult import_mesh_grid(const std::string& path,
                              float voxels_per_unit,
                              const std::string& cache_dir,
                              int threads)
{
    ImportResult out;

    if (!(voxels_per_unit > 0)) {
        out.error = "voxels_per_unit must be positive";
        return out;
    }

    out.sha256 = sha256_file(path);
    if (out.sha256.empty()) {
        out.error = "could not read '" + path + "'";
        return out;
    }

    char res_str[32];
    snprintf(res_str, sizeof(res_str), "%g", voxels_per_unit);
    const std::string key = "mesh:" + out.sha256 + ":" + res_str;

    // Session-level reuse: same content + parameters, same grid
    const uint32_t existing = grid_find(key);
    if (existing) {
        MeshGrid* g = grid_lookup(existing);
        grid_sampled_bounds(g, &out.bounds[0], &out.bounds[1],
                            &out.bounds[2], &out.bounds[3],
                            &out.bounds[4], &out.bounds[5]);
        out.grid_id = existing;
        out.from_cache = true;
        return out;
    }

    // Cross-session reuse: the sidecar cache file
    std::vector<float> data;
    uint32_t dims[3];
    float bounds[6];
    const std::string cfile = cache_dir.empty() ? "" :
            cache_path(cache_dir, out.sha256, voxels_per_unit);
    if (!cfile.empty() && load_cache(cfile, &data, dims, bounds)) {
        out.grid_id = grid_register(std::move(data),
                dims[0], dims[1], dims[2],
                bounds[0], bounds[1], bounds[2],
                bounds[3], bounds[4], bounds[5], key);
        if (out.grid_id) {
            memcpy(out.bounds, bounds, sizeof(bounds));
            memcpy(out.dims, dims, sizeof(dims));
            out.from_cache = true;
            return out;
        }
        data.clear();  // corrupt or inconsistent: fall through and sample
    }

    // Full import: load, build queries, sample
    TriMesh mesh;
    std::string err;
    if (!load_stl(path, &mesh, &err)) {
        out.error = err;
        return out;
    }
    out.tri_count = mesh.tri_count();
    out.stibium_stamp = mesh.header.compare(0, 7, "Stibium") == 0;

    /*  Pad the sampled box by two voxels so every boundary sample is
     *  positive (outside the surface) and the beyond-the-box distance
     *  extension stays continuous. */
    const float vox = 1.f / voxels_per_unit;
    const float pad = 2 * vox;
    for (int a = 0; a < 3; ++a) {
        bounds[a] = mesh.bbox[a] - pad;
        bounds[a + 3] = mesh.bbox[a + 3] + pad;
    }

    for (int a = 0; a < 3; ++a) {
        const float span = bounds[a + 3] - bounds[a];
        uint32_t n = uint32_t(ceilf(span * voxels_per_unit)) + 1;
        dims[a] = n < 2 ? 2 : n;
    }

    const size_t count = size_t(dims[0]) * dims[1] * dims[2];
    if (count > (size_t(1) << 28)) {   // 268M samples = 1 GiB of floats
        char msg[128];
        snprintf(msg, sizeof(msg),
                 "grid would need %zuM samples (max 268M); "
                 "lower voxels_per_unit", count >> 20);
        out.error = msg;
        return out;
    }

    MeshQuery query(mesh);

    data.resize(count);
    const float ox = bounds[0], oy = bounds[1], oz = bounds[2];
    const float sx = (bounds[3] - bounds[0]) / (dims[0] - 1),
                sy = (bounds[4] - bounds[1]) / (dims[1] - 1),
                sz = (bounds[5] - bounds[2]) / (dims[2] - 1);

    // Parallel over z-planes; queries are const and thread-safe
    unsigned nt = threads > 0 ? threads :
            std::thread::hardware_concurrency() / 2;
    if (nt < 1)
        nt = 1;

    std::atomic<uint32_t> next_k(0);
    auto worker = [&]() {
        for (uint32_t k; (k = next_k.fetch_add(1)) < dims[2]; ) {
            const float z = oz + k * sz;
            float* plane = data.data() + size_t(dims[0]) * dims[1] * k;
            for (uint32_t j = 0; j < dims[1]; ++j) {
                const float y = oy + j * sy;
                float* row = plane + size_t(dims[0]) * j;
                for (uint32_t i = 0; i < dims[0]; ++i)
                    row[i] = query.signed_distance(ox + i * sx, y, z);
            }
        }
    };
    std::vector<std::thread> pool;
    for (unsigned i = 1; i < nt; ++i)
        pool.emplace_back(worker);
    worker();
    for (auto& t : pool)
        t.join();

    if (!cfile.empty())
        save_cache(cfile, data, dims, bounds);

    out.grid_id = grid_register(std::move(data),
            dims[0], dims[1], dims[2],
            bounds[0], bounds[1], bounds[2],
            bounds[3], bounds[4], bounds[5], key);
    if (!out.grid_id) {
        out.error = "grid registration failed";
        return out;
    }
    memcpy(out.bounds, bounds, sizeof(bounds));
    memcpy(out.dims, dims, sizeof(dims));
    return out;
}

}  // namespace fab_mesh

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <sys/stat.h>
#include <vector>

#include "fab/mesh/glyph_import.h"
#include "fab/mesh/mesh_import.h"   // sha256_file
#include "fab/tree/grid.h"
#include "fab/tree/eval.h"          // eval_f (Shape field sampling)

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"
#pragma GCC diagnostic ignored "-Wsign-compare"
#endif
#define STB_TRUETYPE_IMPLEMENTATION
#include "../../vendor/stb/stb_truetype.h"
#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif

namespace fab_mesh {
namespace {

struct Seg { float ax, ay, bx, by; };

void flatten_quad(std::vector<Seg>& segs, float x0, float y0,
                  float cx, float cy, float x1, float y1, int K)
{
    float px = x0, py = y0;
    for (int i = 1; i <= K; i++) {
        float t = (float)i / K, u = 1 - t;
        float x = u * u * x0 + 2 * u * t * cx + t * t * x1;
        float y = u * u * y0 + 2 * u * t * cy + t * t * y1;
        segs.push_back({px, py, x, y});
        px = x; py = y;
    }
}

// On-disk grid cache in ~/.stibium/fonts/, same SbGrid1 format mesh import
// uses -- a baked glyph survives across sessions, so only the first ever use
// pays the sample cost.
const char CACHE_MAGIC[8] = {'S', 'b', 'G', 'r', 'i', 'd', '1', 0};

std::string cache_file(const std::string& key)
{
    const char* home = getenv("HOME");
    if (!home) return "";
    std::string base = std::string(home) + "/.stibium";
    mkdir(base.c_str(), 0777);
    std::string dir = base + "/fonts";
    mkdir(dir.c_str(), 0777);
    std::string name;
    for (char c : key) name += std::isalnum((unsigned char)c) ? c : '_';
    return dir + "/" + name + ".grid";
}

bool load_grid(const std::string& path, std::vector<float>* data,
               uint32_t dims[3], float bounds[6])
{
    FILE* f = fopen(path.c_str(), "rb");
    if (!f) return false;
    char magic[8];
    bool ok = fread(magic, 1, 8, f) == 8 &&
              memcmp(magic, CACHE_MAGIC, 8) == 0 &&
              fread(dims, sizeof(uint32_t), 3, f) == 3 &&
              fread(bounds, sizeof(float), 6, f) == 6;
    if (ok) {
        size_t count = (size_t)dims[0] * dims[1] * dims[2];
        data->resize(count);
        ok = count > 0 && fread(data->data(), sizeof(float), count, f) == count;
    }
    fclose(f);
    return ok;
}

void save_grid(const std::string& path, const std::vector<float>& data,
               const uint32_t dims[3], const float bounds[6])
{
    FILE* f = fopen(path.c_str(), "wb");
    if (!f) return;
    fwrite(CACHE_MAGIC, 1, 8, f);
    fwrite(dims, sizeof(uint32_t), 3, f);
    fwrite(bounds, sizeof(float), 6, f);
    fwrite(data.data(), sizeof(float), data.size(), f);
    fclose(f);
}

}  // namespace

GlyphResult bake_glyph_grid(const std::string& font_path, int codepoint,
                            float cap, float thickness, float vpm)
{
    GlyphResult out;

    FILE* f = fopen(font_path.c_str(), "rb");
    if (!f) { out.error = "could not open font '" + font_path + "'"; return out; }
    fseek(f, 0, SEEK_END); long n = ftell(f); fseek(f, 0, SEEK_SET);
    std::vector<unsigned char> buf(n > 0 ? n : 1);
    size_t got = fread(buf.data(), 1, n, f); fclose(f);
    if ((long)got != n) { out.error = "short read on '" + font_path + "'"; return out; }

    stbtt_fontinfo font;
    if (!stbtt_InitFont(&font, buf.data(),
                        stbtt_GetFontOffsetForIndex(buf.data(), 0))) {
        out.error = "stb_truetype: could not init '" + font_path + "'"; return out;
    }

    // Content key: same font+glyph+params reuses the registered grid.
    std::string sha = sha256_file(font_path);
    char pbuf[128];
    snprintf(pbuf, sizeof(pbuf), "glyph:%s:%d:%g:%g:%g",
             sha.c_str(), codepoint, cap, thickness, vpm);
    std::string key = pbuf;
    if (uint32_t existing = grid_find(key)) {
        MeshGrid* g = grid_lookup(existing);
        grid_sampled_bounds(g, &out.bounds[0], &out.bounds[1], &out.bounds[2],
                            &out.bounds[3], &out.bounds[4], &out.bounds[5]);
        out.grid_id = existing; out.from_cache = true; return out;
    }

    // Cross-session reuse: the ~/.stibium/fonts sidecar.
    const std::string cfile = cache_file(key);
    if (!cfile.empty()) {
        std::vector<float> cdata; uint32_t cdims[3]; float cbnds[6];
        if (load_grid(cfile, &cdata, cdims, cbnds)) {
            out.grid_id = grid_register(std::move(cdata), cdims[0], cdims[1],
                    cdims[2], cbnds[0], cbnds[1], cbnds[2],
                    cbnds[3], cbnds[4], cbnds[5], key);
            if (out.grid_id) {
                memcpy(out.bounds, cbnds, sizeof(cbnds));
                memcpy(out.dims, cdims, sizeof(cdims));
                out.from_cache = true; return out;
            }
        }
    }

    // Cap-height scale: top of 'H' in font units maps to `cap` mm.
    int hx0, hy0, hx1, hy1;
    stbtt_GetCodepointBox(&font, 'H', &hx0, &hy0, &hx1, &hy1);
    if (hy1 <= 0) { out.error = "font has no usable 'H' cap box"; return out; }
    float scale = cap / (float)hy1;

    // Outline -> flattened contour segments.
    stbtt_vertex* v; int nv = stbtt_GetCodepointShape(&font, codepoint, &v);
    if (nv <= 0) { out.error = "glyph has no outline"; return out; }
    std::vector<Seg> segs;
    float sx = 0, sy = 0, cxp = 0, cyp = 0; bool open = false;
    for (int i = 0; i < nv; i++) {
        float X = v[i].x * scale, Y = v[i].y * scale;
        if (v[i].type == STBTT_vmove) {
            if (open) segs.push_back({cxp, cyp, sx, sy});
            sx = X; sy = Y; cxp = X; cyp = Y; open = true;
        } else if (v[i].type == STBTT_vline) {
            segs.push_back({cxp, cyp, X, Y}); cxp = X; cyp = Y;
        } else if (v[i].type == STBTT_vcurve) {
            float CX = v[i].cx * scale, CY = v[i].cy * scale;
            flatten_quad(segs, cxp, cyp, CX, CY, X, Y, 8); cxp = X; cyp = Y;
        }
    }
    if (open) segs.push_back({cxp, cyp, sx, sy});
    stbtt_FreeShape(&font, v);
    if (segs.empty()) { out.error = "glyph flattened to no segments"; return out; }

    // Bounds: xy from the outline, z from the extrusion, padded.
    float x0 = 1e9f, y0 = 1e9f, x1 = -1e9f, y1 = -1e9f;
    for (auto& s : segs) {
        x0 = fminf(x0, fminf(s.ax, s.bx)); y0 = fminf(y0, fminf(s.ay, s.by));
        x1 = fmaxf(x1, fmaxf(s.ax, s.bx)); y1 = fmaxf(y1, fmaxf(s.ay, s.by));
    }
    float pad = 2.0f / vpm + 0.05f * cap;
    x0 -= pad; y0 -= pad; x1 += pad; y1 += pad;
    float z0 = -pad, z1 = thickness + pad;

    // Exact 2D signed distance (negative inside) via segment distance + winding.
    auto sdf2d = [&](float px, float py) -> float {
        float d2 = 1e18f; int wn = 0;
        for (auto& s : segs) {
            float ex = s.bx - s.ax, ey = s.by - s.ay;
            float wx = px - s.ax, wy = py - s.ay;
            float t = (ex * wx + ey * wy) / (ex * ex + ey * ey + 1e-12f);
            t = t < 0 ? 0 : (t > 1 ? 1 : t);
            float qx = wx - t * ex, qy = wy - t * ey, dd = qx * qx + qy * qy;
            if (dd < d2) d2 = dd;
            if (s.ay <= py) { if (s.by > py && (ex * wy - ey * wx) > 0) wn++; }
            else            { if (s.by <= py && (ex * wy - ey * wx) < 0) wn--; }
        }
        float d = sqrtf(d2); return wn != 0 ? -d : d;
    };

    uint32_t ni = (uint32_t)((x1 - x0) * vpm) + 1;
    uint32_t nj = (uint32_t)((y1 - y0) * vpm) + 1;
    uint32_t nk = (uint32_t)((z1 - z0) * vpm) + 1;
    if (ni < 2) ni = 2;
    if (nj < 2) nj = 2;
    if (nk < 2) nk = 2;
    std::vector<float> data((size_t)ni * nj * nk);
    float dx = (x1 - x0) / (ni - 1), dy = (y1 - y0) / (nj - 1),
          dz = (z1 - z0) / (nk - 1);

    // Extruded-solid SDF: combine the 2D profile with the z-slab [0, thickness].
    for (uint32_t k = 0; k < nk; k++) {
        float z = z0 + k * dz;
        float zd = fmaxf(-z, z - thickness);   // signed dist to slab, neg inside
        for (uint32_t j = 0; j < nj; j++) {
            float y = y0 + j * dy;
            for (uint32_t i = 0; i < ni; i++) {
                float d2 = sdf2d(x0 + i * dx, y);
                float ox = fmaxf(d2, 0.0f), oz = fmaxf(zd, 0.0f);
                float outside = sqrtf(ox * ox + oz * oz);
                float inside = fminf(fmaxf(d2, zd), 0.0f);
                data[(size_t)i + ni * ((size_t)j + nj * k)] = outside + inside;
            }
        }
    }

    const uint32_t dims[3] = {ni, nj, nk};
    const float bnds[6] = {x0, y0, z0, x1, y1, z1};
    if (!cfile.empty()) save_grid(cfile, data, dims, bnds);

    out.grid_id = grid_register(std::move(data), ni, nj, nk,
                                x0, y0, z0, x1, y1, z1, key);
    if (!out.grid_id) { out.error = "grid_register failed"; return out; }
    out.bounds[0] = x0; out.bounds[1] = y0; out.bounds[2] = z0;
    out.bounds[3] = x1; out.bounds[4] = y1; out.bounds[5] = z1;
    out.dims[0] = ni; out.dims[1] = nj; out.dims[2] = nk;
    return out;
}

// Advance width of a glyph in mm at the given cap height (for string layout).
float glyph_advance(const std::string& font_path, int codepoint, float cap)
{
    FILE* f = fopen(font_path.c_str(), "rb");
    if (!f) return 0.0f;
    fseek(f, 0, SEEK_END); long n = ftell(f); fseek(f, 0, SEEK_SET);
    std::vector<unsigned char> buf(n > 0 ? n : 1);
    size_t got = fread(buf.data(), 1, n, f); fclose(f);
    if ((long)got != n) return 0.0f;
    stbtt_fontinfo font;
    if (!stbtt_InitFont(&font, buf.data(),
                        stbtt_GetFontOffsetForIndex(buf.data(), 0))) return 0.0f;
    int hx0, hy0, hx1, hy1;
    stbtt_GetCodepointBox(&font, 'H', &hx0, &hy0, &hx1, &hy1);
    if (hy1 <= 0) return 0.0f;
    int adv, lsb;
    stbtt_GetCodepointHMetrics(&font, codepoint, &adv, &lsb);
    return adv * (cap / (float)hy1);
}

// --- Redistance: field sign -> true euclidean distance grid (Felzenszwalb) ---

// 1D squared distance transform of f[] (0 at seeds, +INF elsewhere) into d[].
static void edt1d(const float* f, float* d, int* v, float* z, int n)
{
    int k = 0; v[0] = 0; z[0] = -1e20f; z[1] = 1e20f;
    for (int q = 1; q < n; q++) {
        float s;
        while (true) {
            s = ((f[q] + (float)q * q) - (f[v[k]] + (float)v[k] * v[k]))
                / (2.0f * q - 2.0f * v[k]);
            if (s <= z[k]) k--; else break;
        }
        k++; v[k] = q; z[k] = s; z[k + 1] = 1e20f;
    }
    k = 0;
    for (int q = 0; q < n; q++) {
        while (z[k + 1] < (float)q) k++;
        float dq = (float)(q - v[k]);
        d[q] = dq * dq + f[v[k]];
    }
}

// In-place separable 3D squared EDT over an isotropic grid (index units).
static void edt3d(std::vector<float>& g, uint32_t ni, uint32_t nj, uint32_t nk)
{
    uint32_t maxn = ni; if (nj > maxn) maxn = nj; if (nk > maxn) maxn = nk;
    std::vector<float> f(maxn), d(maxn), z(maxn + 1);
    std::vector<int> v(maxn);
    for (uint32_t k = 0; k < nk; k++)
        for (uint32_t j = 0; j < nj; j++) {
            for (uint32_t i = 0; i < ni; i++)
                f[i] = g[(size_t)i + ni * ((size_t)j + nj * k)];
            edt1d(f.data(), d.data(), v.data(), z.data(), ni);
            for (uint32_t i = 0; i < ni; i++)
                g[(size_t)i + ni * ((size_t)j + nj * k)] = d[i];
        }
    for (uint32_t k = 0; k < nk; k++)
        for (uint32_t i = 0; i < ni; i++) {
            for (uint32_t j = 0; j < nj; j++)
                f[j] = g[(size_t)i + ni * ((size_t)j + nj * k)];
            edt1d(f.data(), d.data(), v.data(), z.data(), nj);
            for (uint32_t j = 0; j < nj; j++)
                g[(size_t)i + ni * ((size_t)j + nj * k)] = d[j];
        }
    for (uint32_t j = 0; j < nj; j++)
        for (uint32_t i = 0; i < ni; i++) {
            for (uint32_t k = 0; k < nk; k++)
                f[k] = g[(size_t)i + ni * ((size_t)j + nj * k)];
            edt1d(f.data(), d.data(), v.data(), z.data(), nk);
            for (uint32_t k = 0; k < nk; k++)
                g[(size_t)i + ni * ((size_t)j + nj * k)] = d[k];
        }
}

GlyphResult redistance_grid(::MathTree_* tree,
                            float xmin, float ymin, float zmin,
                            float xmax, float ymax, float zmax,
                            float vpm, const std::string& key)
{
    GlyphResult out;
    if (uint32_t existing = grid_find(key)) {
        MeshGrid* g = grid_lookup(existing);
        grid_sampled_bounds(g, &out.bounds[0], &out.bounds[1], &out.bounds[2],
                            &out.bounds[3], &out.bounds[4], &out.bounds[5]);
        out.grid_id = existing; out.from_cache = true; return out;
    }
    const std::string cfile = cache_file(key);
    if (!cfile.empty()) {
        std::vector<float> cdata; uint32_t cdims[3]; float cbnds[6];
        if (load_grid(cfile, &cdata, cdims, cbnds)) {
            out.grid_id = grid_register(std::move(cdata), cdims[0], cdims[1],
                    cdims[2], cbnds[0], cbnds[1], cbnds[2],
                    cbnds[3], cbnds[4], cbnds[5], key);
            if (out.grid_id) {
                memcpy(out.bounds, cbnds, sizeof(cbnds));
                memcpy(out.dims, cdims, sizeof(cdims));
                out.from_cache = true; return out;
            }
        }
    }

    float pad = 2.0f / vpm;
    xmin -= pad; ymin -= pad; zmin -= pad;
    xmax += pad; ymax += pad; zmax += pad;
    uint32_t ni = (uint32_t)((xmax - xmin) * vpm) + 1;
    uint32_t nj = (uint32_t)((ymax - ymin) * vpm) + 1;
    uint32_t nk = (uint32_t)((zmax - zmin) * vpm) + 1;
    if (ni < 2) ni = 2;
    if (nj < 2) nj = 2;
    if (nk < 2) nk = 2;
    float dx = (xmax - xmin) / (ni - 1), dy = (ymax - ymin) / (nj - 1),
          dz = (zmax - zmin) / (nk - 1);

    const float INF = 1e20f;
    size_t N = (size_t)ni * nj * nk;
    std::vector<float> dout(N), din(N);
    for (uint32_t k = 0; k < nk; k++) {
        float z = zmin + k * dz;
        for (uint32_t j = 0; j < nj; j++) {
            float y = ymin + j * dy;
            for (uint32_t i = 0; i < ni; i++) {
                bool inside = eval_f(tree, xmin + i * dx, y, z) < 0.0f;
                size_t idx = (size_t)i + ni * ((size_t)j + nj * k);
                dout[idx] = inside ? 0.0f : INF;   // seeds = inside cells
                din[idx]  = inside ? INF : 0.0f;    // seeds = outside cells
            }
        }
    }
    edt3d(dout, ni, nj, nk);
    edt3d(din, ni, nj, nk);
    std::vector<float> data(N);
    float cell = 1.0f / vpm;   // isotropic by construction
    for (size_t idx = 0; idx < N; idx++)
        data[idx] = (sqrtf(dout[idx]) - sqrtf(din[idx])) * cell;

    const uint32_t dims[3] = {ni, nj, nk};
    const float bnds[6] = {xmin, ymin, zmin, xmax, ymax, zmax};
    if (!cfile.empty()) save_grid(cfile, data, dims, bnds);
    out.grid_id = grid_register(std::move(data), ni, nj, nk,
                                xmin, ymin, zmin, xmax, ymax, zmax, key);
    if (!out.grid_id) { out.error = "grid_register failed"; return out; }
    memcpy(out.bounds, bnds, sizeof(bnds));
    out.dims[0] = ni; out.dims[1] = nj; out.dims[2] = nk;
    return out;
}

// --- Marching-squares contour of a 2D field (z=0), sub-cell interpolated ----
// Emits the zero-level as line segments with linearly-interpolated endpoints,
// so the boundary is smooth rather than staircased. This is what lets the
// Antimony font (Matt Keeter's fields) bake to clean walls, matching the TTF
// path's exact-distance quality instead of the EDT's grid striping.
static void contour_field(::MathTree_* tree, float x0, float y0,
                          float x1, float y1, float cres,
                          std::vector<Seg>& segs)
{
    int ci = (int)((x1 - x0) * cres) + 1;
    int cj = (int)((y1 - y0) * cres) + 1;
    if (ci < 2) ci = 2;
    if (cj < 2) cj = 2;
    float dx = (x1 - x0) / (ci - 1), dy = (y1 - y0) / (cj - 1);
    std::vector<float> F((size_t)ci * cj);
    for (int j = 0; j < cj; j++)
        for (int i = 0; i < ci; i++)
            F[(size_t)j * ci + i] = eval_f(tree, x0 + i * dx, y0 + j * dy, 0.0f);

    for (int j = 0; j < cj - 1; j++) {
        for (int i = 0; i < ci - 1; i++) {
            float v00 = F[(size_t)j * ci + i],       v10 = F[(size_t)j * ci + i + 1];
            float v01 = F[(size_t)(j + 1) * ci + i], v11 = F[(size_t)(j + 1) * ci + i + 1];
            float xa = x0 + i * dx, xb = x0 + (i + 1) * dx;
            float ya = y0 + j * dy, yb = y0 + (j + 1) * dy;
            float cx[4], cy[4]; int nc = 0;
            auto edge = [&](float va, float ax, float ay,
                            float vb, float bx, float by) {
                if ((va < 0) != (vb < 0)) {
                    float t = va / (va - vb);
                    cx[nc] = ax + t * (bx - ax);
                    cy[nc] = ay + t * (by - ay);
                    nc++;
                }
            };
            edge(v00, xa, ya, v10, xb, ya);  // bottom
            edge(v10, xb, ya, v11, xb, yb);  // right
            edge(v11, xb, yb, v01, xa, yb);  // top
            edge(v01, xa, yb, v00, xa, ya);  // left
            if (nc == 2) {
                segs.push_back({cx[0], cy[0], cx[1], cy[1]});
            } else if (nc == 4) {          // saddle
                segs.push_back({cx[0], cy[0], cx[1], cy[1]});
                segs.push_back({cx[2], cy[2], cx[3], cy[3]});
            }
        }
    }
}

GlyphResult bake_shape_glyph(::MathTree_* tree, float x0, float y0,
                             float x1, float y1, float thickness,
                             float vpm, const std::string& key)
{
    GlyphResult out;
    if (uint32_t existing = grid_find(key)) {
        MeshGrid* g = grid_lookup(existing);
        grid_sampled_bounds(g, &out.bounds[0], &out.bounds[1], &out.bounds[2],
                            &out.bounds[3], &out.bounds[4], &out.bounds[5]);
        out.grid_id = existing; out.from_cache = true; return out;
    }
    const std::string cfile = cache_file(key);
    if (!cfile.empty()) {
        std::vector<float> cdata; uint32_t cdims[3]; float cbnds[6];
        if (load_grid(cfile, &cdata, cdims, cbnds)) {
            out.grid_id = grid_register(std::move(cdata), cdims[0], cdims[1],
                    cdims[2], cbnds[0], cbnds[1], cbnds[2],
                    cbnds[3], cbnds[4], cbnds[5], key);
            if (out.grid_id) {
                memcpy(out.bounds, cbnds, sizeof(cbnds));
                memcpy(out.dims, cdims, sizeof(cdims));
                out.from_cache = true; return out;
            }
        }
    }

    // Contour a touch finer than the bake grid for accurate crossings.
    std::vector<Seg> segs;
    contour_field(tree, x0, y0, x1, y1, vpm * 1.5f, segs);
    if (segs.empty()) { out.error = "shape has no 2D contour"; return out; }

    float pad = 2.0f / vpm;
    x0 -= pad; y0 -= pad; x1 += pad; y1 += pad;
    float z0 = -pad, z1 = thickness + pad;
    uint32_t ni = (uint32_t)((x1 - x0) * vpm) + 1;
    uint32_t nj = (uint32_t)((y1 - y0) * vpm) + 1;
    uint32_t nk = (uint32_t)((z1 - z0) * vpm) + 1;
    if (ni < 2) ni = 2;
    if (nj < 2) nj = 2;
    if (nk < 2) nk = 2;
    float dx = (x1 - x0) / (ni - 1), dy = (y1 - y0) / (nj - 1),
          dz = (z1 - z0) / (nk - 1);

    // Sign is constant in z (the field is 2D): precompute one XY sign plane.
    std::vector<unsigned char> inside((size_t)ni * nj);
    for (uint32_t j = 0; j < nj; j++)
        for (uint32_t i = 0; i < ni; i++)
            inside[(size_t)j * ni + i] =
                eval_f(tree, x0 + i * dx, y0 + j * dy, 0.0f) < 0.0f ? 1 : 0;

    std::vector<float> data((size_t)ni * nj * nk);
    for (uint32_t j = 0; j < nj; j++) {
        float y = y0 + j * dy;
        for (uint32_t i = 0; i < ni; i++) {
            float x = x0 + i * dx;
            float d2 = 1e18f;
            for (auto& s : segs) {
                float ex = s.bx - s.ax, ey = s.by - s.ay;
                float wx = x - s.ax, wy = y - s.ay;
                float t = (ex * wx + ey * wy) / (ex * ex + ey * ey + 1e-12f);
                t = t < 0 ? 0 : (t > 1 ? 1 : t);
                float qx = wx - t * ex, qy = wy - t * ey, dd = qx * qx + qy * qy;
                if (dd < d2) d2 = dd;
            }
            float sd = sqrtf(d2);
            if (inside[(size_t)j * ni + i]) sd = -sd;   // sign from Matt's field
            for (uint32_t k = 0; k < nk; k++) {
                float z = z0 + k * dz;
                float zd = fmaxf(-z, z - thickness);
                float ox = fmaxf(sd, 0.0f), oz = fmaxf(zd, 0.0f);
                float inS = fminf(fmaxf(sd, zd), 0.0f);
                data[(size_t)i + ni * ((size_t)j + nj * k)] =
                    sqrtf(ox * ox + oz * oz) + inS;
            }
        }
    }

    const uint32_t dims[3] = {ni, nj, nk};
    const float bnds[6] = {x0, y0, z0, x1, y1, z1};
    if (!cfile.empty()) save_grid(cfile, data, dims, bnds);
    out.grid_id = grid_register(std::move(data), ni, nj, nk,
                                x0, y0, z0, x1, y1, z1, key);
    if (!out.grid_id) { out.error = "grid_register failed"; return out; }
    memcpy(out.bounds, bnds, sizeof(bnds));
    out.dims[0] = ni; out.dims[1] = nj; out.dims[2] = nk;
    return out;
}

}  // namespace fab_mesh

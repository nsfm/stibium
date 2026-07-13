#include <atomic>
#include <cmath>
#include <cstring>
#include <mutex>
#include <unordered_map>

#include "fab/tree/grid.h"

/*  Sample block edge for the interval-eval min/max tables.  8 keeps
 *  the tables at ~1/256 of the sample data while letting interval
 *  queries touch only a handful of blocks. */
static const uint32_t BLOCK = 8;

struct MeshGrid_ {
    uint32_t ni, nj, nk;
    float xmin, ymin, zmin;
    float xmax, ymax, zmax;
    float inv_dx, inv_dy, inv_dz;

    std::vector<float> d;               // samples, x-major

    uint32_t bi, bj, bk;                // block counts per axis
    std::vector<float> bmin, bmax;      // per-block sample min/max

    std::atomic<int> refs;
    uint32_t id;
    std::string key;
};

////////////////////////////////////////////////////////////////////////////////
// Registry

namespace {

struct Registry {
    std::mutex mtx;
    std::unordered_map<uint32_t, MeshGrid_*> by_id;
    std::unordered_map<std::string, uint32_t> by_key;
    uint32_t next_id = 1;
};

/*  Leaked intentionally: grids referenced by static Python objects
 *  can be released during interpreter teardown, after a normal
 *  static registry would already be destroyed. */
Registry& registry()
{
    static Registry* r = new Registry;
    return *r;
}

size_t grid_bytes(const MeshGrid_* g)
{
    return (g->d.capacity() + g->bmin.capacity() + g->bmax.capacity())
           * sizeof(float);
}

}  // namespace

uint32_t grid_register(std::vector<float>&& samples,
                       uint32_t ni, uint32_t nj, uint32_t nk,
                       float xmin, float ymin, float zmin,
                       float xmax, float ymax, float zmax,
                       const std::string& key)
{
    if (ni < 2 || nj < 2 || nk < 2 ||
        samples.size() != size_t(ni) * nj * nk ||
        !(xmax > xmin) || !(ymax > ymin) || !(zmax > zmin))
    {
        return 0;
    }

    auto g = new MeshGrid_;
    g->ni = ni;  g->nj = nj;  g->nk = nk;
    g->xmin = xmin;  g->ymin = ymin;  g->zmin = zmin;
    g->xmax = xmax;  g->ymax = ymax;  g->zmax = zmax;
    g->inv_dx = (ni - 1) / (xmax - xmin);
    g->inv_dy = (nj - 1) / (ymax - ymin);
    g->inv_dz = (nk - 1) / (zmax - zmin);
    g->d = std::move(samples);
    g->refs = 1;    // the registry's own reference

    // Per-block sample min/max
    g->bi = (ni + BLOCK - 1) / BLOCK;
    g->bj = (nj + BLOCK - 1) / BLOCK;
    g->bk = (nk + BLOCK - 1) / BLOCK;
    const size_t blocks = size_t(g->bi) * g->bj * g->bk;
    g->bmin.assign(blocks, INFINITY);
    g->bmax.assign(blocks, -INFINITY);
    for (uint32_t k = 0; k < nk; ++k) {
        for (uint32_t j = 0; j < nj; ++j) {
            const size_t row = ni * (j + size_t(nj) * k);
            const size_t brow = size_t(g->bi) *
                    (j / BLOCK + size_t(g->bj) * (k / BLOCK));
            for (uint32_t i = 0; i < ni; ++i) {
                const float v = g->d[row + i];
                const size_t b = brow + i / BLOCK;
                if (v < g->bmin[b])  g->bmin[b] = v;
                if (v > g->bmax[b])  g->bmax[b] = v;
            }
        }
    }

    auto& r = registry();
    std::lock_guard<std::mutex> lock(r.mtx);
    g->id = r.next_id++;
    g->key = key;
    r.by_id[g->id] = g;
    if (!key.empty())
        r.by_key[key] = g->id;
    return g->id;
}

uint32_t grid_find(const std::string& key)
{
    auto& r = registry();
    std::lock_guard<std::mutex> lock(r.mtx);
    const auto it = r.by_key.find(key);
    return it == r.by_key.end() ? 0 : it->second;
}

size_t grid_registry_trim()
{
    auto& r = registry();
    std::lock_guard<std::mutex> lock(r.mtx);
    size_t freed = 0;
    for (auto it = r.by_id.begin(); it != r.by_id.end(); ) {
        MeshGrid_* g = it->second;
        if (g->refs.load() == 1) {
            if (!g->key.empty())
                r.by_key.erase(g->key);
            it = r.by_id.erase(it);
            delete g;
            ++freed;
        } else {
            ++it;
        }
    }
    return freed;
}

size_t grid_memory_bytes(uint32_t id)
{
    auto& r = registry();
    std::lock_guard<std::mutex> lock(r.mtx);
    const auto it = r.by_id.find(id);
    return it == r.by_id.end() ? 0 : grid_bytes(it->second);
}

extern "C" {

MeshGrid* grid_lookup(uint32_t id)
{
    auto& r = registry();
    std::lock_guard<std::mutex> lock(r.mtx);
    const auto it = r.by_id.find(id);
    return it == r.by_id.end() ? NULL : it->second;
}

uint32_t grid_id(const MeshGrid* g)
{
    return g->id;
}

void grid_retain(MeshGrid* g)
{
    g->refs.fetch_add(1, std::memory_order_relaxed);
}

void grid_release(MeshGrid* g)
{
    /*  Grids are freed by grid_registry_trim, never here: the
     *  registry's reference keeps the count above zero until the
     *  registry entry itself is dropped inside the trim lock. */
    g->refs.fetch_sub(1, std::memory_order_release);
}

////////////////////////////////////////////////////////////////////////////////
// Evaluation

/*  Signed per-axis excess beyond the sampled box (0 inside). */
static inline float excess(float v, float lo, float hi)
{
    return v < lo ? v - lo : (v > hi ? v - hi : 0);
}

float grid_eval_f(const MeshGrid* g, float x, float y, float z)
{
    const float ex = excess(x, g->xmin, g->xmax),
                ey = excess(y, g->ymin, g->ymax),
                ez = excess(z, g->zmin, g->zmax);

    float u = (x - ex - g->xmin) * g->inv_dx,
          v = (y - ey - g->ymin) * g->inv_dy,
          w = (z - ez - g->zmin) * g->inv_dz;

    int i = (int)u, j = (int)v, k = (int)w;
    if (i < 0)  i = 0;  else if (i > (int)g->ni - 2)  i = g->ni - 2;
    if (j < 0)  j = 0;  else if (j > (int)g->nj - 2)  j = g->nj - 2;
    if (k < 0)  k = 0;  else if (k > (int)g->nk - 2)  k = g->nk - 2;

    const float fu = u - i, fv = v - j, fw = w - k;

    const float* d = g->d.data();
    const size_t s0 = i + g->ni * (j + size_t(g->nj) * k),
                 sj = g->ni,
                 sk = size_t(g->ni) * g->nj;

    const float c00 = d[s0]           + fu * (d[s0 + 1]           - d[s0]),
                c10 = d[s0 + sj]      + fu * (d[s0 + sj + 1]      - d[s0 + sj]),
                c01 = d[s0 + sk]      + fu * (d[s0 + sk + 1]      - d[s0 + sk]),
                c11 = d[s0 + sj + sk] + fu * (d[s0 + sj + sk + 1] - d[s0 + sj + sk]);
    const float c0 = c00 + fv * (c10 - c00),
                c1 = c01 + fv * (c11 - c01);
    const float t = c0 + fw * (c1 - c0);

    const float e2 = ex * ex + ey * ey + ez * ez;
    return e2 > 0 ? t + sqrtf(e2) : t;
}

void grid_eval_g(const MeshGrid* g, float x, float y, float z,
                 float* value, float* dx, float* dy, float* dz)
{
    const float ex = excess(x, g->xmin, g->xmax),
                ey = excess(y, g->ymin, g->ymax),
                ez = excess(z, g->zmin, g->zmax);

    float u = (x - ex - g->xmin) * g->inv_dx,
          v = (y - ey - g->ymin) * g->inv_dy,
          w = (z - ez - g->zmin) * g->inv_dz;

    int i = (int)u, j = (int)v, k = (int)w;
    if (i < 0)  i = 0;  else if (i > (int)g->ni - 2)  i = g->ni - 2;
    if (j < 0)  j = 0;  else if (j > (int)g->nj - 2)  j = g->nj - 2;
    if (k < 0)  k = 0;  else if (k > (int)g->nk - 2)  k = g->nk - 2;

    const float fu = u - i, fv = v - j, fw = w - k;

    const float* d = g->d.data();
    const size_t s0 = i + g->ni * (j + size_t(g->nj) * k),
                 sj = g->ni,
                 sk = size_t(g->ni) * g->nj;

    const float d000 = d[s0],           d100 = d[s0 + 1],
                d010 = d[s0 + sj],      d110 = d[s0 + sj + 1],
                d001 = d[s0 + sk],      d101 = d[s0 + sk + 1],
                d011 = d[s0 + sj + sk], d111 = d[s0 + sj + sk + 1];

    const float c00 = d000 + fu * (d100 - d000),
                c10 = d010 + fu * (d110 - d010),
                c01 = d001 + fu * (d101 - d001),
                c11 = d011 + fu * (d111 - d011);
    const float c0 = c00 + fv * (c10 - c00),
                c1 = c01 + fv * (c11 - c01);
    *value = c0 + fw * (c1 - c0);

    /*  Analytic trilinear partials; a clamped axis contributes no
     *  interpolation slope (the clamped coordinate is constant). */
    const float gu = ex == 0 ?
        ((d100 - d000) * (1 - fv) + (d110 - d010) * fv) * (1 - fw) +
        ((d101 - d001) * (1 - fv) + (d111 - d011) * fv) * fw : 0;
    const float gv = ey == 0 ?
        (c10 - c00) * (1 - fw) + (c11 - c01) * fw : 0;
    const float gw = ez == 0 ? c1 - c0 : 0;

    *dx = gu * g->inv_dx;
    *dy = gv * g->inv_dy;
    *dz = gw * g->inv_dz;

    const float e2 = ex * ex + ey * ey + ez * ez;
    if (e2 > 0)
    {
        const float e = sqrtf(e2);
        *value += e;
        *dx += ex / e;
        *dy += ey / e;
        *dz += ez / e;
    }
}

/*  Sample-index range [lo, hi] covered by a coordinate interval,
 *  clamped to the grid. */
static inline void index_range(float lo, float hi, float origin,
                               float inv_d, uint32_t n,
                               uint32_t* out_lo, uint32_t* out_hi)
{
    int a = (int)floorf((lo - origin) * inv_d);
    int b = (int)floorf((hi - origin) * inv_d) + 1;
    if (a < 0)  a = 0;  else if (a > (int)n - 1)  a = n - 1;
    if (b < 0)  b = 0;  else if (b > (int)n - 1)  b = n - 1;
    *out_lo = a;
    *out_hi = b;
}

Interval grid_eval_i(const MeshGrid* g, Interval X, Interval Y, Interval Z)
{
    uint32_t i0, i1, j0, j1, k0, k1;
    index_range(X.lower, X.upper, g->xmin, g->inv_dx, g->ni, &i0, &i1);
    index_range(Y.lower, Y.upper, g->ymin, g->inv_dy, g->nj, &j0, &j1);
    index_range(Z.lower, Z.upper, g->zmin, g->inv_dz, g->nk, &k0, &k1);

    /*  Whole-block granularity: covering a partial block with the
     *  full block only widens the result (still conservative). */
    float lo = INFINITY, hi = -INFINITY;
    for (uint32_t kb = k0 / BLOCK; kb <= k1 / BLOCK; ++kb) {
        for (uint32_t jb = j0 / BLOCK; jb <= j1 / BLOCK; ++jb) {
            const size_t row = size_t(g->bi) * (jb + size_t(g->bj) * kb);
            for (uint32_t ib = i0 / BLOCK; ib <= i1 / BLOCK; ++ib) {
                const float bl = g->bmin[row + ib],
                            bh = g->bmax[row + ib];
                if (bl < lo)  lo = bl;
                if (bh > hi)  hi = bh;
            }
        }
    }

    /*  Outside-the-box distance term: value(p) = tri(clamp p) + |e(p)|,
     *  so add [min |e|, max |e|] over the query box per-axis. */
    const float exl = X.upper < g->xmin ? g->xmin - X.upper :
                      (X.lower > g->xmax ? X.lower - g->xmax : 0),
                eyl = Y.upper < g->ymin ? g->ymin - Y.upper :
                      (Y.lower > g->ymax ? Y.lower - g->ymax : 0),
                ezl = Z.upper < g->zmin ? g->zmin - Z.upper :
                      (Z.lower > g->zmax ? Z.lower - g->zmax : 0);
    const float exh = fmaxf(fmaxf(g->xmin - X.lower, X.upper - g->xmax), 0),
                eyh = fmaxf(fmaxf(g->ymin - Y.lower, Y.upper - g->ymax), 0),
                ezh = fmaxf(fmaxf(g->zmin - Z.lower, Z.upper - g->zmax), 0);

    Interval out;
    out.lower = lo + sqrtf(exl * exl + eyl * eyl + ezl * ezl);
    out.upper = hi + sqrtf(exh * exh + eyh * eyh + ezh * ezh);
    return out;
}

void grid_sampled_bounds(const MeshGrid* g,
                         float* xmin, float* ymin, float* zmin,
                         float* xmax, float* ymax, float* zmax)
{
    *xmin = g->xmin;  *ymin = g->ymin;  *zmin = g->zmin;
    *xmax = g->xmax;  *ymax = g->ymax;  *zmax = g->zmax;
}

}  // extern "C"

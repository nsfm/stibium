#include <atomic>
#include <cstdlib>
#include <thread>
#include <vector>

#include "fab/tree/render_mt.h"

#include "fab/tree/render.h"
#include "fab/tree/tree.h"

namespace {

int resolve_threads(int threads)
{
    if (const char* env = getenv("STIBIUM_RENDER_THREADS"))
    {
        const int t = atoi(env);
        if (t > 0)
            return t;
    }
    if (threads > 0)
        return threads;

    // Physical cores: this workload is dense FP dependency chains,
    // and hyperthread siblings measurably fight over the FP units
    // (same finding as the contour tracer).
    const unsigned hw = std::thread::hardware_concurrency();
    return hw > 1 ? int(hw / 2) : 1;
}

/*
 *  Splits the region into xy-disjoint chunks: repeatedly bisect the
 *  largest chunk along its longer xy axis.  z is never split, so
 *  every chunk writes its own image cells.
 */
std::vector<Region> split_xy(const Region& r, size_t target)
{
    std::vector<Region> out = {r};
    while (out.size() < target)
    {
        // Find the largest splittable chunk
        size_t best = 0;
        uint64_t best_area = 0;
        for (size_t i = 0; i < out.size(); ++i)
        {
            const uint64_t area = uint64_t(out[i].ni) * out[i].nj;
            if (area > best_area && (out[i].ni > 1 || out[i].nj > 1))
            {
                best_area = area;
                best = i;
            }
        }
        if (best_area == 0)
            break;

        Region a, b;
        if (out[best].ni > out[best].nj)
            bisect_x(out[best], &a, &b);
        else
            bisect_y(out[best], &a, &b);
        out[best] = a;
        out.push_back(b);
    }
    return out;
}

template <typename F>
void run_chunked(MathTree* tree, const Region& r, volatile int* halt,
                 int threads, F&& render_chunk)
{
    threads = resolve_threads(threads);
    threads = int(std::min<uint64_t>(threads,
                                     uint64_t(r.ni) * r.nj / 64 + 1));

    if (threads < 2)
    {
        render_chunk(tree, r);
        return;
    }

    // More chunks than workers, so an empty chunk doesn't idle a core
    const auto chunks = split_xy(r, size_t(threads) * 3);

    // clone_tree writes bookkeeping into the source; clone serially
    std::vector<MathTree*> trees(threads);
    for (auto& t : trees)
        t = clone_tree(tree);

    std::atomic<size_t> next(0);
    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (int t = 0; t < threads; ++t)
    {
        pool.emplace_back([&, t]() {
            size_t c;
            while ((c = next.fetch_add(1)) < chunks.size() && !*halt)
                render_chunk(trees[t], chunks[c]);
        });
    }
    for (auto& th : pool)
        th.join();
    for (auto& t : trees)
        free_tree(t);
}

}  // namespace

void render16_mt(MathTree* tree, Region r,
                 uint16_t** img, volatile int* halt, int threads)
{
    run_chunked(tree, r, halt, threads,
                [=](MathTree* t, const Region& chunk) {
                    render16(t, chunk, img, halt, nullptr);
                });
}

void shaded8_mt(MathTree* tree, Region r,
                uint16_t** depth, uint8_t (**out)[3], volatile int* halt,
                int threads)
{
    run_chunked(tree, r, halt, threads,
                [=](MathTree* t, const Region& chunk) {
                    shaded8(t, chunk, depth, out, halt, nullptr);
                });
}

#include <algorithm>
#include <cstring>
#include <thread>
#include <vector>

#include "fab/tree/parallel_eval.h"

#include "fab/tree/tree.h"
#include "fab/tree/eval.h"
#include "fab/util/region.h"
#include "fab/util/switches.h"

namespace {

/*
 *  Evaluates the field at n points, in chunks of MIN_VOLUME (the
 *  per-node result buffers' capacity).
 */
void eval_span(MathTree* tree,
               const float* xs, const float* ys, const float* zs,
               float* out, size_t n, volatile int* halt,
               std::atomic<uint64_t>* progress)
{
    size_t done = 0;
    while (done < n && !*halt)
    {
        const size_t count = std::min(size_t(MIN_VOLUME) - 1, n - done);
        Region dummy = {};
        dummy.X = const_cast<float*>(xs) + done;
        dummy.Y = const_cast<float*>(ys) + done;
        dummy.Z = const_cast<float*>(zs) + done;
        dummy.voxels = count;

        const float* result = eval_r(tree, dummy);
        memcpy(out + done, result, count * sizeof(float));
        done += count;
        if (progress)
            progress->fetch_add(count, std::memory_order_relaxed);
    }
}

/*
 *  eval_span across a pool of threads, each on its own tree clone.
 */
}  // namespace

void parallel_eval(MathTree* tree,
                   const float* xs, const float* ys, const float* zs,
                   float* out, size_t n, int threads, volatile int* halt,
                   std::atomic<uint64_t>* progress)
{
    if (threads <= 0)
    {
        // Default to physical-core count (hardware_concurrency counts
        // hyperthreads): this workload is dense FP dependency chains,
        // and measured on an 8c/16t machine, 16 threads spent 6x the
        // CPU of 8 threads for 3x worse wall time.
        const unsigned hw = std::thread::hardware_concurrency();
        threads = hw > 1 ? int(hw / 2) : 1;
    }
    threads = int(std::min(size_t(threads), n / MIN_VOLUME + 1));

    if (threads < 2)
    {
        eval_span(tree, xs, ys, zs, out, n, halt, progress);
        return;
    }

    // clone_tree writes bookkeeping into the source, so clone serially
    std::vector<MathTree*> trees(threads);
    for (auto& t : trees)
        t = clone_tree(tree);

    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (int t = 0; t < threads; ++t)
    {
        const size_t lo = n * t / threads;
        const size_t hi = n * (t + 1) / threads;
        pool.emplace_back([=, &trees]() {
            eval_span(trees[t], xs + lo, ys + lo, zs + lo,
                      out + lo, hi - lo, halt, progress);
        });
    }
    for (auto& th : pool)
        th.join();
    for (auto& t : trees)
        free_tree(t);
}



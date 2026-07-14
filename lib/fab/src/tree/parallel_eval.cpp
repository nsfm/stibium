#include <algorithm>
#include <cstring>
#include <thread>
#include <vector>

#include "fab/tree/parallel_eval.h"

#include "fab/tree/tree.h"
#include "fab/tree/tape.h"
#include "fab/util/region.h"
#include "fab/util/switches.h"

namespace {

/*
 *  Evaluates the field at n points, in chunks of MIN_VOLUME (the
 *  per-slot result rows' capacity).
 */
void eval_span(const Tape* tape, TapeCtx* ctx,
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

        const float* result = tape_eval_r(tape, ctx, dummy);
        memcpy(out + done, result, count * sizeof(float));
        done += count;
        if (progress)
            progress->fetch_add(count, std::memory_order_relaxed);
    }
}

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

    // One immutable deck serves every thread; each worker only needs
    // its own workspace (no per-thread tree clones).
    Deck* deck = deck_from_tree(tree);
    const Tape* tape = deck_base(deck);

    if (threads < 2)
    {
        TapeCtx* ctx = tape_ctx_new(deck);
        eval_span(tape, ctx, xs, ys, zs, out, n, halt, progress);
        tape_ctx_free(ctx);
        deck_free(deck);
        return;
    }

    std::vector<std::thread> pool;
    pool.reserve(threads);
    for (int t = 0; t < threads; ++t)
    {
        const size_t lo = n * t / threads;
        const size_t hi = n * (t + 1) / threads;
        pool.emplace_back([=]() {
            TapeCtx* ctx = tape_ctx_new(deck);
            eval_span(tape, ctx, xs + lo, ys + lo, zs + lo,
                      out + lo, hi - lo, halt, progress);
            tape_ctx_free(ctx);
        });
    }
    for (auto& th : pool)
        th.join();
    deck_free(deck);
}

#ifndef PARALLEL_EVAL_H
#define PARALLEL_EVAL_H

#include <atomic>
#include <cstdint>

struct MathTree_;

/*
 *  Evaluates the field at n points (parallel coordinate arrays),
 *  spread across a pool of threads, each working on its own clone
 *  of the tree.  threads <= 0 means physical cores.
 *
 *  progress, if given, is incremented per evaluated point.
 */
void parallel_eval(struct MathTree_* tree,
                   const float* xs, const float* ys, const float* zs,
                   float* out, size_t n, int threads, volatile int* halt,
                   std::atomic<uint64_t>* progress = nullptr);

#endif

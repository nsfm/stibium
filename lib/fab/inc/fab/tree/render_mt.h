#ifndef RENDER_MT_H
#define RENDER_MT_H

#include <cstdint>

#include "fab/util/region.h"
#include "fab/tree/tape.h"

/*
 *  Multithreaded render16 / shaded8: the region is split into
 *  xy-disjoint chunks (never z, so each chunk owns its image cells
 *  and depth stays exact) rendered by a pool of workers sharing the
 *  given immutable deck (compile it once per shape - see
 *  Shape::deck - not per frame), each with its own workspace.
 *
 *  threads <= 0 means physical cores (hardware_concurrency / 2);
 *  the STIBIUM_RENDER_THREADS environment variable overrides.
 */
void render16_mt(const Deck* deck, Region r,
                 uint16_t** img, volatile int* halt, int threads=-1);

void shaded8_mt(const Deck* deck, Region r,
                uint16_t** depth, uint8_t (**out)[3], volatile int* halt,
                int threads=-1);

#endif

#ifndef FAB_GPU_H
#define FAB_GPU_H

#include <stdbool.h>
#include <stdint.h>

#include "fab/util/region.h"
#include "fab/tree/tape.h"

#ifdef __cplusplus
extern "C" {
#endif

/*  GPU depth renderer: the MPR pipeline (interval tile classify ->
 *  on-device tape simplify -> per-pixel march on the shortened tape)
 *  as three GL 4.3 compute passes over the deck's exported bytecode.
 *  See doc/TAPE-DESIGN.md Rounds 5-6.
 *
 *  Pixels are max-merged into img exactly like render16_tape, so
 *  callers can layer multiple shapes.  Returns false - and touches
 *  nothing - when the GPU path can't run this deck (no EGL/GL 4.3
 *  stack, or the tape references host-side OP_GRID payloads);
 *  callers fall back to the CPU path.
 *
 *  Thread-safe via an internal mutex (one GL context, rebound to the
 *  calling thread per render).  Compiled programs are cached per
 *  (deck bytecode, region dims).  */
bool gpu_render16(const Deck* deck, Region r,
                  uint16_t** img, volatile int* halt);

/*  True if a usable EGL + GL 4.3 compute stack was found (the first
 *  call initializes it).  Diagnostic; gpu_render16 already guards. */
bool gpu_available(void);

/*  Renderer name string for diagnostics ("" when unavailable).  */
const char* gpu_renderer_name(void);

#ifdef __cplusplus
}
#endif

#endif

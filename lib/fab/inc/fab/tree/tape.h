#ifndef TAPE_H
#define TAPE_H

#include <stdbool.h>
#include <stdint.h>

#include "fab/util/interval.h"
#include "fab/util/region.h"
#include "fab/tree/math/math_g.h"

#ifdef __cplusplus
extern "C" {
#endif

struct MathTree_;

/*  Immutable evaluation machinery compiled from a MathTree
 *  (see doc/TAPE-DESIGN.md).
 *
 *  A Deck owns the base Tape and the slot layout.  Tapes are
 *  immutable and refcounted: tape_push() specializes a tape for an
 *  interval region (pruning decided min/max branches) and returns a
 *  new reference; tape_release() drops one.  Because tapes never
 *  mutate, one Deck serves any number of threads, each with its own
 *  TapeCtx workspace - no per-thread tree clones.
 */
typedef struct Deck_    Deck;
typedef struct Tape_    Tape;
typedef struct TapeCtx_ TapeCtx;

/** @brief Compiles a deck from a math tree.
    @details The tree is only read (never mutated) and may be freed
    afterwards; the deck retains any OP_GRID payloads it references. */
Deck* deck_from_tree(struct MathTree_* tree);

void deck_free(Deck* deck);

/** @brief The full (unspecialized) tape.  Borrowed reference:
    call tape_retain if it must outlive the deck. */
Tape* deck_base(const Deck* deck);

/** @brief Number of value slots (= workspace entries per eval mode). */
unsigned deck_slots(const Deck* deck);

/** @brief Same bit mask as active_axes(): (x << 2) | (y << 1) | z. */
uint8_t deck_active_axes(const Deck* deck);

/*  Per-thread workspace: slot-indexed result arrays for all four
 *  eval modes plus push scratch.  Slot ids are stable across pushes,
 *  so one ctx serves a whole tape stack. */
TapeCtx* tape_ctx_new(const Deck* deck);
void tape_ctx_free(TapeCtx* ctx);

/*  Evaluation.  Mirrors eval_f/i/r/g: same math_* primitives in the
 *  same order, so results are bit-identical to the MathTree walk.
 *  _r and _g read r.X/Y/Z[0..r.voxels-1]; r.voxels must be
 *  <= MIN_VOLUME for _r and <= MIN_VOLUME/4 for _g (the g rows alias
 *  the r rows, exactly like Results.r today). */
float tape_eval_f(const Tape* tape, TapeCtx* ctx,
                  float x, float y, float z);
Interval tape_eval_i(const Tape* tape, TapeCtx* ctx,
                     Interval X, Interval Y, Interval Z);

/** @brief Batched interval evaluation: classifies `count` boxes in
    ONE pass down the tape (the MPR trick that makes wide tile
    fan-outs pay - tape traversal is amortized and the simple ops
    run across SIMD lanes).  Bounds are exactly equal to `count`
    scalar tape_eval_i calls: hot ops are elementwise mirrors of
    math_i and everything branchy delegates to it per lane.
    count <= TAPE_BATCH.  Does NOT record push records (ClauseIv) -
    ambiguous boxes still take a scalar tape_eval_i before their
    push, exactly like today.  */
#define TAPE_BATCH  64
void tape_eval_i_batch(const Tape* tape, TapeCtx* ctx,
                       const Interval* Xs, const Interval* Ys,
                       const Interval* Zs, Interval* out, int count);
const float* tape_eval_r(const Tape* tape, TapeCtx* ctx, Region r);
const derivative* tape_eval_g(const Tape* tape, TapeCtx* ctx, Region r);

/** @brief Pruning modes for tape_push.
    @details STANDARD drops decided min/max branches (readers are
    remapped to the surviving input - values are exact).  BINARY
    additionally replaces sign-determined subtrees in boolean
    (min/max/neg) context with their interval upper bound, matching
    disable_nodes_binary; only valid when the caller consumes signs,
    not distances (i.e. the raster renderer). */
typedef enum TapePushMode_ {
    TAPE_PUSH_STANDARD,
    TAPE_PUSH_BINARY,
} TapePushMode;

/** @brief Specializes a tape against the interval results of the most
    recent tape_eval_i on this ctx (which must have been over X/Y/Z).
    @returns A new reference: either a shorter tape whose parent is
    the input (bounds recorded for tape_base), or the input itself
    retained when nothing could be pruned.  Pair with tape_release. */
Tape* tape_push(Tape* tape, TapeCtx* ctx,
                Interval X, Interval Y, Interval Z, TapePushMode mode);

Tape* tape_retain(Tape* tape);
void tape_release(Tape* tape);

/** @brief Walks up the parent chain until the point lies inside the
    tape's recorded interval bounds (or the base is reached).
    Borrowed reference - useful when evaluating points that may
    escape the region a tape was pushed for. */
Tape* tape_base_for_point(Tape* tape, float x, float y, float z);

/** @brief Like tape_base_for_point, for an axis-aligned box: walks up
    until the whole box lies inside the tape's recorded bounds. */
Tape* tape_base_for_region(Tape* tape,
                           float xmin, float xmax,
                           float ymin, float ymax,
                           float zmin, float zmax);

/*  Introspection (tests / diagnostics) */
unsigned tape_length(const Tape* tape);
bool tape_is_terminal(const Tape* tape);

/** @brief Dumps clauses (and, with a ctx, the ctx's interval results
    and maybe-NaN taint per slot) to stderr.  Diagnostics only. */
void tape_dump(const Tape* tape, const TapeCtx* ctx);

#ifdef __cplusplus
}
#endif

#endif

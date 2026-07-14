# Tape campaign: the next moves

*Written 2026-07-13, end of the keystone night, specifically to
survive a context compact.  Companion to [TAPE-DESIGN.md](TAPE-DESIGN.md)
(design + results so far) and the Fidget re-recon in [TODO.md](../TODO.md).
Everything here is implementation-ready: anchors, gotchas, and the
numbers to beat are all recorded.*

## State of the world (what's already true)

- Immutable shared tapes replaced disable_nodes everywhere; register-
  allocated decks (pinned consts/axes + linear-scan registers);
  pushes emit OP_COPY, min/max choices + per-clause bounds/taint are
  recorded in `ClauseIv` during `tape_eval_i` (push REQUIRES a prior
  eval_i on the SAME tape).  STANDARD push = exact values (mesher);
  BINARY = sign-preserving only (raster renderer).
- `tape_eval_i_batch` (tape.cpp): 64 boxes per pass, bounds exactly
  equal to scalar (hot ops mirror math_i elementwise, branchy ops
  delegate per lane).  Lives in the r-rows: 64 lo + 64 hi = one
  MIN_VOLUME row, zero allocation.  **Parity rule is load-bearing**:
  fill decisions must never depend on which evaluator ran.
- The depth pass fans 64-wide + batch-classifies (default ON;
  `STIBIUM_TILE_RENDER=0` restores bisection, `=N` tunes).  The
  shaded pass pushes a STANDARD tape per 64px tile (`SHADE_TILE`)
  against the depth buffer's z-range.
- Knobs: `STIBIUM_TAPE_STATS=1` (shrinkage curve + deck/affine
  stats), `STIBIUM_AFFINE=1` (opt-in transform-chain folding),
  `STIBIUM_RENDER_THREADS`, `STIBIUM_FUZZ_COUNT`/`_SEED`.
- Verification kit: `[tape]` suite, `[fuzzer]` (run 100k after ANY
  tape/interval change), `[golden]` dumps (bit-identical expected),
  ground truth = rebuild with `PRUNE 0` in switches.h and cmp
  renders, bench = `python3 $SP/peak.py <cmd>`.  showcase_gear
  renders are compared against unpruned ground truth (re-baselined
  2026-07-13; the OLD render had a 5,842-px pruning artifact).
  Never golden a multithreaded detect-features STL (seam swaps are
  scheduling-dependent, ±2 triangles in ~1M).

Numbers to beat (8c/16t):

| workload | current champ |
|---|---|
| merged Zeiss --render 512px | ~3.3 s (was 3.5) |
| merged Zeiss --render 2048px | **12.36 s** (was 14.97) |
| merged Zeiss --export r7, 1/4 scale | ~19 s |
| merged Zeiss --export r7, full scale | ~310 s / 14.4 GB peak |
| gear.sb --export r60 detect | ~1.6 s |

## 1. Mesher batching (biggest CPU win on the board)

The mesher still classifies octants one scalar `tape_eval_i` at a
time (`Mesher::triangulate_region(const Region&, Tape*)`,
triangulate/mesher.cpp).  Exports are the heavy workload; give the
recursion the same treatment the renderer got:

- Above the packed-block threshold (`!has_data` and
  `voxels >= MIN_VOLUME * 2`, mirroring render.c), split into up to
  `TAPE_BATCH` children (`split()` from region.h, like render.c) and
  `tape_eval_i_batch` them on the CURRENT tape.  Children proven
  empty/full are skipped wholesale; ambiguous children recurse into
  the existing scalar eval_i + push path (which keeps ClauseIv
  coherent for the push - do NOT try to push straight off batch
  results).
- **Gotchas:** (a) progress accounting - skipped children must
  `progress->fetch_add(child.voxels)` exactly like the scalar skip
  branch, or the export progress bar lies; (b) keep the `has_data`
  guard - never fan out or push inside a packed block; (c) the skip
  test is `lower > 0 || upper < 0` (note: mesher skips FULL regions
  too, unlike the renderer - copy the existing condition exactly);
  (d) `load_packed` must still run on the tape that covers the block
  (`block_tape` borrow feeds `eval_zero_crossings`/`get_normals`).
- **Verify:** goldens bit-identical (batch bounds are exact ⇒ same
  skips ⇒ same mesh), then bench gear r60 + quarter-scale merged +
  (patience permitting) full-scale merged.

## 2. Per-tile work queue in render_mt

`run_chunked` (render_mt.cpp) splits xy into `threads * 3` chunks -
far coarser than MPR's tile queue.  A worker that lands on the
detailed corner of the Zeiss idles everyone else.  Replace with a
fine tile queue (e.g. 64-128 px tiles, same `split_xy`, count =
area/tile²), workers pulling via the existing atomic counter.  Tiles
are xy-disjoint so no image contention; keep one `TapeCtx` per
worker (they're the expensive part).  Measure repeat renders in one
process (GUI-shaped load), not one-shot CLI.  Possible refinement:
per-tile decks are unnecessary, but per-tile ROOT pushes (against
the tile's frustum, before descending) would let workers start from
a pre-shrunk tape - test separately.

## 3. Autovectorizer audit (cheap, do first)

Confirm the batch kernels actually vectorize:
`ninja -v` the tape.cpp compile line, re-run with
`-fopt-info-vec-missed` (gcc) and look at `batch_add/sub/mul/min/max`.
If blocked on aliasing, `__restrict__` the row pointers (safe: the
locals-then-store pattern already handles the R==A case, but the
compiler can't see that) or `#pragma GCC ivdep`.  Then `perf record`
a 2048px merged render and confirm where time actually goes now -
candidates: batch kernels, leaf `eval_r`, the qsort, `split()`.

## 4. The dowry (Fidget re-recon, clone at ~/code/fidget, MPL-2.0)

Port order:

1. **Register spilling** (small, unblocks GPU): adopt Load/Store ops
   with LRU eviction from `fidget-core/src/compiler/alloc.rs` +
   `lru.rs` + `op.rs` (`Load(u8, u32)` / `Store(u8, u32)`).  Only
   needed when targeting a bounded register file - the CPU path can
   keep unbounded registers.  Merged Zeiss currently needs 465;
   GPU-class caps are ~128.
2. **GPU renderer** (the big one, now a port not research):
   `fidget-wgpu/src/shaders/*.wgsl` is the complete MPR pipeline -
   study order: `common.wgsl` → `interval_tiles.wgsl` (batch
   classify, our tape_eval_i_batch on device) → `tape_simplify.wgsl`
   (pushing on device - the part we didn't know how to structure) →
   `voxel_tiles.wgsl` + `normals.wgsl` → the sort/merge/repack glue;
   `voxel.rs` shows host-side orchestration.  Our packed-clause
   export should follow `fidget-bytecode` (u32 words, forward eval
   order, `iter_ops` introspection - format explicitly unstable, so
   version-check).  Prototype headless first: GPU-render a deck to
   PNG, cmp against the CPU render before touching the viewport.
   WGSL runs on Vulkan/Metal/DX/WebGPU - this same work feeds the
   WASM web-gallery moonshot.
3. **fidget-solver** (`solve()` in fidget-solver/src/lib.rs):
   constraint solving over tape gradients - the Tier 3
   differentiable-CAD item.  We already have gradient eval; the
   solver is the driver.  Later, after the GPU rung.

## 5. Deferred / parked

- Affine collapse default: flip `STIBIUM_AFFINE` on only when a
  workload shows base-tape evals dominating (the tile renderer's
  coarse passes are the candidate); it changes float rounding, so
  Nate makes the fidelity call.
- Interval-derived bounds (recon item, M) - unrelated to rendering,
  removes the per-shape bounds tax.
- Keeter's 2026-07-03 "Please Steal my Meshing Algorithm Idea"
  (QEF points + Delaunay) - the meshing successor, separate track.

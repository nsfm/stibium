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
| merged Zeiss --render 2048px | **7.6 s** (round 4; was 14.97 pre-tape-batch) |
| merged Zeiss --export r7, 1/4 scale | ~18 s |
| merged Zeiss --export r7, full scale | ~310 s / 14.4 GB peak |
| gear.sb --export r60 detect | ~2.25 s (same-machine A/B; old "1.6" was measured differently) |

## 1. Mesher: batching measured FLAT, dedup cache landed (2026-07-14)

**Octant batching was built, measured, and reverted.**  The spec'd
batch-classify (8 octants per `tape_eval_i_batch` in
`triangulate_region`, exact-parity skips, progress mirrored) was
bit-identical and wall-neutral: gear 2.2 s → 2.2 s, quarter Zeiss
18.0 → 18.1.  Why it can't win what the renderer won: the mesher
recursion tracks the surface, so most children of any recursed
region are ambiguous — they still need their own scalar eval_i for
the push, and the batch pass just duplicates classification.  The
renderer won because view tiles skip in huge swaths.  Don't rebuild
this without new evidence.

**What the export profile actually said** (gdb-sampled, quarter
Zeiss r7): ~30% `interpolate_between` (fixed, below), ~18% scalar
eval_i (irreducible while pushes need ClauseIv), ~14% vertex-intern
hashtable, visible `std::list<InterpolateCommand>` node churn, rest
geometry/IO.  The export lever is the GEOMETRY pipeline, not
evaluation.

**Landed: interp dedup cache.**  `interpolate_between` walked the
whole queue per edge with an Eigen compare (O(edges × queue)) —
now an unordered_map keyed on the bit-packed vertex pair.  TRAP
THAT BIT ONCE ALREADY: `Vec3f` is `Eigen::Vector3d` — doubles, 8
bytes; packing 4-byte words collides and corrupts fan topology
(check_feature spins forever).  Semantics mirrored exactly: -0
normalized to +0, canonical pair order, NaN never enters.
`INTERP_CACHE_DIFF` (compile-time define in mesher.cpp) re-arms a
differential harness that runs old walk + new cache side by side
and aborts on divergence — full suite ran clean under it.
gear r60 detect 2.25 → 2.05 s, quarter Zeiss r7 18.0 → 15.7 s.

**Remaining export candidates, by profile share:** vertex-intern
hash (~14%; FNV + std::unordered_map — try a flat open-addressing
table), `std::list` queue → vector (node allocs visible),
`check_feature`/`get_contour` Eigen-double geometry (Vector3d
throughout the fan logic — precision audit before touching).

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

## 3. Autovectorizer audit — DONE 2026-07-13 (TAPE-DESIGN Round 4)

Found something much bigger than missed vectorization: `min_f`/
`max_f` were libm `fmin`/`fmax` PLT calls — ~40% of render wall
time.  Now bit-exact inline replicas (see Round 4 for the measured
platform semantics: ties incl ±0 → B, NaN-vs-NaN → A).  ivdep on the
batch kernels, `-fno-trapping-math` scoped to SbFab.  2048 px:
11.8 → 7.6 s, everything bit-identical.

Leftovers if profiles ever point here again: `batch_mul` still
scalar (SSE2 can't if-convert the 4-product NaN-select chain; fix =
`-march=x86-64-v2` or explicit reformulation); no `-march` in the
build at all, so everything is SSE2-width — a `-march` bump doubles
vector width but is Nate's portability call.  Profiling recipe that
works here (no perf/valgrind on the box): LD_PRELOAD a
`prctl(PR_SET_PTRACER_ANY)` shim (ptrace_scope=1), then loop
`gdb -batch -ex "thread apply all bt 8" -p PID` — scratchpad
`pmp.sh`/`traceme.c`, aggregate with awk.  Post-round-4 profile
guess for the next lever: leaf `eval_r` transcendentals + the
remaining batch-kernel share; re-profile before believing anything.

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

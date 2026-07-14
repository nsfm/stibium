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

## 2. Per-tile work queue — MEASURED NEGATIVE, dropped (2026-07-14)

Built it (STIBIUM_MT_TILE knob, count = area/tile², same atomic
queue) and swept the grain on the 2048 px merged Zeiss:

| chunks | wall |
|---|---|
| threads×3 = 24 (status quo) | 7.8-8.7 s |
| 64 (256 px tiles) | 8.2 s |
| 256 (128 px tiles) | 11.1-11.4 s |

Monotone: finer = slower, 40% worse at MPR-style grain.  Every
chunk pays its own descent from the base tape, and the top-of-tree
pushes are exactly what fat chunks amortize — the recursion INSIDE
a chunk already is the shared-tape hierarchy that makes fine tile
queues viable in MPR (where parent tiles hand pruned tapes to
children across the whole frame).  The existing atomic-counter
queue over threads×3 chunks sits within noise of optimal on this
workload; imbalance was a theory, the descent tax is a measurement.
Don't rebuild without sharing pushed tapes ACROSS chunk boundaries
(that's the GPU pipeline's breadth-first structure, not a CPU
retrofit).

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

**Stage A LANDED (2026-07-14, TAPE-DESIGN Round 5):**
`tape_export_blob` (u32 bytecode, format v1, OP_GRID refused) +
headless EGL/GL4.3 compute interpreter in tests/gpu.cpp.  Hidden
tags: `SbFabTest "[.gpu]"` (referee - ZERO mismatched pixels vs the
CPU renderer on both the Intel iGPU and the 1650 Ti; prime-run
drives EGL directly) and `"[.gpubench]"` (STIBIUM_GPU_BENCH_N=512
etc).  Brute force on the 1650 Ti is within ~2x of the
fully-armed single CPU thread - the pruning arsenal is the 20-100x,
and porting IT is stage B.  Gotchas already solved in gpu.cpp:
GLSL min()/max() have undefined NaN behavior (hand-mirror math_f),
bake opcode values + blob offsets into generated source, old-catch
hidden tags need the dot (`"[.gpu]"` not `"[gpu]"`).

Port order from here:

1. **Register spilling - NOW THE CRITICAL PATH** (Round 7 measured
   it: the Zeiss renders bit-perfectly via STIBIUM_GPU=1 but takes
   MINUTES vs the CPU's 3.1s at 512px, because its 877 slots make a
   3.5KB per-thread array and occupancy collapses; the 106-slot
   gear is fast).  Adopt Load/Store ops with LRU eviction from
   `fidget-core/src/compiler/alloc.rs` + `lru.rs` + `op.rs`
   (`Load(u8, u32)` / `Store(u8, u32)`).  Only needed for a bounded
   register file - the CPU path keeps unbounded registers.  Target
   ~64-128 registers + a spill buffer in the tape blob.
2. **Stage B - LANDED (2026-07-14, TAPE-DESIGN Round 6)**: interval
   classify + on-device tape simplify + shortened-tape march, one
   subdivision level, in tests/gpu.cpp (`[.gpu2]` referee - zero
   mismatches on both GPUs; `[.gpu2bench]`).  Beats the CPU thread
   at 512³ already.  Known gaps = **stage C**: (a) z-slab / 3D tile
   classification + a second subdivision level (CPU's z-octree
   still wins at 1024³); (b) run on a big deck - needs register
   spilling (item 1) since Zeiss wants 465 slots; (c) per-tile tape
   memory is a fixed n_clauses-sized slot - fidget allocates via
   atomic counters, port that before big decks; (d) normals/shading
   pass; (e) viewport integration (the Qt viewport is already GL).
   Gotchas recorded in gpu.cpp comments: Mesa miscompiles
   interleaved SSBO read/write loops (locals-then-store), GLSL
   min()/max() NaN-undefined, `precise` on interval endpoints,
   conservative-is-sound means exotic ops can just widen to
   [-inf,inf]+taint.  WGSL/GLSL runs everywhere - this same work
   feeds the WASM web-gallery moonshot.
3. **fidget-solver** (`solve()` in fidget-solver/src/lib.rs):
   constraint solving over tape gradients - the Tier 3
   differentiable-CAD item.  Later, after stage B.

## 5. Deferred / parked

- Affine collapse default: flip `STIBIUM_AFFINE` on only when a
  workload shows base-tape evals dominating (the tile renderer's
  coarse passes are the candidate); it changes float rounding, so
  Nate makes the fidelity call.
- Interval-derived bounds (recon item, M) - unrelated to rendering,
  removes the per-shape bounds tax.
- Keeter's 2026-07-03 "Please Steal my Meshing Algorithm Idea"
  (QEF points + Delaunay) - the meshing successor, separate track.

# Shortened-tape evaluation: design

*2026-07-13. The Tier 1 keystone from [LIBFIVE-RECON.md](LIBFIVE-RECON.md) §1/§8:
replace in-place `disable_nodes`/`enable_nodes` pruning with immutable,
shareable, contiguous shortened tapes (libfive `eval/tape.cpp` pattern,
adapted to this kernel's C core and four eval modes).*

## Why (recap)

The current pruning mutates the `MathTree` (swap-to-back + `ustack`
counts) and stores results *inside* nodes (`Results` = `f`, `i`,
`r[MIN_VOLUME]` per node). Consequences:

- Every worker thread needs a full `clone_tree` (triangulate.cpp,
  render_mt.cpp, parallel_eval.cpp all clone serially).
- Pruning at one recursion level can't be handed to sibling subdivisions
  as a compact artifact; each `eval_*` still walks the full leveled node
  array with pointer-chasing per node.
- `disable_node` must `fill_results` so parents read sane stale values —
  a mutation that makes tapes unshareable by construction.

The tape design fixes all three and is the prerequisite for the CPU SIMD
tile viewport and any GPU renderer (MPR-style), which specialize tapes
per tile.

## Core structures (new module: `fab/tree/tape.h` + `src/tree/tape.cpp`)

C++ implementation (house precedent: grid.cpp, contour.cpp), **opaque
types + extern "C" API** so C consumers (render.c) and C++ consumers
(mesher, parallel_eval) both work.

- **Clause** `{op, out, a, b, m, imm, payload}` — one operation, SSA-ish:
  `out/a/b/m` are *slot* indices (m = OP_GRID's z input; payload = its
  MeshGrid*). `imm` backs constant-fill clauses emitted by binary
  pruning (see below).
- **Deck** (per model, immutable once built): `num_slots`, the base
  Tape, constant slot/value pairs, X/Y/Z slot ids, active-axes mask.
  Compiled from a `MathTree` in one level-order pass —
  constants/axes get slots but **no clauses**; the eval driver fills
  their slot rows directly (libfive does the same).
- **Tape** (immutable, refcounted, `parent` strong-ref): contiguous
  `Clause` vector, root slot, optional X/Y/Z bounds (set when pushed
  from an interval eval), `terminal` flag (no min/max left ⇒ push is a
  no-op).
- **TapeCtx** (per thread, one allocation): result arrays indexed by
  slot — `f[num_slots]`, `i[num_slots]`, `r[num_slots * MIN_VOLUME]`
  (g-mode reuses `r` cast to `derivative`, capacity MIN_VOLUME/4,
  exactly like today's `shaded8`) — plus push scratch
  (`disabled[]`, `remap[]`, libfive Deck-style).

Slot ids are stable across pushes: a shortened tape references the same
slots, so no remapping of workspaces ever happens and `getBase`-style
walk-ups are free. Slot count = node count (no register allocation in
v1; a Fidget-style allocator is a later, orthogonal win — noted, not
attempted).

## Evaluation

`tape_eval_f/i/r/g(tape, ctx, ...)` walk the clause vector forward,
calling the **existing `math_f/i/r/g` primitives** — same operations on
same inputs in same order as today's level walk, so results are
bit-identical by construction. `tape_eval_i` leaves per-slot intervals
in `ctx->i`, which is what `tape_push` consumes.

## Push (replaces disable_nodes / disable_nodes_binary / enable_nodes)

`tape_push(tape, ctx, X, Y, Z, mode)` — after an interval eval; walks
clauses backward (reverse-topological = the old top-down level walk):

- min/max branch choice uses the **exact comparisons** from
  `disable_nodes` (MAX: `lower(a) >= upper(b)` keeps a, etc.). The
  losing branch is dropped and readers of the min/max's slot are
  **remapped** to the survivor's slot (libfive remap-chain). No
  `fill_results` needed — parents read the surviving child directly,
  which is the same value min/max would have produced. Values are
  therefore *exact*, not approximated: pruning never changes output.
- **binary mode** (render-only, replaces `disable_nodes_binary`):
  boolean-context propagation (min/max/neg chains from the root) as
  today; a sign-determined subtree root becomes a single `imm` clause
  writing `i.upper` — byte-for-byte what `fill_results(node, i.upper)`
  made parents see. Only render8/16 use it (sign-only consumers);
  the mesher never does (needs true distances).
- OP_GRID clauses are KEEP_ALWAYS (oracle-style) but keep their three
  children live.
- unchanged/terminal ⇒ returns the same tape (ref++), so deep uniform
  regions cost nothing.

Recursion does `push` on the way down, `tape_release` on unwind —
the ustack dance becomes ordinary object lifetime. `tape_base(tape, p)`
walks up parents until the point is inside the recorded bounds
(needed later by DC-style meshers; cheap to include now).

## Consumer migration (each step gated on goldens)

1. **parallel_eval**: shared base tape + per-thread TapeCtx; the serial
   clone loop is deleted outright.
2. **render.c / render_mt**: recursion carries `Tape*`; `#if PRUNE`
   blocks become push/release pairs; region8/16 use `tape_eval_r`,
   shaded8 `tape_eval_g`. Public signatures keep taking `MathTree*`
   (compile a Deck at entry) so app code doesn't change in this pass.
   render_mt shares one Deck; per-thread ctxs; no clones.
3. **mesher/triangulate**: Mesher gets Deck+ctx instead of a cloned
   tree; `triangulate_region` pushes per octree level (the old code
   only pruned at packed-block level — pushing per level is the actual
   point of the change and cannot alter values, only skip dead work).
   `load_packed`/`unload_packed` lose their disable/enable calls.
4. **Delete** `disable_nodes`, `disable_nodes_binary`, `enable_nodes`,
   `ustack`, `clone_tree`, `Node.clone_address`, and the
   `fill_results`-on-disable path. `MathTree` remains as the
   build/parse/serialize representation; `eval.c` remains for the
   low-rate consumers (bounds.cpp, contour.cpp) until they migrate
   with the SIMD-viewport work.

## Verification

- New `[tape]` unit tests: bit-equivalence of `tape_eval_*` vs
  `eval_*` across the parser corpus (points, intervals, regions,
  gradients), push-vs-disable_nodes equivalence on interval splits,
  nested push/release, binary-mode equivalence on sign-determined
  trees, OP_GRID tapes.
- `[golden]` mesh dumps must `cmp` bit-identical across every step.
- `--render` PNGs before/after must `cmp` identical.
- `--resave`/determinism untouched (no serialization changes).
- Bench: `[.bench]` gyroid + zeiss_id02 render/mesh wall time before
  vs after (expect wins from per-level mesher pushes, contiguous
  clause walks, and deleted clone setup; regression = investigate,
  not ship).

## Results (landed 2026-07-13)

All of the above ran clean: `[tape]` equivalence suite (1687
assertions), full SbFabTest (598,759 assertions / 53 cases), all four
golden mesh dumps and all --render PNGs (csg, gear, showcase_gear,
import_bead with OP_GRID, zeiss_id02) byte-identical old-vs-new.

Wall time (8c/16t, examples at stated res):

| workload | before | after |
|---|---|---|
| gyroid mesh r10, detect, 1 thread | 7.69 s | 7.42 s |
| gyroid mesh r10, detect, 8 threads | 2.38 s | 2.15 s |
| zeiss_id02 --render | 3.27 s | 3.00 s |
| gear.sb --export r60, detect | 2.33 s | 1.62 s |
| zeiss_id02_merged --export r7, 1/4 scale | 70-91 s | **17-23 s (~4x)** |
| zeiss_id02_merged --export r7, full scale | **DNF** (killed at 801 s, no output yet) | **310 s, 9.8 GB STL** |

The scaling law is the MPR paper's: pruning wins grow with model
size and min/max density. The gyroid (almost no min/max) barely
moves; the merged Zeiss (dozens of unioned assemblies) gets ~4x and
goes from cannot-complete to done-in-five-minutes at full scale -
the old evaluator re-walked the entire node array at every octree
level, the tape build's stack shrinks as it descends.

Peak RSS flat everywhere (the N per-thread clones each carried full
per-node result buffers, so per-thread workspaces are a wash on
memory until register allocation lands).

## Round 2 (landed 2026-07-13, same night)

The follow-up sweep after the keystone landed: soundness, memory,
and the instrumentation the tile renderer will steer by.

**Pruning fuzzer + maybe-NaN taint.** New `[fuzzer]` test: random
expressions (min/max-heavy, domain-error ops included on purpose) x
random regions, comparing pushed tapes against base pointwise, plus
nested pushes and binary-mode sign checks (`STIBIUM_FUZZ_COUNT` /
`STIBIUM_FUZZ_SEED`; 100k trees = 10.6M assertions green). It found
four real bugs, three of them pre-dating the tape work:

1. **Domain-error bounds poisoning pruning** (ancient, conceptual):
   sqrt/log/asin/acos/0-div of out-of-domain inputs produce garbage
   finite intervals, and pruning trusted them.  Fix: per-clause
   maybe-NaN taint recorded during interval evaluation (libfive's
   `Interval::maybe_nan`, plus infinity rules it lacks: inf-inf,
   0*inf, sin/cos/tan(inf) are NaN factories with clean-looking
   bounds).  Tainted records veto pruning and binary collapse;
   interval values are never altered, so culling is unchanged.
2. **pow_i truncated the exponent to int(B.lower)** (ancient,
   upstream): pow(x, Z) got the bounds of pow(x, 0).  Real/varying
   exponents now get sound corner bounds for positive bases and
   +/-inf + taint otherwise; the integer fast path is byte-identical.
3. **Binary collapse used lower >= 0** (ancient, upstream): a bound
   touching zero isn't one sign class under negation (neg(+0) stays
   >= 0, neg(upper) goes negative), so exactly-zero fields could
   flip fill state.  Now strictly > 0.
4. OP_COPY missing from arity() (ours, same session - the fuzzer
   caught it in nested pushes within minutes).

**The taint fix repaired a visible artifact:** the old build's
showcase_gear render differs from unpruned ground truth by 5,842
pixels (3.5%, involute-gear domain errors driving bogus collapses);
the new render matches ground truth bit-for-bit.  That render is
re-baselined as a correctness fix - every other golden and render
compares identical.

**Register allocation.** Deck compile now linear-scans clause
outputs onto a minimal register file (constants/axes stay pinned;
evaluators are elementwise, so a dying input donates its register
in-place).  Merged Zeiss: 3157 clauses run in 465 registers -
workspaces shrink ~4x (and ~7x on smaller decks).  Two structural
consequences, both MPR-shaped: pushes emit OP_COPY instead of
remapping readers (a survivor's register may be redefined between
def and use), and min/max choices + per-clause bounds/taint are
recorded during the interval pass (`ClauseIv`), because the register
file only holds values until their last read.  Wall-time neutral on
CPU meshing today; it's the precondition for the SIMD tile viewport
and any GPU path (bounded per-thread registers).

**Deck caching.** `Shape::getDeck()` compiles lazily once (thread-
safe, shared across copies); the multithreaded renderers take a
Deck* instead of recompiling the tree every frame.  Headless one-
shots are unchanged; interactive viewport frames stop paying a
per-frame compile.

**Block-escape fix.** `get_normals`' +/-epsilon probes can step just
outside the packed block whose specialized tape is in effect; they
now walk up the tape stack (`tape_base_for_region`) until covered.
Goldens unchanged (the escape has to actually flip a pruned branch
to matter), but the hole is closed.

**Churn control.** Thread-local spare-tape freelist (tape objects
and their clause capacity recycle without locks) and no more pushes
inside packed blocks, where the curve shows shrinkage has flatlined
- push count on the merged-Zeiss export drops 4.7M -> 0.6M.

**Shrinkage curve** (`STIBIUM_TAPE_STATS=1`, merged Zeiss export):

| depth | tapes | avg clauses |
|---|---|---|
| 0 | 1 | 3157 |
| 1 | 57 | 1479 |
| 3 | 1,575 | 591 |
| 6 | 172,019 | **156 (20x)** |
| 9 | 2,318,029 | 200 |

Keeter's MPR shrinkage, reproduced on our worst model; the flat tail
is the near-surface selection effect.  Depth ~6 saturation is the
number the tile renderer's tile size should be tuned against.

**Affine collapse** (`STIBIUM_AFFINE=1`, off by default): folds
add/sub-by-constant and mul-by-constant chains (stacked transforms)
at deck compile - 836 folds on the merged Zeiss, 3157 -> 2914
clauses.  Wall-neutral on exports (hot loops run on pushed tapes
where decided mins already dropped those chains), so it stays
opt-in: reassociation changes float rounding, and bit-identity is
the default until the tile renderer gives it a workload that pays.

Verification: full suite 625,588 assertions green; 100k-tree fuzz
green; all four goldens byte-identical; all renders byte-identical
except showcase_gear, which now matches unpruned ground truth
exactly (see above).

## Round 3 (landed 2026-07-13, still the same night)

**`tape_eval_i_batch`: the vectorized interval evaluator.**  The MPR
tile hierarchy finally pays on CPU.  One pass down the tape
classifies up to 64 boxes at once: tape traversal (the clause fetch,
the switch dispatch) is amortized across the whole batch, and the
hot interval ops - add/sub/mul/min/max, the bulk of any CSG tape -
run as elementwise loops the compiler vectorizes.  Every branchy op
(div/pow/mod/trig/grid) delegates to the scalar math_i function per
lane, so batched bounds are *exactly* the scalar bounds by
construction - fill decisions cannot depend on which evaluator ran.

The layout is the good kind of pun: an r-row is MIN_VOLUME = 128
floats per register, which is exactly 64 lower bounds followed by 64
upper bounds.  The batch workspace is the region workspace; constant
rows (one value replicated across the row) are already correct
interval constants.  Zero new allocation.

The renderer's fan-out (which lost ~13% as a serial probe) now
batch-classifies its 64 children: provably-inside tiles fill,
provably-empty tiles vanish, neither ever re-walks the tape; only
ambiguous children recurse (scalar eval_i + push, as before, so
ClauseIv stays coherent).  Near-to-far ordering feeds the occlusion
cull.  Default ON; `STIBIUM_TILE_RENDER=0` restores binary
bisection, `=N` tunes the fan-out.

| merged Zeiss --render | binary bisection | batched fan-out 64 |
|---|---|---|
| 512 px | 3.52 s | 3.32 s |
| 2048 px | 14.97 s | **12.36 s** |

Pixel-identical at every size; the win grows with resolution (the
~2.5 s Python graph eval is constant, so the render-only fraction
improves ~20%+).  Verified: suite (627,666 assertions), 100k-tree
fuzz with batch-parity checks folded into every iteration (11.4M
assertions), goldens and all renders byte-identical.

Remaining rungs: per-tile worker scheduling (xy-chunking is coarser
than MPR's tile queue), and the GPU compute shader (Tier 3) - which
now needs nothing the kernel doesn't already have: bounded
registers, 8-byte-able clauses, batch classification, copy-based
pushes.

Notes for the next rung (the SIMD tile viewport):

- Pushed tapes are malloc'd per push; a spares freelist (libfive
  `Deck::claim`) is the first thing to try if profiles show churn.
- Slot count = node count. Fidget-style register allocation would
  shrink workspaces dramatically and improve cache behavior.
- Known pre-existing wobble: multithreaded *detect-features* meshing
  is not run-to-run deterministic (each worker's Mesher shares its
  swappable-edge map across the chunks it happens to grab, so seam
  swaps vary with scheduling; ±2 triangles in ~1.1 M observed).
  Confirmed present before the tape change (old build differs from
  itself the same way). Single-threaded output - the golden referee -
  is bit-stable.

## Round 4 (2026-07-13, post-compact): the fmin toll

The autovec audit (TAPE-NEXT §3) went looking for missed
vectorization in the batch kernels and instead found the kernel's
oldest tax.  A gdb-sampled profile of the 2048 px merged-Zeiss
render put ~64% of busy samples in `tape_eval_r` — and most of THOSE
inside `fmin@plt`/`fminf64`.  `min_f`/`max_f` were `fmin(A, B)`:
C's *double* fmin, so every min/max node × every voxel paid
float→double×2, a PLT call into libm, and a truncation back.  A CSG
model is mostly min/max nodes.  This predates Antimony's first
commit.

Fixes, in order of landing:

1. **`#pragma GCC ivdep`** on the batch kernels and fill loops in
   tape.cpp — each lane k touches only indices k and BH+k, so the
   pragma is a statement of fact; it removes gcc's runtime aliasing
   check and versioned scalar twin (~2.5%).
2. **`min_f`/`max_f` written out** (math_f.h) as the bit-exact
   replica of this platform's libm/x86 semantics: `isnan(B) → A`,
   `isnan(A) → B`, else `A < B ? A : B` (ties **including ±0 return
   B**, NaN-vs-NaN returns A — measured, not assumed: glibc's fmin
   here is minsd-based asm, NOT its own C template, which returns +0
   for fmin(+0,-0)).  Parity proven by exhaustive-specials + 100 M
   random-bit-pattern harness against real libm (`-fno-builtin`),
   0 mismatches; sNaN quieting is the one divergence, and arithmetic
   can't produce sNaN.  `mul_i`/`div_i`/`batch_mul`/`min_g`/`max_g`
   fmin chains routed through them (~31% — the headline).
3. **`-fno-trapping-math`** scoped to SbFab — the kernel never reads
   FP exception flags, and the flag is what lets gcc if-convert the
   NaN selects: `min_r`/`max_r` DUAL loops now vectorize (~6%).

| merged Zeiss --render 2048 px | wall |
|---|---|
| round 3 champ (same night, same machine state) | 11.8 s |
| + ivdep | 11.5 s |
| + inline min_f/max_f | 8.1 s |
| + -fno-trapping-math | **7.6 s** |

gear r60 detect-features export: 2.37 s → 2.25 s; quarter-scale
merged r7 export ~18 s (mesher time is octant classification +
feature detection, not min chains — mesher batching stays the next
big export lever).

Verified at every step: full suite (627,666), 100k fuzz (11.4 M),
goldens bit-identical old-vs-new, 2048 px render byte-identical to
the round-3 binary's.  `batch_mul` still doesn't vectorize (the
4-product NaN-select chain defeats SSE2 if-conversion) — it's
inlined now, and batch kernels totalled ~16% of the profile, so it
keeps.  Candidates if it ever matters: `-march=x86-64-v2` (SSE4.1
blends; changes no arithmetic results) or an explicit two-phase
min/NaN-repair formulation.

## Round 5 (2026-07-14, same long night): first light on the GPU

Stage A of the GPU rung (TAPE-NEXT §4) is real: the deck now
serializes to a flat u32 blob (`tape_export_blob`, format v1 -
constants, axis slots, 5-word clauses; refuses OP_GRID), and a
headless EGL + GL 4.3 compute prototype (`tests/gpu.cpp`, hidden
tags `[.gpu]` / `[.gpubench]`) interprets it per pixel-ray and
writes the same heightmap the production renderer does.

The referee: because the renderer's BINARY pushes preserve signs
exactly and a pixel is a pure function of signs (highest k with
field < 0), a bit-exact GPU eval must reproduce the CPU heightmap
bit-for-bit.  It does: **zero mismatched pixels** on every test
model - algebraic CSG, spheres (sqrt), nested min/max, even the
sin/cos sheet - on BOTH the Intel UHD iGPU and the GTX 1650 Ti
(prime-run drives EGL directly).  min/max in GLSL are hand-written
mirrors of math_f's exact semantics (GLSL min() has undefined NaN
behavior); the C opcode values and every blob offset are baked into
the generated source, so there is no protocol to drift.

Brute-force bench (nested-CSG model, 1650 Ti, dispatch+readback,
shader compile excluded as a per-deck cost):

| n | GPU brute force | CPU 1-thread (full tape arsenal) |
|---|---|---|
| 256³ | 14.3 ms | ~5 ms |
| 512³ | 76.2 ms | ~42 ms |

Read that carefully: a laptop GPU evaluating EVERY voxel of the
FULL tape with NO pruning is within ~2x of a CPU thread wielding
everything this document describes.  The pruning arsenal is worth
~20-100x; it currently lives only on the CPU.  Stage B ports it:
fidget-wgpu's pipeline (interval tile classification =
tape_eval_i_batch on device; tape_simplify.wgsl = tape_push on
device) - a port, not research, and the bytecode it consumes is
the blob above.

## Round 6 (2026-07-14, dawn): the arsenal on the device

Stage B: the MPR pipeline itself now runs on the GPU, one
subdivision level, three passes (tests/gpu.cpp, tags `[.gpu2]` /
`[.gpu2bench]` / diagnostics `[.gpudbg]` `[.gpudump]`):

1. **interval_tiles** - one thread per 16px tile evaluates the tape
   in interval arithmetic over the tile's xy-frustum (full z),
   recording per-clause min/max choices under the CPU's exact
   conditions (<=/>= comparisons, operand-taint veto) and
   classifying EMPTY / FULL / AMBIGUOUS.
2. **tape_simplify** - one thread per ambiguous tile replays
   tape_push's backward liveness walk (STANDARD verdicts:
   dead/elide/copy/keep) and emits a compacted per-tile tape.
3. **voxel march** - per pixel, the stage-A interpreter runs on the
   tile's SHORTENED tape; EMPTY/FULL tiles fill without marching.

The soundness frame that makes the GLSL tractable: the GPU's
choices need not MATCH the CPU's - any sound choice set leaves
values on the shortened tape exact, and pixels are functions of
those values' signs.  So hot interval ops are `precise` mirrors of
math_i, and anything exotic (pow/mod/trig widening) returns
[-inf,inf] + taint: conservative, sound, and the referee can't
tell the difference.  Verified: zero mismatched pixels vs the
production CPU renderer, all 8 models including the decided-min
PRUNABLE, on Intel and NVIDIA.

Bench (nested-CSG, 1650 Ti Max-Q, dispatch+readback):

| n | GPU brute | GPU pipeline | CPU 1-thread |
|---|---|---|---|
| 256³ | 14.0 ms | 7.6 ms | ~5.2 ms |
| 512³ | 77.3 ms | **33.1 ms** | 42.3 ms |
| 1024³ | 425.7 ms | 121.8 ms | 75.3 ms |

The pipeline beats brute force everywhere and crosses over the CPU
at 512³ - with ONE subdivision level, no z-culling, and a ~30-clause
model where tape shortening barely matters.  The CPU retakes 1024³
because its octree culls in z; the fix is known (3D tiles / z-slab
classification + a second subdivision level, i.e. the rest of
fidget's pipeline).  Big decks are where the device should shine:
the merged Zeiss drops 3157 -> ~200 clauses under pushing, and the
GPU pays the full-tape price only once per tile instead of per
voxel.

War story for the file: Mesa's Intel compiler miscompiled the
simplify pass - interleaved SSBO reads/writes made the emitted
copy's out-slot word read back 0 (NVIDIA ran the same source
correctly).  Locals-then-store fixed it; pattern kept in the
shader with a comment.  Trust the referee, not the driver.

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

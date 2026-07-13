# libfive vs. Antimony: what Keeter learned, and what Stibium should borrow

*Recon sweep, 2026-07-13. Grounded in a fresh clone of
`github.com/libfive/libfive` (MPL/GPL, 2022-era master), the MPR 2020
SIGGRAPH paper ("Massively Parallel Rendering of Complex Closed-Form
Implicit Surfaces"), Keeter's project/blog writing, and read-only
inspection of this repo. Companion to the upstream issues/PR recon in
TODO.md.*

A framing fact that reshapes the whole analysis: **libfive is not
Keeter's latest thinking — Fidget is.** The lineage is `kokopelli →
Antimony (2013, C++/Qt node graph) → Ao/libfive (2016, Scheme text +
C++ kernel) → Fidget (2022, Rust kernel + JIT)`. libfive is the
*middle* generation; items refined further in Fidget are flagged.
Keeter now calls libfive "40K lines of mostly C++ … extremely
challenging to hack on, even as the original author" (Fidget
"Origins") — so the play for Stibium is *idea reuse and selective code
port*, not adopting libfive wholesale.

---

## 1. Representation & evaluation

**libfive's Tree.** A `Tree` is a reference-counted handle to a
`TreeData`, a `std::variant` (`libfive/include/libfive/tree/data.hpp`)
over `TreeNonaryOp` (X/Y/Z/const), `TreeUnaryOp`, `TreeBinaryOp`,
`TreeConstant`, `TreeOracle`, `TreeRemap`, `TreeApply`, `TreeInvalid`.
Two design points matter:

- **`TreeRemap` is a first-class node.** Coordinate substitution (what
  Stibium does with `part.map(Transform(...))`) is a *node in the
  DAG*, lazily flattened via cached `TREE_FLAG_HAS_REMAP` bits.
  Cleaner than textually rewriting prefix strings inside `Transform`,
  and why libfive's entire transform/loft/revolve/twirl stdlib is
  "just remap."
- **Deduplication + affine collapse** happen in `Tree::optimized()`
  (`src/tree/tree.cpp:612`): common-subexpression elimination via a
  `canonical` map keyed on `[opcode, lhs_id, rhs_id]`/constant value,
  plus collapse of nested affine forms (`AffineMap =
  unordered_map<Tree,float>`, `tree.cpp:626`) — e.g. `(2x + 3x)` →
  `5x`, sums-of-scaled-terms folded into one linear combination.
  Stibium dedups structurally at parse time but has **no
  affine-collapse pass**; long transform chains (each `move`/`scale`
  emits arithmetic) accumulate un-folded, inflating the node array the
  evaluator walks.

**Evaluation.** The tree compiles once into a `Deck`
(`eval/deck.hpp`) producing a `Tape` — a flat `std::vector<Clause>`,
each `Clause` being `{op, id, a, b}` (`eval/clause.hpp`). Four eval
modes over that tape (point / SIMD-array via
`LIBFIVE_EVAL_ARRAY_SIZE` / interval / jacobian-feature), mirroring
Antimony's four backends. **There is no JIT in libfive core** — that
arrived in Fidget (SSA bytecode + hand-written x86-64/AArch64 asm,
"31× on brute-force eval, ~25% once the smart algorithm is in play").
JIT is a Fidget idea, not a libfive one, and a large lift.

**The tape-shortening trick, and where Stibium already stands.**
`Tape::push` (`src/eval/tape.cpp:23`) is the CPU ancestor of the MPR
paper's core optimization. Given an interval region, it walks the tape
and for each `min`/`max` clause asks a `KeepFunction` which branch
survives; the losing subtree is dropped, and a **new, shorter,
immutable tape** is emitted with a `parent` pointer back to the fuller
one. `getBase(point/region)` (`tape.cpp:168`) walks *up* the parent
chain when a point falls outside the current cell (needed in DC where
vertices leave their cell). Key refinements over Antimony:

- Tapes are **immutable and shareable** (`shared_ptr`) → thread-safe
  across the mesher's worker pool with no per-thread tree clone.
- The shortened clause vector is **contiguous and re-emitted**
  (cache-friendly sequential walk), with a `terminal` flag
  short-circuiting when no `min`/`max` remain.

Stibium's evaluator instead does **in-place
`disable_nodes`/`enable_nodes`** on a rank-leveled node array, swapping
disabled nodes to the back of their level and pushing counts onto a
`ustack` for reversal. Same *idea* (interval-driven pruning of
`min`/`max` branches), but in-place mutate+stack is inherently
single-tree/single-thread-per-clone and doesn't leave a contiguous
specialized tape to hand to many child cells. **The single most
backportable-with-real-payoff item** (see §8).

**Interval arithmetic differences.** libfive's `Interval`
(`eval/interval.hpp`) wraps `boost::numeric::interval<float>` **plus a
`maybe_nan` flag** — every op tracks whether the result could be NaN
(domain errors, `0/0`, `∞−∞`), with `min`/`max` propagating NaN
semantics to match `std::min`. More rigorous than `math_i.c`, and it
directly prevents meshing artifacts when a `sqrt`/`log`/`acos` domain
error silently produces garbage bounds. Keeter's caveat (2025-05-14
gradients post): interval arithmetic "forgets correlations" (`x·x`
over `[-1,1]` gives `[-1,1]` not `[0,1]`) and degrades as transforms
stack.

**Backportable to the CPU C core, realistically:** (a) affine-collapse
pass — **yes, M**; (b) `maybe_nan` interval tracking — **yes, S-M**;
(c) immutable shared shortened-tape eval — **yes, M-L, highest
value**; (d) JIT — **no** (Fidget, L); (e) SIMD array eval — Stibium
already parallelizes, marginal.

---

## 2. Meshing

**libfive ships three meshers** (`render/brep/settings.hpp`:
`DUAL_CONTOURING`, `ISO_SIMPLEX`, `HYBRID`) plus an optional `VolTree`
acceleration structure, all over a shared hierarchical octree (`XTree`
base). The flagship is **Manifold Dual Contouring**:

- Adaptive octree subdivided until edge < `min_feature` (default 0.1),
  cells collapsed when combined QEF error < `max_err` (default 1e-8)
  → **adaptive resolution**: flat regions become big cells, detailed
  regions stay fine.
- Each ambiguous leaf places a vertex by minimizing a **QEF** (`AtA`,
  `AtB`, `BtB` in `DCLeaf`, `dc/dc_tree.hpp:84`; writeup at
  mattkeeter.com/projects/qef) over Hermite samples → **sharp edges
  and corners**.
- **Feature rank** (1=face, 2=edge, 3=corner) drives collapse;
  **manifold criterion** from [Ju et al. 2002]/[Gerstner et al. 2000]
  guarantees output **watertight and manifold** even across
  octree-level transitions (`dc_mesher.cpp` handles multi-level shared
  edges explicitly).

**How Stibium compares — closer than you'd think.** Our `Mesher` is
*not* plain marching cubes: it implements feature-sensitive extraction
(Kobbelt et al. 2001, extended marching cubes), and `check_feature()`
already does an SVD least-squares normal-fit distinguishing edge vs.
corner features. That *is* a QEF-flavored sharp-feature recovery, so
Stibium already gets sharp edges. **What Stibium lacks vs. libfive
DC:**

1. **Adaptive resolution.** We mesh a uniform grid; libfive spends
   triangles where curvature demands. Big win for large models with
   local detail.
2. **Manifold/watertight guarantees.** EMC can produce
   cracks/non-manifold edges at resolution boundaries; libfive's
   octree + Ju-manifold-tables are provably watertight+manifold.
3. **Dual topology** (one vertex per cell) gives better-conditioned
   triangles than MC's per-tet cases.

**Honest counterweight (Keeter's own words):** Manifold DC "will not
necessarily preserve thin features; the resulting meshes may include
self-intersections; vertex positioning is susceptible to adversarial
cases" (Fidget meshing notes; see also the 2023-04-23 adversarial-model
post). His **2026-07-03 "Please Steal my Meshing Algorithm Idea"**
post proposes the successor: separate "generate surface points (QEF)"
from "build manifold mesh (Delaunay tetrahedralization)." Don't treat
libfive DC as an endpoint — but it is a clear step up from uniform-grid
EMC.

**License (verified from file headers):** the mesher
(`src/render/brep/**`) is **MPL 2.0** without the
incompatible-with-secondary-licenses marker → GPL-compatible. Literal
copy into Stibium is allowed (retain per-file MPL headers; those files
stay MPL, the combined binary ships GPLv3). Do **not** copy the Guile
`.scm` files or Studio (GPLv2+).

---

## 3. The GPU story (MPR 2020)

MPR (Keeter, ACM ToG 39(4) Art. 141; ref impl `github.com/mkeeter/mpr`,
C++/CUDA) is the tape+interval machinery reprojected onto SIMT:

- **Compile once** (CPU) to an SSA tape: CSE'd DAG → topo sort →
  greedy register allocation → **8-byte clauses** (1B opcode, 1B
  out-slot, two 1B in-slots *or* 1B slot + 4B immediate), 128
  slots/thread.
- **Shallow high-branching hierarchy, not a deep octree.** 2D:
  `64×64 tiles → 8×8 subtiles → per-pixel`; 3D: `64³ → 16³ → 4³ →
  per-voxel`. The 64× branching factor makes each subdivided tile
  exactly two 32-lane warps → no thread divergence.
- **Interval pass per tile** classifies inside/outside/ambiguous and
  records which branch each `min`/`max` took (2 bits/choice).
- **Tape shortening** (one backward pass): 6056-clause model → 356
  clauses at 64² tiles (17×), 28 clauses at 8² subtiles (216×).
- **3D output = heightmap + normals** (atomic-max depth, occlusion
  early-out), normals via forward-mode autodiff over the pruned tape.
- Interpreter is 19× slower per-op than a model-specific compiled
  kernel, but pruning wins overall: beats the compiled kernel above
  ~1536², <8 ms/frame at 4096² on a 1080 Ti.

**Minimal path for Stibium, cheapest-first:**

1. **CPU SIMD tile renderer** (no GPU): tile the viewport,
   interval-cull tiles, specialize the tape per ambiguous tile (§1c),
   AVX float-eval survivors. SIMD lanes replace warp lanes; the thread
   pool replaces the GPU grid. Pan/zoom is a transform matrix on
   X/Y/Z — rebuild the tape only on model edit. Effort M, captures
   most of the interactive-latency win.
2. **Then GPU** if needed: compute shader (Vulkan/`QRhi`) running the
   8-byte-clause interpreter + 64×-tile scheme. Effort L.

Gotcha Keeter flags: smooth-blend-heavy models prune poorly (interval
bounds stay ambiguous) — worst-case frames on any hardware.

---

## 4. Modeling stdlib gaps (all MPL, all "borrow the math")

Every libfive stdlib function is pure arithmetic + remap, bounds-free
(`libfive/stdlib/stdlib_impl.cpp`) — trivially expressible as Stibium
prefix strings. Confirmed missing from `py/fab/shapes.py`:

- **Exact-SDF primitives** (`box_exact`, `rectangle_exact`,
  `rounded_box`, `stdlib_impl.cpp:206`): true Euclidean distance
  `min(0,max(dx,dy,dz)) + |max(d,0)|`. Our `rectangle`/`cube` use the
  cheap max-form, which is **not** a real distance field — so
  `offset`/`shell`/blends misbehave near corners. Highest-leverage
  single port in the stdlib.
- **`blend_expt`/`blend_expt_unit`** — log-sum-exp smooth union
  (`-log(exp(-m·a)+exp(-m·b))/m`, `stdlib_impl.cpp:79`) and
  **`blend_difference`** (smooth subtraction; ours is union-only).
  `OP_LOG`/`OP_EXP` already landed, so these are ~3 lines each.
- **`half_space(norm, point)`** = `dot(p−point, norm)`: arbitrary cut
  planes, exact polygons/triangles.
- **`gyroid`/TPMS** (`stdlib_impl.cpp:318`): `sin·cos` sum wrapped in
  `shell` — ~5 lines, big lattice/infill value.
- `clearance(a,b,o)` (press-fit gaps), `loft_between` (skewed
  z-interpolated remap between profiles), `symmetric_x/y/z`
  (abs-remap mirror).

**Not worth porting:** libfive's `array_*` are naive O(n) unions —
Stibium's `mod`-based Repeat family is strictly stronger (O(1),
unbounded). We're ahead there.

---

## 5. Bounds handling

libfive went **bounds-free by construction** — a `Tree` has no bbox;
bounds are supplied only at mesh/render time. Keeter first tried
interval-derived automatic bounds (2016-04-10, 2017-09-29 posts), then
removed the `findBounds` API "due to unpredictability."
Antimony/Stibium thread explicit `Bounds` through every shape — and the
tax shows: `offset` "assumes a linear distance field," `revolve_y`
hand-recomputes extents, `attract`/`repel` hand-pad by `r/e`, `invert`
needs bounds for its complement cube. Every primitive author must
reason about bounds; several are wrong for non-exact fields.

**Recommendation:** adopt the middle path Keeter used first —
**interval-derived bounds** (contract a seed region with the interval
evaluator we already have) to remove per-shape bounds threading (M).
Don't chase full bounds-free (L, fights the 2D/3D split and export
code).

---

## 6. Oracles vs. OP_GRID

libfive's **Oracle** (`oracle/oracle.hpp`): an opaque callable in the
tree implementing `evalInterval/Point/Array/Derivs/Features`,
`checkAmbiguous`, and crucially `push(Tape::Type)` returning an
`OracleContext` — the oracle **participates in tape-shortening**,
spatially specializing itself per region.

Our `OP_GRID` is a different and in some ways nicer design: pure
serializable data, composes with transforms like any node, no callback
machinery. What to steal from the oracle contract when import needs
grow:

1. **Interval tightness**: an oracle doing native interval arithmetic
   over its source (mesh BVH sub-tree, heightmap tiles) beats bounding
   a trilinear cell → better pruning. For **heightmap import**, prefer
   an oracle-style payload over resampling to a 3D grid.
2. **Derivatives**: trilinear normals are C0-discontinuous at cell
   boundaries → faceted shading, confused feature detection. An
   explicit deriv path fixes import render quality.
3. **Self-specialization** (`OracleContext.push`): let a huge-source
   payload carry its own spatial index instead of one flat grid.

Keep `OP_GRID` as the default; add a "smart payload" path if/when
large-source import demands it.

---

## 7. Where Stibium is already ahead

- **The node-graph product itself.** Keeter's Graphene retrospective
  frames abandoning the graph UI as a solo-maintainer bandwidth
  decision, not a verdict — "I've struggled to maintain a kernel and
  full-fledged CAD UI simultaneously." A modernized node-graph
  Antimony is the tool he wished someone else would build on a good
  kernel. (His two graph-engine wishes: correct-by-construction name
  scoping, and collapsing the node/subgraph distinction.)
- **O(1) domain repetition** beats stdlib `array_*`.
- **Mesh import** — libfive has no built-in mesh import.
- **2D vector export pipeline** (SVG/DXF corner recovery + DP
  simplification) — nothing comparable in libfive core.
- **Headless CLI verbs, deterministic saves, analytics** — product
  features a library doesn't have.
- **Color preservation + chamfer/fillet families.**

The honest asymmetry: **libfive is ahead on the kernel; Stibium is
ahead on the product.** Exactly the split Keeter predicted.

---

## 8. Verdict — ranked borrowable shortlist

Effort S = days, M = week-two, L = month+.

| # | Item | Source | Effort | Impact |
|---|------|--------|--------|--------|
| 1 | Exact-SDF primitives + log-sum-exp blends (`blend_expt_unit`, `blend_difference`) | stdlib | **S** | **High** — fixes offset/shell/fillet quality; blends are ~3 lines now that log/exp landed |
| 2 | `half_space`, `gyroid`/TPMS, `clearance`, `symmetric_*`, `loft_between` | stdlib | **S** | Med-High |
| 3 | Affine-collapse + CSE optimization pass | `tree.cpp:612` | **M** | Med-High — shrinks transform-chain bloat before every eval |
| 4 | `maybe_nan` interval tracking | `eval/interval.hpp` | **S-M** | Med — kills domain-error meshing artifacts |
| 5 | Immutable shared shortened-tape eval (replaces disable_nodes) | `eval/tape.cpp` | **M-L** | **High** — the keystone; thread-safe specialization, enables 6 and 10 |
| 6 | CPU SIMD tile viewport renderer | MPR ideas | **M** | High — interactive latency without GPU dependency |
| 7 | Adaptive manifold DC mesher (port MPL code or reimplement) | `render/brep/dc/**` | **L** | High — adaptive resolution + watertight guarantee |
| 8 | Interval-derived bounds | ao bbox posts | **M** | Med — removes the per-shape bounds tax |
| 9 | Oracle-style smart-payload import path | `oracle/oracle.hpp` | **M** | Med, situational |
| 10 | GPU compute-shader viewport (after 6) | MPR / `mkeeter/mpr` | **L** | High but large |

**Strategic note:** 1–4 are landable now. 5 is the keystone — the CPU
realization of what Keeter calls the trick "that makes large
expressions feasible" — and unblocks 6 and 10. 7 is the marquee
quality upgrade. Meta-lesson from three rewrites: the durable value he
kept is all kernel-level; the graph product layer is ours and he
wished someone would build it. **Borrow his kernel, keep our
product.**

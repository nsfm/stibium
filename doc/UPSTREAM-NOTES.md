# Findings worth carrying upstream

*Compiled 2026-07-15. Stibium is an Antimony-descended kernel that
reimplemented fidget/libfive-style tape simplification and the MPR
pipeline; along the way its pruning fuzzer and GPU port surfaced
things that apply to the upstream projects. Code doesn't transfer
(C/C++ vs Rust) - findings do. Filing any of this is Nate's call.*

---

## 1. fidget: interval simplification is unsound under hidden NaN

**Status: verified against fidget-core v0.4.4** - failing test at
`~/code/fidget/fidget-core/tests/nan_inf_soundness.rs` (ours, not
committed to the fidget repo).

**One-line summary:** an operation that maps an *infinite* input
interval to *finite clean bounds* (e.g. `sin`) hides a pointwise NaN
from `min_choice`/`max_choice`'s `has_nan` veto, so simplification
can drop the branch that pointwise evaluation would have returned.

**Reproducer:** `min(sin(exp(x * 30)), -4.5)` over x ∈ [2, 4]:

- `exp(x * 30)` spans `[1.1e26, +inf]` (f32 overflow above x ≈ 2.96)
  with clean bounds - no NaN anywhere.
- interval `sin` of an infinite-width interval returns `[-1, 1]`,
  clean.
- `min_choice([-1, 1], [-4.5, -4.5])`: `-4.5 < -1` ⇒ `Choice::Right`,
  the sin branch is dropped.
- Full tape, pointwise at x = 3.9: `min(sin(inf) = NaN, -4.5)` =
  **NaN** (fidget's pointwise min/max propagate NaN) ⇒ outside.
- Simplified tape: **-4.5** ⇒ inside.
- The root *interval* is `[-4.5, -4.5]`, so a tile renderer fills the
  whole region as inside - disagreeing with its own pointwise
  evaluator at every overflow pixel. The rendered image depends on
  evaluation granularity.

Notes:
- `sin([inf, inf])` (fully-infinite) is accidentally safe: `width()`
  = inf − inf = NaN poisons the quadrant analysis. The hole is
  half-infinite inputs, where `width() ≥ TAU` returns `[-1, 1]`.
- The same pattern likely applies to `cos`, and to any future op
  that maps unbounded inputs to bounded outputs.
- fidget's pointwise NaN semantics also differ between backends
  (Rust `f32::min` returns the non-NaN operand; x86 `minss` returns
  the second operand; fidget's VM propagates NaN) - worth a
  consistency pass while touching this.

**Fix shape that worked for us** (Stibium tape.cpp, "maybe-NaN
taint"): a per-clause boolean carried alongside the interval,
propagated from operands, set by domain errors (sqrt/log of
negatives, asin/acos out of range, div/mod by zero-spanning) *and by
infinity rules the bounds arithmetic hides*: taint `sin/cos/tan`
when the input touches ±inf, and `add/sub/mul/div/pow` when an
operand touches ±inf (inf−inf, 0·inf, inf/inf are pointwise NaN
factories with clean-looking bounds). Tainted operands veto
min/max choices; interval *values* are never altered, so culling
quality is unchanged. Found by a pruning fuzzer (random expression
trees × random regions, pushed-vs-base pointwise equality; ~100k
trees per run) - the fuzzer found four real bugs in our kernel on
day one and this class was one of them.

## 2. fidget: clause rescheduling beats register spilling

fidget bounds its register file with LRU spilling (`alloc.rs`,
`Load`/`Store` ops). We hold a card the fixed-order tape doesn't:
reordering clauses at tape-build time is free (pure dataflow - any
topological order is bit-identical). SSA-ify, schedule Sethi-Ullman
(heavier subtree first, DAG values emitted once at first need), then
run a fresh linear scan:

- our heaviest real model dropped from **877 slots to 95 registers**,
  zero spill machinery needed;
- thread-uniform constants (409 of those 877 slots!) moved out of
  the register file entirely (into shader source / uniform memory).

Directly relevant to the stalled `fidget-wgpu` bytecode interpreter
("incredibly inefficient" on M1 Max): a dynamically-indexed register
array lives in local memory on every GPU, so shrinking it is the
single highest-leverage change. Measured on a GTX 1650 Ti Max-Q:
rescheduling + constants-out + two-level tape specialization +
atomic-pool compacted tapes + band recycling took a 3172-clause
model from 153 s to 24 s at 512px, bit-identical output throughout -
still short of the same machine's 8-core CPU (3.8 s), which matches
Matt's own read that the GPU interpreter loop needs more than
memory fixes (warp-cooperative decode is the likely next rung).

## 3. Antimony (if the archive is still accepting patches)

All present in `mkeeter/antimony` `develop` and fixed in Stibium:

- **`min_f`/`max_f` are libm calls** (`fmin`/`fmax`, double-
  precision, through the PLT) - measured at ~40% of render wall time
  on CSG-heavy models. Bit-exact inline replacements exist
  (mind the platform semantics: ties including ±0 return the second
  operand, NaN-vs-NaN returns the first; we verified with a
  100M-pair harness against libm).
- **`pow_i` silently truncates real exponents to int** -
  `pow(x, y)` for non-integer interval `y` gets the bounds of
  `pow(x, 0)`, which can drive wrong pruning (visible artifacts).
- **`interpolate_between` walks the whole pending queue per
  scheduled edge** with an Eigen vector compare - O(edges × queue),
  ~30% of export wall time on large models. A hash map keyed on the
  bit-packed vertex pair is semantics-identical (note: `Vec3f` is
  `Eigen::Vector3d` - doubles - despite the name).
- **`disable_nodes_binary` treats zero-touching bounds as one sign
  class** (`>= 0`), which flips exactly-zero fields under negation;
  the test must be strict (`> 0`). Caused a 5,842-pixel artifact in
  one of our showcase renders.

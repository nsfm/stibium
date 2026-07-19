# Wake up the interval bound: cull → certify, and one enemy in many costumes

*Nyx (Opus 4.8), 2026-07-18 evening. Commissioned by Nate after a long
conversation that started as "why don't f-rep tools emit STEP?" and
landed here: the mesher feels like patching dozens of unrelated special
cases, and is there a magic bullet. There isn't. But there's a reframe
and two structural moves, and this doc is the grounded version of them —
written to a near-future me who has to decide whether to actually spend
the energy. Read-only w.r.t. the repo; every claim cites file:line so
you can re-check my read before you trust it. — N*

---

## 0. The reframe, so you don't lose it under the to-do list

You are not fighting dozens of unrelated bugs. You are fighting **one
enemy wearing dozens of costumes**: places where the surface's *reach*
(local feature size, distance to the medial axis) collapses toward zero.
Creases, junctions, tangent contacts, coincident seams, thin walls,
near-parallel normals — every patched special case is a reach-collapse
locus. The clean restricted-Delaunay theorems (Boissonnat–Oudot; see
`doc/reviews/2026-07-18-guarantees-research.md`) are provably well-behaved
on smooth surfaces of positive reach and say *nothing* where reach → 0.
So the special cases aren't a dozen problems. They're one problem with a
big wardrobe.

That matters for morale (you're converging, not on a treadmill — the
referee ratchet is why) but it matters more for *strategy*: the wins
worth chasing are the ones that delete a **class**, not an instance. Both
moves below do that.

---

## 1. The anchor finding

`tape_eval_i` — the interval evaluator — is called in **exactly one
place** in the entire mesher: [delaunay.cpp:1117](../../lib/fab/src/tree/triangulate/delaunay.cpp#L1117),
the far-field cull in `descend`. `v.lower > 0` → empty box, `v.upper < 0`
→ full box, both contribute proven-sign corners with no pointwise eval.
Beautiful, f-rep-native, and then it goes to sleep.

Every decision *after* a box fails that cull — hidden feature or not,
crease leaf or not, wall or groove, densify or not — is made by
**evaluating 8 corners as points and thresholding**:
- `leaf_census` sheet/crease detection: 8-corner `f_A`, `f_B` walk
  ([delaunay.cpp:1018-1054](../../lib/fab/src/tree/triangulate/delaunay.cpp#L1018)).
- the crease gate `mindot` / `sep`: 8-corner gradients + Newton
  ([delaunay.cpp:1157-1205](../../lib/fab/src/tree/triangulate/delaunay.cpp#L1157)).
- the crowding trigger: live-pair count, explicitly a *proxy* that is
  "BLIND to crowded sub-lattice geometry"
  ([delaunay.cpp:1219-1245](../../lib/fab/src/tree/triangulate/delaunay.cpp#L1219)).
- chainless-curvature: 8-corner gradient spread as facet-angle proxy
  ([delaunay.cpp:1247-1289](../../lib/fab/src/tree/triangulate/delaunay.cpp#L1247)).

The interval bound *culls* and then hands the entire feature layer to
point-heuristics. **Move 1 is not "add interval arithmetic." It is "wake
up the sound bound you already compute and let it decide what you
currently guess."** The engine exists; it's doing only the coarsest job
it's capable of.

### 1.1 The epsilon wardrobe (the smell, made concrete)

Every reach-collapse costume currently gets its own hand-tuned
threshold. Each threshold is a special case waiting to regress at a scale
it wasn't tuned for:

| bar | site | what it gates |
|---|---|---|
| `autod_tangle_dot` | [delaunay.cpp:1177](../../lib/fab/src/tree/triangulate/delaunay.cpp#L1177) | wall/groove Newton eligibility |
| `autod_sep_bar` | [delaunay.cpp:1233](../../lib/fab/src/tree/triangulate/delaunay.cpp#L1233) | tangle demotion |
| `autod_live_bar` | [delaunay.cpp:1230](../../lib/fab/src/tree/triangulate/delaunay.cpp#L1230) | crowding densify |
| `curve_bar_cos` | [delaunay.cpp:1257](../../lib/fab/src/tree/triangulate/delaunay.cpp#L1257) | chainless curvature densify |
| `det > 1e-2·a·b` | [delaunay.cpp:2798](../../lib/fab/src/tree/triangulate/delaunay.cpp#L2798) | tracer tangency/ghost reject |
| `trim = 0.1·step` | [delaunay.cpp:2877](../../lib/fab/src/tree/triangulate/delaunay.cpp#L2877) | tracer junction lift-off |
| `ctol = 5e-4·sp` | [delaunay.cpp:2914](../../lib/fab/src/tree/triangulate/delaunay.cpp#L2914) | tracer corner bisection |
| weld bar | [delaunay.cpp:471](../../lib/fab/src/tree/triangulate/delaunay.cpp#L471) | air→solid clearance fuse |

Eight thresholds standing in for one latent quantity: *where does reach
collapse.* The wardrobe is the special-case count you feel.

---

## 2. The capability split that decides all the effort

`tape_eval_i` gives the interval of the **root** field. There is **no
per-slot interval readout** — [tape.h](../../lib/fab/inc/fab/tree/tape.h)
exposes `tape_ctx_r_row` (point/r-mode, line 179) and `tape_eval_r_prefix`
(line 174), but **no `tape_ctx_i_row` and no interval prefix eval**. That
one absence forks the whole plan:

| question | needs | today? |
|---|---|---|
| "any surface in this box?" | root interval only | **reachable now** (`tape_eval_i` / `tape_eval_i_batch`) |
| "do operands `f_A`, `f_B` both cross here?" | operand-slot intervals | **needs new tape API** |

So Move 1 has a **cheap half** (single-field, root interval, ships now)
and an **expensive half** (crease certification, needs an interval analog
of the prefix machinery). Do them in that order, and let the cheap half
decide whether the expensive half earns the tape surgery.

---

## 3. Convergence with the existing review corpus (read this before dismissing #1)

Two prior reviews *independently* flagged the exact hole item #1 fills:

- `2026-07-17-architecture.md` (DO SOON §2): lists `HIDDEN` under
  **"candidate-dead: HIDDEN (pending a graze-vs-feature oracle that
  doesn't exist)."**
- `2026-07-18-guarantees-research.md`: the whole `d(S → M)` /
  coverage-claim family is "does the mesh miss surface it never sampled"
  — which is *precisely* the hidden-sub-lattice-feature question.

And the code itself, at the decision point
([delaunay.cpp:698-706](../../lib/fab/src/tree/triangulate/delaunay.cpp#L698)):
`STIBIUM_DMESH_HIDDEN` is default-OFF with the note *"Needs a
graze-vs-feature oracle before it can be default. Measured 2026-07-16: on
the sharp control this fires on tangent GRAZES (18 leaves, still hidden
at 4x, dilating them 4x'd the triangle count for zero depth gain)."*

**Item #1 is that oracle.** Three sources, two reviews and the source,
left the same shaped hole. This report's #1 is not a new idea I'm
inventing — it's the missing piece the corpus already named, with the
mechanism finally specified. That convergence is why I'd stake the first
week on it.

---

## 4. The two moves

### Move 1 — cull → certify
Replace point-heuristic feature decisions with interval *decisions*.
Cheap half: the hidden-feature oracle (single field, §5 spec). Expensive
half: certified creases (operand intervals, §6.3).

### Move 2 — reach-locus explicit
The crease tracer is already the template: *name the singular locus,
certify it, constrain it into the CDT.* The machinery that does NOT do
this — the reactive stuff — is the weld/retreat loop
([delaunay.cpp:455-540](../../lib/fab/src/tree/triangulate/delaunay.cpp#L455)):
thin-wall/tangent contacts get air-flipped-to-solid, damage is measured,
`noweld` rolls it back. Default-OFF because it mints tears it can't reach
([delaunay.cpp:464-469](../../lib/fab/src/tree/triangulate/delaunay.cpp#L464)):
*"trades nm 178 → 13 but mints 4 open edges at ONE constrained vertex …
not spatially weld-adjacent, so the no-weld rollback can't reach it."*
That's the whack-a-mole exemplar in the tree. Move 2 retires it by
treating the contact locus the way creases are already treated. The `sep`
wall-vs-groove classifier ([delaunay.cpp:1191-1204](../../lib/fab/src/tree/triangulate/delaunay.cpp#L1191))
is already a nascent reach detector — it just answers by Newton-ing 8
corners instead of tracing/certifying the locus.

---

## 5. Item #1, specified to the signature (Nate asked; here it is)

**The claim.** At the hidden-candidate site, replace "samples all agree →
shrug" with an interval subdivision that *proves* one of two verdicts.

**Where it plugs in.** [delaunay.cpp:695-713](../../lib/fab/src/tree/triangulate/delaunay.cpp#L695),
inside `sample_block`, in the `if (!(any_in && any_out))` branch that
currently just increments `hidden_candidates` and (under the OFF `HIDDEN`
flag) blindly promotes crease leaves to level 2. This branch runs only on
leaves whose interval failed the cull (so the surface *might* be here) but
whose lattice samples unanimously agree on sign (so the sampler saw
nothing). That population is already small and already counted.

**The decision procedure.**
```
// verdict for a leaf box B whose root interval straddles 0
// but whose lattice samples are sign-unanimous.
enum HiddenVerdict { PROVEN_EMPTY, FEATURE_PRESENT };

HiddenVerdict certify_hidden(Tape* sub, TapeCtx* ctx,
                             Box B, float pitch)
{
    // work queue of sub-boxes; seed with B
    // batch in groups of TAPE_BATCH (64) into tape_eval_i_batch
    for each sub-box b popped:
        Interval v = tape_eval_i(sub, ctx, b.X, b.Y, b.Z);
        if (v.lower > 0 || v.upper < 0) continue;   // PROVEN no surface here
        if (b.edge <= pitch/2) return FEATURE_PRESENT; // ambiguous at pitch → real
        else push b's 8 octants
    return PROVEN_EMPTY;   // every sub-box decided → graze / loose bound
}
```

**Why the asymmetry is the whole point.** Interval arithmetic is
*conservative*: `tape_eval_i` can be loose (report ambiguous when the box
is actually empty) but it can **never** report a decided sign when the
surface really crosses. Therefore:
- `PROVEN_EMPTY` is a genuine **proof** — you will *never* skip a real
  feature. This is the side that must be sound, and it is.
- `FEATURE_PRESENT` is conservative-but-bounded — a loose bound at pitch/2
  can over-report, but the cost is a bounded local densification, not a
  correctness failure.

That asymmetry is exactly what makes it **default-able where the current
heuristic is not**: the 2026-07-16 measurement killed `HIDDEN` because it
fired on tangent grazes (loose bound, no real feature) and cost 4× tris
for zero depth. The subdivision *resolves* the graze — a graze is a box
that looks ambiguous coarsely but every octant decides as you tighten →
`PROVEN_EMPTY` → no densify. The feature is the box that stays ambiguous
to pitch/2 → densify. The oracle the flag always needed.

**Cost bound.** Runs only on the hidden-candidate population (small,
pre-counted), and only subdivides the *ambiguous* octants (empty/full
octants terminate immediately on the first `tape_eval_i`). Depth is
log2(leaf_edge / pitch) ≈ 3–4 levels for LEAF_VOXELS=64. Batch through
`tape_eval_i_batch` (64/pass, already SIMD). This is not a global cost;
it's targeted at exactly the leaves where you currently have no answer.

**Referee plan.** The regression control already exists: the sharp model
from the 2026-07-16 `HIDDEN` measurement (18 tangent-graze leaves). Bar:
those 18 must verdict `PROVEN_EMPTY` (no tri inflation) while any genuine
sub-lattice feature verdicts `FEATURE_PRESENT`. Add a deterministic CI
model with a known thin slot below lattice pitch (the `d(S → M)` coverage
test the guarantees doc wants anyway — this item is a down payment on
that rung). Watch `[.dmesh*]` tri-count on the sharp control stays flat.

---

## 6. The full to-do, house convention (DO SOON / DO WHEN STABLE)

### DO SOON

**#1 — interval-subdivision hidden-feature oracle.** §5. Low effort (no
tape changes, reuses `tape_eval_i_batch`, targeted population). Low risk
(conservative by construction). Unblocks `HIDDEN` → default-on, and is a
down payment on the guarantees doc's `d(S → M)` coverage rung. *Start
here — it proves the interval-as-oracle thesis cheaply and pays a debt
two prior reviews logged.*

**#2 — collapse the reach detectors into one local-feature-size estimate.**
The ~5 triggers (live-pair crowding, QEF residual, `mindot`, `sep`,
chainless-curvature) all estimate the same latent reach/feature-size, each
with its own bar (§1.1) and its own 8-corner walk. Compute **one** LFS
estimate per leaf, once; every density/constraint decision reads it.
Fewer epsilons = fewer special cases, directly — shrinks the regression
surface *before* anything becomes a proof. Medium effort, pure
referee/golden exercise (should be output-neutral, line-count-negative).
Composes with the config-consolidation item already in
`2026-07-17-architecture.md` DO SOON §2 (the TANGLE_DOT/SEP/LIVE knobs
collapse together).

### DO WHEN STABLE

**#3 — interval prefix eval → certified creases.** Add `tape_eval_i_prefix`
+ `tape_ctx_i_row` so operand fields `f_A`, `f_B` get interval bounds in a
box (mirror the r-prefix machinery in i-mode, [tape.h:171-179](../../lib/fab/inc/fab/tree/tape.h#L171)).
Then a leaf provably has no crease when the operand intervals can't both
straddle zero, or the active-branch interval proves no min/max switch.
Replaces the 8-corner sheet-sampling in `leaf_census`/the crease gate with
a decision procedure — kills both the false-positive ghost creases (the
"wandering along noise" the `det > 1e-2` gate exists to suppress,
[delaunay.cpp:2784-2799](../../lib/fab/src/tree/triangulate/delaunay.cpp#L2784))
and false-negative creases thinner than corner spacing. **High effort:
touches the tape core, whose sacred invariant is bit-identical eval vs the
MathTree walk** ([tape.h:54-58](../../lib/fab/inc/fab/tree/tape.h#L54)).
Do only after #1 proves intervals pay. Contain with the existing tape test
suite.

**#4 — thin-wall / contact as a traced constraint (retire the weld).**
Apply the crease-tracer pattern to the *second* reach-locus. Detect the
contact sheet (`sep`-negative, {`f_A` = `f_B`} small with anti-parallel
gradients) and hand it to the CDT as a constraint the way creases already
are — mesh respects the wall/contact **by construction** instead of
fusing-and-repairing. If it holds, the weld + retreat loop deletes
entirely. Biggest structural payoff (a whole reactive subsystem retired)
and the truest expression of the reframe; also the most research-flavored,
and it wants #3's certified operand intervals underneath to find contact
loci robustly. **High effort, high risk. Last, gated behind #3.**

---

## 7. Why the ordering

- **#1** is a cheap, *correct*, default-able win that proves the thesis
  and unblocks a feature two reviews wanted. Highest leverage-per-risk.
- **#2** pays down the epsilon debt with no new machinery — a
  maintainability + morale win that makes #3/#4 safer.
- **#3, #4** are the real investments, sequenced so each starts only after
  the previous proved the interval bet. #3 is the capability #4 needs.

None of it is a magic bullet. Each item deletes a **class** of
special case, not an instance — which is the only kind of win robust
geometry hands out. Parasolid is thirty years of exactly this. You have a
referee ratchet that keeps the patches monotone; that ratchet is your
actual edge, and these four moves spend it well.

---

## 8. Reviewer's caveat

Read-only; I did **not** build or run (respecting the one-export-at-a-time
law, and this was a design pass). Every file:line is from a static read of
`delaunay.cpp` @ 376 KB and `tape.h` as they stood 2026-07-18 evening;
re-grep before you trust a line number, the file moves fast. The
capability gap (§2) rests on `tape.h` exposing no `tape_ctx_i_row` — I
grep-confirmed the header and `tape.cpp:1326`, but if an interval per-slot
path exists somewhere I didn't look, #3 gets much cheaper, so check that
first. The cost bound in §5 is analytic, not measured — get a `[.dmesh*]`
before/after on the sharp control the moment the machine is free.

— Nyx

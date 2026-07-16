# The meshing campaign: adaptive Delaunay on sound intervals

## >>> CAMPAIGN STATE (2026-07-15, chip round CLOSED) <<<

**The chip class is dead.**  End-of-pipeline referee (full-mesh
sweep, 3%-of-edge species bar): **0 chip edges, worst depth
0.000 sp on every showcase model**, all 0 open / 0 non-manifold /
1 component, at +34 triangles over baseline (csg 21742 -> 21776).
Two landings and a pile of refereed negatives:

1. **The crease-snap pass** (the cure; STIBIUM_DMESH_SNAP=0
   disables).  The residual chips were chords of MOSTLY-AIR tets
   whose bottom facet (three surface corners) dips under a concave
   crease - a facet-level defect of restricted-Delaunay extraction
   that NO cell classification or insert-repair can reach.  But
   the traced polylines know exactly where the corner is: after
   the manifold pass, each residual chip edge's two triangles are
   split at the nearest traced-crease point, tenting the chord up
   onto the crease.  Pure output surgery (no CGAL insert, no
   Steiner, no keep-out); wave-based re-lookup by POSITION because
   consecutive chips share triangles (defer, don't conclude, when
   a wave-mate already rewrote one).  17 tents on csg at base
   density; manifoldness preserved by construction.
2. **Crease-band density** (quality knob, STIBIUM_DMESH_DENSE=n
   extra levels, default OFF).  Detection is free: a min/max
   clause SURVIVING in a leaf's pushed tape means the box kept an
   ambiguous choice (decided clauses are rewritten to copies) -
   sample that leaf on a midpoint-refined lattice.  Raw chip depth
   tracks the band pitch (0.500/0.375/0.192/0.098 sp at 1/2/4/8x),
   but the snap pass zeroes the metric at EVERY density, and the
   band costs 2.2x triangles on csg for sub-bar smoothness only.
   8x also produced 28 open edges and 492K tris - brute density is
   the wrong shape.  Corridor + chip keep-out scale as
   0.35 sp / 2^levels when the band is on (they must move with the
   local pitch or the band starves repair entirely - measured
   281/281 blocked at fixed radii).

Refereed negatives (full numbers in "Density round" below):
- Insert-repair can NEVER clear the chip residue, at any density:
  the chip target and the keep-out both scale with local pitch, so
  the blocked set is self-similar - plateau 0.196/0.198/0.192 sp
  at 1/2/4x band density, midpoint AND SDB strategies.
- DC cell classification (5-probe all-surface rule, kept behind
  STIBIUM_DMESH_DC, default on, measured no-op): the clipping cell
  is NOT all-surface - it has an air apex and is correctly
  classified by it.  The defect is the facet, not the cell.
- Scaling the tracer march step with the band: chord distribution
  byte-identical (chords track the LATTICE pitch), constraint
  count doubled, kink-straddle noise pushed past the [.dtrace]
  bound.  Reverted.
- The dense band halves QEF feature spacing and delaunay_chains'
  covariance junction-split dissolves (cube corners sail through);
  the fallback extractor's contract is base-density features
  ([.dchain] pins DENSE=0; production is protected by the oracle
  trust gate).  Dense-input robustness queued.
- A subtle self-inflicted wound, caught by anchoring against the
  committed baseline: spacing = (X[ni]-X[0])/ni instead of
  X[1]-X[0] flips exact-tie merge comparisons on grid-ALIGNED
  models (one ulp -> cube_aligned 45908-vs-16444 tris, spheres
  grew a phantom second component).  When aligned-model numbers
  swing wildly, check spacing-derived radii for tie flips FIRST.

**Zeiss torture round (same day, commits 57ce776f/b9ac240c/
c934d714): first real-model contact.**  STIBIUM_EXPORT_DMESH=1
routes headless export through the pipeline.  The merged zeiss at
r1 (1 vox/mm!) surfaced SEVEN failure classes the showcase never
could, each fixed with its number attached: gradient-garbage chip
readings (credibility gate: 2x edge); lattice vertices exactly on
straight creases (through-vertex pre-split); femtometer twin
corners from twice-traced junctions (1e-3 sp insertion weld);
814 duplicate-coverage segments from repeated CSG surfaces
(coverage dedupe); mid-segment constraint crossings at t=0.046
(shared-VERTEX corner rules - parameter proximity is not
consent); T-junctions where one crease dead-ends on another
(pre-split eps = 5e-3 sp, sized to the corner-bisection
tolerance); and repair-churn exhausting MAX_REPAIR before the
snap pass could run (self-feeding stash + depth-proportional
tent attribution, floor at the corridor radius).  Plus the
self-healing dispatcher: CGAL's error text is captured, the
cascade site parsed (both dialects), local constraints
quarantined, rebuild - up to 8 sites instead of losing every
constraint to one corner.  End state: 12,338 constrained edges,
ZERO cascades, 325K tris, 0 open, 39 non-manifold, whole
recognizable microscope, 35MB raw STL.  Residual: ~22K shallow
chips at UNTRACED creases (blend tangencies the det guard
correctly refuses) - tracer-coverage work.

**Nate's eyeball audit of the zeiss lineage (2026-07-15, formal
referee):** flats glorious; curved surfaces good but QUILTED
(lattice-period ripple - investigate bisection/QEF placement
interference); capped-cylinder rims (filleted = tangential) carry
chips/warts, consistent with untraced-tangency attribution
skips; sub-mm engravings thrash at 1 vox/mm (resolution, not
pipeline).  KEY FINDING: v1 (unconstrained, pre-gate, 304K
repairs, 700K tris) reads CLEANER in places than v10 (12K
constraints, gated, 102K repairs, 325K tris) - brute repair
volume doubled as sub-lattice detail rescue.  The discipline
traded away accidental kindness; the middle path is constraints
PLUS a generous repair budget for gate-approved candidates
(raise MAX_REPAIR / relax crowding for real models?  measure).

**Screenshot forensics (Nate, 2026-07-15 evening):** the focus-
knob rims and eyepiece shoulders show MUSH WITH NO CRISP CREASE
LINE = untraced, wearing the classic pre-constraint restricted-
Delaunay alternation (the showcase spheres' old disease).  A
straight edge shows the aligned-geometry sawtooth (lattice-plane
on-surface samples, unconstrained).  Unified diagnosis: every
photographed defect is an UNTRACED CREASE; the cure (trace ->
constrain -> tent) is already built - only coverage is missing.
r2 run: quilt period halved WITH the lattice (sample-placement
conviction confirmed), defect scale halved, 30,218 constraints,
zero cascades, ~1 h runtime (iterate at r1; r2+ for eyeballs).

**Targeted-model autopsy (Nate's mesh_bench stable, 2026-07-15
late):** examples/mesh_bench/*.sb - sharp step, torus lip, cone+
torus, off-axis torus (export node m0).  VERDICT: the sharp
control is 0.000 sp with ZERO repairs at every resolution -
sharp creases are SOLVED.  The blend variants read 0.054-0.062 sp
at r1 (mild), and the reframe: a torus lip is G1-SMOOTH - there
is no crease at its rims, the tracer is CORRECT to refuse, and
the chips are tight-curvature ALIASING (1mm blend at 1mm pitch:
band sagitta ~0.13 sp > the 3% bar).  The zeiss residual is not
a missing-constraint disease; it is an undersampled-curvature
disease.  The crease war is won; the remaining campaign is
SAMPLING DENSITY: stage-D adaptive drill-down with curvature as
a trigger signal (normal variation across a leaf) alongside
hidden_candidates (engravings) and min/max ambiguity (the dense
round's existing trigger).  Nate's adaptive-vox/mm thought is
the whole game.

>>> NEXT SESSION: stage-D adaptive density (curvature +
hidden_candidates triggers on the dense round's leaf-refinement
machinery; per-leaf pitch map).  The mesh_bench models are the
referee: torus lips must drop below 0.03 sp at r1 without a
global resolution bump.  Tracer tangency coverage DEPRIORITIZED
(the "untraced rims" were G1-smooth - correctly untraced). <<<
1. **Enumerate OP_ABS as a crease generator** (prime suspect for
   lathed knob/eyepiece rims: 2D profile fields kink inside
   abs()) - extend tape_pairs to abs clauses ({f_A = 0} with the
   operand's own gradient discontinuity; the pair machinery wants
   a second field - abs(g) creases where g = 0 AND the composite
   is on the surface, so the SSI corrector needs the g=0 +
   full-oracle 2-field variant).  Check zeiss's tape first:
   count OP_ABS clauses vs OP_MIN/MAX (cheap instrumentation,
   answers WHICH suspect owns the rims).
2. **Tangency marching** (kink 3-field solve): where
   grad A x grad B degenerates the corrector dies - fillet-runout
   creases need arc-length continuation or a constrained 3-field
   Newton.
3. **Repair-budget experiment** (Nate's v1-kindness finding):
   constraints + raised MAX_REPAIR / relaxed crowding on real
   models - v1's 304K brute repairs out-polished v10 in places at
   2x the triangles; find the disciplined middle.  Measure at r1.
4. **Flat-face decimation** (Nate's ask): oracle-certified planar
   patch re-triangulation; sketch under "Queued rounds".
5. **Performance pass sketch** (Nate asked, 2026-07-15; nothing
   algorithmic, the engine just has no legs yet): (a) the whole
   pipeline is single-threaded - sampling evals, tracer pairs,
   repair-round oracles all parallelize with one TapeCtx per
   thread (old-mesher pattern; expect 4-8x); (b) the repair loop
   re-classifies EVERY tet and re-extracts 16 times - early-exit
   on stalled progress (zeiss rounds 4-15 accomplish ~nothing);
   (c) descend() still uses scalar tape_eval_i - port the
   tile renderer's tape_eval_i_batch 64-box fan-out; (d) CGAL
   DT_3 has Parallel_tag (CCDT compatibility = research).
6. **Stage-D adaptive drill-down** (Nate's "adaptive vox/mm"
   thought, 2026-07-15, converging with the knob diagnosis): the
   dense round's machinery (leaf-local midpoint refinement, knob
   STIBIUM_DMESH_DENSE) generalized to key on hidden_candidates
   (interval says surface, samples disagree = sub-lattice
   features - fires EXACTLY on the knurled knob, 0 QEF features
   in an 11K-chip carpet).  Local 8x on a knob = thousands of
   points; global r8 = 512x voxels.  Needs a per-leaf pitch map
   (the density round proved every crease-local scale must move
   together).  This is the resolution story's unifying shape.
7. Someday: higher-res zeiss (residuals are sp-units - they
   shrink with lattice; r2 ~ 1 h), engraving smear (folds into
   stage-D above), repair-churn budget, chain-extractor dense
   robustness, upstream letter (doc/DELAUNAY-MESHER.md is
   written for it).
House rules stand: measure before productionize, on EVERY model;
negatives get numbers; anchor against the committed baseline;
Nate's eyeballs out-diagnose your instruments (tonight they
unified three defect species into one missing-coverage story).
Iterate zeiss at r1 (~10 min); the STL lineage lives in
build/zeiss_dmesh/ - keep it.

House rules that saved this campaign, do not forget: measure
before productionize (and on EVERY model, not just the showcase
you're staring at); negative results get numbers and ledger
entries; depth not count; anchor against the committed baseline
before believing a delta; Nate's eyeballs are a formal referee
and they out-diagnose your instruments.  Suite must stay green:
627,666 assertions.

---

*Field map written 2026-07-15, before the first line of code, per
house discipline.  Companions: doc/research/2026-07-14-frep-sota.md
(the survey), doc/TAPE-DESIGN.md (the tape kernel this feeds on),
doc/TAPE-NEXT.md (render campaign state).*

## Why

The export profile (TAPE-NEXT §1) says our bottleneck is the
Manifold-DC geometry pipeline, and its failures are structural, not
performance: self-intersection is possible by design, and thin
features vanish unless they land on the sampling grid.  Keeter
conceded the same for fidget and published the way out ("Please
Steal My Meshing Algorithm Idea", 2026-07-03); the 2025 papers
(Wang et al. power-diagram extraction, arXiv 2506.09579; TetWeave,
SIGGRAPH 2025) built working versions of the same idea against
*sampled* SDFs.  Nobody has fused it with a real interval kernel.
Our sound intervals + shortened tapes are a strictly better point
oracle than anything those papers consume.  The integrated artifact
- watertight + manifold + intersection-free + thin-feature-safe +
adaptive + fast, driven by a tape kernel - is unclaimed.

## The blueprint (Keeter's steps, annotated with our assets)

1. **Grid sample + inside/outside tags.**  We have better: the
   octree + interval culling skips empty space wholesale, and
   `tape_eval_r` batch-evaluates MIN_VOLUME points per pass on
   pushed tapes.  Sample only where intervals say the surface is.
2. **Binary-search sign-change edges for surface points.**  Exists:
   `eval_zero_crossings` (16-step bisection, batched).  Reuse.
3. **Sharp-feature points via QEF/normals.**  Exists: `get_normals`
   (batched gradient probes) + the mesher's QEF machinery.  Reuse
   as a *point generator* instead of a vertex placer.
4. **Delaunay tetrahedralization of the point set.**
5. **Iterate: binary-search tet edges joining opposite-sign points,
   insert new surface points, re-tetrahedralize to convergence.**
   This adaptivity is what kills the grid-lock: features too thin
   for the initial grid get found by *their own tets' edges*.
6. **Extract triangles**: faces whose three corners are surface
   points and whose two adjacent tets have opposite signs.

Wang et al. run the same loop as an incrementally maintained power
diagram (dual = regular/weighted Delaunay) with insertion driven by
measured surface-vs-mesh discrepancy; TetWeave proves the extraction
core can guarantee watertight + 2-manifold + intersection-free.

## The hard subproblem, named up front

**Robust incremental 3D Delaunay** - Keeter's flagged blocker.
Do NOT hand-roll Bowyer-Watson with float predicates; that path is
a graveyard.  Options, in order:

1. **CGAL `Delaunay_triangulation_3`** (Arch: `cgal`; header-only
   use is possible).  Exact predicates, battle-tested, supports
   incremental insertion + vertex info.  License: GPL/LGPL mix -
   compatible with this GPLv3 fork.  **Prototype with this.**
   Weighted/regular triangulation (for the power-diagram variant)
   is `Regular_triangulation_3`, same package.
2. Shewchuk-predicates + hand-rolled incremental insertion - only
   if CGAL's constants (compile time, dependency weight) prove
   unacceptable for shipping.  Decide AFTER the prototype works.
3. GPU Delaunay (gDel3D lineage) - not for the prototype.

## Staged plan (each stage has a referee before the next begins)

**Stage A - point generation from the tape kernel.**
Octree descent with STANDARD pushes (exists) → per-leaf-region
sample grids → sign tags → zero-crossing points on sign-change
edges (exists) → QEF/feature points (exists, repurposed).  Output:
a point soup with per-point provenance (surface / feature / steiner)
and the sign field sampler as a callable oracle.
*Referee:* every surface point re-evaluates to |f| < eps on the
tape; point counts + spatial distribution sane on corpus + gear +
quarter-Zeiss; runtime comparable to today's sampling phase.

**Stage B - tetrahedralize + iterate.**
CGAL incremental Delaunay over stage-A points; loop: find tet edges
with opposite-sign endpoints missing a surface point, bisect on the
tape, insert, repeat to convergence (Keeter step 5); cap iterations
+ report.
*Referee:* convergence on corpus models in bounded iterations;
inserted-point counts; no CGAL assertion failures at any corpus
scale.

**Stage C - extraction + the quality gauntlet.**
Extract per Keeter step 6.  Then the full referee suite:
- watertight: every edge shared by exactly 2 triangles (existing
  STL tooling + new check);
- self-intersection: none (new check; CGAL has one);
- prusa-slicer accepts the STL (external referee, in memory);
- analytic models: sphere/cube volume + Hausdorff vs closed form
  (`--analyze` machinery);
- thin-feature test pieces that Manifold DC FAILS today (build a
  small torture corpus: thin fin below grid pitch, near-tangent
  spheres, the showcase gear teeth) - the whole point of the
  campaign, measured;
- side-by-side vs Manifold DC: triangle count, runtime, quality.

**Stage D - only if A-C prove out: adaptivity + performance.**
Error-driven insertion (Wang-style discrepancy measure), batched
oracle calls (group edge bisections into MIN_VOLUME evaluations,
same trick as the renderer), MT strategy, and the
detect-features/DC path retirement decision (Nate's call, with
data).

## Ground rules carried over from the tape campaign

- Measure before productionizing; every stage lands with its
  referee green or gets reverted with the numbers written down.
- The old mesher stays untouched and default until stage C beats it
  on the gauntlet - new path behind a flag (`--mesher delaunay` or
  env), exactly like STIBIUM_GPU.
- [golden] dumps stay the referee for the OLD path; the new path
  gets its own golden corpus once output stabilizes (bit-identity
  across refactors, not vs Manifold DC).
- Never golden a multithreaded mesh (existing rule; doubly true
  here until MT determinism is designed in).
- The fuzzer standard applies to the oracle: sign tags and
  bisection points must be pushed-tape-exact (they are, by
  construction, if we only ever evaluate on STANDARD-pushed tapes).

## Status (2026-07-15, one session in)

**Stages A-C: LANDED and green** (commit 343ae540 + follow-up;
tests `[.dmesh]`, `[.dmeshBC]`, `[.dmeshVS]`):

- Stage A surface points re-evaluate to |f| < 4e-7; interval culls
  feed proven-sign far-field corners for free.
- First extraction ever: **0 open edges, 0 non-manifold edges on
  every model**; sphere volume 0.21% from analytic, cube 1.8% low
  (rounded edges - feature points still pending).
- The refinement mechanism is real: strip ALL surface points from
  the soup and it rebuilds the sphere alone (3,170 insertions, one
  round, closed, 0.35% volume error).
- The thin-plate torture piece: **0 surface points (invisible to
  point sampling - ours OR Manifold DC's) but 16 hidden-feature
  candidate blocks flagged by the intervals.**  The stage-D
  drill-down trigger works; no sampled-SDF pipeline can see this.
- Head-to-head at 64^3 (STLs in build/dmesh_*.stl vs
  build/dcref_*.stl): ~3x FEWER triangles than Manifold DC with
  guarantees DC cannot make, at ~45x the runtime (sphere: 818 ms /
  14,936 tris vs 17.7 ms / 44,784 tris).
- The 45x is profiled, not guessed: CGAL exact-predicate fallbacks
  (Mpzf) firing on our maximally-cospherical lattice samples, plus
  full-triangulation edge iteration per refinement round.  Stage-D
  leads, in order: interior sample thinning to a near-surface shell
  + sparse far field, incremental edge scanning (only cells touched
  since last round), and batched-oracle sign evaluation for
  all-surface cells.  **Jitter: tried and REVERTED** - deterministic
  sub-cell jitter measured slower AND broke manifoldness (cube grew
  110 non-manifold edges: jittered samples landing within epsilon
  of flat faces make degenerate slivers next to their own bisected
  surface points).  Any future jitter needs a minimum-distance-to-
  surface guard (the interval oracle can provide one).

**Feature points: LANDED** (Keeter step 3).  Candidate cells
(>= 3 crossing edges) probe normals at their crossings (batched
central differences); a normal spread past ~25 degrees marks a
crease and a singular-value-clamped SVD QEF places the feature
point, rejected if it escapes its cell.  Verified: the sphere
grows ZERO feature points (no false creases on smooth surfaces);
the cube grows 320 sitting on the surface to 2e-7 and its volume
error drops 1.76% -> 0.62%; union creases get points to |f| <=
1.3e-3.  Nate's eyeball findings (chamfered cube edges, blisters
at edges, rough union seams) were all this one missing organ.

**The eyeball round** (Nate reviewing STLs drove four fixes):
sawtooth cube edges -> feature vertices now REPLACE their cells'
crossings (DC semantics; suppression gaps get rebuilt by the
refinement loop, which thereby became load-bearing); half-sharp
cube -> the test cube was grid-aligned and coincident inserts were
clobbering sign-witness info (guards added; both aligned and
de-aligned cubes ship in the showcase); remaining aligned sawtooth
-> samples with f == 0 exactly now enter the triangulation AS
surface vertices.  Verdict: both cubes 'PERFECT'.

**The one open quality problem: curved creases.**  Symptoms, all
one neighborhood: ~20 gentle warts on the spheres-union
intersection (96^3), 'chipped' csg crease, csg 91 sharp folds
(dot < -0.2) + 3 non-manifold edges at 48^3 - the csg crease is a
sphere meeting a GRID-ALIGNED cut plane (z=0.2), both hard cases
stacked.  A fold/wart detector now runs in [.dmeshVS] and
(graduated sharp/medium/gentle) in [.dmeshSTL].  Reading the 96^3
numbers taught the next spec: ANGLE ALONE CANNOT TELL A CREASE
FROM A WART - the cube's 444 medium folds are its twelve intended
edges (~37 segments each), the union's 84 mix legitimate crease
segments with Nate's ~20 warts, and even the visually-perfect
sphere carries 128 sub-visual sliver folds.  The next detector
validates folds against the surface itself: a wart is a fold whose
edge midpoint is off-surface (an |f| probe - the oracle gives it
for free); a crease is a fold that belongs.  Design direction: feature points along
a curved crease need to act as an ordered CHAIN (crease polyline),
and refinement inserts near a chain should harmonize with it
(snap-to-crease or keep-out) instead of fighting it.

**The wart hunt, round two** (error-driven repair landed; two
approaches measured and rejected): the repair loop (fold edges
with off-surface midpoints -> Newton-project -> insert -> repeat)
KILLED the csg model's non-manifold edges and proved by exclusion
that the union's residual ~20 warts are ON-surface crease
crowding - refinement rebuilds landing a hair from the chain.
Chain mediation (outlier feature points pulled to their neighbors'
segment) found nothing to mediate; feature keep-out (drop rebuilds
near the chain) broke closure - 79 open edges - because those
rebuilds perform topological separation.  Closure outranks
cosmetics; keep-out reverted.  The open design: separation points
near a crease need to exist AND cooperate with the chain
(candidate: snap the rebuild onto the crease direction rather than
the segment's f-zero, keeping separation while joining the chain).

**The chord theorem** (Nate's diagnosis, proven by two failed
cures): the union warts are triangles with one vertex on the seam
and one on EACH sphere - the sphere-to-sphere edge is a chord
crossing the crease, every vertex perfectly on-surface, the chord
itself sagging off it.  Gate-free error-driven repair (check ALL
edge midpoints, press down the deviant) fixed the CUBE to
exactly-analytic volume (1.728) but CANNOT converge on shallow
creases: a crease-crossing chord is SELF-SIMILAR - deviation and
length shrink together under midpoint splitting, so the ratio
never passes a relative tolerance; crease-seeking splits (seed at
the |f| peak along the chord, the branch switch) reduce but do not
eliminate the residue (60 -> 19 pinches at cap).  Both variants
trade ~20 gentle warts for non-manifold pinches: reverted, closure
outranks cosmetics.  THE PRACTICAL CURE
(landed - Nate refused to let it die): crease-seek + a CROWDING
GUARD - no repair insert within a quarter edge of an existing
vertex (CGAL nearest_vertex).  The guard converts self-similar
churn into bounded convergence: the union's warts press down in 5
rounds, 28 repairs, ZERO topology damage; csg improves 244 to 187
sharp folds with zero pinches (was 3,377 inserts and 60 pinches
unguarded).  Gate-free crease-seek is now the DEFAULT
(STIBIUM_DMESH_REPAIR=1 for the conservative fold-gated mode).
The dial sweep (Nate's second wind): tolerance 5% / guard 0.25 =
0 pinches, shallow repair; 2% / 0.15 = pinches return (tighter
tolerance floods the crease, needing MORE guard, not less);
**3% / 0.25 / 16 rounds = the sweet spot**, now baked in: the
union presses to 36 repairs with ZERO pinches and ZERO sharp
folds; csg converges at 13 rounds keeping exactly 2 pinches at its
sharp point - structurally beyond greedy repair, the concrete
CDT_3 case.
The COMPLETE cure for csg's residue remains crease topology - and
CGAL 6.2 ships Conforming_constrained_Delaunay_triangulation_3,
so the feature polyline can be genuinely constrained edges.  Next
design round, foundation confirmed present.

Stage D remains: the crease-polyline constraint design (above),
hidden-candidate drill-down, error-driven adaptive insertion,
performance (CGAL predicate fallbacks - jitter tried and reverted,
see above; next: shell-thinning + incremental edge scans), MT, and
the app-facing flag.  Also noted: the delaunay mesher already rides
STANDARD-pushed tapes for descent/sampling; bisection batches still
use the base tape (bit-identical values; switching to covering
pushed tapes is a pure speed move when performance round opens).

## The constrained-crease round: LANDED (2026-07-15, same night)

**Status: default-on** (`STIBIUM_DMESH_CCDT=0` reverts to the plain
DT path).  Every showcase model now meets or beats the DT baseline:
sphere/cube identical (cube carries 428 constrained edges at zero
cost), cube_aligned medium folds 840 -> 432, spheres CLEAN WITH
ZERO REPAIRS (DT needed 36 - the union crease is now separated by
construction, the chord theorem case closed structurally), csg
gated back to DT semantics (2 pinches) pending the junction round.
Full suite green, 627,666 assertions.

**Construction facts** (hard-won, do not rediscover):
- The documented CCDT_3 wrapper only builds from a finished PLC;
  incremental `insert()` / `insert_constrained_edge()` /
  `restore_Delaunay()` live on
  `Conforming_constrained_Delaunay_triangulation_3_impl<T_3>`,
  which IS the triangulation.  T_3 must be a
  `Delaunay_triangulation_3` (the impl uses its conflict-region
  machinery; plain `Triangulation_3` lacks `Conflict_tester_3` -
  the docs' claim is about the accessor type).  Delaunay under the
  hood also keeps `nearest_vertex()` for the crowding guard.
- Info bases STACK: both CCDT bases take their base class as the
  second template parameter.  `with_info` does NOT initialize the
  info, and the machinery's Steiner vertices arrive without one:
  sweep `is_Steiner_vertex_on_edge()` -> info 0 (they sit on the
  crease = on the surface) after EVERY insert batch.
- Steiner-coincides-with-vertex THROWS (no reuse/retyping);
  `set_segment_vertex_epsilon` is only a validator.  Exactly
  on-segment vertices are incorporated into the polyline by the
  segment traverser at insert_constrained_edge time.

**The measured pinch mechanism** (provenance dump, every single
non-manifold edge of the first runs): a separator insert (refine
round, repair round, or the QEF duplicate clump) landing 0.15-0.25
cells from a chain/Steiner vertex -> sliver tets -> 4-triangle
edge.  The five landed counter-measures, in causal order:
1. **Oracle-refereed segments**: |f|/|grad| at 1/4, 1/2, 3/4 vs 3%
   of length; mislinked shortcut segments never become law.
2. **The trust gate**: >10% rejected segments = untrusted chains =
   no constraints at all for that model (csg: 35% rejected - two
   creases crossing at the sharp point mislink the extractor; a
   half-trusted polyline measured WORSE than none: 7 pinches).
3. **Crease-band drop** (0.35 cells around accepted segments):
   bisected surface points AND non-chain feature duplicates (26-51%
   QEF clump on curves!) are redundant with the polyline; each
   duplicate near a segment encroaches it and forces a Steiner at
   its projection = guaranteed near-coincidence.  Steiner count on
   spheres fell 131 -> 65, aligned cube 492 -> 0.  Sign witnesses
   are never dropped.
4. **Repair keep-out** (0.75 cells): repair pressing points onto a
   constrained edge destroys it and the re-conform lands a sliver
   (fold-gated mode escalated to a 1792-repair war with 716
   Steiner).  The constraint owns the crease; repair owns the
   smooth field.
5. Refinement inserts stay unguarded (closure outranks cosmetics);
   `STIBIUM_DMESH_SLIDE` (default off) can slide near-coincident
   separators away, but stalls on short band edges (slid point
   exits the conflict zone -> holes at cap).  Rejected for now.

**Two dead ends, measured and reverted**: (a) shadowing band edges
+ mixed-cell-oracle signs + emit-and-project extraction - closed
but the oracle sign alternation around band edges makes its own
pinches (2-6) plus 19 sharp folds from skinny projected triangles;
(b) fold-gated repair under constraints - see the war above.
Mixed-sign cells now ask the centroid oracle (kept: strictly more
correct at halt/cap exits than first-witness-wins).

**Same-night follow-ups (both landed, suite green)**:
- **The noise floor** (1edb78f4): csg's 96 "rejected" segments were
  never shortcuts - all on the plane circle, failing a purely
  relative 3% bar because their QEF ENDPOINTS carry ~1e-3 placement
  noise and the densest segments are only 0.5-0.9 cells long.
  Tolerance is now max(3% len, 5% spacing); csg accepts all its
  segments, passes the gate, meshes fully constrained with ZERO
  off-surface warts and 7 sub-visual pinches (Nate accepts pinches
  when visual quality holds).  `STIBIUM_DMESH_SEG_DEBUG=1` dumps
  rejected segments - USE IT before believing a rejection theory.
- **The junction-split walk** (4889d596): covariance classification
  (eigenvalue ratio > 0.2 over the 2-cell neighborhood) cuts
  junction reps out before linking and reattaches them as SHARED
  chain endpoints.  Cube: exactly 8 corners classify, 12 edge
  chains share them (previously two of three edges per corner
  stopped one rep short).  LIMIT: csg's seam-meets-circle crossings
  are TANGENTIAL - locally parallel curves look like a line to
  covariance, so no junction fires; its walk still passes through.
  Coverage 97.5% regardless.

**Diagnostics live**: `STIBIUM_DMESH_NM_DEBUG=1` prints every
non-manifold edge with vertex provenance (sample/bisect/feature/
refine/repair + steiner + n-constraints) - it convicted every pinch
mechanism tonight; `constrained`/`steiner` counts are in DMesh and
[.dmeshSTL] output.  CCDT timing note: sphere 953 -> 1612 ms at
64^3 in [.dmeshVS] - the machinery (time stamps, hierarchy) costs
~70%; part of the performance round.

**Research round (2026-07-15, two opus agents; full reports in
doc/research/2026-07-15-*.md)**, the load-bearing findings:
1. **Our extraction rule is restricted-Delaunay class and NOT
   manifold-by-construction.**  Manifoldness is conditional on the
   Topological Ball Property (Edelsbrunner-Shah 1997); a pinch IS a
   TBP violation (a Voronoi edge stabbing the surface twice).
   TetWeave/Wang never pinch because they use MARCHING-TET
   interpolation (manifold by construction) - which cannot pass
   through feature vertices, i.e. cannot do sharp creases.  Our
   rule is the price of sharpness.
2. **The principled pinch cure for our exact case is Manifold-DC
   style vertex duplication** (Schaefer/Ju/Warren TVCG 2007): when
   one vertex is forced to serve two surface sheets, emit one
   vertex per sheet and rewire - near-coincident but topologically
   distinct.  Literature explicitly warns AGAINST snapping/welding
   in the triangulation (we measured that too: every snap variant
   made sliver pairs).  Post-extraction weld is the accepted lossy
   backstop for artifact pinches.
3. **Feature protection has two schools**: weighted protecting
   balls (Cheng-Dey-Ramos; CGAL Mesh_3) need a REGULAR
   triangulation - does not compose with constrained CDT; the
   constrained school is Shewchuk diametral-sphere reject-and-split
   - our crowding guard/keep-out/band-drop are ad-hoc versions of
   its insertion-radius bound.
4. **The CSG tree is the unclaimed lever** (next campaign,
   'crease tracing'): creases are pairwise primitive intersections
   f_A = f_B = 0; trace with predictor t = grad_A x grad_B +
   2-eq Newton corrector; junctions detected ANALYTICALLY (third
   field crossing zero, or |grad_A x grad_B| -> 0 = TANGENCY -
   exactly csg's sharp point, terminate and share the node).
   Trimming = evaluate the full CSG f.  Gotchas mapped: smooth-min
   blends (skip - no crease), nested subtrees (recurse; new seams
   only leaf-vs-leaf), coincident faces (exclude).  This replaces
   normal-spread QEF heuristics with exact geometry and makes
   junction topology bookkeeping instead of inference.

## The manifold pass + crease tracer (2026-07-15, later the same day)

Both research recommendations LANDED (commits 8577cb88, 085664ac);
suite green throughout.

**Manifold pass** (STIBIUM_DMESH_MANIFOLD, default on): Manifold-DC
style vertex duplication.  Around a pinch edge, the two extracted
facets bounding the same INSIDE run of the cell ring are one wedge
of the solid = one sheet; facets pair by inside-run, and any vertex
whose facet fan is disconnected splits into one coincident output
vertex per sheet.  Zero geometric change; 2-manifold at every
vertex by construction.  Gotchas paid for: fan adjacency must read
a SNAPSHOT of the triangle list (per-vertex rewrites perturb later
edge keys - measured as 6 open edges); unrecoverable ring pairings
must KEEP the pinch (union everything), never tear.  DMesh
.split_verts reports.

**The crease tracer** (STIBIUM_DMESH_TRACE, default on;
delaunay_trace): the campaign's crown.  Tape APIs (tape.h):
tape_pairs enumerates min/max clauses; tape_eval_r_prefix +
tape_ctx_r_row read a pair's operand values f_A, f_B with all
upstream transforms applied - prefix [0, clause) evaluation works
because linear-scan register allocation keeps operand slots live
until their reader executes.  March {f_A=0, f_B=0} with the SSI
predictor-corrector (tangent = grad A x grad B by central
differences, minimal-norm 2-eq Newton corrector), seeds = QEF
feature points, trim by the full oracle (|f|/|grad f| > 0.05 sp =
a third branch owns the surface = the junction).  Junction
endpoints bisect onto the corner with a TIGHT tolerance (5e-4 sp -
bisecting against the marching trim leaves endpoints trim-short;
measured worst |f| exactly 0.05 sp before the fix), then cluster
(0.75 sp) to a shared Newton-projected corner vertex, so incident
curves carry IDENTICAL coords and the triangulation's coincidence
guard makes them one constrained vertex.

Failure modes measured and guarded (all three bit on csg):
- tangency stall: the corrector pulls every predictor step back to
  the same point, emitting float-noise micro-segments (three
  "points" 1e-15 apart threw CGAL's vertex-on-constraint error) ->
  minimum-progress check per step;
- duplicate traces: seeds converge onto an already-traced curve
  from beyond the consumption radius (Newton moves seeds up to 1.5
  cells) -> reject by converged-seed-position against same-pair
  polylines;
- anything else the conforming machinery hates -> try/catch around
  the constrained path, falling back to the plain DT mesher
  (degrade, never crash).

Traced chains REPLACE the QEF chain extractor when present; QEF
features stay as plain surface points (the crease band drop
retires the redundant ones).  Referee [.dtrace] (cube 4 pillar
edges + 2 closed face loops = all 12; union 1 loop at |f|=1e-7;
csg seam arc + closed cut-boundary loop traced through its
branch-switch kinks, worst |f| 7.8e-4 from kink-straddling central
differences - under the 1.3e-3 QEF noise).  Showcase 96^3: EVERY
model 0 open / 0 non-manifold, 0 repairs, 0 splits, csg fully
constrained (483 edges) with no rescue machinery at all.
STIBIUM_DMESH_TRACE_DEBUG=1 dumps per-polyline stats.

Remaining tracer niceties (not blockers): kink-point refinement
(solve the 3-field corner system instead of accepting loose Newton
at branch switches), seeds independent of QEF features (interval-
guided seeding would catch creases QEF misses entirely, e.g. below
lattice pitch), OP_ABS as a crease generator (|x| kinks at x=0;
enumerate ABS clauses the same way), and blend-aware skipping once
smooth-min ops exist.

**Nate's eyeball verdict on the traced build** (2026-07-15):
sphere "transcendent" (1 MB, Blender-regular topology, vs 3.5 MB
optimized / 144 MB raw marching-cubes); the two-sphere weld crease
on spheres AND csg "perfect - no visual disruption whatsoever."
**CHIP ROUND 2 (a524e2f9): two structural cures measured DEAD -
read before re-attempting.**  (1) Crease-edge exposure (flip the
ring cell with the deepest wrong-side centroid where a crease
edge's cell ring has zero sign transitions): only 3-12 flat rings
exist vs 17 chip chords, and ALL are razor wedges with no
legitimately-inside cell - zero flips ever fired.  The buried-edge
theory is a minority mechanism.  (2) Opening the keep-out annulus
(0.10 segment shield + 0.3 sp crease-vertex-only crowding floor;
justified because blocked repairs all land 0.2-0.3 sp out and
sliver pairs are a VERTEX-distance phenomenon): 36 more repairs
land, detection plateaus at ~20 - fresh shallow chips replace
pressed ones one-for-one.  Greedy repair is at its asymptote for
concave corners.  THE REMAINING DESIGNED CURE, for a fresh round:
extraction-level - the literature's answer is the local
Topological-Ball-Property check (Dey-Levine DelIso): a facet
whose dual Voronoi edge crosses the surface more than once is the
divot; refine it by inserting its surface-Delaunay-ball CENTER (a
point ON the surface, oracle-projected), which is a principled
insert location unlike our chord-derived targets.  Alternatively
DC-semantics facet substitution at the wedge.  Both need the
divot FACETS identified first (sharp fold edges with f<0
midpoints already locate them - reuse the wart machinery's
classification, not its repair).

**CHIP ROUND 3 (b0769526): TBP/SDB refinement built, measured
EQUAL, defaulted off.**  The DelIso cure is implemented and works
mechanically (tri_cells at extraction, kernel circumcenters -
dt.dual() needs the Delaunay cell base our info stacks lack, dual-
segment oracle sampling, SDB-center inserts bypassing the segment
keep-out) - but the new WORST-CHIP-DEPTH metric (the honest one;
the count is polluted by an evergreen marginal-chord population,
128 even on the perfect sphere) shows an identical trajectory to
plain midpoint repair: 0.500 -> 0.196 sp plateau on csg at 96^3.
Three repair strategies now hit the same plateau; it is set by the
crowding guard + sampling density, NOT target choice.  Conclusion
for the next attempt: the ~0.2 sp residue needs either (a) finer
LOCAL density near concave creases (adaptive band refinement of
the initial lattice - the sampler's octree could descend one extra
level in crease bands), or (b) the DC-semantics extraction change.
Repair itself is exhausted - three strategies, one plateau, stop
pulling that lever.  STIBIUM_DMESH_SDB=1 re-arms the SDB path;
depth prints in STIBIUM_DMESH_CHIP_DEBUG.

**CONCAVE CHIPS: round landed same day** (918f3c89).  Mechanism
confirmed by count (STIBIUM_DMESH_CHIP_DEBUG): 87 of csg's 90
repair candidates were material-crossing chords (f < 0 at
midpoint) and the keep-out blocked ALL of them (candidate
mechanism 2 - repair starvation).  Fix: chips are a distinct
repair species - midpoint-projection onto the nearest face
(crease-seek would target the corner the constraint owns) with a
0.35-cell keep-out.  Sweep: 0.75 starves (0 fixed), 0.15
escalates (near-crease inserts spawn new sliver chords, 17 -> 31,
pinch splits appear), 0.35 presses 87 -> 17 monotonically, zero
topology damage.  Residue = the shortest corner-hugging chords;
the structural cure is extraction-level (facets forced to use
crease edges at the wedge - DC semantics), a designed future
round if eyeballs demand it.

**The standalone technical record lives in doc/DELAUNAY-MESHER.md**
(how it works, sources, uniqueness claim, everything tried and
rejected, measured state).  Keep it current when rounds land.

## The original round map (2026-07-15, kept for the record)

**Primitive confirmed**: CGAL 6.2's
`Conforming_constrained_Delaunay_triangulation_3` exposes
`insert_constrained_edge(va, vb)` returning a
`Constrained_polyline_id`, with conforming Steiner insertion.
Crease chain segments as constrained edges make crossing chords
geometrically impossible - the complete cure for csg's sharp-point
residue and the last two pinches.

Plan, referee-gated as always:

1. **Chain extraction** (CGAL-independent, buildable first):
   order feature points into polylines by nearest-neighbor walking
   with the collinearity test from mediation (corners break
   chains).  Referee on known topology: the de-aligned cube yields
   12 open chains (edges) meeting at 8 corners; the union yields
   ONE closed loop; csg one loop + the plane-circle.  Chain count,
   closure, and length statistics are all checkable.
2. **CDT_3 integration**: swap/augment the DT with the CCDT_3
   class (check vertex-info support in its Tds; the extraction and
   refinement loops port over).  Insert chains as constrained
   edges BEFORE the refinement loop; Steiner points the CCDT adds
   along constraints must be projected onto the crease (they land
   on the straight segment, not the curve - project via the
   crease-seek sampler).
3. **Referee**: csg sharp point clean (0 pinches, folds -> crease-
   legitimate only); union loop constrained (warts structurally
   impossible); cubes unchanged; suite green; STLs to Nate -
   glorious, for Keeter.

Also banked from the wart campaign: the mental model, stated
plainly (Nate's formulation, verified): f-rep -> point soup ->
Delaunay -> COMPARE the mesh against the f-rep via the oracle
(|f|/|grad| at edge midpoints = exact local deviation) -> insert
on-surface points where it sags -> repeat.  Two nested loops share
the oracle: topological (inside/outside separation) then geometric
(surface fidelity).  Vertices are never moved, only added; the
mesh converges by densification-where-needed.

## Queued rounds

- **Flat-face decimation** (Nate, 2026-07-15): a final output pass
  that detects planar regions and re-triangulates them - massive
  tri reduction on the cube and mechanical/prismatic parts with
  zero shape change.  Sketch: group triangles by normal clustering
  seeded from the oracle (a patch is flat iff every interior
  vertex's |f| stays ~0 under the patch plane - the oracle can
  CERTIFY flatness, unlike normal-only heuristics, which is how
  gentle curves survive); patch boundary = crease chains + patch
  rim; re-triangulate via 2D CDT in the plane.  Isolated,
  approachable, high yield.  NOT YET STARTED.

## Open questions for Nate

- CGAL as a dependency (prototype: yes; shipping: revisit).
- Where the new mesher lives: `lib/fab/src/mesh/delaunay/`?
- Does detect-features (Kobbelt) stay as an option forever, or is
  feature capture via step-3 points expected to replace it?

**Repair-budget experiment (2026-07-15 night, REFUTED with
numbers):** zeiss d1 at STIBIUM_DMESH_ROUNDS=48 vs 16: 3.5x the
repairs, 2.4x the triangles, worst depth IDENTICAL (0.599 sp,
same point), and non-manifold edges exploded 147 -> 1747 (pinch
splits 3.8K -> 34K).  The 16-round cap is a safety rail, not a
budget shortage - extra rounds churn self-similar sliver chips
into topology damage (the csg keep-out escalation at model
scale).  The v1-kindness hypothesis dies with it.  The rough-
interface residual is the CROWDING interaction: stacked features
1-2 cells apart whose corridors/keep-outs overlap and starve the
strip between.  Nate's boss models (1mm cylinder on flat, 0.5/2mm
variants) are the referee for that fight.

**Boss-model bisection (Nate's emboss/engrave stable, 2026-07-15
night):** cylinders stacked/sunken at 0.5/1/2mm steps, r1.
CROWDING THEORY DEAD: 1mm and 2mm steps (emboss AND engrave) mesh
PERFECT - 0.000 sp, ~5K tris, zero repairs.  The law is sharp:
features >= 1 lattice cell are solved; BELOW one cell they do not
fail gracefully - repair churn sprays 68K inserts and 28x
triangles for a metrically-tiny finish (the zeiss rough
interfaces and smeared inlay in one sentence).  The dense band
RESCUES sub-lattice features: engrave 0.5mm at DENSE=1 = 16
chips / 0.000 sp / 16K tris; emboss 0.5mm needs DENSE=2 (58K
tris; lattice-phase asymmetry, dig later).  Grand unification:
ambiguity-triggered local density cures BOTH residual diseases
(blend aliasing + sub-lattice steps).  Remaining stage-D work:
(1) automatic per-leaf level selection (feature-scale detection
instead of the global env knob), (2) churn suppression when a
feature stays unresolvable (stall detection, stop early).
Repair-budget theory also refuted same night (48 rounds = 12x
topology damage, zero depth gain - 16 is a safety rail).

**Performance round 1 (2026-07-15 night): profile then aim.**
STIBIUM_DMESH_TIME profile of zeiss d1 convicted the TRACER (957 s
of 1387 s - 69%; the last suspect on the intuition list, which is
why the house rule applies to speed).  The seed gate
(tape_eval_f_pairs: one recording walk per seed captures every
pair's operand values in flight; Newton only where both fields
< 8 sp) cut it to 9.57 s - 100x - with byte-identical output
(12,235 constraints, same tris/nm/tents).  Zeiss d1: 23 min ->
7 min.  Remaining pie: extract+repair 272 s (64%; early-exit +
incremental re-extraction is the designed next round), insert
85 s, snap 26 s.  Gate prunes 16M exhaustive (pair, seed) Newton
attempts to 257K candidates (1.6%).

**Perf round 2 (same night): repair stall exit.**  Two flat
rounds = churn ahead, stop (round >= 3; stash keeps chips for the
snap).  extract+repair 272 s -> 34 s; zeiss d1 TOTAL 23 min ->
2.7 min across the evening (seed gate + stall exit).  Bonus: the
churn rounds were the damage - nm 147 -> 50, tris 1.04M -> 407K
(600K churn vertices), depth unchanged, constraints identical.
AWAITING Nate's eyeball on v2fast vs the reviewed d1 (v1-kindness
caveat: 17.5K repairs vs 333K).  Next perf target: insert points
(84 s, now 52% of the pipeline) - CGAL single-threaded insertion;
then threading the eval side.

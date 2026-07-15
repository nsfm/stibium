# The meshing campaign: adaptive Delaunay on sound intervals

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
outranks cosmetics.  THE cure is crease topology: extract the
feature chain as an ordered polyline and make its segments
CONSTRAINED edges of the triangulation (or retriangulate the
crease band against the polyline), so no chord can cross.  That is
next session's field-mapped design.

Stage D remains: the crease-polyline constraint design (above),
hidden-candidate drill-down, error-driven adaptive insertion,
performance (CGAL predicate fallbacks - jitter tried and reverted,
see above; next: shell-thinning + incremental edge scans), MT, and
the app-facing flag.  Also noted: the delaunay mesher already rides
STANDARD-pushed tapes for descent/sampling; bisection batches still
use the base tape (bit-identical values; switching to covering
pushed tapes is a pure speed move when performance round opens).

## Open questions for Nate

- CGAL as a dependency (prototype: yes; shipping: revisit).
- Where the new mesher lives: `lib/fab/src/mesh/delaunay/`?
- Does detect-features (Kobbelt) stay as an option forever, or is
  feature capture via step-3 points expected to replace it?

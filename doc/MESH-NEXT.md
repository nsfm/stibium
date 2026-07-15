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

**The junction round (next design work, csg is the referee)**: the
chain extractor treats crease CROSSINGS (union seam meets plane
circle at csg's sharp point) as ordinary degree-2 walking and
mislinks ~35% of segments.  Needed: junction detection (degree>=3
or crossing-aware merge), chains split at junctions, junction
vertices shared between chains as polyline endpoints.  Then csg's
gate opens and its sharp point gets the structural cure the union
already enjoys.  Diagnostics live: `STIBIUM_DMESH_NM_DEBUG=1`
prints every non-manifold edge with vertex provenance
(sample/bisect/feature/refine/repair + steiner + n-constraints);
`constrained`/`steiner` counts are in DMesh and the [.dmeshSTL]
output.  CCDT timing note: sphere 953 -> 1612 ms at 64^3 in
[.dmeshVS] - the machinery (time stamps, hierarchy) costs ~70%;
part of the performance round.

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

## Open questions for Nate

- CGAL as a dependency (prototype: yes; shipping: revisit).
- Where the new mesher lives: `lib/fab/src/mesh/delaunay/`?
- Does detect-features (Kobbelt) stay as an option forever, or is
  feature capture via step-3 points expected to replace it?

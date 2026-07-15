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

## Open questions for Nate

- CGAL as a dependency (prototype: yes; shipping: revisit).
- Where the new mesher lives: `lib/fab/src/mesh/delaunay/`?
- Does detect-features (Kobbelt) stay as an option forever, or is
  feature capture via step-3 points expected to replace it?

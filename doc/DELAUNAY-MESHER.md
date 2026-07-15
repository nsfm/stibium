# The Stibium Adaptive-Delaunay Mesher

*Technical record, written 2026-07-15 at the close of the campaign's
third day.  Companion documents: doc/MESH-NEXT.md (the living field
map and status ledger), doc/TAPE-DESIGN.md (the tape kernel this
feeds on), doc/research/2026-07-15-junction-extraction.md and
doc/research/2026-07-15-pinch-manifoldness.md (commissioned
literature surveys, full text).*

## What it is

A meshing pipeline for f-rep (function representation) solids that
produces watertight, 2-manifold, adaptively-sized triangle meshes
with **exactly sharp feature creases** - crease polylines recovered
as analytic curves from the model's own CSG expression, not
inferred from samples.  At 96^3 effective resolution the showcase
models mesh with zero open edges, zero non-manifold edges, zero
repair interventions, at roughly one third the triangle count of
the dual-contouring reference at visibly higher quality (the sphere:
1 MB vs 3.5 MB optimized / 144 MB raw marching-cubes, with
regular near-equilateral triangulation).

## Why this combination is new

Each ingredient exists in the literature.  The fusion does not,
because it depends on capabilities that co-exist only in this
kernel:

1. **A sound interval oracle.**  Stibium's evaluator computes
   interval bounds that are *conservative under IEEE semantics*
   (including NaN-taint tracking), so "this box is provably
   empty/full" is a theorem, not a heuristic.  Sampling happens
   only where the surface may be; culled boxes contribute proven
   sign witnesses for free.  The 2025 sampled-SDF pipelines (Wang
   et al., TetWeave) consume point samples of unknown soundness by
   construction.
2. **Immutable pruned tapes.**  Evaluation runs on refcounted,
   region-specialized instruction tapes (doc/TAPE-DESIGN.md).
   Beyond speed, the tape is *introspectable* - which enables:
3. **Exact crease curves from the expression.**  Every min/max
   clause in the tape is a crease generator: the surface kinks
   exactly where its two operands are equal.  A prefix evaluation
   of the tape - guaranteed valid because linear-scan register
   allocation keeps operand slots live until their reading clause
   executes - exposes the operand fields f_A, f_B *with every
   upstream coordinate transform already applied*.  The crease is
   then the solution curve {f_A = 0, f_B = 0}, marched with the
   classical surface-surface-intersection predictor-corrector to
   float precision.  Junctions are not inferred: the full oracle
   trims each curve to the boundary, and the trim boundary IS the
   junction.  No sampled-geometry pipeline can do this; it
   requires owning the compiler.
4. **Conforming constrained Delaunay.**  The traced polylines are
   inserted as constrained edges (CGAL 6.2's
   Conforming_constrained_Delaunay_triangulation_3), making
   crease-violating chords structurally impossible rather than
   iteratively repaired.

The result is, to our knowledge, the first f-rep mesher that is
simultaneously watertight-by-construction, manifold (guaranteed by
an explicit output pass), adaptive, thin-feature-aware (interval
flagged), and *exactly* sharp along CSG creases.

## The pipeline

**Stage A - point soup** (`delaunay_sample`).  Octree descent with
STANDARD tape pushes; interval-culled boxes contribute corner
samples with proven signs; leaf blocks sample a corner lattice;
sign-change lattice edges are bisected (16 rounds, batched
MIN_VOLUME evaluations) into surface points.  Interval-suspicious
blocks where every sample agrees in sign are counted as
hidden-feature candidates (thin features below the lattice pitch -
invisible to any sampling-only method, visible to sound
intervals).

**Feature points** (Kobbelt-style, retained as seeds).  Cells with
>= 3 crossings probe normals by batched central differences; a
normal spread past ~25 degrees marks a crease; a singular-value-
clamped QEF places the feature point; crossings in feature cells
are suppressed (dual-contouring semantics).  Since the tracer
landed, these serve mainly as tracer *seeds* and as fallback
chains for non-min/max creases.

**The crease tracer** (`delaunay_trace`).  For each min/max clause
(via `tape_pairs`): Newton-converge seeds onto {f_A = 0, f_B = 0}
(minimal-norm two-equation steps from 7-tap central differences);
march along t = grad f_A x grad f_B at half-cell steps, re-
correcting each step; trim by the full oracle (|f|/|grad f| beyond
5% of a cell means a third branch owns the surface past here);
bisect the endpoint tightly onto the corner; cluster endpoints
within 0.75 cells into shared, Newton-projected corner vertices.
Closed loops detect by returning to the start.  Guards: minimum
progress per step (tangency stall), duplicate-trace rejection by
converged seed position, degenerate-segment suppression.

**Stage B - triangulate + separate** (`delaunay_mesh_soup`).  All
signed samples and surface points enter a conforming constrained
Delaunay triangulation (vertex/cell info stacked under the CCDT
bases; sign witnesses guarded against coincident-insert
clobbering).  Chain segments are oracle-refereed (|f|/|grad| at
1/4, 1/2, 3/4 against max(3% length, 5% cell) - endpoints carry
noise, so short segments need the absolute floor), then inserted
as constrained edges; a >10% rejection rate distrusts the chains
entirely (fall back unconstrained).  Refinement: any tet edge
joining an inside vertex directly to an outside vertex gets a
bisected surface point until none remain.  Redundant surface
vertices within 0.35 cells of a constrained segment are dropped
(they pair with chain/Steiner vertices into pinch slivers);
machinery Steiner vertices are swept to surface-info after every
insert batch (their info arrives uninitialized; they sit on the
crease).

**Stage C - extract + repair.**  Facets between opposite-signed
cells whose three corners are surface vertices, oriented by the
outside cell; mixed-sign cells ask the oracle at their centroid.
Error-driven repair (after Wang et al.): every edge midpoint is
tested against the oracle; off-surface warts split at the |f| peak
along the chord (the crease-seek - a plain midpoint split of a
crease-crossing chord is self-similar and never converges), with
a crowding guard (no insert within a quarter edge of an existing
vertex) and a crease keep-out (0.75 cells - the constraint owns
the crease).  Dials 3% / 0.25 / 16 were swept, not guessed.

**The manifold pass.**  Restricted-Delaunay-style extraction is
not manifold by construction (its manifoldness is conditional on
the Topological Ball Property; a pinch is a TBP violation).  The
tet complex knows the truth: around a pinch edge, the two facets
bounding the same inside-run of the cell ring are one wedge of the
solid.  Facets pair by inside-run; any vertex whose facet fan is
disconnected splits into one coincident output vertex per sheet
(Manifold-DC style).  Zero geometric change; 2-manifold at every
vertex by construction; unrecoverable pairings keep the pinch
rather than tear a hole.

## Sources

Blueprint and immediate lineage:
- M. Keeter, "Please Steal My Meshing Algorithm Idea," 2026-07-03
  (the staged plan this campaign implements and extends).
- Wang et al., "Power Diagram Enhanced Adaptive Isosurface
  Extraction from Signed Distance Fields," arXiv:2506.09579, 2025
  (error-driven insertion; the repair loop's ancestor).
- Binninger, Wiersma et al., "TetWeave," ACM TOG (SIGGRAPH) 2025
  (extraction-guarantee framing).
- Dey & Levine, "Delaunay Meshing of Isosurfaces," SMI 2007 (the
  restricted-Delaunay extraction class we use).
- Edelsbrunner & Shah, "Triangulating Topological Spaces," IJCGA
  1997 (the Topological Ball Property - why pinches happen).

Sharp features and crease curves:
- Kobbelt, Botsch, Schwanecke, Seidel, "Feature Sensitive Surface
  Extraction from Volume Data," SIGGRAPH 2001; Ju, Losasso,
  Schaefer, Warren, "Dual Contouring of Hermite Data," SIGGRAPH
  2002 (QEF feature points, crossing suppression).
- Schaefer, Ju, Warren, "Manifold Dual Contouring," TVCG 2007 (the
  vertex-duplication cure implemented in the manifold pass).
- Bajaj, Hoffmann, Lynch, Hopcroft, "Tracing Surface
  Intersections," CAGD 1988; Patrikalakis & Maekawa, "Shape
  Interrogation for CAD/CAM," 2002 (the SSI predictor-corrector
  the tracer implements).
- Gumhold, Wang, McLeish, 3DPVT 2001; Pauly, Keiser, Gross,
  Eurographics 2003 (covariance junction classification - the
  chain extractor's junction-split walk).
- Shewchuk, "Delaunay Refinement Algorithms," Comput. Geom. 2002
  (insertion-radius reasoning behind the crowding guard);
  "Constrained Delaunay Tetrahedralizations and Provably Good
  Boundary Recovery," IMR 2002.
- CGAL 6.2, Conforming_constrained_Delaunay_triangulation_3
  (Rineau; GPL) - the constrained triangulation engine.

Full annotated bibliographies: doc/research/2026-07-15-*.md.

## What was tried and rejected (measured, not guessed)

The campaign's discipline: every stage lands with its referee
green or gets reverted with the numbers written down.  The corpses
with lessons attached:

- **Deterministic sub-cell jitter** (against CGAL cospherical
  slowdowns): slower AND broke manifoldness (110 non-manifold
  edges - jittered samples make degenerate slivers against their
  own bisected surface points).
- **Feature keep-out for refinement rebuilds**: 79 open edges -
  those rebuilds perform topological separation; closure outranks
  cosmetics.
- **Gate-free midpoint repair**: exact cube volume, but
  crease-crossing chords are SELF-SIMILAR under midpoint splits
  (deviation and length shrink together; a relative tolerance
  never fires) - 3,807 inserts, 60 pinches.  The chord theorem.
- **Tighter repair tolerance with weaker crowding guard**: pinches
  return; a tighter tolerance floods the crease and needs a
  STRONGER guard, not weaker.
- **Half-trusted chain constraints** (csg pre-tracer): worse than
  none (7 pinches vs 3) - a wrong constraint is FORCED geometry.
- **Edge shadowing + mixed-cell oracle + emit-and-project
  extraction** (avoiding separator inserts near constraints):
  closed, but oracle sign alternation around band edges makes its
  own pinches plus skinny-triangle folds.
- **Separator slide-away**: stalls on short crease-band edges (the
  slid point exits the conflict zone; holes at cap).
- **Purely relative segment-referee tolerance**: starves short
  segments - csg's densest circle segments all failed at 3.0-4.1%
  against a 3% bar because their QEF endpoints carry ~1e-3
  placement noise.  Absolute floor added.
- **Covariance junction detection on tangential crossings**: two
  locally-parallel curves look like one line; provably blind
  there.  The tracer's trim-boundary junctions replaced it.
- **Reading min/max operands from a full tape evaluation**:
  impossible - register allocation reuses slots.  The prefix-
  evaluation liveness argument is what makes the tracer cheap.

And the tracer's own measured failure modes, each now guarded:
tangency stalls emitting float-noise micro-segments (three
"points" 1e-15 apart), loose seeds re-tracing the same curve from
beyond the consumption radius, junction endpoints stopping
trim-short of corners (worst |f| exactly equal to the trim
tolerance - the number confessed its cause).

## Current measured state (96^3 showcase, 2026-07-15)

| model        | tris   | constrained | open | non-manifold | repairs | splits |
|--------------|--------|-------------|------|--------------|---------|--------|
| sphere       | 33,768 | 0 (smooth)  | 0    | 0            | 0       | 0      |
| cube         | 18,240 | 892 traced  | 0    | 0            | 0       | 0      |
| cube_aligned | 16,444 | 878 traced  | 0    | 0            | 0       | 0      |
| csg          | 21,624 | 483 traced  | 0    | 0            | 0       | 0      |
| spheres      | 38,258 | 301 traced  | 0    | 0            | 0       | 0      |

Traced-curve referee ([.dtrace], 64^3): cube = 4 pillar edges + 2
closed face loops (all 12 edges, shared corners); union = one loop
at |f| = 1.2e-7; csg = trimmed seam arc + the closed cut-face
boundary traced through its branch-switch kinks (worst |f| 7.8e-4,
under the 1.3e-3 QEF noise).  Suite green: 627,666 assertions.

Human referee (Nate, prusa view): sphere "transcendent" -
regular Blender-like topology at 1 MB vs 3.5 MB optimized
marching-cubes; the two-sphere weld crease "perfect - no visual
disruption whatsoever."

## Knobs

- `STIBIUM_DMESH_TRACE` (default 1): the crease tracer.
- `STIBIUM_DMESH_CCDT` (default 1): constrained triangulation
  (0 = plain Delaunay reference path).
- `STIBIUM_DMESH_MANIFOLD` (default 1): the manifold output pass.
- `STIBIUM_DMESH_REPAIR` (default 2 = crease-seek; 1 = fold-gated).
- `STIBIUM_DMESH_N`: showcase resolution.
- Diagnostics: `STIBIUM_DMESH_TRACE_DEBUG`, `STIBIUM_DMESH_NM_DEBUG`
  (non-manifold provenance), `STIBIUM_DMESH_SEG_DEBUG` (segment
  referee rejections).

## Open problems

- **Concave chips - largely pressed** (round landed 2026-07-15,
  commit 918f3c89): material-crossing chords are a distinct repair
  species (f < 0 at the sagging midpoint).  The repair keep-out
  was starving them (87 of 90 candidates blocked); they now take
  midpoint-projection onto the nearest face with a 0.35-cell
  keep-out (swept: 0.75 starves, 0.15 escalates).  csg presses
  87 -> 17 residual shallow chips with zero topology damage.  The
  residue: the shortest chords hugging the corner, whose projected
  repairs land inside even the sliver radius - structural cure
  would be extraction-level (facets forced to use crease edges,
  DC-semantics at the wedge), not repair-level.
- **Kink corners**: branch-switch points inside a traced loop get
  loose Newton convergence (gradient flip straddles the kink);
  the principled fix is a 3-field corner solve.
- **OP_ABS creases**: |x| kinks at x = 0 - same generator pattern,
  not yet enumerated.
- **Seeding without QEF**: interval-guided seeds would catch
  creases the lattice misses entirely (thin features).
- **Blends**: once smooth-min ops exist, the tree knows a blend
  has no crease - skip its pair (sampled normals can't know this).
- **Performance**: CGAL exact-predicate fallbacks on cospherical
  lattices; CCDT machinery overhead (~70% on the 64^3 sphere);
  shell-thinning, incremental edge scans, MT - all mapped in
  MESH-NEXT, none attempted yet.
- **Scale**: nothing beyond 96^3 / five showcase models has been
  meshed.  The Zeiss torture corpus awaits.

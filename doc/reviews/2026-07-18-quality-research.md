# Mesh-quality research: the r2 surface-artifact classes

Reviewer pass, 2026-07-18. Scope: the OPEN quality classes behind
Nate's "quite a lot of surface artifacts" on `zeiss_stibnite_r2.3mf`
(all defects < 0.1 mm, visually noticeable). Ranked proposals with
literature, expected impact per class, effort, and referee plans.

Ground truth for this pass (bino r1, THIN+gates, measured today,
reproduces the ledger's `p3-thin-gated` row exactly):

- Repair plateau: worst chip 0.342 -> 0.356 -> **0.242 -> 0.242 sp**
  (flat by round 2). Chip COUNT climbs 2,902 -> 5,203 across the same
  rounds. Depth flat, count rising = self-similar chip regeneration,
  the classic signature of midpoint/bisection refinement with no
  curvature-anchored target.
- Snap pass: 795 edges tented onto traced creases + 1,747 onto the
  surface (SNAP_SURF), final worst **0.170 sp**. But it SKIPS 4,408:
  **993 FAR** (no crease, surface projection failed the acap gate),
  **733 damage-refused**, **2,680 churn-refused**.
- FAR clusters band tightly: 6 top bins all at z=72.6, x in
  [-12.5, +12.5], y in [57, 63] - the eyepiece/barrel tangent-blend
  collar. 104 / 96 / 75 / 63 / 51 / 42 chips per bin, ~430 in the top
  six alone.

Read: the *plateau itself is not the visible defect* - the snap pass
is designed to consume it. The visible artifacts are the chips snap
**cannot** cure: the 993 FAR untraced-blend divots (one banded
feature) plus the ~3,400 tents the damage/churn referees veto. So
Class 1 (plateau) and Class 2 (FAR blends) are the same fight from two
ends, and the highest-leverage move attacks snap's skip list, not the
repair loop.

---

## Where Stibnite already sits in the literature

The pipeline's steps 5-8 are the restricted-Delaunay refinement
lineage: Boissonnat & Oudot's provably-good surface sampling
(restricted Delaunay of an epsilon-sample w.r.t. local feature size is
topologically and geometrically faithful) [BO05], realized for
isosurfaces without a full 3D structure by DelIso [DL08], extended to
sharp/non-smooth input by Cheng-Dey-Ramos-Ray's PSC refinement and
DelPSC [CDRR07, CDR07].

The CDRR move is **protecting balls**: weight the non-smooth curves and
corners so they appear as edges of the (weighted) restricted Delaunay
*before* volume refinement runs. **Stibnite already does the strong
form of this** - it traces the crease curves analytically and inserts
them as hard CGAL constraints (step 4-5). Where CDRR protect an
*undersampled approximation* of the feature with balls, Stibnite pins
the *exact* curve. That is why sharp creases are clean and why the open
defects are exactly the features the tracer CANNOT reach: tangent
blends whose seam is a double root (Class 2), and smooth high-curvature
regions with no crease clause at all (Class 3 cones).

The theory also tells us why the plateau is density-invariant. At a
concave crease the local feature size collapses (the medial axis
touches the surface), so a *uniform* bisection can never reach the
epsilon-sample density BO05 requires - it just mints smaller
self-similar chips, which is precisely the 0.196/0.198/0.192 sp
plateau the code already documents at 1x/2x/4x. Convergence needs a
sizing field that shrinks toward the feature, or the feature pinned
directly. Stibnite pins traced creases; the residue is the un-pinnable
features.

---

## Ranked proposals

### P1 - Hessian curvature-valley seeds (Class 2 FAR + Class 3a cones)

**The unifying win.** The FAR clusters are untraced tangent-blend
seams. The min/max clause tracer can't march them because
{f_A = f_B, f = 0} has a double root at G1 contact (graveyard, correct
on paper). But the seam is still a real *geometric* feature: the
concave valley of a fillet blend is a sharp negative extremum of the
minimum principal curvature, and curvature stays finite and
well-defined straight through G1 contact even though the clause algebra
degenerates. So extract the seam as a **ridge/valley line** instead of
a clause intersection [OBS04].

Mechanics, all on the existing field oracle:
- The shape operator is the tangent-plane restriction of the Hessian
  of f, normalized by |grad f|. Stibnite already takes central
  differences of f everywhere (repair, snap, QEF); the Hessian is one
  more difference stencil (the 7-point stencil already sampled in the
  repair gate gives the diagonal; add mixed terms -> 9-ish extra evals
  per candidate leaf, batched through `eval_points_mt`).
- Principal curvatures = eigenvalues of the 2x2 shape operator.
  Valley = locus where kappa_min is a negative extremum along its own
  curvature direction (the OBS04 extremality condition, third-order but
  computable by differencing kappa along the eigendirection).
- Emit valley points as **tracer seeds**, march them with the existing
  predictor-corrector (the corrector already Newton-projects onto
  f = 0; add a second constraint keeping the point at the curvature
  extremum - a 2-equation corrector Stibnite's SSI machinery is shaped
  for). Output: valley polylines in the SAME `tchains` structure the
  snap pass already consumes.

Then the 993 FAR chips stop being FAR - they snap onto a valley line
instead of failing the acap gate. This mirrors the v25/v26 "law
delivery + local-pitch snap" load-bearing pair: the seam doesn't need
to be a hard *constraint* (that risks the CCDT crossing/fallback
breakage the graveyard warns about at SPREAD_DOT widening), it only
needs to be a snap *target*, which is pure output surgery.

Why not the alternatives: **medial-axis extraction** of the concave
seam is the same feature from the dual side, but it needs a global
medial computation and the graveyard already burned two days on fitted-
geometry forensics ("the FIELD is the only oracle") - the Hessian
route stays on the oracle. **Curvature-line-aligned anisotropic
meshing** [ACSYD03] would cure blends and cones together but is a
ground-up remesher, exactly the ambitious-refactor shape the graveyard
is full of.

- Impact: **highest.** Directly eliminates the dominant visible class
  (the 993-chip banded collar) and, with a curvature bar, doubles as
  the cone-quilting density signal (P3). Expected FAR skip 993 ->
  low-hundreds; worst residual chip below 0.170 sp.
- Effort: **1-2 days.** New seed source + 2-constraint corrector, but
  reuses tracer march and snap consumer end-to-end. Risk: valley
  extraction on near-flat regions is noisy (third-order derivatives) -
  gate emission on a curvature-magnitude floor so it only fires where
  a real blend lives.
- Referee: bino r1, `_CHIP_DEBUG`. Watch `SNAP: ... N far`
  (993 -> target < 200) and `FINAL mesh worst depth` (0.170 -> lower).
  `_DUMP_CHAINS` to eyeball the new valley polylines land on the
  z=72.6 collar. Then lamp r1 for the cone crossover (P3).

### P2 - Snap-coverage recovery (Class 1 residual)

The repair plateau is fine - it is the feedstock for snap. The visible
residue is snap's skip list: 993 FAR (P1 handles these), 733
damage-refused, 2,680 churn-refused. The churn/damage referees were
added for real reasons (unrefereed tents minted the 0.309 sp chip and
180-degree folds), so this is *recalibration*, not removal.

Two levers, both low-risk because they only change *which* residual
chips get a tent, never introduce an unrefereed one:

1. **FAR surface-projection acap.** Today a FAR chip snaps onto the
   Newton-projected surface point only if the projection travels
   < acap = min(local_pitch, max(2.5*depth, 0.35*pitch)). The 993
   skips are projections that overshot that cap - i.e. the chord
   midpoint's gradient pointed at a surface point further than 2.5x
   the divot depth, usually because the midpoint sits in the concave
   throat where central-difference |grad| is under-read (the same
   1/cos(half-angle) inflation the credibility gate already models at
   2x). Widen the FAR-only cap toward the wedge-inflation bound and
   re-project from the *refined* gradient. Pairs with P1 (valley snap
   supersedes most FAR chips; this catches the rest).
2. **Churn-bar per-region.** 2,680 churn refusals (apex < 0.02 sp
   above the split planes) are near-flat sliver bumps - correctly
   refused on flats, but on a *curved* blend a 0.02 sp tent is real
   geometry, not z-fighting. Make the churn floor curvature-aware:
   keep 0.02 sp on low-curvature facets, relax it where the Hessian
   (P1's stencil) shows real bending. Same stencil, no new pass.

- Impact: **medium.** Recovers a chunk of the 3,400 refused tents that
  are genuine cures mis-classified as noise. Cosmetic-to-moderate;
  compounds with P1.
- Effort: **half-day** if P1's Hessian stencil already exists (P2.2
  depends on it); P2.1 alone is hours.
- Referee: bino r1 `_CHIP_DEBUG` skip breakdown (`damage` / `churn` /
  `far` counts) + `tools/fightpix.py` on the flats (the churn relax
  must NOT raise z-fight pixels - that is the guardrail the gate
  protects). Screws tilt.py must hold (churn tents papering pinches).

### P3 - Curvature-adaptive density trigger (Class 3a cone quilting)

Cone quilting (lamp) is the density campaign's "chainless-curvature"
item, and the literature is unambiguous on the mechanism: autodense
triggers on QEF residual, which is **near-zero on a smooth developable
cone** (the clamped solve is an exact plane fit there) - so the trigger
is blind exactly where the cone needs density [DL08-lineage; the code
comment "the QEF residual cannot see it" is correct]. Curvature-
adaptive sizing is the standard cure: derive the target edge length
from max principal-curvature magnitude and an approximation tolerance
[TVCG22, ACSYD03].

Concretely: reuse P1's Hessian to compute max |kappa| per crease-
suspect leaf, and add a density trigger keyed on it (a cone reads high
mean curvature along its circular direction while a bench sphere's is
bounded - calibrate the bar so bench referees stay at 4-5 live pairs
and lamp cones cross it). This is the CHAINLESS complement to the
existing live-pair (chain-crowding) trigger: dense where curvature is
high AND no traced chain passes.

Caveat from Stibnite's own record: density is self-similar for the
plateau, so denser-alone won't zero the quilt - it resolves the
faceting finer but the residual valley chips still need P1's valley
snap. Density + valley-snap is the pair, same as v25/v26. A *proper*
cure (curvature-line-aligned anisotropic sampling, ACSYD03) is the
heavy option; not recommended against the graveyard's refactor
tombstones.

- Impact: **medium**, cosmetic (visible faceting on cones, no
  correctness stake). r2 already halves it per the doc; this is the r1
  and coarse-feature win.
- Effort: **half-day** given P1's stencil (trigger + bar calibration).
- Referee: lamp r1 (~2 min), `tools/roughness.py` dihedral clusters on
  the cone faces (count must drop) + eye. Bench models must not gain
  triangles (bar calibration guard).

### P4 - Hidden-thin graze oracle (Class 3b thin cone caps)

The knob's missing sub-lattice caps are the `_HIDDEN` trigger's target,
but it's parked because it fires on tangent grazes (18 false leaves on
the sharp control, still hidden at 4x). It needs a graze-vs-feature
oracle. This is an interval-arithmetic discriminator, not a literature
gap: a hidden leaf where the interval said "maybe" but all corner
samples agree is either (a) a real thin cap - the field dips to near
|f|=0 somewhere interior but the lattice straddled it - or (b) a graze
- the interval barely kissed zero and |f| stays bounded away.

Discriminator: inside a flagged hidden leaf, probe a finer sub-lattice
(one extra level, batched) and take min |f| / |grad|. A real cap reads
min-distance below a fraction of the cell pitch (surface is *inside*,
just missed); a graze reads min-distance bounded away from zero. Gate
the hidden dilation on that min-distance test.

- Impact: **low-medium**, narrow feature class (sub-lattice caps on
  thin tapered geometry). Matters for the knob and any < 1-cell detail.
- Effort: **half-day** (one sub-lattice probe + bar).
- Referee: knob r1 (~1 min), `_THIN_DEBUG="x,y,z,r"` at a cap witness
  (must flip KEEP), visual material presence at the cap. Sharp control
  must NOT gain the 18 graze leaves (the guard the oracle exists for).

### P-NOT - Sliver exudation (named in the brief; ruled out)

The mission listed sliver exudation. It does not apply. Exudation
[CDEFT00] is a *weighted-Delaunay tetrahedron*-quality method: it
perturbs vertex weights to pump degenerate tets (slivers) out of the
3D triangulation. Stibnite's defects are **surface-facet** defects
(chips/divots on the extracted boundary), not tet slivers, and the
RAZOR/Z-FIGHT arc already established that the tet-complex disjointness
guarantee makes raw extraction fold-free (strict-interior fightpix
reads zero). There is no sliver population to exude, and adding a
weighted-DT perturbation pass would jeopardize the exact-predicate
constraint conformance that the crease law depends on. Filing here so
it is not re-derived.

---

## Class 4 - Are the r2 artifacts new classes or the same at half scale?

**Same classes, half scale - reasoned from the doc + today's bino
data, no new zeiss run needed** (and the full model must never run).

The r2 economics (doc, measured on bino): time x3.4, tris x3.2, and
*every sp-scaled defect halves in absolute size*. The FAR-cluster band
and cone quilting are sp-scaled - so at r2 each divot is ~half as deep
but the finer tessellation puts **~3x more edges** across the same
blend, so the chip COUNT in the band roughly triples while each drops
under the eye's per-divot threshold. The visible result of "many
sub-0.1 mm divots tiling one collar" is exactly a patch that reads
rough without any single defect being large - which is precisely
"quite a lot of surface artifacts, all below 0.1 mm but visually
noticeable." That is the **P1 FAR class at r2**, not a new class.

The one class that does NOT shrink is the coincident fold/scar class
(zero-thickness / z-fighting geometry - resolution-immune, doc-
confirmed: bino r2 was "absolutely fucking beautiful, only the
fold/scar class survives"). Halving everything else *relatively
promotes* it: at r2 the coincident class is the largest absolute
defect left. So the r2 quality program is two-front:

1. **FAR-cluster valley snap (P1)** - the dominant *textured* artifact
   at r2, tripled in count. Highest leverage.
2. **The coincident fold/scar class** - already in flight as the SLICER
   HYGIENE arc (SEAL + WINDING); the bino residual-336-open anatomy is
   the open hunt there. Not a mesher-geometry fix, a topology-hygiene
   one.

No evidence points to a *new* r2-only class. Recommend confirming
cheaply on bino r2 (198 s, known-safe) with `_CHIP_DEBUG` after P1
lands - the FAR skip count at r2 is the number that proves the band is
the same feature at half scale.

---

## Referee plan summary

| Proposal | Dish | Instrument | Number that must move |
|----------|------|-----------|----------------------|
| P1 valley seeds | bino r1; then lamp r1 | `_CHIP_DEBUG`, `_DUMP_CHAINS` | SNAP `far` 993 -> <200; final worst 0.170 sp -> lower; valley chains land on z=72.6 collar |
| P2 snap coverage | bino r1 | `_CHIP_DEBUG`, `fightpix.py`, `tilt.py` | `far`+`churn`+`damage` skips down; fightpix flats not up; screws tilt held |
| P3 curvature density | lamp r1 | `roughness.py`, bench tri-count | cone dihedral clusters down; bench tris unchanged |
| P4 hidden oracle | knob r1 | `_THIN_DEBUG`, sharp control | cap witness flips KEEP; sharp control gains 0 graze leaves |
| P4 confirm | bino r2 (198 s, safe) | `_CHIP_DEBUG` | FAR skip count at r2 == same band, halved depth |

## Impact / effort ranking

1. **P1 Hessian valley seeds** - highest impact (dominant r1 *and* r2
   visible class), 1-2 days, reuses tracer + snap.
2. **P2 snap coverage** - medium impact, half-day, compounds with P1,
   low risk (referee-gated by construction).
3. **P3 curvature density** - medium/cosmetic, half-day on P1's
   stencil, cone-only.
4. **P4 hidden oracle** - low-medium/narrow, half-day, knob-only.

P1's Hessian stencil is the shared dependency under P2.2, P3, and P4 -
build it once, four classes get cheaper. That is the sequencing
recommendation: land the curvature stencil first, then P1 -> P3 -> P2
-> P4.

---

### References

- [BO05] Boissonnat, Oudot. *Provably good sampling and meshing of
  surfaces.* Graphical Models, 2005. (restricted-Delaunay epsilon-sample
  fidelity; Lipschitz sizing.)
  https://link.springer.com/article/10.1007/s00454-008-9109-3
- [DL08] Dey, Levine. *Delaunay meshing of isosurfaces.* The Visual
  Computer, 2008 (DelIso). https://link.springer.com/article/10.1007/s00371-008-0224-1
- [CDRR07] Cheng, Dey, Ramos, Ray. *Delaunay refinement for piecewise
  smooth complexes.* SODA/SoCG 2007 (protecting balls).
  https://cse.hkust.edu.hk/~scheng/pub/soda2007a-psc.pdf
- [CDR07] Cheng, Dey, Ramos. *DelPSC: a Delaunay mesher for piecewise
  smooth complexes.* / *A practical Delaunay meshing algorithm for a
  large class of domains.* IMR 2007.
  https://link.springer.com/chapter/10.1007/978-3-540-75103-8_27
- [OBS04] Ohtake, Belyaev, Seidel. *Ridge-valley lines on meshes via
  implicit surface fitting.* ACM ToG (SIGGRAPH) 2004.
  https://dl.acm.org/doi/10.1145/1015706.1015768
- [ENG16] Engwirda, Ivers. *Off-centre Steiner points for
  Delaunay-refinement on curved surfaces.* CAD 2016 / *Size-optimal
  Steiner points* arXiv:1501.04002 (curvature-adaptive off-center
  Steiner placement - the anti-self-similarity insert rule, relevant
  if a refinement-time cure is ever preferred over snap).
  https://arxiv.org/pdf/1501.04002
- [CDEFT00] Cheng, Dey, Edelsbrunner, Facello, Teng. *Sliver
  exudation.* J. ACM 2000. (Ruled out - see P-NOT.)
- [ACSYD03] Alliez, Cohen-Steiner, Yvinec, Desbrun et al.
  *Anisotropic / curvature-line meshing.* (Heavy alternative for
  blends+cones; not recommended vs the graveyard.)
- [TVCG22] *Adaptively isotropic remeshing based on curvature smoothed
  field.* IEEE TVCG 2022 (curvature sizing field).
  https://ieeexplore.ieee.org/document/9978684/
- [Recent] *Better sampling bounds for restricted Delaunay
  triangulations.* arXiv:2603.19826 (epsilon <= 0.3245 suffices).

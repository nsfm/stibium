# From instruments to certificates: a roadmap to engineering-grade guarantees

*Opus research agent, 2026-07-18. Commissioned to chart the path from
Stibnite's current measured instruments to formal/certifiable output
claims — the product differentiator in the DRAFT briefing
(doc/DRAFT-stibnite-briefing.md, Claim 4 and Product thesis §2:
"this mesh deviates ≤ X mm from the model; watertight, proven; zero
repairs required"). Read-only w.r.t. the repo. Sources cited inline.*

---

## 0. The one idea that organizes everything: two-sided vs one-sided

"Deviation ≤ X mm" sounds like a single number. It is not. The
Hausdorff distance between the output mesh **M** and the true surface
**S = {f = 0}** is the *max of two one-sided distances*:

- **d(M → S)** — every point of the mesh is within δ of the surface.
  "No triangle sits far from the model." Catches *bulges, overshoot,
  spurious skin*.
- **d(S → M)** — every point of the surface is within δ of the mesh.
  "No part of the model is missing from the mesh." Catches *dropped
  thin walls, missed holes, unmeshed pockets*.

These have wildly different costs. **d(M → S) is facet-local and
cheap** — you interrogate the field on each triangle you already have.
**d(S → M) requires reasoning about surface you did *not* mesh** — it
is a coverage/topology claim and is fundamentally harder. Every rung
below is really a statement about which of these two it proves, and
over what set it proves it. The single most important honesty
discipline for marketing: **never say "Hausdorff ≤ X" when you have
only proved the one-sided M → S bound.** Say what you proved.

A parallel decomposition governs the other two guarantee families:

- **Watertight / manifold** — combinatorial, already counted *exactly*
  (0 open edges under bit-exact vertex identity; MESH-WAR). Solved.
- **Intersection-free** — geometric, true by construction of the raw
  extraction, *breakable by post-passes*. Verifiable cheaply (§3).
- **Isotopy to S** — the topological d(S → M) claim in full strength.
  Research-grade for sharp CSG (§Rung 5).

---

## 1. What is proven today, and what is only sampled

| property | current status | proof strength |
|---|---|---|
| watertight (0 open edges) | counted on bit-exact welded ids | **exact / proven** (combinatorial) |
| worst surface deviation | `field_probe`/repair-loop metric: \|f\|/\|∇f\| at edge **midpoints**, 7-tap finite-difference gradient (`delaunay.cpp:5288`, `6798`) | **sampled** — a point-estimate, not an enclosure |
| intersection-free | strict-interior `fightpix` reads 0 on raw meshes | **instrumented, not proven**; post-passes can break it |
| topology / no missing features | eye referee + repair loop | **not claimed** |

The deviation metric is better than the mission brief assumed: it
already divides by the gradient (`dist = |f| / |grad|`,
`delaunay.cpp:5288`), so it reports a *length in model units*, not raw
\|f\|. But it has two weaknesses that keep it a measurement, not a
certificate:

1. **It samples, it does not enclose.** The metric probes the field at
   edge midpoints only. The worst point on a facet is generically *not*
   a midpoint. The reported "worst chip" is a lower bound on the true
   worst deviation, dressed as if it were the worst.
2. **The gradient is a finite-difference point estimate, not a proven
   lower bound.** `|grad|` comes from a 7-tap central difference at one
   point (`h=0.001`). A rigorous distance certificate needs a proven
   *lower* bound on \|∇f\| *over the whole facet* (dividing by too-large
   a gradient understates the true gap). A point sample of the gradient
   cannot deliver that.

**Two foundational gaps the code read exposed (both gate Rung 1 —
correct the briefing on the first).** The briefing and
doc/research/2026-07-14-frep-sota.md call the kernel's intervals
"sound." The meshing interval path does **not** currently look sound:

- **(P0a) Intervals appear unsound on the CPU meshing path.** The
  interval op library `lib/fab/src/tree/math/math_i.c` computes
  endpoints in round-to-nearest with **no directed rounding and no ULP
  widening** (`add_i` is literally `{lo+lo, hi+hi}`; no
  `nextafter`/`fesetround`/`fenv`). `Interval` (`inc/fab/util/interval.h`)
  is a bare `{float lower, upper}` with no rounding metadata. A
  certificate built on these bounds is *almost* proven but can be
  violated by accumulated round-off — not good enough to put "proven"
  in front of a customer. **Action item: locate where the "sound
  intervals / recent commit" claim actually lives** (a GPU/WGSL path? a
  wrapper that widens?) and, if the CPU `math_i.c` path is the one the
  mesher uses, add outward-rounding / conservative ULP widening before
  any deviation number is called "certified." This is the prerequisite
  for Rungs 1, 3, and 4.
- **(P0b) No interval-gradient evaluator exists.** The AD `derivative`
  type (`inc/fab/tree/math/math_g.h`: `{float v, dx, dy, dz}`, driven
  by `tape_eval_g`) is **pointwise forward-mode AD**, not interval-AD.
  There is no combined type carrying an interval value *and* interval
  partials. Rung 1 needs exactly that (an `interval_derivative` type +
  a `math_ig` op set + `tape_eval_ig`) to get a proven \|∇f\| lower
  bound over a box. This is real new code, not a reuse of an existing
  evaluator — the cost note in Rung 1 accounts for it.

There is also a procedural gap the correctness review flags
(doc/reviews/2026-07-17-correctness.md #1): the §3 post-passes and the
**export-time meshoptimizer QEM** (`app/export/export_mesh.cpp`,
`STIBIUM_DMESH_SIMPLIFY`) run *after* the metrics are computed, and
`DMesh` (`delaunay.h:162`) doesn't even persist a worst-deviation
float — the number lives only in a debug print. Today's reported
deviation/open/nm counts describe the *pre-decimation* mesh, not the
shipped STL. Any certificate must be a struct field on `DMesh`,
computed on the *final* triangle set, last thing before write-out.

---

## 2. Rung 1 — Per-facet certified one-sided deviation bound  ★ do this first

**The claim it earns:** *"Every triangle of this mesh is certified, by
interval arithmetic on the exact model, to lie within δ mm of the
mathematical surface. δ is reported per export."* This is Claim 4 of
the briefing turned from *sampled* into *proven*, at trivial cost. It
is the highest-ROI item on the board.

### Math

For a facet T with vertices in a box B ⊇ T, using a *sound* interval
evaluator (P0a) and an interval-gradient evaluator (P0b):

1. **Field enclosure.** Interval-evaluate f over B (via `tape_eval_i`,
   `tape.cpp:771`, once it is outward-rounded) → F = [f_lo, f_hi]. Let
   v(T) = max(\|f_lo\|, \|f_hi\|) — a proven upper bound on \|f\|
   anywhere on T.
2. **Gradient lower bound.** Interval-evaluate ∇f over B (requires the
   new interval-AD path, P0b — the existing `tape_eval_g` is pointwise).
   Each component ∂f/∂x_i encloses to [a_i, b_i]; its distance from
   zero is c_i = max(0, a_i, −b_i). Then
   **m(T) = √(Σ c_i²)** is a proven lower bound on \|∇f\| over B.
3. **Distance certificate.** Standard first-order lemma: if
   \|f(x)\| = v and \|∇f\| ≥ m > 0 on the closed ball of radius v/m
   around x, then S has a point within distance v/m of x (gradient
   flow reaches the zero set before leaving the ball). Applied over
   the whole facet: **δ(T) = v(T) / m(T)** is a proven bound on
   d(T → S), *provided* B is dilated by δ(T) so the ball stays inside
   the evaluated box. One fixed-point pass settles this: compute a
   provisional δ on a box dilated by v/m_prelim, re-evaluate, accept
   when self-consistent.

Report **δ* = max_T δ(T)** as the certified per-export bound. This is
a rigorous d(M → S) ≤ δ* statement — the one-sided direction, over the
mesh you actually emit.

### Implementation shape

Prerequisites P0a (sound outward-rounded intervals) and P0b
(interval-gradient evaluator) first. Then a post-pass after the *final*
simplify (including the export-time QEM), before STL/3MF write:
- add a `float worst_deviation` (+ unresolved-facet count) field to
  `DMesh` (`delaunay.h:162`) so the number survives to the STATS
  dialog the integration sprint plans (MESH-WAR);
- iterate the final facets; build the dilated AABB (or better, an
  oriented box or the barycentric-Bernstein enclosure, §6, for slivers);
- one interval-f eval + one interval-∇f eval per facet, fixed-point on
  the dilation;
- reduce to the max; emit δ* alongside the existing deviation line;
- flag any facet with m(T) = 0 as **"unresolved"** and either
  (a) subdivide its box a few times and retry, or (b) report the count
  of unresolved facets honestly (see failure modes). Never silently
  drop them.

### Compute cost

The per-facet pass itself is cheap; the build cost is the two
prerequisites. From the TIME=2 profile (MESH-WAR, bino r2): field eval
runs ~0.7 µs/point-eval single-threaded. A sound interval eval is ~2×,
an interval-gradient eval ~3–4× a point eval; call it ~10 point-eval-
equivalents per facet for f + ∇f + the fixed-point retry. At 540K
facets (zeiss) that is a few seconds single-threaded, and the pass is
embarrassingly parallel (per-facet independent; the kernel already has
per-thread `TapeCtx` — P5) → **sub-second** wall. The real work is
**P0a** (outward-round `math_i.c` — small, localized, but must be done
carefully and re-validated against the 627K-assertion suite) and
**P0b** (a new interval-AD evaluator — a `math_ig` op library mirroring
`math_g.c`, ~a few days). Medium build, trivial run. Still the highest
ROI on the board because it directly earns the headline product claim.

### Failure modes (name them in the product, do not hide them)

- **Gradient vanishing → uncertifiable facets.** Where ∇f can be zero
  inside B, m = 0 and δ = ∞. This is *not* the sharp crease case: a
  true CSG `min`/`max` of unit-gradient SDFs keeps \|∇f\| ≈ 1 on each
  side of a crease (the gradient *kinks* but does not vanish), so
  hard creases certify fine. Genuine vanishing happens at: smooth-min
  **fillets/blends** (a whole low-gradient region), **medial-axis-like
  interior points** that a fat facet's box can reach, and
  **poorly-normalized primitives**. Remedy ladder: dilate less / use a
  tighter enclosure (§6) → bisect the box → report unresolved. A small
  unresolved-facet count with its own worst-case fallback number is an
  honest certificate; a hidden one is not.
- **Sliver-facet box slack.** A near-degenerate triangle's axis-
  aligned box is far larger than the triangle → loose v(T) → loose δ.
  Oriented boxes or Bernstein-over-the-triangle (§6) fix this; it is
  the main tightness lever.
- **This proves M → S only.** It says nothing about missing surface.
  A mesh that dropped an entire thin wall can still pass Rung 1 with a
  tiny δ*, because every triangle it *did* emit is close to S. That is
  why Rung 3 exists and why the marketing line must not overreach.

**Honest-marketing wording at this rung:** *"Certified maximum surface
deviation: δ* mm — every facet proven within δ* of the model by
interval arithmetic (N facets unresolved, bounded by δ_fallback)."*
NOT "Hausdorff ≤ δ*".

---

## 3. Rung 2 — Intersection-free certificate  ★ do this second

**The claim it earns:** *"Zero self-intersections, verified by
exhaustive triangle-triangle sweep"* (verification posture) and,
optionally, *"intersection-free by construction"* (stronger, harder).
Turns the briefing's "instrumented, not proven" into a real number.

### Which post-passes preserve it, which break it

The raw extraction (`delaunay_mesh_soup`, `delaunay.cpp:5087`) is a
subcomplex of an interior-disjoint tet complex → **intersection-free by
construction** — asserted in code (`delaunay.cpp:5090` "facets of a tet
complex cannot self-intersect", `delaunay.h:14`) and corroborated by
MESH-WAR (strict-interior `fightpix` reads 0 on every un-decimated
mesh; the earlier "extraction mints folds" story was an instrument
artifact). *Note the guarantee is a property of the CGAL triangulation,
not independently checked in code.* It is then *spent* by mutating
post-passes (order in `PROG_STAGES[]`, `delaunay.cpp:282`):

| post-pass | file:fn | vertex action | can self-intersect? |
|---|---|---|---|
| `repair_winding` | `delaunay.cpp:7680` | reorders winding only | **no** — safe |
| `seal_seams` (ε-weld + twin annihilation) | `delaunay.cpp:7540` | merges coincident, deletes degenerate/twins | **no** at ε = 3e-4 sp (merging coincident points / deleting zero-area tris cannot create a crossing; only a *large* ε could pull a vertex through a facet) |
| pinch-split (manifold split) | `delaunay.cpp:5800` | duplicates coincident vertex | **no** — separation, not overlap |
| snap tents (crease-snap) | `delaunay.cpp:6179` | **adds** apex vertex, retriangulates | **yes** — a tent apex can pierce a neighbor (MESH-WAR: suspected in bino's residual fights) |
| `weld_slivers` / `flip_slivers` | `delaunay.cpp:7227` / `7392` | merges/moves/flips sliver fans | **yes** |
| `decimate_flats` | `delaunay.cpp:7063` | **collapses** flat-fan vertices (≤30 passes) | **yes** — MESH-WAR caught it minting 522 fights via the double-remap bug; the in-plane orientation guard cannot see a fold through a *non-adjacent* triangle |
| export QEM (`simplifyMesh`/meshoptimizer) | `app/export/export_mesh.cpp` | true quadric collapse, *outside the fab lib* | **yes** — same family; tiles flats with skinny tris (MESH-WAR post-QEM note); **no guard at all** |

So: **winding, seal, pinch-split preserve it; snap tents, the sliver
passes, decimate_flats, and export QEM can break it.** These are exactly
the passes the correctness review flags as running without a
link-condition check (#1), and the export QEM runs in the app layer
where no fab-side guard reaches it.

### Two postures, cheapest first

**(a) Verify — global triangle-triangle sweep, per export.** Broadphase
+ narrow-phase. **No tri-tri test exists anywhere in the repo today**
(grep confirms only prose asserting disjointness), so this is new — but
the broad-phase building block is already in-tree: `MeshQuery`
(`lib/fab/inc/fab/mesh/mesh_query.h`) is an **AABB tree** currently used
only to turn *imported* meshes into SDFs (`mesh_import.cpp:294`); reuse
its tree for the self-intersection query instead of importing CGAL PMP.
- *Broadphase:* the `MeshQuery` AABB tree (or a uniform spatial hash at
  ~mean-edge-length). CGAL's `Polygon_mesh_processing::self_intersections`
  is the reference algorithm to mirror. Expected O(n log n).
- *Narrow-phase:* Möller or Devillers triangle-triangle, **skipping
  shared-vertex/shared-edge adjacent pairs** and using exact/adaptive
  predicates for the coplanar case — this is the exact false-positive
  that produced the `fightpix` instrument artifact, so the narrow
  phase must be built on robust predicates (Shewchuk; the
  exact-predicates-for-mesh-CSG reference, arXiv 2405.12949).
- *Cost:* ~a few seconds at 540K tris, parallelizable — comparable to
  or cheaper than Rung 1. This is the belt-and-suspenders certificate:
  "0 self-intersections, verified."

**(b) Guarantee — restrict the mutating passes to provably-safe moves.**
Each candidate collapse or tent is accepted only if a *local* triangle-
triangle test against the affected 1-/2-ring passes, on top of the
link-condition (topological safety). If raw is intersection-free and
every accepted move preserves the invariant locally, the output is
intersection-free *by construction* — no global sweep needed. This is
the natural extension of MESH-WAR's existing per-tent DAMAGE/CHURN
referees, which today probe *field depth* but **not triangle
crossings**; adding the local geometric guard closes the actual hole.

**Recommendation:** do both — local guards make the passes safe
(construction), the global sweep is the cheap per-export receipt
(verification). Ship the verification number; earn the construction
claim over time.

**Honest-marketing wording:** *"Self-intersections: 0 (verified by
exhaustive triangle-triangle sweep)."* Only after the local guards are
proven-complete: *"intersection-free by construction."*

---

## 4. Rung 3 — Topology audit: "no missing/spurious feature above scale h"  ★ cheap, high-trust

**The claim it earns:** the *honest, affordable* half of d(S → M):
*"No connected component of the model larger than scale h was omitted,
and no output component lacks a corresponding field zero-crossing —
verified by interval sweep of the octree."* This is deliberately
*weaker than isotopy* and much cheaper, and it catches the failures
customers actually fear (dropped holes, missing thin walls, spurious
blobs) without claiming a topology theorem.

### Math / implementation

The octree descent already interval-evaluates f per cell for culling.
Reuse it as a coverage audit at leaf scale h:
- **Missing surface:** any leaf whose interval F straddles 0 (surface
  provably passes through) but which contains **no mesh facet** → a
  candidate omission. Flag and count.
- **Spurious surface:** any mesh component in a region where interval F
  is sign-definite (no zero can exist) → provably spurious.
- The residual after both checks is a proven statement: *at resolution
  h, every field zero-crossing cell owns a facet and every facet owns
  a zero-crossing.* Genus/handle agreement (full Euler-characteristic
  or Morse critical-point matching) is a further, harder step — defer.

### Cost

Near-free: the interval-f-per-leaf is already computed during descent;
the audit is a second bucketed pass associating facets to leaves
(spatial hash). Sub-second.

### Why this is the right product move (not full isotopy)

Full certified isotopy is Plantinga–Vegter (below) and it **does not
terminate on sharp CSG**. Rung 3 gives 80% of the customer-facing trust
("nothing went missing at print scale") for 1% of the cost, and it is
*honest*: it states the scale h at which it is verified, and it does
not pretend to have proved topological equivalence.

**Honest-marketing wording:** *"Feature-completeness verified to scale
h mm: no model feature larger than h omitted, no spurious component
(interval-verified against the exact model)."* NOT "topologically
correct."

---

## 5. Rung 4 — Two-sided Hausdorff (upgrade Rungs 1+3 to a real δ)

**The claim it earns:** the genuine article — *"Hausdorff distance to
the model ≤ δ mm, both directions, proven."*

### Math

Add the missing S → M direction. For every octree leaf that S crosses,
prove that a mesh facet passes within δ of the surface patch inside
that cell. The surface patch is enclosed by interval-Newton /
**interval Krawczyk** (arXiv 2602.07718, Feb 2026 — the most on-point
recent result: certifies existence/uniqueness of a zero in a box and
bounds the Hausdorff distance from an approximation to the enclosed
surface). Combined with Rung 1's M → S bound, the max is a certified
two-sided Hausdorff bound over the smooth part of the model.

### Cost & risk

Medium. Interval-Krawczyk per surface-crossing leaf is more expensive
than plain interval f, but the leaf set is the surface band, not the
volume. Real but not research-hard **on smooth regions**. The catch is
the same crease problem as Rung 5: Krawczyk needs ∇f ≠ 0 and a regular
zero, so the *crease curves themselves* must be certified separately
(next rung). Practical form: two-sided δ on the smooth patches +
one-sided δ + a certified crease-curve deviation, reported as a
composite.

**Honest-marketing wording:** *"Hausdorff ≤ δ mm on smooth regions;
crease curves within δ_c; both interval-proven."*

---

## 6. Cross-cutting enabler — tighter enclosures (affine / Bernstein)

Rungs 1 and 4 are only as good as the enclosure tightness. Axis-
aligned interval arithmetic over a fat box is loose, which inflates δ*
and manufactures unresolved facets. Two levers, both flagged as
unexploited in doc/research/2026-07-14-frep-sota.md §2/§4:

- **Bernstein-form enclosure over the triangle's barycentric domain.**
  The convex-hull property of the Bernstein basis bounds f over the
  *simplex itself* (not its AABB), giving much tighter v(T) for slivers
  — "similar performance to affine arithmetic" per the enclosure
  literature. Directly attacks the sliver-box-slack failure mode.
- **Affine arithmetic** (revised/reduced): tracks correlations between
  subexpressions, tighter than IA on deep correlated tapes (common in
  CSG), at higher per-op cost. Tightens both f and ∇f enclosures →
  smaller certified δ, fewer unresolved facets, shallower subdivision
  near blends. Nobody in the fidget lineage runs AA in a JIT — it is
  also a rendering-culling win, so it compounds.

These are not their own rung; they are the multiplier that makes Rungs
1 and 4 report *small* numbers instead of merely *true* ones.

---

## 7. Rung 5 — Certified isotopy (the research horizon)

**The claim it would earn:** *"Provably topologically identical to the
mathematical model"* — the full d(S → M) topological guarantee.

**Why it is hard here specifically.** Classic Plantinga–Vegter (octree
+ interval f + interval ∇f small-normal-variation test) proves isotopy
and a Hausdorff bound — **but assumes 0 is a regular value, i.e. ∇f ≠ 0
everywhere** (confirmed across the PV literature and the complexity
analyses, Springer DCG 2022; arXiv 2004.06879). **CSG models are all
creases**; the small-normal-variation test never passes on a cell
containing a crease, so PV *subdivides forever* — non-terminating on
exactly this workload. The interval-Krawczyk surface paper and the
singular-algebraic-surface work (arXiv 0903.3524) both explicitly fail
at ∇f = 0 and creases.

**The route that could work — and why Stibnite is uniquely placed.**
Certify *piecewise*: the smooth patches by PV/Krawczyk, and the crease
curves *separately* by certifying the traced feature polylines. This is
where **Stibnite's tape-derived crease law becomes a certification
asset, not just a meshing asset** (briefing Claim 2): the mesher
already identifies each crease {f_A = f_B, f = 0} from surviving
min/max clauses and marches it with a field-verified predictor-
corrector. Certifying that traced curve (interval enclosure of the
1-D system along each segment) turns the sharp-feature machinery into
the missing crease-certification ingredient that generic PV lacks. No
one in the fast-f-rep world has fused certified topology with a
production JIT tape (frep-sota §3); doing it *crease-aware* is the
defensible, genuinely-open niche — but it is quarters, not weeks.

**Honest-marketing wording:** do **not** claim isotopy until this
lands. Until then Rung 3 ("no missing feature above scale h") is the
honest topology statement.

---

## 8. Prior art on shipped certification (what to claim novel, precisely)

| system | what it guarantees | proven per-output deviation bound? |
|---|---|---|
| **CGAL Mesh_3 / Surface_mesher** | homeomorphism + Hausdorff/Fréchet bound — but **only** when the radius criterion < ε·local-feature-size with ε < 0.16, for **smooth** surfaces; *termination* guaranteed, *bound* holds conditionally | **No** — guarantee-by-theorem-under-hypotheses, not a number measured/reported per export |
| **nTop** | watertight + manifold + self-intersection-free **by construction** ("error-proof meshes," "no repair"); increasingly *avoids* meshing by slicing the implicit directly (Magics reads nTop implicit natively) | **No** — construction guarantee, no published proven deviation number vs the model |
| **TetWeave** (SIGGRAPH 2025) | watertight + 2-manifold + intersection-free by construction (Marching Tets on Delaunay) | **No** — differentiable-opt method on sampled/directional SDF; no exact-model certificate |
| **Metrology / inspection** (ZEISS INSPECT/GOM, PolyWorks, Geomagic Control X, Hexagon PC-DMIS) | per-part deviation **reports**: colorized deviation maps, tolerance pass/fail, PDF | **Partial & different** — compares a *scan* of a physical part to CAD, and the deviation is *sampled* point-to-surface, not a *proven enclosure*; done at inspection time, not generation time |
| **Slicers / repair** (Prusa, Netfabb, Magics) | auto-repair ~80–90% of STL defects heuristically | The ecosystem *assumes meshes are broken* — "zero repairs required" is a real differentiator against it |

**The defensible novelty (calibrated).** No shipping mesher reports a
*proven, per-export, model-exact* deviation bound (a general web survey
confirms certified deviation bounds "are not commonly highlighted in
standard mesh-generation documentation"). The competitive shape:

- nTop says *"trust our construction / don't mesh at all."*
- Metrology software says *"we'll measure the physical part's deviation
  after the fact, sampled, for a five-figure license."*
- **Stibnite can say *"here is a proven deviation number for this mesh
  against the exact model, computed at generation time, for free."***
  That is the metrology *deliverable* (a deviation certificate
  customers already know how to read) with a *provenance* nobody else
  has (interval-proven, model-exact, at export). That is the real
  Claim-4 wedge, and Rung 1 alone earns the honest version of it.

---

## 9. Ranked roadmap (ROI = value × 1/cost × honesty)

| # | rung | proves | depends on | cost | horizon | claim earned |
|---|---|---|---|---|---|---|
| **P0** | sound intervals (a) + interval-AD (b) | *prerequisite for 1,3,4* | — | small+medium build | **days** | (nothing customer-facing; makes "proven" honest) |
| **1** | Per-facet deviation certificate | d(M→S) ≤ δ* over the mesh | P0a+P0b | sub-second run | **days** (after P0) | "every facet proven within δ mm of the model" |
| **2** | Intersection-free certificate | 0 self-intersections | **none** | ~seconds sweep + local guards | **days–1 wk** | "0 self-intersections, verified" → later "by construction" |
| **3** | Topology audit (scale h) | no missing/spurious feature > h | P0a | near-free (reuse octree) | **days** (after P0) | "no feature above h mm omitted" |
| 4 | Two-sided Hausdorff | d(M↔S) ≤ δ on smooth part | P0 + Krawczyk | medium (Krawczyk on band) | **weeks** | "Hausdorff ≤ δ, smooth regions" |
| — | tighter enclosures (Bernstein/affine) | *tightens 1 & 4* | P0 | medium | weeks | (smaller δ, fewer unresolved) |
| 5 | Certified crease-aware isotopy | full d(S→M) topology | P0 + research | high, research | **quarters** | "provably identical topology" |

**Sequencing logic.** **Rung 2 is prerequisite-free** — it is pure
geometry on the output triangles (robust tri-tri predicates), depends on
neither interval fix, and reuses the in-tree `MeshQuery` AABB tree. It
is the cheapest *fully-proven* rung; ship it first while P0 is built.
**P0 unlocks Rungs 1 and 3**, which are then days each: they reuse the
octree descent and tet disjointness the pipeline already has, but they
sit on the soundness fix, so the honest "proven" word waits for P0.
Ship 1+2+3 as a bundle: a per-export **certificate block** (δ*,
0 self-intersections, feature-complete to h, watertight exact) as a
`DMesh` field rendered into the STATS dialog the integration sprint
already plans (MESH-WAR). That bundle *is* the "engineering-grade
artifact." Rung 4 upgrades δ* from one-sided to two-sided when a
customer needs the true Hausdorff word. Rung 5 is the moat, not the
MVP — pursue it because the tape-derived crease law makes Stibnite
uniquely able to, but do not gate the product on it.

**Two hygiene fixes that gate any certificate (both already flagged in
doc/reviews/2026-07-17-correctness.md):** compute all certificates on
the **final** post-simplify triangle set (today's metrics precede
decimation, so they describe a mesh you don't ship), and add the
link-condition + local-intersection guards to decimation/snap so the
mutating passes stop silently spending the guarantees the raw
extraction hands them.

---

## Key sources

- Plantinga & Vegter, "Isotopic meshing of implicit surfaces," *Visual
  Computer* 2007 — graphics.stanford.edu/courses/cs164-10-spring/Handouts/paper_vegter-isotopic.pdf;
  Plantinga PhD thesis, cs.nyu.edu/~exact/pap/mesh/vegter/plantingaThesis.pdf
- On PV complexity & the ∇f≠0 / regular-value assumption — *Discrete &
  Comput. Geom.* 2022, link.springer.com/article/10.1007/s00454-022-00403-x;
  arXiv 2004.06879
- Burr & Byrd, "Certified simultaneous isotopic approximation of curves
  via subdivision," arXiv 2302.04908 / 2407.16911
- **"Certified surface approximations using the interval Krawczyk
  test," arXiv 2602.07718 (Feb 2026)** — isotopy + Hausdorff via
  Krawczyk; fails at ∇f = 0 / creases
- Ambient isotopic meshing of implicit algebraic surfaces with
  singularities, arXiv 0903.3524
- CGAL 3D Surface Mesh Generation (ε < 0.16·lfs guarantee) —
  doc.cgal.org/latest/Surface_mesher/index.html; CGAL Mesh_3 —
  doc.cgal.org/latest/Mesh_3/index.html
- CGAL AABB_tree & Polygon_mesh_processing self-intersections —
  doc.cgal.org/latest/AABB_tree/index.html
- Exact predicates & combinatorics for mesh CSG (robust
  triangle-triangle), arXiv 2405.12949
- Binninger et al., TetWeave, TOG/SIGGRAPH 2025 — arXiv 2505.04590;
  alexandrebinninger.com/TetWeave (watertight/manifold/
  intersection-free by construction)
- nTop, mesh-from-implicit "error-proof" / implicit-native slicing —
  ntop.com/resources/blog/meshing-in-fea-cfd-manufacturing/;
  support.ntop.com (Methods for creating surface meshes)
- Metrology deviation reporting — ZEISS INSPECT/GOM, PolyWorks
  Inspector, Geomagic Control X, Hexagon PC-DMIS (per-part sampled
  deviation vs CAD, colorized maps + tolerance PDFs)
- Bernstein convex-hull enclosure over simplices —
  interval.louisiana.edu/reliable-computing-journal (Bernstein basis);
  onlinelibrary.wiley.com/doi/10.1155/2022/9156188
- Affine arithmetic for implicit surfaces (tighter than IA) — Knoll et
  al., *CGF* 2009, doi 10.1111/j.1467-8659.2008.01189.x; revised affine
  arithmetic ray-tracing
- Project internals: doc/DRAFT-stibnite-briefing.md,
  doc/MESH-WAR.md, doc/research/2026-07-14-frep-sota.md,
  doc/reviews/2026-07-17-correctness.md
- Code (verified by read, 2026-07-18): the dmesh pipeline
  lib/fab/src/tree/triangulate/delaunay.cpp (extraction 5087,
  deviation metric 5288/6798, post-passes 5800/6179/7063/7227/7392/
  7540/7680); interval library lib/fab/src/tree/math/math_i.c +
  inc/fab/util/interval.h (naive endpoints — P0a); AD path
  inc/fab/tree/math/math_g.h + tape_eval_g (pointwise — P0b);
  tape interval eval lib/fab/src/tree/tape.cpp:771; AABB tree
  lib/fab/inc/fab/mesh/mesh_query.h; export QEM
  app/export/export_mesh.cpp

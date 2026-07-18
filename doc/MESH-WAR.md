# MESH-WAR: the battle doc

Post-compact entry point.  Read THIS, not the full MESH-NEXT
history (that file is the archive/lab notebook; this is the war
room).  Keep it current: prune won battles, add new fronts.

MISSION: features >= ~0.1 mm retained BY DEFAULT (traditional
3D-print-ready); stretch = finer, for resin.  Watertight is law
(0 GEOMETRIC open edges - counted on welded ids).  Nate's eyes
are the FORMAL referee and sit in the loop BEFORE any default
changes.  One full-model export at a time, ALWAYS (the laptop
incident).  Damage is measured, not predicted.  Count is a lying
metric; use depth + renders.  When the mesh looks wrong, first
ask what the model actually IS.

## ARMORY

Referee models (fast -> slow):
- examples/torture/zeiss_id02_plinth.sb   (~40 s; junction/air-chord dish)
- examples/torture/zeiss_id02_screws.sb   (~30 s; tight-feature dish; ALL
  surfaces axis-aligned -> any tilted tri is error BY DEFINITION)
- examples/torture/zeiss_id02_knob.sb     (~1 min; thin tapered cones)
- examples/torture/zeiss_id02_lamp.sb     (~2 min; cone quilting)
- examples/torture/zeiss_id02_bino.sb     (~3-5 min; the classic referee)
- examples/torture/zeiss_id02_merged.sb   (m20, r1: ~25 min, ~16 GB peak.
  ONE AT A TIME.  Wrap in scratchpad peak.py for wall+RSS.)
- examples/mesh_bench/*.sb (m0, seconds; sharp-crease truth)
- build/zeiss_old_algo_ref.stl (old mesher @7vox, simplified; noisy ref)

Harness: ninja -C build; build/lib/fab/SbFabTest (suite 627,666)
+ hidden tags [.dmesh] [.dtrace] [.dchain] after EVERY change.
Export: STIBIUM_EXPORT_DMESH=1 app/stibium --export X.stl
--resolution 1 <model>.  STL lineage + renders in
build/zeiss_dmesh/ - KEEP ALL, Nate reviews rev-by-rev; ledger
every handed-over rev in doc/MESH-REVISIONS.md.

Instruments (env): _TIME (phase profile), _CHIP_DEBUG (depth
referee + repair counters), _CENSUS[=path] (leaf populations),
_FPROBE="x,y,z;..." (field oracle at points - ends is-it-real
arguments), _PROBE="x,y,z" (constraint forensics), _SEG_DEBUG
(referee rejects), _NM_DEBUG (pinch anatomy), _TRACE_DEBUG,
_THIN_DEBUG="x,y,z,r" (witness keep/drop classes near a point),
REFINE line under _TIME (rounds/inserts/leftover i/o edges),
_DUMP_CHAINS=x.stl (chains as tube STL), _STAGES=prefix
(formation film).  Analysis scripts (scratchpad; rebuild from
here if lost - all ~50-100 lines numpy): stlview.py (z-buffer
region render), roughness.py (dihedral clusters), openedges_
exact.py (independent watertight check), teeth/tilt/densdiff
(anatomy censuses - tilt.py is the screws error-by-definition
benchmark), fightpix.py (z-fighting pixel census).  All in
tools/ in-repo.  Reviews: doc/reviews/2026-07-17-*.md
(correctness / performance / architecture / quadric).

Knob roster: ~26 STIBIUM_DMESH_* - defaults are all measured-in;
see the architecture review for the product/frozen/dead tiers.

## CURRENT STATE (2026-07-18 early am)

Best zeiss: autod31 (727K tris / 0 open / 0.465* / 42 min pre-
femto-guard).  With THIN + gates the zeiss projection is SINGLE-
DIGIT minutes - the validation run is the next big dish.
Bino: THIN+gates 59 s (Nate: identical to base-gated), r2 in
198 s / 1.77 GB (impossible pre-THIN).  Screws v20-class needs
MAX=3 (queue: referee level-3 into the default).  THIN DEFAULT
ON since 2026-07-18 (Nate: "I couldn't spot any visible
difference... Worth it"; the 0.097-vs-0.170 depth delta did
not register to the eye referee).  QEM dial verdict:
STIBIUM_DMESH_SIMPLIFY=0.01 acceptable, 0.02 too far, 0.05 =
octagon bores.  Analysis scripts live IN-REPO at tools/
(Nate's call, 2026-07-18): fightpix.py (z-fighting pixel
census - THE separator for coincident-geometry classes),
razors/overlaps/sheets/tilt/openedges_exact/peak.py.
R2 economics (bino, measured): time x3.4, tris x3.2, RSS
1.77 GB; every sp-scaled defect HALVES in absolute size -
Nate on r2 bino: "absolutely fucking beautiful", only the
fold/scar class survives (resolution-immune, coincident).
Zeiss r2 projection ~30-40 min / ~8-10 GB - the post-fold-cure
validation dish.

## FRONT 1: PERFORMANCE (Nate's vote - before the next big run)

Bino cost accounting (v29 log): insert samples 97 s (56%!),
extract+repair 41 s, refinement 17 s, sample 14 s, tracer 10 s.
What cursed us (in cost order):
1. LIVE TRIGGER OVERSHOOT: flags 1,653 of 1,663 bino leaves ->
   level 2 nearly BLANKET (zeiss: 6,927/7,703).  Point count
   drives insert time AND the 16 GB.  The bar (16) is calibrated
   to bench models, not to real-model live distributions.
2. Point count generally: every signed sample enters the DT;
   deep-interior witnesses contribute nothing (perf review #1).
3. Strips attempt-0 re-run: one extra full stage-A when strips
   promote (every model with joints).
4. Small change: dup guard O(polys x pts), trace-local extra
   march steps, weld/flip/decimate/recount (~4 s bino total).
QEM/decimation shrink the OUTPUT, not the runtime - the curse is
upstream point count.

TIME=2 PROFILE, r2 bino (2026-07-18, the perf-war map; TIME=2 =
sub-stage lines + per-repair-round anatomy + EVAL tally):
- EVAL grand total 73 s / 105M pts / 36.5K calls = ONE THIRD of
  the 221 s wall is single-threaded field evaluation.  The P5
  threading prize in one line.
- extract+repair 97 s = 5 rounds x ~20 s (signs+extract 7->12 s
  growing, detect+insert ~10 s each); depth plateaus by round 2
  - rounds 3-4 are ~40 s of heat.  Stall-exit tuning + P4
  incremental detection both aim here.
- insert points 51 s (one-by-one CCDT surface/Steiner inserts -
  the batching investigation).
- fix stages 5.3 s total (weld 1.3 / flip 0.9 / decimate 2.6 /
  recount 0.6) - the specialized stages are CHEAP; suspicion of
  them was misplaced.
- Instruments: STIBIUM_DMESH_TIME=2, WITNESS counter (review #3:
  173 sign-witness overwrites per bino - real, consequence
  unknown, fix waits for a correlated defect).

PERF ROUND 2 - LANDED 2026-07-18 (3f4da55e): P5 parallel descend
+ range-parallel bisect (THREADS env, default hw-2, deterministic
merge), surface-insert spatial sort + hints, repair stall exit
(STALL env, default 1, round-2 floor).  r2 bino 221 -> 165 s;
r1 bino 49 s (170.6 at the night's start = 3.5x); depths/water-
tightness identical, screws tilt 2.84 -> 3.12% (stall trade, eye
referee pending).  RSS 1.78 -> 1.93 GB (worker transients).
INCREMENT 2 (same day): eval_points_mt - range-parallel base-tape
batch evals, MIN_VOLUME-aligned (BIT-identical output, verified),
at the 14 big sites.  r2 bino extract+repair 54 -> 32 s, total
165 -> 141 s.  REMAINING ROCKS, now genuinely CGAL-bound:
insert points 45 s (CCDT bookkeeping - P6 research), insert
samples 28 s, tracer 18 s (march batches too small for mt;
chain-parallel is the shape), refinement 15 s, referee 9 s
(non-eval dominated).  STALL knob surfaced per Nate (default 1,
=2 restores patience; tilt 2.84 -> 3.12% trade accepted
"sub-slicing scale").

Battle plan (each is A/B-able on bino in minutes, referee-gated):
- P1 LIVE-BAR CALIBRATION - RUN 2026-07-17 night.  Census: bino
  live median 91 (bench 4-5; live is partly CSG-clause count, not
  pure crowding), so bar 16 = blanket by construction.  LIVE=96
  DEAD (flood+cores rebuild the blanket, 0 refund).  LIVE=128:
  1.3x only, and with THIN in place the blanket's cost mostly
  vanishes - VERDICT: keep LIVE=16, the witnesses were the cost,
  not the surface sampling.
- P2 SAMPLE THINNING - LANDED 2026-07-17 night, default OFF
  (STIBIUM_DMESH_THIN=rings, -1 off): band + shell + block
  corners + exact-zeros enter the DT, rest dropped (bino: 93%).
  Bino 170.6 -> 56.3 s (insert 97 -> 8 s), screws@MAX3 12.1 ->
  6.2 s; tris/constraints/watertight/tilt all held on both.
  ONE OPEN DEFECT gates the default: bino worst chip 0.092 ->
  0.309 sp at ONE step corner (-8.31, 57.37, 74.09).  Cause
  NON-local: site witnesses all KEPT (THINDBG), shell 1/2/4
  byte-identical chip, refinement converges (1 round, 0 leftover
  i/o edges), AND separator counts match (31.5K thin vs 31.7K
  base - the chord-separator theory died too, 5th of the night).
  CONVICTED: the SNAP PASS minted it POST-repair (both repair
  loops plateau 0.242 at the same site; snap tents baseline down
  to 0.092 but one wrong-attribution tent under THIN builds
  0.309; SNAP=0 kills it; snap also deepened screws 0.125 ->
  0.177).  CURE LANDED: per-tent referees - DAMAGE (probe the 4
  new tent edges, refuse tents out-defecting the chip they cure)
  + CHURN (refuse apexes < 0.02 sp above the split triangles'
  planes).  Gated results: bino THIN 0.309 -> 0.170 (187 damage
  + 2,868 churn refused), bino base 0.092 -> 0.097 (held),
  screws 0.177 -> 0.125.  Cost: index nm rises (243 -> ~352,
  both paths - churn tents were papering over pinch sites);
  geometric 0-open holds everywhere.
- P3 STRIPS WITHOUT RE-RUN (half-day): promote from the pass-2
  soup incrementally instead of full stage-A round 2.
- P4 INCREMENTAL REFINEMENT (half-day, medium): only re-scan
  edges near fresh inserts (48-round full sweeps today).
- P5 EVAL-SIDE THREADING (day): octree subtrees independent;
  tape_ctx per thread is designed-for.  ~8x on sample phases.
- P6 CGAL PARALLEL SPLIT (days, research): bulk points in a
  Parallel_tag plain DT + CCDT only for crease law.  The only
  route to parallel insert - CCDT TDS is not thread-safe.

## MARATHON 2 QUEUE (2026-07-18, five-agent fleet synthesis)

Reports: doc/reviews/2026-07-18-{quality-research, correctness-
round2, bino-hygiene-autopsy, perf-round3-design, guarantees-
research}.md.  Fleet's own fixes already landed: TRACE=0 unbraced-
if bug, fresh dialog-knob reads, weld_slivers fan-claim.

TIER A - cheap certain wins (each <= half-day):
1. Decimation ORIENTATION-VETO on collapse (autopsy: reversed
   150 -> ~6; DECIMATE=0 proves decimation mints them).
2. NON-MANIFOLD EDGE-SPLIT at the tail (autopsy: bino's 336
   prusa-opens are 318 real pinch seams, net-winding 0, ZERO
   geometric holes; split -> manifold=yes, +~600 verts, geometry
   byte-identical).  The "0 auto-repairs" finish line.
3. STRIPS RE-RUN KILL (perf 1a: reuse soup.tchains - creases are
   density-invariant; ~28 s/run hidden tax, log swallows it).
4. HASH-GRID coincidence guard for surface inserts (perf 2a:
   10-18 s at r2, bit-identical).
5. Minors: cancelled export writes partial file; SEAL/WINDING
   decouple from the DECIMATE master switch.

TIER B - quality marathon core (1-2 days each):
6. HESSIAN CURVATURE-VALLEY SEEDS as SOFT snap targets (quality
   P1): the 993 FAR untraced-blend divots march as min-curvature
   valleys [Ohtake-Belyaev-Seidel] where the clause system has
   its double root.  Build the Hessian stencil once; then
7. chainless-curvature density trigger (cone quilting) + per-
   curvature churn-gate recalibration (quality P2/P3).
8. OPEN-EDGE PREDICTOR for retreat (perf 1b, gate on a stability
   instrument first): doomed attempts skip the repair tail,
   ~3 min off the full assembly.
9. CHAIN-PARALLEL TRACER (perf 3a: consumed[] is an optimization
   not correctness; per-thread tracers + post-merge dedup).

TIER C - the guarantees ladder (the differentiator):
10. RUNG 2 FIRST (prerequisite-free): tri-tri sweep via the
    in-tree MeshQuery AABB + guards on mutating passes ->
    "0 self-intersections, VERIFIED" per export.
11. P0a outward-rounded intervals + P0b interval-AD evaluator ->
    RUNG 1: per-facet certified deviation (max|f|/min|grad f|
    over facet boxes) -> "every facet proven within X mm" - the
    claim NOBODY ships (nTop/TetWeave: by-construction only;
    metrology: sampled scan-vs-CAD only).  NOTE: kernel
    intervals are decision-sound (NaN taint, fuzzer-proven) but
    NOT outward-rounded enclosures - P0a is real work.
12. RUNG 3 topology audit ("no missing feature above scale h,
    interval-swept") - the honest cheap cousin of certified
    isotopy (full P-V is research-grade, dies on creases).

PARKED/DOA (fleet-verified, do not revisit blind): DT-reuse
across retreat attempts (CCDT bookkeeping risk), SEAL epsilon
raise (fuses REAL 0.012-0.5 mm walls), WINDING gate loosening
(mints backwards edges, measured), "r2 will dilute razors"
(counts triple, resolution-immune), sliver exudation (no tets in
the output).

## OPEN DEFECT CLASSES (quality queue, post-perf)

- RAZOR/Z-FIGHT SCARS - ***SOLVED*** (2026-07-18, three rounds,
  ten+ theories).  FINAL STORY: raw extraction is INNOCENT (the
  tet-complex disjointness guarantee is perfect - strict-interior
  fightpix reads ZERO on every nodec mesh; the earlier
  "extraction mints 3,600 folds / decimation janitors" story was
  an INSTRUMENT ARTIFACT: boundary-tolerant rasterization counted
  shared-edge pixels of coplanar NEIGHBORS as fights.  Verify the
  instrument - one full-precision anatomy dump exposed it).
  The trio bisection convicted DECIMATE_FLATS SOLO (522 fights;
  weld_slivers solo 0).  THE HOLE: collapses decided on ORIGINAL
  positions, applied at pass end - a triangle with two vertices
  collapsing in one pass gets double-remapped into a
  configuration neither orientation guard evaluated.  THE FIX:
  independent-set rule (an accepted collapse claims every fan
  vertex for the pass).  RESULTS: screws 571 -> 0/0/0 all axes
  (tris/depth unchanged); bino 443 -> 59 (razors 1,019 -> 823;
  residue = small next hunt, snap tents suspected), time/depth/
  watertight held everywhere.  Solo gates (WELDSLIV / FLIPSLIV /
  DECFLATS) + DEGEN pass + histograms remain as instruments.
- CLOSE-PERIMETER can of worms (bino): air-chords + pockmarks on
  close-together perimeters.  The strips/level-3/weld composite
  class - reopen AFTER perf (Nate's call).  Assets: NM_DEBUG
  anatomy, strips machinery, WELD (default off, 4-open tear
  autopsy pending at (-9.07, 47.64, 51.94) w/ weld=0.1).
- QEM PHASE 2 (vertex_lock on constrained creases + LockBorder):
  the octagon insurance; ~3-5 h; design in doc/reviews/ quadric.
- POST-QEM CLEANUP: meshopt legally tiles flats with edge-to-edge
  skinny tris (the "scratch" chords, top-edge-to-bottom-edge).
  Fix lives in the EXPORT path after simplifyMesh.
- CONE QUILTING (lamp): r2 cures it; signal = chainless curvature
  (Nate's design seed).  Density campaign trigger 2 of 3.
- THIN TAPERED CONES (knob): missing material at sub-lattice
  caps; hidden-thin trigger (density campaign 3 of 3).
- SLICER HYGIENE - ROUND 1 LANDED 2026-07-18 late night.
  Referee: prusa-slicer --info (the formal harness judge).
  Passes at the pipeline tail: SEAL (epsilon-weld at 3e-4 sp +
  degenerate cull + opposite-twin ANNIHILATION - sealed pinch
  seams mint coincident zero-thickness walls that cancel in
  pairs) and WINDING (flood-fill orientation across 2-facet
  edges, ONE field vote per component, gated to components >= 8
  tris with >50% decisive vote - ungated fragment flips MINTED
  backwards edges, measured).  Env: STIBIUM_DMESH_SEAL (sp
  units, 0 off), STIBIUM_DMESH_WINDING=0 off.
  SCORE: screws 56 open/228 rev -> MANIFOLD=YES, 0 open, 76
  rev.  Bino 528 open/14 rev -> 336 open/150 rev (nm 504 ->
  318).  OPEN QUESTIONS: bino's residual 336-open anatomy
  (T-junctions? seams past 3e-4 eps?), and reversed-count
  attribution (p6 stall-exit-era baseline never prusa'd -
  the 14 -> 150 delta may predate the seal).  Screws-class =
  DONE; bino-class = next hunt with the same referee.
- Regression tests for the eyeball-only bug classes (architecture
  review #6: two-pass liveness, ship-best, geometric-vs-index
  opens, spacing tie-flip, trust-gate cliff).
- INTEGRATION SPRINT (2026-07-18, Nate's go): the mesher is
  named STIBNITE (antimony -> Stibium -> its razor-sharp crystal
  form; env prefix stays STIBIUM_DMESH_* for compat).  UI:
  Stibnite/Classic dropdown, quality Standard(r1)/Beautiful(r2)/
  Extreme(r3+, untested), simplify mm spin; advanced: density
  cap, decimation, snap, stall patience, auto-density (the
  sparkle reveal); threads env-only (auto hw-1).  Post-export
  STATS DIALOG from DMesh diagnostics + phase times.  Classic
  mesher kept (low-RAM constituency + release insurance).
- MANPAGE: ship stibnite(7) covering all five knob tiers
  (visible/advanced/frozen/experimental/instruments) with
  defaults, ranges, and graveyard citations - MESH-WAR is the
  source text.

## GRAVEYARD (tried + rolled back; do NOT re-derive)

- Seam tracing {fA=fB, f=0}: double root at G1 contact, Newton
  cannot march it.  Refuted on paper 2026-07-16.
- Live-count as TANGLE gate: collars and damage sites both read
  55-201.  (As DENSITY trigger it landed - different use.)
- Anti-parallel mindot tangle gate alone: grooves read anti-
  parallel too; only the SIGNED wall gap separates - and even
  that overlaps.  Final answer: retreat rollback (measure damage).
- Blanket level 2/3: nm 241 / 7x nm.  Coverage must be CONTIGUOUS
  (flood+cores); isolated dense islands scored WORSE than d0.
- Straight-to-flood retreat demotion: 8x pitch cliff mid-band
  minted 26 opens from 0.  Graduated (one level per conviction).
- Phantom QEF rejection: punched gaps in fallback chains
  ([.dchain]).  Project-don't-reject + flag stage-D pre-projection.
- Per-chain trust conviction: strips oracle-verified law (bino
  0 -> 14 open).  Model-global nuke on traced chains: torched all
  law over local rejects (autod20).  Gate is FALLBACK-ONLY.
- Repair keep-out: starved groove repairs (thorn crowns); default
  0 since 2026-07-17; crowding guard covers slivers.
- Law-blind crowding dial (CROWD_LAW): plateau - fresh self-
  similar chips replace pressed ones.  Parked at 1.0.
- Repair budget 48 rounds: 12x topology damage, zero depth gain.
- SPREAD_DOT widening (0.95/0.97): poisons suppression + fallback
  radius graph ([.dchain] csg).  Answer = tracer-only shallow
  channel (landed).
- GLOBAL march-step halving: byte-identical chips, [.dtrace]
  noise.  LOCAL band stride + stride-scaled trim: landed.
- Strip level-from-gap: 895/1077 gaps < quarter-cell (real
  near-tangent fits) - refunds only 9%.  Dial parked at 0 (=3).
- Clearance weld default-on: 4-open tear at ONE constrained
  vertex; rollback couldn't reach.  Default off; autopsy queued.
- DC bias inside crease bands: roofed junction air wedges (the
  air-chord conviction).  Suppressed within DC_BAND of law.
- Hidden trigger (_HIDDEN): fires on tangent grazes; needs a
  graze-vs-feature oracle first.  Opt-in, parked.
- Fitted-geometry forensics (my ellipse models): lied twice in
  one day.  The FIELD is the only geometry oracle - FPROBE first.
- LIVE bar raise (96): flood+cores rebuild the blanket from the
  residual trigger's own flags - raising the live bar cannot
  shrink coverage while core dilation spreads every @2 one ring.
  128 refunds 25% but THIN makes the point moot (witnesses, not
  lattice coverage, were the cost).
- Thinning-chip local theories, four buried in one night: local
  witness drop (THINDBG: all kept), shell thickness (1/2/4
  byte-identical), refinement stall, round cap (converges 1
  round, 0 leftover).  The chip is a NON-local tessellation
  shift.  Instruments, not theories.

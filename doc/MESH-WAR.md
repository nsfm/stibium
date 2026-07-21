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

## THE SEED-POACHING FIX (2026-07-20 - the dense-detail class ROOT CAUSE)

Nate's ladder models (examples/detail_spacing_test.sb + _cyl.sb -
PERMANENT BENCH: groove pairs at gaps 8/4/2/1/0.5 sp, flat and
cylindrical) reduced the class to a 30-second repro, and the
autopsy chain ran: defects radiate from junction corners whose
exact vertices are missing -> the incoming chains were never
traced -> their seeds exist but were CONSUMED by an earlier
pair.  MECHANISM (opus agent, verified): with two walls 0.5 sp
apart, the earlier pair's correct() Newton-drags the later
pair's seeds across the strip onto its own already-traced
crease (the moved-away gate rightly allows 1.5 sp), the dup
guard brands them duplicates, and consumed[] is GLOBAL - the
owning pair finds its seed pool looted and traces nothing.
Autodense-dependent because dense seeding covers the whole
wall.  FIX: consume only GENUINE duplicates - the seed's
ORIGINAL position must lie within the dup radius of the matched
curve; drag-in poachers are skipped without consuming.
RESULTS: flat ladder 0 chips / 0 airs / FACE 0.000% (tilted
slivers 8 -> 2); cylinder rims at gap 1 and 0.5: 190 defects
each -> 4 shallow 0.002 sp whispers - THE DEFECT-VS-GAP CURVE
IS FLAT; screws at DEFAULT config: constraints 1,199 -> 1,674
(+40% law), nm -> 0, repairs halved, tilt 2.75%/0.581 deg =
best ever measured, at HALF the level-3 grant's triangles -
"dense feature compliance separate from triangle density",
delivered.  Suite 9/9.  RESIDUAL (root-caused, separate
follow-up): the plate-outline chain marches THROUGH strip
corners without an exact vertex (10.9116 -> 11.0366) and the
through-vertex pre-split does not fire on it - pins both strip
corners 0.0366 off (pre-existing, also on never-poached c5);
plus shallow straight-line chip rows on screws slot profiles
(<= 0.063 sp).  Instruments minted for this hunt (all
probe-gated envs): STIBIUM_DMESH_MARCH_PROBE (march exit
reasons), STIBIUM_DMESH_SEED_PROBE (SEEDINV inventory + SEEDDBG
per-pair verdicts), scratchpad ladder_referee.py /
cyl_referee.py (one-command acceptance: corner exactness +
tilt census + defect-vs-gap curve).  Next: re-referee bino +
full zeiss with the fix (expect the dense-detail class to
shrink model-wide), then the pre-split follow-up.
POST-FIX EYE SCORES (Nate): flat 99.25 -> 99.75%, screws 90 ->
95%, cylinder 97 -> 99.9% ("basically perfect... huge
success").  RESIDUES IDENTIFIED: (1) flat's two remaining tris
both fan onto the TWIN displaced plate-outline vertices
(10.5366 / 11.0366 - the 0.0366 pre-split residual, one per
strip side, mirror-symmetric exactly as Nate's eyes read it).
(2) screws' three spikes radiating into the socket hollow =
fans at three azimuth gaps in the hollow-rim circle law
(z=67.25): surviving-but-untraced features sit exactly there
(PROBE: features-without-polylines = seeds refused; window
census 518 correct-fail / 362 consumed / 112 moved), tris
roof +0.063 sp into the void, MIRRORED on both screws.  NOT a
regression: the same site read 0.197-0.206 sp pre-fix (its
worst faces) - the fix improved it 3x and unmasked it.  Open
question: why three ~120-degree gaps - rim-junction breaks or
another consumption flavor.
TWIN-CORNER RESCUE LANDED (2026-07-20, the pre-split follow-up,
second opus agent): at crowded junctions the wall crease
DETOURS on its final step (last segment (11,-10,2.375) ->
(10.75,-10,2.5) instead of marching to the corner), so no chain
vertex is planted; the exact QEF corner feature exists but
sole_owner corridor-drops it (one visible owner = "redundant"),
and the through-vertex pre-split has nothing to split on -
the outline chain keeps its 0.0366-off vertex.  FIX: keep a
QEF FEATURE that lies on an accepted constraint's INTERIOR
(perp < 1e-3 sp, endpoints > 5e-3 sp clear = genuinely
stranded) - it becomes exactly the through-vertex the
pre-split wants.  Two variants measured-rejected first
(exempt-all flooded the pre-split 32 -> 2,289).  RESULTS:
***THE FLAT LADDER IS PERFECT*** - all six grooves EXACT,
0 tilted, 0 chips, 0 airs, referee VERDICT PASS; cylinder
unchanged; screws tilt 2.71%/0.580 deg (best ever), repairs/
snaps down, constrained 1,674 -> 1,780 (+106 = the slot-
profile chip-row class receiving its stranded-corner law -
the predicted follow-up dividend, deviation accepted).
Rescue counters print under SEG_DEBUG.  REMAINING on the
focus models: ONLY the three socket-rim spikes (hunt #2).

## PRODUCT TIERS + TEXT CAMPAIGN BRIEF (Nate, 2026-07-20)

STAGING (Nate's framing, the campaign's success criteria):
r1 = "suitable for FDM prints" (the clean bench is this tier);
r2 = "suitable for resin SLA" (the text gauntlet's target tier);
r3 = exotic (defer).  The text model (examples/torture/
engraved_text.sb - "stibnite mesher!", the word only 3.8mm
wide, m-legs 0.15mm) is NOT expected to pass at r1; at r2 it
must.  TEXT BASELINE: r1 tilt 20.5%/5.47deg, r2 15.1%/3.11deg
(zero angled faces expected - tilt is error-by-definition);
r2 tilted area by band: walls 6.4 > floors 3.6 > top 1.4 (the
"noise inside engraved areas" per Nate's eyes; edges already
crisp).  OPENING THESIS: engraved strokes are THIN AIR - the
material-only witness gate and the 0.1mm weld doctrine
(correct for assembly contacts) collide with intentional
engraving; the mesher must learn "two parts touching" (fuse)
vs "a letter someone carved" (resolve).  Hidden-cell witnesses
minted only 32 bisections on the whole model - the air class
is deliberately unserved today.

BACKLOG EXPERIMENT (Nate): scale model 2x -> mesh at r1 ->
scale mesh 0.5x, vs native r2.  Prediction from architecture:
nearly identical (resolution enters only via region voxel
counts; all downstream tolerances are sp-relative) EXCEPT
knobs expressed in absolute mm (clearance weld 0.1mm, simplify
deviation).  Any visible difference = a dial that is
accidentally resolution-dependent and should be re-expressed
in sp.  A diagnostic wearing a hack's clothes.

## FIFTH FIX - CORNER-SUPPORT KEEP (2026-07-20 afternoon, a00e4ce2)

Nate's eye sweep of the four-fix bench caught a cylinder
regression (crowded-groove inside corners, 0.177 sp air roofs)
- one-env bisected to the SOLID graze branch stripping corner
support from buried wall-plane sheets.  Fifth agent, fifth
kill: the sheet-normal axis probe separates thin plates (air
along the normal at <=2 cells: demote) from concave corners in
thick solid (air only diagonal: un-demote).  Cylinder restored
to the clean baseline, screws byte-identical, ladder PASS,
suite green.  Bench state: flat 100 percent (Nate), cylinder
clean, screws 99 percent, bino "feels better" at HALF the
triangles (law doubled to 23.8K - the compliance thesis paying
at model scale).  Eye-refuted hypothesis worth keeping: the
screws underside tents are real thin-plate surface, NOT graze
collateral - measured, not assumed.

## FOCUS-BENCH ENDGAME (2026-07-20 late - four fixes, one residue)

The dense-detail campaign closed FOUR mechanisms in one day, each
with its own commit and referee battery: seed poaching (04531b9a),
twin-corner rescue (228c9d82), tangent-graze demotion (3d8d7607),
deterministic razor tie-break (fdec9167).  BENCH STATE: flat
ladder VERDICT PASS (all corners EXACT, zero tilted/chips/airs -
Nate: visually 100%); cylinder ladder flat at every gap; screws
FACE 0.000% over bar / tilt 1.91%/0.438 / repairs 12 / snapped 0
/ 0 open / 0 nm (Nate: "stunning... above 99%").  Suite green
throughout.

ONE RESIDUE, FULLY MAPPED (the ~4-tri underside dimple at
(6.92, -16.5, 64.2), worst 0.025 sp, Nate-visible): an
OFF-SURFACE RIDER - a neighbouring clause pair's curve runs
~0.009 sp inside the trim margin (0.0125 at stride 0.125) below
the true underside edge, traces FIRST by clause order, and
consumes the seeds; the exact plane-pair that owns the edge
never traces.  The poaching pattern BELOW the dup radius: two
curves 0.01 apart cannot be separated by the 0.1 sp identity
story.  FIX DIRECTIONS for a fresh session (touches the
tracer's oldest invariants - do not rush): (a) absolute fdist
cap ~2-5e-3 sp on accepted march points (csg's measured real-
crease convergence is 1.9e-3; the stride-scaled trim admits
riders in dense bands), (b) surface-distance priority in curve
competition (prefer the pair whose curve has smaller fdist when
curves collide within the dup radius).  Evidence anchors: chain
segs at z=64.240021 (FPROBE: f=0.036 there, |g|~4), the dip
vertex (6.8577, -16.5, 64.2006), site improved 0.069 -> 0.025
by graze demotion.  Eye-calibration note: Nate SEES coherent
0.025 sp dimples on flat ground - FACE may want an area-
coherence term (cluster of shallow same-sign faces), not a
lower bar.

NEXT after the rider: bino + full zeiss generalization sweep
(all four fixes are model-agnostic; expect the dense-detail
class to shrink everywhere), then NATE'S TEXT GAUNTLET
(engraved text on flat + cylindrical surfaces - the graduation
exam).

## CURRENT STATE (2026-07-20, post metric-audit)

THE RANKING METRIC IS NOW FACE (area-weighted centroid
deviation, dmesh_face_sweep + the --facedev CLI verb): Nate's
audit caught worst-depth improving with no visual change -
it is an L-infinity EDGE-midpoint statistic, blind to
triangles whose edges hug the surface while their INTERIOR
spans air/material.  FACE is calibrated against his verdicts
(screws base "hideous" = 0.612% of area deviant vs grant
"clean" = 0.024% - a 25x separation worst-depth scored a tie).
Worst-depth remains the REPAIR referee only.  --facedev sweeps
any past STL/3MF against the tape; every 3MF now carries
provenance metadata (config + env + stats), STLs a header
stamp.

DEFAULT ZEISS (the shipping config): 577K tris / 0 open /
worst 0.217 / FACE 0.022% deviant / 14.7 min / 3.87 GB -
statistically EQUAL on FACE to the 1.03M grant build (0.020%)
at 56% of the mass; Nate's slicer bar (prusa complains at 1M)
and his "600k at 99% is the right move" verdict decided the
default.  CROWD_MAX=3 is the documented DETAIL DIAL (screws-
class FACE 0.233 -> 0.024%, eye-visible) for part-scale
exports.  Hidden-cell witnesses stay DEFAULT ON: ~free on
zeiss (577K vs 582K), the screws depth cure, and on bino
their 2.2x mass buys the improvement Nate's eyes verified
(chip removal, air-chords reduced; FACE 0.044 -> 0.038%,
mean 0.0009 -> 0.0006 sp; the mass is real thin-wall sheets
resolving - a worthiness discriminator for sub-pitch walls is
an open question).  NEXT BOSS (found by FACE, invisible to
the old metric): the shared worst face at (-5.26, 53.5, 68.0),
0.39 sp, IDENTICAL in curve_r1 and night_r1 - resolution-
immune, class unknown.

(2026-07-18 state, kept for context:) autod31 was 727K / 0
open / 0.465* / 42 min pre-femto-guard.
Bino: THIN+gates 59 s (Nate: identical to base-gated), r2 in
198 s / 1.77 GB (impossible pre-THIN).  THIN DEFAULT
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

PRUSA INTERNALS VERDICT (2026-07-18, doc/reviews/2026-07-18-
prusa-internals.md - the case cracked): prusa runs TWO judges by
format.  STL -> admesh bit-exact weld + orientation flood (the
150 "reversed" = pure flood artifact, 0 on 3MF; the weld re-fuses
our coincident splits).  3MF -> index-preserving, no weld, no
flood - it TRUSTS our topology (bino: 292/150 STL vs 73/0 3MF,
and its residual nm there is honest self-contact).  REFEREE
PROTOCOL CHANGE: slicer-hygiene checks run on 3MF.  Path to
manifold=yes/0/0: ship 3MF (app default already!), then
SEAM-CLOSURE-CONSISTENT pinch-split on the indexed rep (sheet
labels walked along seam chains so sheets re-close; naive split
in 3MF minted 2,202 boundary opens - measured).  EPSILON
SEPARATION FORMALLY DEAD (nearby-weld ladder climbs ~0.02 mm,
fuses real walls).  Fan-claim guard STAYS (prusa's 3MF manifold
test is topological only).

SEAM-CLOSURE PINCH-SPLIT ***LANDED*** (2026-07-18 night, the
manifold=yes finish line).  Three-part design in
dmesh_split_pinches (public, export-tail):
1. MATERIAL-WEDGE PAIRING: at each nm edge, radial-sort the
   incident faces around the edge axis; a b->a-wound face has
   material on its CCW side, an a->b face on its clockwise side,
   so each material wedge runs dir- -> next dir+ CCW (zero-width
   wedges skipped - doubled-sheet kisses put faces at coincident
   angles and a noise pairing glues the sheets).  The pairs
   COUNT as fan connectivity in the vertex union-find, so sheets
   re-close by construction: 0 boundary edges minted anywhere
   (naive barrier split: 2,202).
2. FROZEN-SNAPSHOT ROUNDS: all reads from a per-round copy of
   the index array; the live-read pattern let later vertices see
   rewritten ids as unknown edges -> spurious barriers.
3. MIDPOINT SUBDIVISION for stuck seams (bino: 157/274 edges
   whose sheets reconnect through BOTH endpoint fans - no vertex
   partition can separate them): cut the edge at its collinear
   midpoint (zero geometric deviation); the fresh vertex is
   entangled with nothing but the seam and splits next round.
   Endpoints stay fused = bowtie vertices, legal manifold-test
   topology.
POSITION: export tail, STRICTLY after QEM - meshopt tears
coincident copies into real boundary holes if the split runs
first, and its collapses mint fresh pinches of their own (bino
274 -> 325 nm).  Tail position cures both.  SCORE (3MF, r1):
screws no/8-open -> MANIFOLD=YES 0/0 parts=1, hist pure
{2: 8448}; bino no/73-open/274 nm -> MANIFOLD=YES, 0 boundary,
1 residual nm (one non-alternating 3.8-degree kiss fan);
bino+QEM0.01 no/80-open/325 nm -> MANIFOLD=YES, 0 boundary, 10
residual nm.  Geometric welded 0/0/0 on all finals (odd-
incidence census); suite 627,666 green.  Instrument:
tools/idxedges.py (index-edge incidence histogram on 3MF - the
hygiene referee).  DEFAULT ON (Nate's call, 2026-07-18 night;
STIBIUM_DMESH_PINCHSPLIT=0 disables).

TIER A SESSION VERDICTS (2026-07-18 day, autonomous run - the
plan met reality and reality won several rounds):
- PINCH-SPLIT: built, works (index-nm 274 -> 0), but coincident
  copies make prusa WORSE (it re-welds them; screws manifold
  yes -> no) and the epsilon nudge opened 2,198 geometric edges
  (per-vertex components disagree along seam curves).  OPT-IN
  (PINCHSPLIT=1) until seam-consistent sheet labeling is
  designed.  KEY REFRAME: prusa's residual opens/reversed on
  pinch models measure THROUGH its own exact-position weld -
  our geometric counter + bit-exact welder + fightpix are the
  geometry truth; prusa-clean is a separate (epsilon-
  separation) engineering goal.
- DECIMATION WINDING-VETO: built, then refuted - the autopsy's
  attribution dissolved under solo gates (prusa reversed is
  mostly its pinch-walk artifact; weld owns ~4-6 real ones).
  The WELD variant tripled the razor census (627 -> 1,744):
  removed.  Decimate veto OPT-IN (DECVETO=1).
- TRACE CARRY: built, law-refuted (re-trace on the promoted
  soup mints MORE law; carried = -4.5% constraints).  OPT-IN
  (TRACE_CARRY=1).  LEAD: the carried mesh read worst depth
  0.097 vs 0.170 with LESS law - is some law near (-2.7, 37.1,
  61.7) hurting depth?  Fresh-eyes candidate.
- HASH-GRID GUARD: kept (semantics-exact, ~1 s - the CCDT
  insert wall is bookkeeping, P6 confirmed).
- OPEN TRADE for Nate: the weld FAN-CLAIM (fleet correctness
  guard) costs screws prusa-manifold (yes -> no, 68 opens) via
  collapse reordering - structural-guard-vs-slicer-cleanliness
  needs an eyes decision.
- NATE'S LEAD (2026-07-18): bino is mirror-symmetric but ONE
  side carries significantly more flaws - order-dependence
  somewhere (insertion order / tie-breaks / float sign under
  mirroring).  Investigation seed.

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

QUALITY ROUND 1 - CURVATURE STENCIL + CHURN CARVE-OUT LANDED
(2026-07-18 night).  The stencil: curvature_probe() - 19-point
second-difference Hessian on the field oracle, principal
curvatures + kappa_min eigendirection.  PITCH LESSON (measured):
h=0.01 sp reads pure static (h^2 denominator at float-eval
noise; 0 of 993 valley attempts survived); h=0.1 sp is the
working pitch.
ANATOMY CORRECTION (SNAP_DUMP instrument, bino r1): the quality
report's FAR model was WRONG for bino - the 993 "FAR" chips are
sub-floor depth noise (median 0.006, max exactly 0.0100 sp = the
floor boundary; surf projection almost never overshoots).  The
DEEP skipped populations are CHURN (2,123 of 2,680 >= 0.01 sp,
bunched at 0.012-0.019 RIGHT UNDER the 0.02 bar - a gate
artifact: on a curved blend a chord's sagitta L^2*kappa/8 reads
0.01-0.03 sp, real bend geometry) and damage (672 of 733).
THE CARVE-OUT: churn refusals re-judged by measured curvature -
allow the tent when apex height <= 2x the chord sagitta at the
probed kappa (flats read kappa~0 and still refuse; 0.005 sp
absolute floor stays; saved tents still face the damage
referee).  STIBIUM_DMESH_CHURN_CURVE=0 disables.
BINO SCORE: worst depth 0.170 -> **0.093 sp** (UNDER the 0.1
product bar; better than the 0.097 no-THIN baseline - THIN's
whole quality cost recovered).  Snapped 795 -> 1,495 (875 churn
refusals saved, 173 pushed to damage - the referee chain works).
GUARDRAILS: fightpix +z 17 -> 14 (down!), +x/+y flat 0/1
unchanged, screws fightpix 0, screws tilt 4.03 -> 4.20% area
(+0.17pp, stall-trade scale), screws worst 0.125 held, suite
627,666 green.
VALLEY SNAP: implemented (valley_project - surface-project,
kappa_min eigendirection, cross-valley parabola to the OBS04
extremum, iterate x3; polyline-attribution acap; damage/churn
refereed; STIBIUM_DMESH_SNAP_VALLEY=0 disables) but IDLE at bino
r1 - no deep FAR population exists there.  The dish is zeiss
(FAR divots 0.01-0.59 sp per the surf-floor measurement).
Instruments: STIBIUM_DMESH_SNAP_DUMP=path (per-chip outcome +
depth), VALLEY counter line under CHIP_DEBUG.
ZEISS r1 VALIDATION (curve-r1, 13.9 min / 3.52 GB): 3,297 churn
refusals saved by curvature (collar-texture cure at scale);
VALLEY verdict: 24 attempts, 1 genuine deep-FAR valley cure, 20
floor-rejects - SNAP_SURF + the carve-out own the deep
population now; the valley path stays as the safety net for
true blend-throat divots (costs nothing when idle).  Worst
0.364 sp UNCHANGED at (4.15, 15.94, 89.01) - not a churn-class
chip; separate hunt.  AND: the exported 3MF (pinch split now
default) reads incidence {2: 809676} PURE - **prusa
MANIFOLD=YES on the FULL ZEISS**, the whole hygiene arc closed
end-to-end.

QUALITY ROUND 2 - THE WIREFRAME SPLOTCH CASE ***CRACKED***
(2026-07-18 late night; Nate's lead: quilted low-tri polygonal
patches, random on cylinder walls).  Instrument:
tools/areamap.py (per-pixel front-triangle area heatmap; red =
big).  TWO mechanisms separated by A/B, THREE theories buried:
- Big irregular plates: decimate_flats merging TRUE FLATS
  (dihedral autopsy 0.000 deg, axis/45-deg normals) - legal and
  correct, NOT the splotch.  A curvature VETO now caps plate
  growth on curved surface anyway (batched Hessian at each
  candidate; refuse when fan-span sagitta kappa*L^2/8 beats
  STIBIUM_DMESH_DECCURVE*mesh_edge, default 0.05; flats merge
  freely).  Bino vetoes 10.2K, tris +2K.
- THE SPLOTCH: TANGLE-DEMOTION SPECKLE.  Bino is blanket
  level-2 (1,562/1,665 @2) but 404 leaves tangle-demote to
  level 1 wherever the per-leaf sep reading (min over corner
  pairs - noisy) flickers across the 0.7-cell bar on thin
  walls: isolated 2x-coarse islands = the polygonal splotches.
  CURE: DESPECKLE - a real tangle band is contiguous (its
  neighbors demote too, majority filter leaves it alone); only
  isolated level-1 islands with >= 5 touching level-2 neighbors
  re-promote.  Bino: 102 islands -> map fully uniform @2,
  splotches GONE in areamap, worst 0.093 held, 0 open, +2K
  tris.  Buried theories: live=0 pockets (bino has exactly 2
  smooth leaves), leaf-corner mindot trigger (thin walls put
  both sheets in every leaf, mindot pegs -1 - measured on
  lamp).
- CHAINLESS-CURVATURE TRIGGER (density campaign 2 of 3) landed
  with the honest signal: PROJECTED Hessian probe at the leaf
  (leaf_curve_theta - project center to surface on the BASE
  tape, kappa*cell = facet angle/cell).  Band bar-to-25deg
  (STIBIUM_DMESH_CURVEBAR, default 3 deg/cell; ceiling =
  crease grade).  Bino 158 flags / lamp 84 - mostly subsumed
  by blanket coverage today; the trigger matters where coverage
  is sparse.  Smooth-pocket fill also in (live=0 leaves
  recorded + filled; rare - bino 2).
- LAMP CONE VERDICT: already blanket level-1/2 at r1 - the
  residual quilt needs LEVEL 3 (the queued MAX=3 referee), not
  a new trigger.
- Dense-knob smoke test isolated from the trigger (CURVEBAR=0
  around the DENSE=1 smoke; bench pitch legitimately reads
  several deg/cell).
- SMOOTH-BRANCH CURVE TRIGGER - ***QUILT SOLVED*** (2026-07-20,
  Nate's V-patch coordinate (-3.775, 32.14, 40.98)): the
  chainless-curvature trigger only probed CREASE-suspect
  leaves, but a tilted cylinder wall is a live = 0 SMOOTH
  expanse - and the smooth-pocket fill needs >= 3 touching
  promoted CREASE leaves, which a big smooth region's interior
  never has.  DENS autopsy: 9 leaves at the patch, all smooth,
  7 at level 0 - the V is exactly the set beyond fill reach.
  Fix: the same leaf_curve_theta probe in the smooth branch
  (same band, same magnitude ladder; the projection clamp
  already refuses out-of-reach leaves).  Zeiss referee: the V
  region 125 -> 876 tris, max tri area 0.813 -> 0.061 (13x),
  worst/opens/tris/nm all held.  Instrument minted: the DENS
  census rows (C/S tagged) - built for this hunt, kept for the
  next one.
- HIDDEN-CELL WITNESSES - density campaign 3 of 3 LANDED
  (2026-07-20, the dense-details autopsy; Nate's #1 class:
  "when details are densely packed, chipping and air chords
  abound").  AUTOPSY CHAIN, each step measured: all screws-head
  defects cluster at z 64-67 -> 1,023 of 1,067 sit within
  0.25 sp of a traced chain (NOT seed starvation - law is
  complete and the mesh disobeys it at close range) -> the
  worst-site corridor (two law rings 0.25 sp apart) has 7
  mid-row vertices where ~40 stand on each ring -> corridor
  drop, decimation, thinning all exonerated (one-variable) ->
  THINDBG ground truth: 136/136 lattice points OUT across the
  whole band even at level 3 - the rim is thinner than the
  LOCAL pitch, the lattice cannot see the material, and the
  surface hangs from law chains alone.  The leaf-level hidden
  oracle can't reach it (it requires an ALL-quiet leaf; the rim
  hides inside a mixed leaf).  THE RUNG: cell-level certify
  witnesses - in pass-2 dense leaves (level >= 2),
  sign-unanimous cells within 2 cells of a crossing cell are
  batch-interval-checked; straddlers run certify_hidden with a
  PROVEN-BOX COLLECTOR (cap 4), and each proven center pairs
  with 4 tetrahedral cell corners as guaranteed sign-change
  segments into the ordinary bisection queue.  Proof-carrying
  witnesses: minted only on interval PROOF of opposite-sign
  volume, never on grazes.  MATERIAL-ONLY GATE (bino referee,
  measured): the symmetric version certified 8,225 cells /
  127K bisections on bino and dragged worst 0.091 -> 0.372 -
  proven-AIR cells inside solid are near-tangent assembly gaps,
  the class the clearance weld fuses on purpose (sheet-
  separation verdict); resolving them re-mints the pinches.
  Hidden MATERIAL in air is the rim/thin-wall class and cures.
  RESULTS: screws worst 0.177 -> 0.061 sp (3x, well under the
  0.1 product bar; 29 cells, 464 witnesses; tilt 2.99% /
  0.662 deg - the best of the whole density round); bino
  neutral (224 cells, worst 0.096, nm slightly down, tris
  flat); knob/suite untouched.  STIBIUM_DMESH_HIDDEN_CELLS=0
  disables.  Instruments minted: STIBIUM_DMESH_CHIP_DUMP=path
  (every defective edge midpoint, BOTH signs - the
  localization referee), STIBIUM_DMESH_DROP_PROBE="x,y,z,r"
  (corridor-drop survivor census).
- LEVEL-3 CROWDING GRANT, TWO-SIGNAL RULE - ***LANDED DEFAULT
  ON*** (2026-07-20, the MAX=3 referee round; three configs
  measured + Nate's eyes on every leg).  VERDICTS: (1) blanket
  MAX=3 refuted AGAIN (bino: residual formula blankets 1,555 of
  1,665 leaves @3 -> worst 0.093 -> 0.372, nm 6.4x - same
  failure as 2026-07-17, now with the mechanism named).  (2)
  live-magnitude targeting refuted: bino's blend bands honestly
  read live >= 64 in 1,288 leaves - crowding alone re-blankets.
  (3) the MEASURED discriminator is signal COINCIDENCE: screws'
  grant leaves are 92% residual-HOT (sub-lattice slots show
  both extreme crowding AND QEF residual), bino's false grants
  86% residual-SILENT (assembly collars crowd the clause count
  without sub-lattice geometry).  RULE AS SHIPPED: level 3 iff
  live >= 4x bar AND flag_leaf fired on the leaf (Collector::
  resid_hot) AND not tangle; decided POST-survey in
  delaunay_sample (the in-descent gate cannot see residuals);
  ceiling crowd_max_level() (STIBIUM_DMESH_CROWD_MAX, default
  3) - the residual formula stays behind AUTODENSE_MAX=2.
  RESULTS: screws tilt 4.03 -> 3.13% of area, weighted 0.947 ->
  0.711 deg, 0 geometric opens (Nate: base "hideous", all
  level-3 variants visually identical - the 0.125 -> 0.177
  depth delta does not register); bino 664 @3, worst 0.091
  (BEST EVER, beats 0.093), Nate-verified visual improvement
  (chip removal at sharp edges, air-chord severity down; chips
  remain near certain sharp-edge classes - a named residual
  class, not a regression); zeiss 1,472 @3 of 7,882 (19%,
  targeted), WORST 0.364 -> 0.217 (the eternal (4.15, 15.94,
  89.01) site moved; new worst at (-28.65, -23.56, 17.84)),
  1.02M tris / 26.5 min / 6.99 GB (1.75x/1.9x/2x - RAM well
  inside the laptop envelope; Nate 2026-07-15: speed is not a
  concern for final meshes).  COST NOTE: the bill is real -
  CROWD_MAX=2 reverts to the pre-round default in one env.
  Instruments minted this round: DENS rows in the census dump
  (final per-leaf density map, C/S tagged - the density-shadow
  referee), tools-side join scripts in the session scratchpad.
  Screws' index nm 1 -> 123 is kiss-fan class: openedges_exact
  reads 0 geometric opens; pinch split owns it at export.
  Determinism note (Nate's eye): the two screw instances now
  fail IDENTICALLY (mirrored defects) - residual error is pure
  deterministic geometry response now, no randomness left in
  the failure class.

INTERVAL ARC OPENED (2026-07-19, strategy doc: doc/reviews/
2026-07-18-interval-certify-reach-locus.md - "cull -> certify",
one enemy = reach collapse; Nate's go).  MOVE-1 CHEAP HALF
LANDED: certify_hidden - the graze-vs-feature oracle the HIDDEN
flag always needed (named missing by two prior reviews + the
code comment).  Interval subdivision on the leaf's pushed tape
via tape_eval_i_batch (batch mode records no push state = safe
mid-pipeline).  VERDICT RULE, refined by the first referee
round: "ambiguous at pitch/2" MISFIRES on internal tangencies
(the stacked_cylinders 18-leaf control stays ambiguous at every
depth - material kissing material never decides) - so
FEATURE_PRESENT iff subdivision PROVES a sub-box of the sign
OPPOSITE the lattice's unanimous read (a hidden feature is,
by definition, missed opposite-sign volume).  Tangency shells,
loose bounds, and sub-pitch/8 unrecoverables resolve EMPTY.
SCORE: stacked_cylinders 18 -> 0 features/18 proven grazes (the
2026-07-16 killer population, dismissed by proof); off-axis lip
144 -> 3 honest features; knob 51 -> 0 (caps taper below
pitch/8 - consistent with "still hidden at 4x"; tip cure is
explicit geometry, not density); bino 1,188 -> 5 features.
STIBIUM_DMESH_HIDDEN: 2 = oracle (DEFAULT ON), 1 = legacy
blind, 0 = off.  Instrument (Nate's wish): STIBIUM_DMESH_
DUMP_HIDDEN=prefix writes verdict boxes as STLs
(prefix_feature.stl / prefix_graze.stl - viewable beside the
mesh like the chain tubes).  Down payment on the guarantees
d(S->M) coverage rung.  Next in the arc: #2 collapse the reach
detectors into one LFS estimate (plan as a refereed
recalibration campaign, NOT a neutral refactor); #3 interval
prefix eval (tape surgery) gated on the arc proving out; #4
contact-as-constraint last.

#4a CONTACT TRACER LANDED (2026-07-19, Nate's go - the snap-
target on-ramp): trace_contact_chains - seeds from contact
verdicts, pulled onto the curvature-extremum line by
valley_project, marched predictor-corrector (predictor
n x t_min, corrector = valley_project; the clause system
double-roots at contact but curvature stays finite).  PROVEN:
off-axis rim traced as 10 chains/5,200 pts (laps its loops -
closure detection TODO), bino 28 chains/492 pts from 1,077
contact boxes, ~5 s cost.  Chains join the snap law;
STIBIUM_DMESH_DUMP_CONTACT=path dumps them as fat tubes.
DELIVERY REFEREE'D HONEST: the off-axis air-chord class reads
SECOND-ORDER field depth (a tangency chord hovers ON the
surface while crossing the seam LINE) - below every snap
floor.  Wart admission built (warts near contact chains enter
the snap sweep) and correctly admits ZERO there.  VERTEX RAIL
built (chain points as plain DT vertices,
STIBIUM_DMESH_CONTACT_POINTS=1, OPT-IN: near-site skinny 1.2 ->
2.0%, eye verdict pending).  CONCLUSION: output surgery cannot
own the tangency class - the seam must be IN the triangulation.
Full #4 (contact chains as CCDT constraints) is the cure, now
DE-RISKED: the marcher is battle-tested, only constraint
insertion remains, gated behind #3 per the strategy doc.  Ctx
lesson: collector ctxs die at the descent merge - post-merge
passes make their own.

TRACER ROUND 2 (2026-07-19, d8365c8f - Nate's tube inspection
cracked it): the traced ring was a BURIED ZERO-SHEET - CSG-
coincident faces leave internal f=0 surfaces (the cylinders'
interface disk + the tilted cylinder's wall persisting inside
the union) that samplers never see and intervals honestly
cannot decide.  FIXES: (1) SIGN-CROSS filter - real surface
reads opposite signs at +/-0.25sp along the normal, buried
sheets touch zero one-sided; applied per seed + per march step
(off-axis: 106/122 seeds provably buried, bino: 278/1,077).
The contact dump is now ALSO a CSG-coincidence hygiene
instrument.  (2) RIDGE MODE in valley_project (kappa_max
extremum lines) - the defect rim is a fading CONVEX crease
with no valley to find; contact traces try valley then ridge
at a 10x lower floor.  Verified: off-axis 2 lip ridges, bino
43 chains.  REMAINING GAP: graze-azimuth seeds Newton-project
onto the adjacent flat, outside the extremum capture radius -
next increment is ENDPOINT SEEDING (the crease tracer stops
exactly where the rim fades; its open chain endpoints are
perfect ridge-mode seeds - continuation, not discovery).

TRACER ROUND 3 - THE PRIMITIVE VERDICT (2026-07-19, Nate's
semicircle question + the on-axis control): ridge-as-PEAK is
the WRONG detector for fillet boundaries.  The on-axis lip
traces ZERO chains (28/28 vproj-rejected): both surfaces are
constant-curvature plateaus (torus meridians included) and the
junction is a curvature STEP the stencil smears into a
monotonic ramp - no local max exists ANYWHERE on the honest
geometry.  The off-axis semicircle was tilt-created kappa
peaks, not the junction.  Nate's instinct ("shouldn't this be
a ridge no matter how it's tilted?") was right twice: the
eigen-ordering also swaps with tilt (kappa_max direction
rotates 90 deg where 1/R_local crosses 1/r_tube).  NEXT
PRIMITIVE (fresh session): a curvature-STEP detector - march
the extremum of d(kappa)/ds across the line (the ramp's
steepest point), sampled from the same stencil; cross
direction from the MARCH FRAME (n x dir), never from eigen
ordering.  Note: the MESH is unaffected throughout - contact
chains feed only the snap pass and the tangency chips are
sub-floor; the tracer is an instrument still earning its
delivery mechanism (full #4 constraints).

TSEED STEP-TRACING LANDED - #4a ROUND 4 (2026-07-19, both
pre-agreed acceptance tests PASS): trace_step_seams runs after
delaunay_trace, seeds step-marches from soup.tseeds (the
shallow channel), claim map preloaded from tchains + contact
chains so crease law and step law never double-cover.  Seed
ladder per unclaimed tseed: newton7 -> mc_tgrad (wake bar
0.005/sp - blend interiors are mean-curvature-flat and reject
free) -> step_project -> re-check claim at the CENTERED point
(first seed on a ring traces it, every later seed centers onto
the claimed line and stops - 1,492/1,496 on-axis seeds dedupe
this way) -> sign-cross -> march both ways.  LOOP CLOSURE built
(march_xline returns closed when it re-enters 0.75 step of the
start after real progress; the round-2 lap TODO is done).
ACCEPTANCE MEASURED: on-axis lip = 1 CLOSED ring, 165 pts,
24/24 azimuth bins at z=16/r~13.3 (the z=15/r=14 ring is
CREASE-LAW-OWNED - the c0-rim clause pair is a real tape crease
even though the union is G1; correct dedup, both rings covered).
Off-axis = 1 CLOSED ring, 168 pts, 24/24 bins after untilt
(r 13.32-13.53, z 15.81-16.01) - ORIENTATION INDEPENDENCE IS
NOW AN EXECUTABLE REFEREE AND IT PASSES.  Referees: bino
IDENTICAL to splotch-r1 (133,316 tris / 0 open / 287 nm /
12,109 cstr / worst 0.093) with 0 false step seams (909 seeds:
904 crease-claimed, 4 step-rejected, 1 short); screws worst
0.125 held, pass inert (no tseeds on thread geometry); suite
9/9 green.  REFACTOR: tracer probes are shared statics now
(cc_key/cc_claimed, probe_normal7, probe_real_surface, newton7,
mc_tgrad, march_xline) - the mc-gradient block was duplicated
twice inside trace_contact_chains and is one function.
STIBIUM_DMESH_STEP_TRACE=0 disables.  Exhibits for Nate:
build/zeiss_dmesh/{onaxis,offaxis}_tseed_{ring,mesh}.stl.

OPEN - THE SPIKE-RUNG QUESTION (opt-in STIBIUM_DMESH_TSEED_
RIDGE=1; three discriminators measured-refuted 2026-07-19,
read before resuming): the cone-torus lip carries a true ~13
deg crossing ring at z=15/r=10.5 (the torus outer equator
meets the CLIPPED cone's rim - note the model is union(torus,
c1, c0 INTERSECT cone): the cone stops at z=15, there is no
upper crossing; a chord-geometry derivation missing the clip
cost an hour).  NOBODY owns that ring: crease law traces only
z=16/r=9.5 + the sharp rims (fine census), and step_project
rightly peak-rejects its seeds (spike, not ramp - the same
strictness that keeps bino clean).  The ridge rung traces it
full-azimuth (24/24) but ALSO chains the smeared-step SHOULDER
0.5 sp inside the off-axis fillet.  REFUTED SEPARATORS, with
mechanisms: (1) PROMINENCE (bars 0.02/0.05/0.1 per sp each
fail one side; RIDGE_PROM knob kept, default 0) - at sp ~ tube
radius the smeared kink height theta/2h ~ 1.1/sp coincides
with the fillet plateau 1/sp.  (2) FOLD-EXCESS (normals +/-d
minus curvature-explained rotation) - total normal turn is
CONSERVED under smearing, so every step "folds"; worse, a
kappa station inside the smear zone reads the spike as plateau
and explains the kink away (cone-side kappa read 0.907 vs
honest ~0); at d=0.4 the shoulder out-folded the kink 0.34 vs
0.13 rad.  Function deleted; do not rebuild it.  (3) SCALE
DIVERGENCE (kmax at h=0.05 vs h=0.1; kink ~ theta/2h diverges,
C1 step saturates - clean theory, right primitive) - reads
1.001 EVERYWHERE because the ridge projector converges ~0.5 sp
AWAY from the kink onto smooth torus (candidates z~15.475 vs
ring z=15.0), and even 0.02-sp station sweeps along kmax-dir
from there never touch the kink.  ROOT CAUSE TO HUNT NEXT: why
valley_project ridge mode is displaced half a cell from a real
13-deg kink (candidate suspects: eigen-direction rotation near
the kink, the 3-iteration parabola walk, project_points_impl
pulling along the surface).  Instrument first: dump the kmax
profile along the meridian through (10.5, z 14.5..16) before
touching the projector.  The scale-divergence referee is
WIRED (STIBIUM_DMESH_SCALE_DIV, default 1.4, station sweep
0.02 sp) and becomes correct the moment candidates land on the
kink.  Damage context: the unowned ring costs 0.008 sp today
(sub-floor) - this rung is about zeiss-class shallow-crossing
coverage, not the bench.

STEP DETECTOR CORE BUILT (2026-07-19 late, 2bfbcd04): hg_from_
vals + kappa_along (raw Hessian shared), step_project (five
kappa stations across the march frame, strict d(kappa)/ds peak,
prominence bar STIBIUM_DMESH_STEP_K default 0.05/sp), tracer
step mode with mean-curvature-gradient seeding (basis-free).
SEED-SOURCE VERDICT (the on-axis control, third referee round
in a row to teach something): contact boxes are the WRONG DOOR
for fillet junctions - the oracle flags tangencies and hidden
things, but a fillet boundary is ORDINARY SAMPLED SURFACE and
never becomes a hidden candidate.  Buried-ring seeds Newton-
project straight up to featureless flats (all 28 land on z=15,
gradient exactly 0; walk-toward-box is the zero vector after
tangent projection).  NEXT: seed step-marches from the SHALLOW
channel (soup.tseeds) - on-surface QEF points in exactly the
sub-feature-grade curvature band a fillet boundary occupies;
they already exist and are already positioned.  Bino referee:
step mode mints no false seams (9 valley chains only) - the
strict rules hold.
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

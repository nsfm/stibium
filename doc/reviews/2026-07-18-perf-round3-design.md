# Perf round 3 — design proposals (2026-07-18)

Read-only design study. Nothing here is landed. Each proposal carries
expected savings, implementation effort, risk, and an A/B measurement
plan gated on the open-edge / watertight referee.

## Grounding measurement (this machine, 16 cores → 15 threads)

Fresh r2 bino profile, `STIBIUM_DMESH_TIME=2` (HEAD b5a8c32b):

| phase | s | notes |
|-------|----|-------|
| sample+bisect+QEF (attempt 0) | 9.85 | strips re-run — see below |
| sample+bisect+QEF (attempt 1) | 10.16 | second full stage-A |
| crease tracer | 18.03 | **only attempt 1 measured; attempt 0's ~18 s is swallowed by the strips `continue`** |
| segment referee | 8.53 | |
| insert samples | 28.94 | 928K sign witnesses, CCDT one-by-one |
| **insert points** | **46.76** | 883K surface pts, one-by-one + per-pt nearest_vertex |
| constraints+restore | 0.37 | |
| refinement | 14.58 | 1 round, 73K inserted |
| extract+repair | 30.72 | r0/r1/r2 signs+extract = 6.7/8.1/9.3 s |
| manifold pass | 2.38 | |
| snap pass | 4.93 | |
| fix stages | 6.35 | weld/flip/decimate/seal/winding/recount |
| **B+C total** | **138.38** | |
| **EVAL grand total** | **33.17** | 18% of wall — the P5/eval war is won; the rest is CGAL bookkeeping |

Bino does **not** retreat (it watertights), so the ~185 s wall here is a
*single* B+C plus the strips re-sample. The full-assembly 13.1 min with 3
retreat attempts is ≈ 3× the B+C tail plus 3× sample+trace.

Key measured facts driving the design:

1. **The strips re-run is an always-on tax the profile UNDER-reports.**
   Every jointed model hits the `attempt == 0 && promote.empty()` strips
   block (delaunay.cpp:8051). It fills `promote` and `continue`s
   (:8124) — re-running `delaunay_sample` AND `delaunay_trace` from
   scratch. Because the `continue` jumps over `pt.mark("crease tracer")`
   (:8127), attempt 0's tracer (~18 s) is invisible in the log. Real
   cost on bino ≈ 10 s (re-sample) + 18 s (re-trace) ≈ **28 s wasted per
   run**, and it scales with the model.

2. **Eval is no longer the bottleneck** (33/185 s). Insert points,
   insert samples, tracer, refinement, extract are all CGAL-TDS
   single-threaded bookkeeping — the exact rocks the task names.

---

## ROCK 1 — the retreat re-run tax

The retreat loop (`delaunay_mesh`, :8018) re-runs the ENTIRE
sample→trace→insert→refine→extract→repair pipeline per attempt. The
`demote`/`promote`/`noweld` maps that change between attempts touch only
**tens of leaves**, yet everything rebuilds globally. Three independent
levers, in ascending effort:

### 1a. Kill the strips re-trace (and ideally the re-sample)

**What.** The strips promotion needs `soup.tchains` (traced chains) and
`soup.dense_boxes` — both already exist after attempt 0's trace. The
crease GEOMETRY is invariant to `promote` (density does not move
creases). So re-running the tracer on the strips `continue` is pure
waste. The *only* thing that legitimately changes is the tracer's
**local march step** inside newly-promoted strip bands (v25/v26: law
delivery at local pitch is load-bearing).

**Design (staged):**
- *Minimum:* on the strips re-run, skip `delaunay_trace` and reuse the
  existing `soup.tchains`/`soup.tclosed`. The strip bands are promoted
  from level ≤2; the base-pitch chains already cover them. Re-mark only
  if `TRACE_LOCAL` law-delivery in the promoted bands is measurably
  needed (A/B the joint depth). **Saves ~18 s/run** on bino, more on
  zeiss.
- *Better:* fold strip detection into the FIRST sample's stage-D so
  there is no re-sample either. Strip pairs are computed from
  `tchains`; the dense_boxes they promote could be requested in the same
  drill-down round that already runs (`want_dense`). Removes the whole
  attempt-0 restart. **Saves ~28 s/run** on bino.

**Effort.** Minimum: ~half day (guard the trace call + thread the cached
soup through the `continue`). Better: ~1 day (strip detection must move
before the drill-down survey consumes `crease_leaves`).

**Risk.** LOW-MEDIUM. The only real risk is under-delivering local-pitch
law in promoted strip bands (regresses joint snap tents — the v25/v26
class). Mitigation: keep re-tracing ONLY the strip-band creases (a tiny
subset), not the whole model. Gate on joint worst-depth (bino collar
site) + watertightness.

**Measure.** Bino + plinth (both jointed) with TIME=2; compare summed
wall, `worst` depth, open/nm counts vs HEAD. The strip promotion count
and joint-region depth are the referee.

### 1b. Light open-edge predictor (the big lever on the full assembly)

**What.** The retreat only needs `open_edges` to decide demotion. It
currently pays the FULL tail — 16 repair rounds + snap + decimate + seal
+ winding + recount — on every *doomed* attempt before throwing it away.
The open-edge count is a property of the extracted tet-complex boundary;
repair/snap/decimate target wart/chip DEPTH, not boundary topology.

**Design.** Split B+C into a **predictor** and a **finisher**:
- Predictor = insert samples + insert points + constraints + refinement
  + ONE extract pass + `recount_quality` (open/nm only). Skip the repair
  rounds 2..16, the snap pass, and all fix stages.
- Retreat loop runs the predictor per attempt. Once an attempt is clean
  (or is the ship-best), run the **finisher** (full repair tail) exactly
  ONCE on that attempt's DT.

**Expected savings.** Per doomed attempt: skip ~30 s repair + 5 s snap +
6 s fix + refinement heat ≈ **~40 s/attempt saved** at r2 bino scale;
the full assembly's 2 extra passes (~6 min today) would drop toward
~2–3 min. This is the single largest wall lever named in the brief.

**Effort.** MEDIUM (~1–1.5 days). Requires refactoring `mesh_impl` so
the repair loop is a callable stage that can run 0 times, keeping the DT
+ `prov` map + `snap_*` stashes alive between predictor and finisher.

**Risk.** MEDIUM. **The load-bearing assumption is that `open_edges`
after ONE extract equals `open_edges` after the full repair tail.** Must
be verified before building: repair inserts CAN split cells near the
boundary and change the geometric welded-id open count (the
`weld_ids`-based count is what the retreat reads, :8183). If repair
moves the count, the predictor needs the minimum repair rounds that
stabilize it (ledger says depth plateaus by round 2 — check whether
*open_edges* plateaus even earlier, likely round 0-1). Instrument first:
add an `open_edges`-per-repair-round print under TIME=2 and confirm on
plinth/bino/screws that the count is flat after round 0.

**Measure.** Add the per-round open-edge instrument (½ day, standalone).
Run plinth + screws (screws densifies, closest small proxy for the
retreat dynamics). Confirm the predictor's round-0 open count matches
the final count on every model. Then A/B a synthetic demote (force a
retreat via `STIBIUM_DMESH_SEP`) and compare wall.

### 1c. Incremental re-sample of changed leaves only

**What.** `delaunay_sample` re-runs the whole octree sample+bisect+QEF
(~10–20 s) each attempt though only demoted/promoted leaves change their
point emission. Cache the base sample set keyed by leaf; on retry, only
re-emit the leaves whose level changed.

**Expected savings.** ~8–18 s/attempt (the sample phase). Smaller than
1b but composes with it.

**Effort.** MEDIUM-HIGH. The sampler is a streaming octree pass; keying
per-leaf output and splicing changed leaves is a real refactor.

**Risk.** MEDIUM. THIN/weld state and the flood+cores dilation are
global-ish; a spliced set must reproduce the contiguity laws exactly.
Bit-identical verification (git-stash A/B) is mandatory.

**Verdict.** Do 1a + 1b first; 1c only if the sample phase is still
material after the predictor lands.

### NOT recommended: DT reuse / local re-mesh across attempts

Reusing attempt N's CCDT with localized vertex `remove()`+re-insert
(brief's option) is HIGH risk / HIGH effort: the CCDT constraint
bookkeeping around demoted leaves (which sit *near creases* by
construction — that's where open edges appear) is exactly the machinery
that is not safely mutable. Local patch re-meshing + stitching is worse.
The predictor (1b) captures most of the win without touching CCDT
internals. Park these.

---

## ROCK 2 — insert points (46.76 s)

883K surface points inserted one-by-one into the CCDT (:4262). Two costs:
the CGAL point-location+insert per point, and — CCDT-only — a
`nearest_vertex` coincidence query per point (:4281) on top of it.

### 2a. Replace the per-point nearest_vertex guard with a hash grid (near-term)

**What.** The coincidence guard (:4277-4305) exists to catch
grid-aligned points converging within `1e-3 sp` of an existing vertex
(reuse the vertex instead of minting a femto-segment). It calls
`dt.nearest_vertex(sp3, shint)` for EVERY surface point — a second walk
per insert. Replace with an O(1) spatial hash of already-placed vertices
at `1e-3 sp` cell size; only when the hash reports a candidate do we pay
a distance check. No DT query in the common (no-coincidence) case.

**Expected savings.** The guard is roughly half the per-point work on
grid-aligned models. Plausibly **10–18 s** of the 46.76 s. Cheapest win
on this rock.

**Effort.** LOW (~half day). A `std::unordered_map<cellkey, vertex>` fed
as points are placed; mirrors the existing `vgrid` bin-grid pattern used
15 lines down for through-vertex splits (:4351).

**Risk.** LOW. Pure optimization of an existing guard; output must stay
bit-identical (verify with [golden] + git-stash A/B). The 1e-3 sp weld
threshold is preserved exactly.

### 2b. CCDT range-insert investigation (near-term, may be a dead end)

**What.** The plain-DT path bulk-inserts via `dt.insert(begin,end)`
(:4218) — CGAL's internal spatial-sort + BRIO, far cheaper than the
CCDT one-by-one loop. Check whether the CCDT `_impl` exposes a range
insert that still lets us stamp the `CORNER` vertex type + info. If it
does, `insert samples` (28.94 s) and possibly `insert points` collapse
toward plain-DT bulk cost.

**Effort.** LOW to investigate (read CGAL 6.2 CCDT_3_impl headers),
unknown to implement.

**Risk.** LOW to try. Likely blocked: the CORNER tagging + info
assignment is why the code went one-by-one in the first place. Worth 2
hours to confirm/refute before committing to 2c.

### 2c. Parallel plain-DT + CCDT-only-for-creases (P6, the real prize)

**What.** The brief's P6. Insert the bulk sign-witnesses and non-crease
surface points into a **concurrent** `Delaunay_triangulation_3`
(`Parallel_tag` TDS + `Spatial_lock_grid_3` + range insert), and build
the CCDT only over the crease neighborhood. Confirmed feasible by CGAL
docs: `Parallel_tag` on the TDS, a lock grid sized to the region bbox,
range insert with automatic spatial sort, **requires TBB**.

**Feasibility check (done):** TBB is present system-wide
(`libtbb.so.12`) but NOT currently linked into the fab build — the CGAL
link line needs `CGAL::TBB_support` / `find_package(TBB)` added.
`Triangulation_vertex_base_with_info_3<int8_t,K>` composes with
`Parallel_tag`.

**Expected savings.** insert samples (28.94 s) + the bulk portion of
insert points → potentially 3–5× on the bulk phase with 15 lanes; the
crease neighborhood (small) stays sequential CCDT. Could take the two
insert phases from ~76 s toward ~20–30 s. Largest structural win on the
whole B+C.

**Effort.** HIGH / RESEARCH (days). The hard part is not the parallel DT
— it is the **hand-off**: the crease points and constraints must land in
a CCDT that shares vertex identity with the bulk DT. Two viable shapes:
(i) build the whole thing as CCDT but bulk-insert non-crease points via a
concurrent scratch DT then transfer — CGAL has no clean transfer, likely
a rebuild; (ii) restrict the constrained triangulation to a crease-band
sub-region and glue. Both are genuine research.

**Risk.** HIGH. Determinism (already ±2 tris under MT, Nate accepts),
the lock-grid granularity tuning, TBB as a new hard dependency, and the
identity hand-off. Prototype the plain concurrent DT in isolation first
(measure bulk-insert speedup on the 928K witness set) BEFORE attempting
the CCDT glue — if the isolated speedup disappoints, P6 is not worth the
glue risk.

**Sequencing.** 2a now (cheap, certain). 2b as a 2-hour probe. 2c only
after 2a/2b and only if the isolated concurrent-DT prototype shows the
expected bulk speedup.

---

## ROCK 3 — tracer (18 s), refinement (14.6 s), insert samples (28.9 s), extract (30.7 s)

### 3a. Chain-parallel tracer

**What.** `delaunay_trace` (:3124) loops over `npairs` min/max clauses,
each over its seed list. Parallelize across clauses (pairs).

**Data-hazard assessment (the brief's ask):**
- `tr` (CreaseTracer) carries per-pair mutable state (`tr.tp`, :3129) →
  each thread needs its OWN CreaseTracer copy (cheap: it holds deck/ctx
  pointers + scalars + the band map). The `ctx` (TapeCtx) is NOT
  thread-safe → one `tape_ctx_new` per thread (the established
  `eval_points_mt` pattern).
- `polys`/`closed`/`poly_pair` — appended. Give each thread local
  buffers, concatenate after the join.
- **`consumed[]` (crossing-suppression) is the one real hazard.** It is
  shared write state that prevents a seed from re-seeding another chain.
  But it is an *optimization, not correctness*: worst case a seed near
  two clauses' creases gets traced twice → a duplicate polyline. The
  within-clause dup guard (:3171, `poly_pair[c] != tp.clause` skip)
  already tolerates same-clause dups; cross-clause dups are new but rare
  (a seed sits on one crease). Handle by a **post-merge dedup**
  (distance-to-segment, the same 0.1 sp identity the dup guard uses) or
  by keeping `consumed[]` as an atomic-flag best-effort (races only
  cause redundant work, never corruption).
- **Junction sharing (:3315) is already a sequential post-pass** over
  all endpoints — leave it after the merge. No hazard.

**Expected savings.** Tracer 18 s → ~3–5 s with 15 lanes if clause work
is balanced. Clauses are unbalanced (a few big creases dominate), so
realistic ~3–4×. Combined with 1a (which removes the *second* 18 s trace
entirely), the tracer stops being a rock.

**Effort.** MEDIUM (~1 day). Per-thread tracer + buffers + merge dedup.

**Risk.** LOW-MEDIUM. MT non-determinism already accepted (memory:
±2 tris). The merge dedup must not eat real neighbor rings — reuse the
measured 0.1 sp point-to-segment identity (the additive-joint autopsy
radius), NOT a point-to-vertex radius (the buried ring-eater class).
Gate on chain count + `DUMP_CHAINS` overlay + bino joint depth.

### 3b. Incremental refinement (P4)

**What.** Each refinement round (:4467) sweeps ALL finite edges to find
i/o edges. After the first round, only edges near freshly-inserted
points can newly become i/o. Track fresh-insert cells and scan only
their incident edges.

**Expected savings.** On bino refinement is 1 round (14.58 s) so limited
here; on zeiss (more rounds) the full-sweep-×-N is the cost. Modest
model-dependent win.

**Effort.** MEDIUM. Needs a fresh-vertex → incident-edge frontier and
proof the incremental set is complete (a missed i/o edge = a chip).

**Risk.** MEDIUM. Correctness-sensitive (missed edges chord). Gate hard
on open/nm parity + worst depth vs full-sweep on all referees. Lower
priority than 3a/1b.

### 3c. Incremental extraction

**What.** The repair loop re-classifies ALL finite cells every round
(:4699) though only cells touched by the round's inserts change sign.
r0/r1/r2 signs+extract = 6.7/8.1/9.3 s and GROWING with cell count.
Restrict re-classification to cells adjacent to fresh inserts.

**Expected savings.** Rounds 1..N drop toward the fresh-cell fraction;
extract+repair 30.7 s → plausibly ~15 s. Composes with 1b (predictor
runs one extract; finisher benefits from incremental rounds).

**Effort.** MEDIUM-HIGH. Same completeness-proof burden as 3b, plus the
DC/all-surface oracle bookkeeping is per-cell.

**Risk.** MEDIUM-HIGH. Defer until 1b lands (1b removes most repair
rounds from doomed attempts, which is the cheaper way to the same wall).

---

## Sequencing & risk summary

| # | proposal | expected saving | effort | risk |
|---|----------|-----------------|--------|------|
| 1a | kill strips re-trace/re-sample | ~28 s/run, every jointed model | ½–1 d | low-med |
| 2a | hash-grid coincidence guard | 10–18 s (insert points) | ½ d | low |
| 1b | light open-edge predictor | ~40 s/doomed attempt; ~3 min on full assembly | 1–1.5 d | med (verify open-count stability first) |
| 3a | chain-parallel tracer | 18 s → ~4 s | 1 d | low-med |
| 2b | CCDT range-insert probe | maybe collapses insert samples | 2 h probe | low |
| 3c | incremental extraction | 30 s → ~15 s | 1.5–2 d | med-high |
| 3b | incremental refinement | model-dependent | 1 d | med |
| 1c | incremental re-sample | 8–18 s/attempt | 2 d | med |
| 2c | P6 parallel-split DT | insert 76 s → ~25 s | days/research | high |

**Recommended order:** 2a (cheap, certain, bit-identical) → 1a (biggest
always-on tax, under-reported) → the 1b open-count instrument (½ day,
standalone, de-risks the predictor) → 1b predictor → 3a tracer. Then
re-profile; 2b probe and 2c prototype gate the parallel-insert research;
3b/3c/1c are second-wave once the predictor reshapes the repair cost.

**Referee gates (every proposal):** geometric open-edge count (welded
ids) = 0, nm within noise, `worst` depth held at each referee's site,
and — for 1a/3a — the joint-region depth (bino collar, plinth teeth) and
`DUMP_CHAINS` overlay. Bit-identical [golden] + git-stash A/B for 2a and
2b. All A/B-able on plinth/screws/bino in minutes; NEVER the full zeiss.

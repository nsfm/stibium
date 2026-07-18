# Opus correctness review round 2 - 2026-07-18 (read-only; nothing fixed)

Scope: the ~15 commits since `5ffae8db` (sample thinning) through `b5a8c32b`
(slicer hygiene). Diff base `5ffae8db^`. Targets: thinning, snap referees,
DEGEN/PANCAKE, decimate independent-set, P5 threading, seal/winding, progress
sink, export integration. Verified with instruments where cheap
(SbFabTest 627,666 green; plinth A/B export).

## Findings (ranked)

### 1. Dialog knobs SNAP / DECIMATE / STALL are read through `static` getenv - a 2nd interactive export in one session silently reuses the FIRST export's values (MEDIUM; certain; moderately reachable)
`export_mesh.cpp:187-195` `qputenv`s the advanced dialog knobs into the process
env before each `delaunay_mesh` call. But three of them are read in the mesher
through function-local `static const char*` caches that initialize exactly once:
`STIBIUM_DMESH_SNAP` (delaunay.cpp:6194), `STIBIUM_DMESH_STALL` (5359),
`STIBIUM_DMESH_DECIMATE` (8136, 8257). After the first export those pointers are
pinned; `qputenv`/`setenv` on the second export does not update them.
`STIBIUM_DMESH_AUTODENSE` / `AUTODENSE_MAX` (901, 833) are read fresh each call,
so they DO update - the result is a confusing *partial* application where some
advanced settings take effect on re-export and others don't.
Verified the C mechanism directly (scratchpad `envtest.cpp`): a `static` getenv
returns `1` even after `setenv("FOO","0")`; a fresh `getenv` returns `0`. The
first export is always correct; the headless/harness path never `qputenv`s, so
referee reproductions are unaffected - this is GUI-only.
- **Scenario**: open a model, export with Decimate=on, then export again in the
  same session with Decimate=off (or Snap off, or Stall patience changed). The
  second STL ships with the first export's decimation/snap/stall behaviour.
- **Test**: two sequential dialog exports of screws, toggling Decimate between
  them; the second export's triangle count matches the first's decimated count
  instead of the raw count. Or add a one-line `fprintf` at 8136/6194 printing
  the read value and watch it not change on the 2nd export.

### 2. `STIBIUM_DMESH_TRACE=0` no longer disables tracing - a progress-stage insert detached `delaunay_trace` from its guard (MEDIUM; certain; default output unaffected)
delaunay.cpp:8028-8031:
```
static const char* tr_env = getenv("STIBIUM_DMESH_TRACE");
if (!tr_env || atoi(tr_env) != 0)
    prog_stage(3);
    delaunay_trace(deck, r, &soup, halt);   // NOT under the if
```
`git blame`: `delaunay_trace` was originally the guarded statement
(`ade2b437a`); commit `7ef2762a` ("Progress reporting for app integration")
inserted `prog_stage(3);` between the `if` and the call, so the brace-less `if`
now guards only the progress update and tracing runs unconditionally. Default
behaviour (TRACE unset -> tracing on) is unchanged and correct, so no mesh is
corrupted. The damage is: (a) the disable knob is dead - any A/B or perf
measurement taken believing tracing was off is silently wrong; (b) with TRACE=0
the "tracing features" progress stage is skipped while the time is still spent
there, so the bar mis-maps. The misleading indentation will also trap the next
reader.
- **Test**: `STIBIUM_DMESH_TRACE=0 STIBIUM_DMESH_TIME=1 ... --export` on bino;
  the `tracer` phase line still shows its ~10-18 s instead of ~0.

### 3. `weld_slivers` keeps the decide-on-original / apply-at-pass-end double-remap hazard that `decimate_flats` was fixed for (LOW-MEDIUM; hazard certain, benign in current measurements)
The fold cure (`5422f...`, independent-set rule) was applied to `decimate_flats`
only. Its sibling `weld_slivers` (delaunay.cpp:7289-7377) decides collapses on
original positions and applies the `remap[]` at pass end (7367), but marks only
`touched[drop]` and `touched[keep]` (7361-7362) - NOT the whole fan like
decimate does (7193-7195). Reachable case: triangle `t = {drop1, drop2, x}`.
`drop1` collapses first from some other incident triangle (marking only
drop1/keep1); `drop2` is left un-touched, so a different incident triangle can
collapse `drop2` too. Both are accepted; at apply time `t` becomes
`{keep1, keep2, x}` - a configuration NEITHER collapse's orientation guard
(7332-7356) nor `link_ok` evaluated (each guard substituted only its own drop and
used the other endpoint's original position). This is exactly the fold class
decimate had. It is empirically benign today (MESH-WAR: "weld solo 0 fights") -
the motion is a razor short edge and the guards catch the common cases - but the
guard-blindness is real and could mint a non-manifold edge or a flipped sub-bar
sliver on clustered-razor geometry; `recount_quality` runs last so it would at
least be counted, not hidden.
- **Test**: WELDSLIV on, FLIPSLIV=DECFLATS=0, independent nm/open sweep
  (`tools/openedges_exact.py`) on a razor-dense model; compare against a
  full-fan-claim variant of weld_slivers.

### 4. Cancelled Stibnite export still writes the output file (LOW; certain)
`export_mesh.cpp:307` guards `simplifyMesh` with `!halt`, but the
`save_3mf_indexed` / `save_stl_indexed` calls at 322-328 are not halt-guarded.
On dialog cancel, `delaunay_mesh` returns early on `halt`, `_stats` stays empty
(so the report dialog is correctly suppressed), but `async()` still writes the
partial/empty `verts`+`indices` to the chosen filename. The classic path has the
same shape, so this is not new, but the Stibnite branch inherits it.
- **Test**: start a long bino export, cancel; the target file exists with 0/near-0
  triangles.

### 5. `seal_seams` + `repair_winding` are gated behind the `STIBIUM_DMESH_DECIMATE` master switch (LOW; knob coupling)
delaunay.cpp:8137 and 8258 wrap the entire tail block - weld/flip/decimate AND
seal/winding/recount - in `if (!de || atoi(de) != 0)`. So
`STIBIUM_DMESH_DECIMATE=0` silently also disables slicer-hygiene seal, winding,
and the post-mutation recount, even though those have their own SEAL / WINDING
knobs. Defaults are on so nothing ships wrong, but the coupling is surprising and
means "decimate off" ships an un-sealed, un-recounted (stale-ledger) mesh.

### 6. Trivia
- `DMeshProgress` is a process-global singleton never reset at mesher entry
  (first `prog_stage` is `prog_stage(1)` at 2385). The export mirror thread
  spins up (export_mesh.cpp:214) before `delaunay_mesh`, so a 2nd export briefly
  reads the prior run's `overall=1.0` and flashes 100% before stage 1 resets it.
  Cosmetic.
- `run()`'s `!halt` guard on the report dialog (export_mesh.cpp:108) is dead:
  `runAsync` resets `halt=0` (export_worker.cpp:76) before returning, so `halt`
  is always 0 here. Harmless - the `_stats.isEmpty()` check already suppresses
  the dialog on cancel.

## Verified sound (with evidence)

- **P5 threading determinism / bit-identity (target 5)**: plinth exported at
  `STIBIUM_DMESH_THREADS=1` and `=8` produced **byte-identical STL** (md5
  `91c7034a...`, 13466 tris, 0 open/0 nm both). The parallel path is exercised
  (collect_tasks depth-2 -> up to 64 tasks; thousands of edges > bisect's 4096
  gate). This covers: `eval_points_mt` MIN_VOLUME alignment (tape eval is
  batch-position-invariant in practice), `bisect_edges` range-parallel (ranges
  are NOT MV-aligned, but the result is still bit-identical - the "each edge
  depends on nothing but itself" reasoning is incomplete as stated but the
  conclusion holds empirically), parallel `descend` + fixed-order
  `merge_collector`, and the crossing-index remap (h_of inversion is total
  because add_edge keeps edge_index dense; feature-cell pts remap correctly).
  Determinism is per-thread-count: a fixed THREADS gives a fixed mesh.
- **Sample thinning (target 1)**: THIN=-1 is a clean bypass
  (`shell < 0 || keep[p] || vals[p]==0` short-circuits before any `keep[]`
  access). Every crossing edge's two endpoints are always kept (band-marking
  loop at 604-605 uses the same forward-neighbour stencil as the add_edge loop
  at 685-688, and add_edge only fires on a sign change, which only exists in a
  mixed block where the band loop ran). Keep is a UNION across blocks via
  add_sample's coord_hash dedup; the parallel path preserves the same set
  (per-block keep is identical, merge dedups shared faces). Shell dilation
  truncates at block boundaries but the boundary crossing is re-seen and
  re-dilated by the neighbour block, so the band is never under-witnessed.
- **Snap tent referees + wave loop (target 2)**: stale deferral is correct - a
  `pem` entry whose triangle was rewritten this wave is skipped as `stale` and
  the edge defers (`continue`, not `done`) rather than concluding; new tent
  triangles absent from this wave's `pem` are covered because the touched
  original still sits in `pem` for the moved edge and trips `stale`. Termination
  is guaranteed: `progress` is set only when a tent is actually applied and each
  applied snap sets `done[s]=1` permanently, so not-done strictly decreases each
  progressing wave. Damage and churn referees gate both the polyline and the
  surf-projection tents. C2/D2 point into `out->verts` but are consumed before
  the `verts.push_back` that could reallocate; A2/B2 are copied locals.
- **DEGEN/PANCAKE (target 3)**: sign convention confirmed (`info==-1` inside,
  `+1`/infinite outside, 4806/5115), so `want = votes>0 ? 1 : -1` breaks ties to
  -1 = INSIDE exactly as the comment claims. Fixpoint is Gauss-Seidel over a
  deterministic cell order, capped at 8 rounds (guaranteed termination). Crucially
  the pass only RE-LABELS cells; the extracted surface is the boundary of a
  cell-set, so it stays closed/watertight regardless of labels - DEGEN cannot
  open a hole. Pancake discriminator uses the field gradient (7-pt stencil,
  eval_points_mt) to keep DC wedges and disenfranchise flats; sound.
- **decimate_flats independent-set (target 4): the fan-claim IS complete.** A
  triangle is rewritten only if one of its vertices was collapsed, i.e. it lies
  in that vertex's fan; a triangle shared by two collapse fans would have both
  collapse sources as its vertices, but the first collapse marks every fan vertex
  (including the second source) touched, so the second collapse is refused. Hence
  no triangle is ever double-remapped, and the orientation guard + `link_ok`
  (evaluated on pass-start `inc`, valid because accepted collapses are
  triangle-disjoint) see exactly the applied configuration. No triangle outside
  both fans is ever touched.
- **seal_seams budget math (target 6)**: `par[k]` = net signed winding over
  coincident copies of a vertex triple; `left=par` is the decrementing budget.
  The loop keeps a majority-orientation copy while `budget!=0`, decrementing by
  `p` - keeping exactly `|net|` copies (there are always >= |net| majority
  copies). The REMOVED set therefore has net winding 0, and since all coincident
  copies share the same three edges, each of those edges' directed balance
  changes by 0 - so annihilation provably preserves geometric edge balance
  (watertightness). `repair_winding` flood-fills coherent orientation per
  orientable component across clean 2-facet edges (pinch/boundary edges are
  barriers), then flips a whole component only on a decisive field vote
  (>=8 tris, |vote|>0.5*mag); a whole-component flip preserves edge multiplicity,
  so open/nm counts are unchanged and `recount_quality` (which runs AFTER seal
  and winding, in both the convergent 8133 and exhaustion 8254 paths) reports the
  shipped mesh - closing round-1 finding #1.
- **Progress sink thread-safety (target 7)**: all `DMeshProgress` fields and the
  export mirror's `progress_*` are `std::atomic`; `prog_*` is called only from
  the main mesher thread (eval_points_mt workers touch only their eval range), so
  the only cross-thread reader is the UI/mirror and relaxed atomics are adequate.
  The prog_emit throttle race is benign (at worst an extra PROGRESS line) and
  documented. No real bug.
- **Export mirror lifecycle (target 8)**: `mesh_done=true; mirror.join()` runs
  unconditionally after the synchronous mesher returns (export_mesh.cpp:251-252),
  on every path including cancel; `runAsync` waits on the future before returning,
  establishing happens-before for the GUI-thread `_stats` read.

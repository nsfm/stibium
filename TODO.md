# TODO — the campaign

Shipped work moves out to [CHANGELOG.md](CHANGELOG.md); this file is only
what's still ahead.

## Tier 1

- **Redistance primitive + true morphology (printability filter,
  min-wall check).** Two nested offsets do NOT implement open/close:
  our offset is field-minus-constant, so they cancel textually
  (discovered 2026-07-13 when a min-wall check built that way passed
  everything; that node was deleted before shipping - Checks has no
  min-wall until this lands). The real fix: `redistance(shape, res)` samples
  the sign on a grid, runs an exact euclidean distance transform
  (Felzenszwalb, O(n)), and registers the result as a grid field -
  the OP_GRID machinery from mesh import already handles the rest.
  With true distances, open/close become exact-to-voxel, unlocking
  Check: Min Wall (erode by t/2, see what vanished) and the "show me
  what my printer will actually produce" preview. Also generally
  fixes offset/shell on non-exact fields.

- **CPU SIMD tile viewport renderer (MPR ideas, now unblocked).**
  Shipped 2026-07-13: shortened-tape evaluation (the libfive
  keystone) landed - immutable shared tapes replaced disable_nodes /
  clone_tree throughout the kernel (see doc/TAPE-DESIGN.md +
  CHANGELOG). The next rung is the tile renderer it was for: tile
  the viewport, interval-cull tiles, specialize the tape per
  ambiguous tile, evaluate survivors with the array evaluator
  (thread pool stands in for the GPU grid; SIMD lanes for warp
  lanes). Round 2 (2026-07-13) already banked the prerequisites:
  register allocation (bounded workspaces), spare-tape freelist,
  per-shape deck caching, the shrinkage curve (saturates ~depth 6 on
  the merged Zeiss - tune tile size against it), and an opt-in
  affine-collapse pass to enable when this lands.
  SHIPPED 2026-07-13, three pieces: (a) the shaded pass pushes per
  64px tile against the depth buffer's z-range, so gradients run on
  pruned tapes (~30% off the shading portion); (b)
  `tape_eval_i_batch` - batched interval evaluation, 64 boxes per
  pass down the tape, bounds exactly equal to scalar (hot ops are
  elementwise math_i mirrors, branchy ops delegate per lane), living
  in the existing r-rows (64 lo + 64 hi = MIN_VOLUME floats, zero
  new allocation); (c) the depth pass fans ambiguous regions into 64
  children and batch-classifies them - decided tiles fill/vanish
  without re-walking the tape, only ambiguous ones recurse.  Default
  ON (STIBIUM_TILE_RENDER=0 restores bisection, =N tunes): merged
  Zeiss at 2048px renders ~20% faster, win grows with resolution,
  pixel-identical everywhere.  Note the fan-out was a LOSS until the
  batch evaluator landed - serial re-specialization was already
  near-optimal; amortized traversal is what pays.  Still open here:
  per-tile worker scheduling (current xy-chunking is coarser than
  MPR's), and the GPU compute-shader rung (Tier 3) which reuses all
  of this machinery.

- **Parser hardening (upstream #198).** A malformed math expression
  hits a lemon assert and aborts the whole app. We OWN this parser
  (extended it twice) - malformed input must error, not core-dump.
  Cheapest safety win on the board, and the gate on trusting foreign
  imports and agent-authored strings.

## Tier 2

- **Mesh import, round two.** The landed core is STL-only, dense
  grids, blocking sampling. Follow-ups in rough priority order:
  3MF import (needs zip reading - reuse zlib, small parser); an
  import _dialog_ (file browser + resolution suggestion from mesh
  bbox/feature stats, memory estimate readout - today you type into
  datum fields); progress reporting during sampling (big scans block
  the graph thread silently); narrow-band grids for memory (dense
  512³ = 512 MB; band + far-field fallback could be ~10×); script
  print()/advisory surfacing in headless verbs (--validate swallows
  script stdout, so the ouroboros advisory only shows in the GUI);
  a File > "copy mesh into project" helper for the unsaved-project
  nag; consider exposing `fab.shapes.import_mesh` bounds padding.

- Shipped 2026-07-13: light-sweep GIF, startup splash, Ctrl-C
  graceful close, and the animation export dialog (size/frames/AA)
  all landed - see CHANGELOG.

- **Render polish, next round.** Supersampled AA (default 2x, --aa)
  and turntable/wiggle/stereo verbs shipped 2026-07-13; remaining:
  edge-aware AA refinement (only where depth/normals step),
  backgrounds beyond transparent/dark (solid color picker,
  gradients), exact-viewport-crop (match pan/zoom framing), GUI menu
  entries for the animation verbs, section-sweep animations,
  geometric-diff heatmaps, batch gallery renders.

- **Deep-zoom render quality (micron-scale features as a feature).**
  Now that extreme zoom no longer clips (28d365ad), surfaces get
  visually harsh up close: capped depth-voxel budget quantizes the
  z-buffer, and gradient evaluation epsilons aren't scaled to the
  zoom. Approach thoughtfully: adaptive z-budget near surfaces,
  eval epsilon proportional to feature scale, possibly float64 field
  eval at extreme magnification. F-rep has no meshes to run out of -
  if deep zoom is nailed, modeling micron-scale detail in
  large-format 2D work (photolithography) becomes a genuine
  differentiator no mesh-based tool can follow.

- **Print-centric checks.** Minimum-wall-thickness detection from the
  field itself (walls thinner than nozzle flagged before slicing;
  waits on the redistance primitive to be exact); overhang-angle
  visualization is plausible too (gradient vs build direction). Ends
  guess-and-check on functional parts.

- **Audit the smooth-CSG lineup.** Smooth Union / Smooth Difference
  nodes now exist (log-sum-exp `blend_expt_unit` / `blend_difference`
  from the raid) alongside blend/shell/offset/morph. Remaining:
  double-check blend's falloff quality against the standard
  polynomial/exponential smin family, add a smooth-intersection
  variant if it's the gap, and make sure all of them land in the
  wiki showcase.

- **Viewport measurement probe.** Click two points: distance, plus
  the field value under the cursor (= exact distance to surface,
  it's an SDF!). The field is a built-in ruler; expose it.

- **Multi-shape export.** One click → N files (print plates, assemblies).
  The export hook currently hard-rejects multiple export tasks per script.

- **Node editor QoL, round two.** Fuzzy-search add menu shipped
  (Tab / double-click palette). Still wanted: a minimap for big
  graphs, and canvas annotations (sticky notes + named zones behind
  nodes, persisted in the .sb JSON like node positions).

- **Expose rim/AO strengths in the render UI.** The key-light
  direction already has a draggable trackball gizmo (shipped); rim
  and ambient-occlusion strengths are still config-file-only. Give
  them the same real UI.

- **Smart resolution suggestion.** The dialogs now floor at sane
  defaults (7 for mesh, 60 for vector/raster) and show the implied
  minimum feature size (~2 voxels); the real version analyzes the
  model - find the smallest feature via interval probing and suggest
  a resolution that resolves it at the target tolerance (0.1mm print
  vs micron litho).

- **Enum/choice datums with dropdowns.** Nodes that take one of a known
  set (dovetail male/female, teardrop axis X/Y/Z, screw-size presets)
  should render a combo box instead of a free-text field - e.g.
  `input('kind', str, 'male', options=['male','female'])` flowing through
  to a QComboBox in the datum row.

- **UI glow-up, remaining milestones.** Tiers A+B shipped (palette,
  dark theme, wire/port/node facelift, enhanced render mode, zoom
  LOD, fuzzy palette). Full prioritized plan in
  [doc/UI-MODERNIZATION.md](doc/UI-MODERNIZATION.md): viewport shader
  modernization, icons/fonts, and the Qt6 render-path polish - with
  the QML question already settled (framework verdict: keep Qt,
  reject QML rewrite; base Qt6 port itself is done).

- **Autosave v2.** Rotating timestamped backups instead of
  overwrite-in-place; crash-recovery prompt on next launch; interval and
  retention in a real preferences dialog (the fork's founding feature,
  un-hacked).

- **Preferences dialog.** Autosave settings, default export deviation,
  default resolution heuristic, viewport colors. The app currently has
  essentially no user configuration surface.

- **Reusable node groups / user parts library.** Save a subgraph cluster
  (e.g. a measurement block + subsystem) as a named, reusable node in a
  user library folder. Foundation for a community parts ecosystem, and
  the highest-demand feature theme on the upstream tracker
  (#217/#68/#25/#127/#22).

- **Geometric diff, follow-ups.** The `--diff` verb shipped (JSON
  volumes + gray/red/green composite render). Still wanted: a
  sub-surface heatmap of |f_a - f_b| magnitude, and a `--diff`
  exit-code mode for CI gates ("fail if anything outside this region
  changed").

- **Assertion nodes, follow-ups.** The Checks category shipped (Fits
  Box / Volume / Clearance over `fab.shapes.measure()`; red in canvas,
  --validate exits nonzero, wired into CI via showcase_gear's spec).
  Min Wall waits on the redistance primitive (Tier 1). Future:
  Check: COM ("stands upright"), per-check resolution presets, a
  summary panel.

- **Projected footprint node (project_z).** F(x,y) = min over z of
  f(x,y,z): the true shadow of a 3D part as a first-class 2D shape -
  straight into the SVG/DXF pipeline for baseplates, gaskets, and
  laser-cut cradles matched to 3D parts.

- **Color propagation fix.** Color currently vanishes when a colored shape
  merges with any uncolored node. Rule: CSG results inherit color from
  colored operand(s) (union of colored+uncolored keeps the color;
  colored+colored could blend or keep-first). Prereq for color export.

- **Color/multi-material 3MF export.** The color nodes already exist and
  die in the viewport; 3MF carries material/color regions and
  multi-material printers are mainstream. ~70% built and doesn't know it.

- **Analytics panel, round two.** The engine, `--analyze` verb, and a
  live Render ▸ Analytics overlay (volume / COM crosshair / tight
  bounds) all shipped. Still wanted: mass-per-material once color
  propagation lands, a stands-upright check, projected footprint,
  and octree interval pruning to make integration near-instant (the
  contour tracer's trick, in 3D).

- **Procedural noise opcodes.** Perlin/simplex in the C evaluator →
  knurling, stipple, woodgrain-as-a-node (with seed input so "every
  render unique" becomes reproducible-unique). Tactile surfaces on
  functional prints.

- **Lithophane / image-displacement node.** Image import's sibling:
  thickness-modulated surface from a photo. (The archive is full of
  camera and darkroom gear; this audience overlap is not a coincidence.)

- **Slicer modifier-volume export.** Export designated shapes as 3MF
  modifier meshes (PrusaSlicer supports them): "this region 100% infill,
  this one 0.1mm layers" - print-tuning designed in the graph.

- **Inline wire math / cross-node datum references** (upstream
  #26/#66/#24): "clearance = other_node.dim - 0.2". Core parametric
  workflow, recurring ask.

## Tier 3

- **GPU field evaluation for the viewport.** CPU raymarching holds back
  complex models (Zeiss ID02: dozens of assemblies). Fidget proves the
  JIT/wide-evaluation approach; upstream has an abandoned `gl-render`
  branch to study. Months, not days - but it makes the tool feel current.
  The shortened-tape prerequisite shipped 2026-07-13; the CPU SIMD tile
  viewport (now Tier 1) is the stepping stone.
- **Adaptive meshing (libfive-style manifold dual contouring).** Fewer,
  feature-aligned triangles at the source. MPL-2.0 core, license-safe to
  adapt. The meshoptimizer post-pass covers most of the value meanwhile.
- **Image import node.** PNG/heightmap → interpolated field (inverse of
  the existing heightmap export). Engraving, litho, textures.
- **Exact layer preview / field-to-toolpath.** A slicer computes planar
  contours per layer; the field gives exact contours for free. Start with
  "show me layer N" in the viewport; unhinged endgame is math → gcode
  with no mesh in between.
- **Proven clearances (interval-verified fit).** eval_i gives
  guaranteed bounds, not samples: subdivide until the interval proves
  min(gap) ≥ tolerance over the whole domain - a mathematical proof
  that parts mate, that threads clear, that nothing collides. No
  sampling tool can promise this; interval arithmetic can. CAD with
  certificates.
- **Shape morphing over time.** morph.node exists for static blends;
  drive its weight from the kinematic-scrubber timeline for tweened
  geometry animation, and print physical morph sequences.
- **Differentiable-CAD optimization.** Fields are closed-form and exactly
  differentiable: "minimize material s.t. wall ≥ 0.8mm", "solve for the
  parameter where these parts stop colliding." Gradient descent over the
  graph's free datums.
- **Kinematic scrubber.** A driver datum on a timeline slider - assemblies
  articulate in the viewport (the Zeiss focus mechanism actually racking).
  Also: motion GIFs for the gallery.
- **WASM web viewer.** Field evaluation compiled to WASM (Fidget proves
  it) → the GitHub gallery becomes interactive: visitors orbit models in
  the browser instead of squinting at JPGs.
- **Convex hull** (upstream #79/#134). Recurring OpenSCAD-migrant ask,
  but no exact closed-form in f-rep; genuinely a research item.

_The thesis binding all of these: everyone else converts to triangles and
then fights the triangles. Stay math longer - slice, optimize, texture,
weigh, and animate the field, upstream of any mesh._

## Agent-friendly modeling (LLM/automation as a first-class user)

Antimony already speaks Python-into-math; make the loop closable without
a GUI so agentic tools can contribute to modeling cleanly. Full design
plan in [doc/AGENT-SURFACE.md](doc/AGENT-SURFACE.md). Shipped so far: the
read/verify loop (headless verbs, `--describe-nodes` machine-readable
node reference, structured errors on stderr/exit) and deterministic
protocol-7 serialization. Remaining write surface:

- **Programmatic .sb authoring.** Deterministic serialization is done
  (protocol 7, byte-identical round-trips, meaningful diffs). Still
  wanted: document the schema + connection encoding and provide a
  tiny writer library so agents can emit .sb files directly.
- **`fab` as an installable pure-Python package.** The shape library
  already imports standalone with a small Shape/Transform shim (proven
  during chamfer testing). Publish it so scripts/agents can compose
  shapes and emit .sb or math strings outside the app.
- **Live-session control (MCP server, stretch).** The file-watcher
  reload covers agent-edits-on-disk; a socket API into the running
  session (load/save graph, list/add/edit nodes, read errors, eval
  bounds, render viewport) is the full answer - any agent-capable
  editor could then drive a live Stibium session interactively.

## Distribution / project health

- freedesktop thumbnailer: .thumbnailer file + shared-mime-info XML
  for .sb in deploy/ - file managers preview models via `--render`
  (the headless renderer that powers it is shipped). Cheap and
  delightful.
- Examples refresh, continued: showcase_gear.sb landed (Parts +
  Repeat + CSG, protocol-7-native, in CI). Still wanted: ISO-thread
  assembly, chamfered-CSG piece, litho-style 2D mask for the vector
  pipeline.
- GitHub Actions CI (build on push), AppImage/Flatpak releases - "runs on
  other people's machines" is a feature. Windows build + package
  managers are the single biggest raw-demand cluster upstream
  (#71 + 10 others); distribution, not code.
- Getting-started wiki (GitHub Pages); node reference partially generated
  from the self-describing .node files (the `--describe-nodes` JSON is
  the seed). Nothing comparable exists.
- Batch shrink tool for existing STL archives (prototype works; decide
  home + overwrite-vs-suffix policy).
- Port more libfive stdlib (MPL-2.0-compatible), continued: gyroid,
  Schwarz P/D, rounded primitives, elongate, bend, twirl, half_space
  all landed; still wanted are more TPMS variants and text-on-path for
  the vector font.

## Small bugs / cleanups

- **spur_gear.node UI loop fails under --validate**: "Cannot declare
  multiple UI elements on same line (without unique keys)" from the
  wireframe-in-a-loop. Check whether the GUI path hits it too (red
  script box?) and either key the loop calls or teach the hook about
  loops.

- **Root-cause the Python-error crash (the autosave origin story).** The
  live code-checker segfaulted randomly mid-typing; the workaround
  (fork-era, see `lib/graph/src/util.cpp` getPyError - the line-number
  extraction is commented out and pinned to 0) disables error line
  numbers entirely. Fix properly with defensive PyObject extraction so
  script errors report real line numbers again without the crash.

- **Infrequent crash when deleting groups of nodes at once.** Suspect
  the multi-delete undo path (`app/undo/undo_delete_multi.cpp`) or
  dangling proxies during batch removal. Needs a repro harness /
  ASan session like the detect-features hunt.

## Recon digests

Condensed pointers to the full recon writeups; only the still-actionable
items remain here (shipped actions - the .vm importer, the stdlib raid,
protocol 7, mesh import, license declaration - moved to CHANGELOG).

### Upstream issues + PRs (172 issues, 69 PRs scanned)

~70% was dead build/packaging noise. Live signal:

- **Crash/geometry audit targets** (repro files on the issues; none
  reproduce in daily use, so these are audit targets against our
  rewrites): #198 parser core-dump (promoted to Tier 1); #20/#29/#177
  interaction segfaults (handle drag, view rotate, wire dropped into
  empty canvas - canvas overhaul may have fixed; test); #174 union of
  exactly-coincident faces leaves a gap in exported meshes; #200
  subtracted extruded text flattens in STL export (probably the dead
  heightmap path); #135 datum formatLink assert on cross-type
  connections.
- **Feature demand not yet covered:** node-reuse cluster and graph
  organization (#223/#221/#226/#73) → the Reusable-node-groups Tier 2
  item; inline wire math (Tier 2); #165 partial torus/revolve wedge -
  now a one-node quick win since half_space landed (wedge
  intersection); #222 scalar math nodes (add/multiply/number) - the
  Math category still has none, confirm and add.
- **Fragile-area map for future upgrades:** Boost::Python ABI matching,
  lemon version sensitivity (yy_find_shift_action signature), Qt
  deprecations, C-locale-after-QApplication (we already do this).
- PR #228 autosave (unmerged) folds into the Autosave-v2 Tier 2 item.

### libfive — see [doc/LIBFIVE-RECON.md](doc/LIBFIVE-RECON.md)

Keeter's post-Antimony decade as a ranked steal-list. Headline: libfive
is ahead on the kernel, we're ahead on the product, and his own
retrospective says the node-graph tool is the thing he wished someone
else would build. Shipped from the raid: exact-SDF box/rect primitives,
log-sum-exp smooth blends + smooth difference (log/exp opcodes),
half_space, gyroid/TPMS, clearance, mirrors. Still open:

- **Shortened-tape evaluation** - SHIPPED 2026-07-13 (Tier 1 → done;
  doc/TAPE-DESIGN.md).
- **maybe_nan interval tracking** - SHIPPED 2026-07-13 (per-clause
  taint in the tape evaluator; found by the pruning fuzzer, fixed a
  visible showcase_gear artifact - see TAPE-DESIGN "Round 2").
- **affine-collapse pass** - SHIPPED 2026-07-13 as opt-in
  STIBIUM_AFFINE=1 (wall-neutral on exports; flip on when the tile
  renderer stresses base-tape evals). **interval-derived bounds** (M)
  still open.
- **Adaptive manifold DC mesher** (MPL 2.0, literal port permitted) -
  the marquee item, tracked in Tier 3.
- **Steal Keeter's Meshing Algorithm Idea.** Per his 2026-07-03 post
  literally titled "Please Steal my Meshing Algorithm Idea": decouple
  surface-point generation (QEF sampling) from manifold mesh
  construction (Delaunay tetrahedralization), sidestepping manifold
  DC's thin-feature and self-intersection warts. He asked. We're
  considering it. The audacity of it all.

### Foreign project import — see [doc/FOREIGN-IMPORT.md](doc/FOREIGN-IMPORT.md)

Tier A shipped (`--import-vm` converts Fidget tapes, deterministic + in
CI; prospero renders). Remaining build order:

- **[S-M] libfive `.io` import** via its own tree-dump API (their
  interpreter does the hard part - never parse Scheme).
- **[M] straight-line transpile to parametric script nodes** once the
  stdlib raid fully lands, so imports come in editable rather than as
  one frozen field node.

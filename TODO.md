# TODO — the campaign

Running list of upgrades for this fork. If it goes far enough: upstream the
good parts to Keeter, or first-class the fork as **Stibium** (the .sb files
were named for the element symbol all along).

## Tier 1 — would die on these hills

- **DXF export for 2D shapes.** SVG shipped (see Done); DXF is the
  remaining format for laser-cutter toolchains — LWPOLYLINE entities
  from the same `contour_field` output.
- **Cross-section preview.** Draggable slice plane in the viewport showing
  part interiors (walls, voids, clearances). Nearly free in f-rep — render
  the field on one plane. Formalizes the accidental z-clip slicing the 2D
  view path already does (see task.cpp screen-space z clipping).
- **Mesh import as a distance field.** STL/3MF in → sampled voxel SDF node
  (trilinear interpolation). Design around existing parts: PCB models,
  scans, vendor STEP-derived meshes. Upstream's most-requested feature
  (mkeeter/antimony#153), never landed.

## Tier 2 — strong wants

- **Print-centric checks.** Minimum-wall-thickness detection from the
  field itself (walls thinner than nozzle flagged before slicing);
  overhang-angle visualization is plausible too (gradient vs build
  direction). Ends guess-and-check on functional parts.
- **Multi-shape export.** One click → N files (print plates, assemblies).
  The export hook currently hard-rejects multiple export tasks per script.
- **Node editor QoL.** Fuzzy-search add menu (type "cyl"), minimap for big
  graphs, canvas annotations (sticky notes + named zones behind nodes,
  persisted in the .sb JSON like node positions).
- **Floating node labels for the zoomed-out canvas.** Card labels
  currently scale with the scene; replace/augment with labels rendered
  at constant screen size floating near their nodes, nudging each other
  to avoid overlap, and prioritizing custom-named nodes over
  default-named ones (default names match `[a-z]\d+`). "Web-map style"
  chrome per Nate. Complements the LOD cards.
- **Light-direction UI for the Enhanced render mode.** The key light is
  configurable today (`render/key_light = "x,y,z"` in the config file);
  give it real UI - a drag gizmo in the viewport or at least a dialog.
  Consider exposing rim/AO strengths the same way.
- **Recent-files menu.** Config now remembers the last directory; grow
  it into File > Recent with a small MRU list.
- **Enum/choice datums with dropdowns.** Nodes that take one of a known
  set (dovetail male/female, teardrop axis X/Y/Z, screw-size presets)
  should render a combo box instead of a free-text field — e.g.
  `input('kind', str, 'male', options=['male','female'])` flowing through
  to a QComboBox in the datum row.
- **UI glow-up.** Full prioritized plan in
  [doc/UI-MODERNIZATION.md](doc/UI-MODERNIZATION.md): palette refresh,
  app-wide dark theme, wire/port/node facelift, viewport shader
  modernization, zoom LOD, icons/fonts — with the Qt6 port as its own
  later milestone (framework verdict: keep Qt, reject QML rewrite).
- **Autosave v2.** Rotating timestamped backups instead of
  overwrite-in-place; crash-recovery prompt on next launch; interval and
  retention in a real preferences dialog (the fork's founding feature,
  un-hacked).
- **Preferences dialog.** Autosave settings, default export deviation,
  default resolution heuristic, viewport colors. The app currently has
  essentially no user configuration surface.
- **Reusable node groups / user parts library.** Save a subgraph cluster
  (e.g. a measurement block + subsystem) as a named, reusable node in a
  user library folder. Foundation for a community parts ecosystem.

## Tier 2.5 — adopted from dreaming mode (concrete, just not scheduled)

- **Color propagation fix.** Color currently vanishes when a colored shape
  merges with any uncolored node. Rule: CSG results inherit color from
  colored operand(s) (union of colored+uncolored keeps the color;
  colored+colored could blend or keep-first). Prereq for color export.
- **Color/multi-material 3MF export.** The color nodes already exist and
  die in the viewport; 3MF carries material/color regions and
  multi-material printers are mainstream. ~70% built and doesn't know it.
- **Parametric ISO thread nodes.** Helix = twisted profile; f-rep handles
  it natively. Internal/external, standard pitches, printable clearance
  presets. Every functional-parts tool hits this wall; Stibium doesn't.
- **Analytics panel.** Field integrals via Monte Carlo: volume, mass per
  material, center of gravity, stands-upright check, bounding box &
  projected footprint. Ends guess-and-check counterweighting.
- **Procedural noise opcodes.** Perlin/simplex in the C evaluator →
  knurling, stipple, woodgrain-as-a-node (with seed input so "every
  render unique" becomes reproducible-unique). Tactile surfaces on
  functional prints.
- **Lithophane / image-displacement node.** Image import's sibling:
  thickness-modulated surface from a photo. (The archive is full of
  camera and darkroom gear; this audience overlap is not a coincidence.)
- **Slicer modifier-volume export.** Export designated shapes as 3MF
  modifier meshes (PrusaSlicer supports them): "this region 100% infill,
  this one 0.1mm layers" — print-tuning designed in the graph.

## Tier 3 — moonshots

- **GPU field evaluation for the viewport.** CPU raymarching holds back
  complex models (Zeiss ID02: dozens of assemblies). Fidget proves the
  JIT/wide-evaluation approach; upstream has an abandoned `gl-render`
  branch to study. Months, not days — but it makes the tool feel current.
- **Adaptive meshing (libfive-style manifold dual contouring).** Fewer,
  feature-aligned triangles at the source. MPL-2.0 core, license-safe to
  adapt. The meshoptimizer post-pass covers most of the value meanwhile.
- **Image import node.** PNG/heightmap → interpolated field (inverse of
  the existing heightmap export). Engraving, litho, textures.
- **Exact layer preview / field-to-toolpath.** A slicer computes planar
  contours per layer; the field gives exact contours for free. Start with
  "show me layer N" in the viewport; unhinged endgame is math → gcode
  with no mesh in between.
- **Differentiable-CAD optimization.** Fields are closed-form and exactly
  differentiable: "minimize material s.t. wall ≥ 0.8mm", "solve for the
  parameter where these parts stop colliding." Gradient descent over the
  graph's free datums.
- **Kinematic scrubber.** A driver datum on a timeline slider — assemblies
  articulate in the viewport (the Zeiss focus mechanism actually racking).
  Also: motion GIFs for the gallery.
- **WASM web viewer.** Field evaluation compiled to WASM (Fidget proves
  it) → the GitHub gallery becomes interactive: visitors orbit models in
  the browser instead of squinting at JPGs.

*The thesis binding all of these: everyone else converts to triangles and
then fights the triangles. Stay math longer — slice, optimize, texture,
weigh, and animate the field, upstream of any mesh.*

## Agent-friendly modeling (LLM/automation as a first-class user)

Antimony already speaks Python-into-math; make the loop closable without
a GUI so agentic tools can contribute to modeling cleanly:

- **Headless verbs.** `stibium export model.sb out.3mf`, `stibium render
  model.sb -o view.png`, `stibium validate model.sb` (parse + evaluate +
  report errors as structured output). Render = the agent's eyes; export
  = hands; validate = conscience. (Same plumbing as the CLI renderer.)
- **`fab` as an installable pure-Python package.** The shape library
  already imports standalone with a small Shape/Transform shim (proven
  during chamfer testing). Publish it so scripts/agents can compose
  shapes and emit .sb or math strings outside the app.
- **Programmatic .sb authoring.** The format is JSON; document the
  schema (protocol 6) + connection encoding, provide a tiny writer
  library, and keep serialization deterministic so agent edits diff
  cleanly next to human edits in git.
- **Machine-readable node reference.** Generate JSON (name, inputs,
  types, defaults, category, doc line) from the self-describing .node
  files — same generator feeds the wiki for humans.
- **Structured errors.** Script/node errors surfaced on stdout/exit
  codes in the headless verbs, not just red text in a GUI datum.
- **MCP server (stretch).** load/save graph, list/add/edit nodes, eval
  bounds, render viewport → any agent-capable editor can drive a live
  Stibium session interactively.

## Distribution / project health

- GitHub Actions CI (build on push), AppImage/Flatpak releases — "runs on
  other people's machines" is a feature.
- Deterministic .sb serialization (stable key/node ordering) so models
  diff cleanly in git.
- Getting-started wiki (GitHub Pages); node reference partially generated
  from the self-describing .node files. Nothing comparable exists.
- Headless CLI renderer (`antimony-render model.sb -o front.png`) —
  render core is pure C, no Qt/GL. Unlocks the antimony-models gallery,
  wiki illustrations, and turntable GIFs. 2D path is first-class.
- Batch shrink tool for existing STL archives (prototype works; decide
  home + overwrite-vs-suffix policy).
- Port more libfive stdlib (MPL-2.0-compatible): gyroid/TPMS, rounded
  primitives, elongate, bend, twirl. Text-on-path for the vector font.

## Small bugs / cleanups

- **Root-cause the Python-error crash (the autosave origin story).** The
  live code-checker segfaulted randomly mid-typing; the workaround
  (fork-era, see `lib/graph/src/util.cpp` getPyError — the line-number
  extraction is commented out and pinned to 0) disables error line
  numbers entirely. Fix properly with defensive PyObject extraction so
  script errors report real line numbers again without the crash.
- **Infrequent crash when deleting groups of nodes at once.** Suspect
  the multi-delete undo path (`app/undo/undo_delete_multi.cpp`) or
  dangling proxies during batch removal. Needs a repro harness /
  ASan session like the detect-features hunt.


## Done
- 2026-07-12 — SVG export for 2D shapes (`export.svg` hook + GUI
  dialog). Marching squares over one field plane (`contour_field` in
  lib/fab): exact edge-indexed chaining, saddle cells resolved by
  center evaluation, padded border so clipped shapes close, and 2D
  feature detection — tangent-line intersection recovers sharp corners
  ~5× finer than the grid (<0.002 units on a 0.01 grid). Output in mm,
  evenodd fill, holes correct. Validated against analytic areas (<1%)
  and corner positions; 5 new CTest cases. Slices 3D shapes at their
  z midpoint. Export dialog gained "Contouring..." phase.
- 2026-07-12 — parallel meshing (`triangulate_indexed_mt`): chunked
  region → worker pool on cloned trees → seam-exact merge (re-intern)
  → global dedup/prune. Provably identical to serial with detect off;
  invariant-verified with detect on. Gyroid 2.57M tris: 11.7s → 2.8s
  (16 threads). Day total: 16.9s/723 MiB → 2.8s/335 MiB. Plus real
  export progress bar (exact voxel accounting via atomic, 20 Hz poll).
- 2026-07-12 — 3MF export (new default; STL stays). Minimal streaming
  ZIP writer over zlib (already a dep via libpng, no new vendoring) +
  3MF core-spec model XML in lib/fab/formats/threemf. Extension-driven
  in ExportMeshWorker, dialog defaults to 3MF; new `export.mesh` script
  hook (export.stl kept as alias). Simplify now compacts orphaned
  verts. Verified: python zipfile CRC/topology checks + PrusaSlicer
  (manifold=yes, correct volume, matches STL of same mesh); ~6× smaller
  than STL. No color/multi-material yet (color propagation is the
  prereq, see Tier 2.5).
- 2026-07-12 — indexed mesh storage in the mesher: interned vertices +
  3×uint32 triangles (tombstone erases, sort-based dedup/prune), new
  `triangulate_indexed`/`get_mesh` API, STL export end-to-end indexed
  (no meshopt re-weld; `save_stl_indexed` streams to disk; also fixed
  the sizeof(float) count write). ~4.6× less peak RSS, 1.47× faster
  with detect-features (723→159 MiB, 16.9→11.5 s @ 2.57M tris); 43M-tri
  export ~12→~3 GB. Output bit-identical (golden-dump verified); new
  tests/mesher.cpp property + equivalence suite. (all four eval backends, prefix +
  infix syntax, CTest coverage). Repeat node category: infinite/finite/
  mirror/polar domain repetition + self-similar scale recursion
  (repeat_scale, unlimited depth) + Iterate Scaled. All O(1) or
  documented; numerically verified against the real C engine. NOTE:
  scale-repeat needs factor > shape's radial span for visible gaps
  (else gapless solid → blank); r0 auto-derives from bounds.
- 2026-07-12 — Qt6 port (Qt 6.11, C++17, per-monitor HiDPI). Framework
  verdict: keep Qt, reject QML rewrite. deploy/ + binary rename deferred.
- 2026-07-12 — UI glow-up (Tiers A+B): warm amber palette + app-wide
  Fusion theme, gradient wires / round ports / shadowed nodes,
  type-tinted + operator-tinted headers, enhanced render mode (ambient/
  rim/AO/gamma) with draggable light gizmo, zoom LOD name-cards, fuzzy
  add-node palette (Tab/double-click), sticky wires, eased hovers,
  startup graph-focus, config system (~/.config/Stibium).
- 2026-07-12 — cleanup: removed dead update checker (dropped QtNetwork),
  wired SbGraphTest/SbFabTest into CTest, modernized BUILDING.md.
- 2026-07-12 — Phase 1 node campaign: 46 new nodes / 60+ shape functions
  (primitives, deforms, 2D stroke kit, functional-parts kit with ISO
  tables, ISO threads + involute gears). All numerically verified; new
  Parts/ node category. See doc/LIBRARY-ROADMAP.md.

- 2026-07-12 — meshoptimizer simplification in STL export (weld +
  error-bounded decimation, "max deviation" in the dialog, `simplify=`
  kwarg for scripted exports). ~94-99% smaller files.
- 2026-07-12 — chamfer + fillet union/intersection/difference nodes
  (six), math verified to machine precision.
- 2026-07-12 — fixed SIGFPE on thin models (region.c div-by-zero), parser
  new/free mismatch, mesher hardening (contour guard, normals-overflow
  early-return, safe erase idioms). Detect-features large-model "crash"
  root-caused as memory exhaustion (~275 B/tri), not a logic bug.
- 2026-07-12 — detect features enabled by default, "(experimental)"
  label dropped.

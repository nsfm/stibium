# Changes since Antimony

Everything Stibium has added or fixed relative to upstream Antimony
(`mkeeter/antimony`, `develop`). Unversioned until the first tagged
release; newest work at the top of each section.

## Geometry & export

- **The fmin toll removed** (doc/TAPE-DESIGN.md "Round 4"): since
  Antimony's first commit, every min/max node evaluated at every
  voxel paid a call into libm - `min_f` was C's double-precision
  `fmin`, meaning two float→double converts, a PLT call, and a
  truncation back, per lane. A CSG model's tape is mostly min/max. A
  sampled profile put ~40% of render wall time inside those calls.
  `min_f`/`max_f` are now bit-exact inline replicas of the
  platform's fmin/fmax semantics (proven by a 100M-pair parity
  harness against libm, then goldens + renders byte-identical), the
  batch kernels carry `ivdep`, and the fab kernel builds with
  `-fno-trapping-math` so the NaN-aware selects vectorize. Merged
  Zeiss at 2048px: 11.8s → **7.6s** on the same machine, same night,
  pixels identical. Exports gain ~5%.
- **The MPR tile renderer** (doc/TAPE-DESIGN.md "Round 3"): interval
  evaluation now runs in batches - one pass down the tape classifies
  64 view tiles at once, with the hot ops vectorized and every batch
  bound exactly equal to its scalar counterpart. Ambiguous-region
  subdivision fans out 64 ways into a single batched classification
  instead of bisecting in two, and the shaded pass pushes a pruned
  tape per 64px tile. Together: the merged Zeiss at 2048px renders
  ~20% faster, the win grows with resolution, and every image is
  pixel-identical (STIBIUM_TILE_RENDER=0 restores the old
  bisection). This is Keeter's 2020 MPR architecture running on CPU
  - the same machinery the eventual GPU renderer will reuse.
- **Tape round two: sound pruning, register allocation, deck
  caching** (doc/TAPE-DESIGN.md "Round 2"). A new pruning fuzzer
  (random expression trees x random regions, pushed-vs-base
  pointwise; 100k trees / 10.6M assertions in one run) found and
  fixed four real bugs - three predating the tape work entirely:
  domain-error intervals (sqrt/log/asin of out-of-range inputs)
  could drive unsound pruning decisions (fixed with libfive-style
  maybe-NaN taint, plus infinity rules libfive lacks); `pow_i`
  silently truncated real exponents to integers (pow(x, Z) got the
  bounds of pow(x, 0)); and the renderer's sign-collapse treated
  zero-touching bounds as one sign class, flipping exactly-zero
  fields under negation. The taint fix visibly repaired
  showcase_gear: the old render differs from unpruned ground truth
  by 5,842 pixels; the new one matches it bit-for-bit. Also landed:
  linear-scan register allocation (merged Zeiss: 3157 clauses run in
  465 registers, ~4x smaller per-thread workspaces - the
  precondition for tile/GPU rendering), per-shape deck caching (the
  viewport stops recompiling math every frame), a spare-tape
  freelist + no pushes inside packed blocks (4.7M -> 0.6M pushes per
  big export), mesher normal probes that escape their packed block
  now walk up to a covering tape, `STIBIUM_TAPE_STATS=1` shrinkage
  instrumentation (the MPR clauses-per-level curve: 20x by depth 6
  on the merged Zeiss), and an opt-in `STIBIUM_AFFINE=1` transform-
  chain collapse pass (bit-identity stays the default).
- **Shortened-tape evaluation (the libfive keystone).** The kernel's
  spatial pruning was rebuilt around immutable, refcounted, shareable
  tapes (doc/TAPE-DESIGN.md; the idea from libfive's `eval/tape.cpp`
  and the MPR paper, reimplemented for this kernel's four eval
  modes). A `Deck` compiles a MathTree once into a flat clause tape;
  interval evaluation then *pushes* shorter specialized tapes per
  region - decided min/max branches drop out by remapping readers to
  the surviving input, so pruned values are exact, not approximate.
  Replaces the mutate-in-place `disable_nodes`/`enable_nodes`/`ustack`
  machinery and the per-thread `clone_tree` copies (renderer, mesher,
  and parallel evaluator now share one deck; each worker brings only
  a slot-indexed workspace). The mesher also specializes per octree
  level instead of once per packed block. Verified bit-identical:
  eval equivalence unit tests across all four modes ("[tape]"),
  golden mesh dumps, and every --render PNG (zeiss included) compare
  byte-for-byte against the pre-tape build. Speedups scale with
  model size and CSG density (the MPR paper's law): gyroid ~10%,
  gear export ~1.4x, and the merged Zeiss ID02 export ~4x - with the
  full-scale Zeiss going from never-completes (killed at 13 min,
  nothing written) to done in 5 minutes. Memory flat. This is the
  Tier 1 keystone that unlocks the CPU SIMD tile viewport and an
  MPR-style GPU renderer.
- **Antialiased renders + turntables, wigglegrams, and stereo pairs.**
  All rendered images (--render, --diff composites, thumbnails, the
  GUI image export) now supersample 2x and smooth-downscale (--aa N
  to adjust, 1 disables). Three new headless verbs render motion:
  --turntable FILE.gif (looping rotation, --frames), --wiggle
  FILE.gif (two-frame depth wobble), and --stereo FILE.png
  (side-by-side parallel-view pair) - all with rotation-stable
  circumsphere framing so the model doesn't swim between frames.
  The GIF writer is a from-scratch GIF89a/LZW encoder (like the 3MF
  ZIP writer: no new dependencies), verified frame-exact against an
  external decoder. All three also live in the viewport's Render
  menu, rendered in the background with per-frame progress, spinning
  the model about world z as seen from your current view angle. A
  fourth, --lightsweep / Render > Light sweep, holds the model still
  and circles the key light instead, so shadows tell the form. The
  Render-menu animations get a size/frames/antialiasing dialog, and
  the still image export gains an antialiasing control too.
- **Check nodes: unit tests for geometry** (new Checks category).
  Check: Fits Box, Check: Volume, and Check: Clearance measure the
  actual field (grid integration; tight bounds accurate to a cell)
  and raise with a precise message when the spec is violated - the
  node goes red in the canvas, and `--validate` exits nonzero, so a
  broken spec fails CI exactly like a broken test. Checks pass the
  shape through, so they wire inline ahead of an export node (a
  failing check blocks the export). Built on a new `fab.shapes.
  measure(shape)` primitive exposing volume/area, center of mass,
  tight bounds, and sample counts to every node script.
  `examples/showcase_gear.sb` now carries its own spec.
- **Geometric diff** (`stibium a.sb --diff b.sb [--render out.png]`):
  set algebra on the distance fields splits two models into
  unchanged / removed / added regions - printed as integrated
  volumes (JSON, with an exact-identity shortcut) and optionally
  rendered as a composite (gray / red / green). Meshes can't diff
  geometry cheaply; fields literally subtract. Review a CAD change
  like a code change.
- **Per-shape color in rendered images**: `--render` and the image
  export now render each terminal shape separately and composite by
  depth, so `set_color` colors survive into thumbnails and exports
  (previously everything flattened into one warm-gray union).
- **Fidget `.vm` import** (`stibium model.vm --import-vm out.sb`):
  translates math tapes from Matt Keeter's successor project
  (github.com/mkeeter/fidget) into Stibium projects - one frozen
  shape node holding the full field. Direct op mapping plus algebraic
  rewrites (ceil/round/recip); atan2 via the V2-in-V1 grammar embed;
  discontinuous ops (compare/and/or/not) rejected with a clear error;
  expression-depth guard against the parser stack. Bounds are
  auto-discovered by two-pass interval octree descent (scale-free,
  ~0.2% tightness), with 2D tapes (no var-z) detected and emitted as
  planar shapes. Deterministic output, in CI. Fidget's flagship
  7,866-op `prospero.vm` converts in under a second and renders
  pixel-perfect. New stdlib (from the libfive raid, see
  doc/LIBFIVE-RECON.md): log-sum-exp smooth union/difference
  (`blend_expt_unit`, `blend_difference` + Smooth Union/Difference
  nodes), exact-field `rectangle_exact` and `rounded_box` (+ nodes),
  and a `gyroid()` function (fixing the Gyroid node's period being
  ~23x off).
- **Mesh import** (upstream's most-requested feature, mkeeter/antimony#153):
  the new Import > Import Mesh node turns an STL into a solid distance
  field that composes with the whole vocabulary — CSG, transforms,
  deforms, analytics, re-export. The mesh is sampled once onto a
  signed distance grid (exact BVH distance; inside/outside from the
  generalized winding number, so non-watertight scans still classify)
  and evaluates as a new ternary `OP_GRID` opcode: trilinear point and
  gradient lookups, conservative block-min/max intervals (all interval
  pruning keeps working), heap payload refcounted and shared across
  thread clones. Paths are project-relative (portable projects, no
  baked absolute paths); processed grids cache in a regenerable
  gitignored `.stibium-cache/` beside the project keyed by
  sha256+resolution, so re-opens don't re-sample. An optional `sha256`
  input pins the source file and fails loudly if it changes; importing
  a Stibium-stamped export gets an advisory to open the source `.sb`
  instead. STL exports now carry a Stibium stamp in their 80-byte
  header to power that advisory. Example: `examples/import_bead.sb`.
- **DXF export**: the same contours as R12 closed POLYLINEs (y-up,
  mm) — the flavor laser-cutter and CAM toolchains read. One
  extension-driven vector exporter covers both (`export.vector`).
- **Contour tracing skips empty space**: interval evaluation proves
  regions empty or solid before any sampling, so sparse high-res
  vector exports (photolithography masks) drop from seconds to
  milliseconds of evaluation.
- **SVG export**: 2D shapes trace to resolution-independent vector
  paths (`export.svg`, sharp corners recovered by feature detection to
  ~5× finer than the sampling grid, holes and clipped shapes handled,
  mm units). Replaces giant raster exports for photolithography masks
  and laser cutting. 3D shapes export their z-midpoint cross-section.
  Contour evaluation runs on all physical cores (~5-8× faster), and
  paths are Douglas-Peucker simplified to a max-deviation bound
  (default: quarter cell; `simplify=` kwarg) — corners provably
  survive, redundant chords don't (~4-5× smaller files).
- **Parallel meshing**: exports fan out across all CPU cores (chunked
  region, per-thread tree clones, seam-exact merge). 2.57M-triangle
  gyroid: 11.7 s → 2.8 s on 16 threads; combined with the indexed
  mesher, the full day's delta is 16.9 s / 723 MiB → 2.8 s / 335 MiB.
  Provably identical output to serial meshing with feature detection
  off; watertight-verified with it on.
- **Real export progress**: the exporting dialog's bar now tracks
  actual meshing progress (exact voxel accounting, including
  fast-forward over empty space) instead of spinning indeterminately.
- **3MF export** (the new default; STL stays for compatibility). The
  export dialog and scripted exports pick the format by extension; a
  new `export.mesh` hook is the documented spelling (`export.stl`
  remains as an alias). Written by a minimal streaming ZIP writer over
  zlib - no new dependencies - and verified against PrusaSlicer
  (manifold, correct volume) and python zipfile CRC/topology checks.
  Typically ~6× smaller than the same mesh as STL.
- **Indexed mesh storage in the mesher**: triangles are 12-byte vertex
  index triples over an interned vertex table instead of a linked list
  of double-precision corners with tree-based end-of-run dedup. ~4.6×
  less peak memory and ~1.5× faster on the feature-detection path
  (gyroid, 2.6M tris: 723 → 159 MiB, 16.9 → 11.5 s); the 43M-triangle
  export that needed ~12 GB now fits in ~3. STL export consumes the
  indexed mesh directly (no soup expansion, no meshopt re-weld) via
  `save_stl_indexed`. Output verified bit-identical to the old mesher;
  new mesher test suite (topology, analytic volume/area, sharp-corner
  reconstruction, golden dumps, indexed==soup equivalence).
- **`mod` and `floor` opcodes** in the math engine (prefix `M`/`F`,
  infix `mod()`/`floor()`), with float, interval, gradient, and region
  backends and CTest coverage.
- **`log` opcode** (prefix `l`, infix `log()`) - unlocks log-space
  domain repetition: **Repeat Scale**, infinite self-similar recursion
  about a point at O(1) cost, plus **Iterate Scaled** automating finite
  copy-translate-scale chains.
- **Domain repetition nodes (new Repeat category)**: infinite repeat
  along X/Y/Z, XY grids, mirrored repeat, exact finite repeat, and
  polar repeat - all at O(1) field cost regardless of copy count
  (the Array nodes union N copies at O(N) and crash on large counts).

- **Mesh simplification on STL export**: triangle soup is welded and
  decimated to a user-set "max deviation" (meshoptimizer, vendored).
  Typical exports shrink 94–99% with sub-print-tolerance deviation.
  Scripted exports accept `export.stl(shape, simplify=...)`.
- **Feature detection enabled by default** (no longer "experimental").
  Its large-model crash was root-caused to memory exhaustion, not logic
  (verified to 43M triangles under sanitizers).
- **Chamfer + fillet CSG** (union/intersection/difference × 2): 45°
  chamfers and rounded fillets with parametric radius.

## Node library (Phase 1 campaign, 46 nodes / 60+ functions)

- **Parts category (new)**: ISO metric threaded rod + tapped hole
  (provably mating, exact 60° flank), involute spur gears, racks,
  teardrop holes, clearance/counterbore/countersink holes with M2–M12
  tables, hex nut traps, heat-set bosses, PCB standoffs, snap clips,
  lid tongue-and-groove, dovetails, hex vent grilles.
- **2D stroke kit**: segment, polyline, arbitrary-vertex polygon,
  quadratic bezier; star, pie, arc, annulus, trapezoid, rhombus,
  oriented box, vesica, crescent, cross, ellipse.
- **3D primitives**: capsule, capped cone, rounded cylinder, hex/tri
  prisms, cut sphere (dome), slab, half-space, exact-SDF box.
- **Deforms**: bend_x/y, twirl_x/y/z (fixing a per-axis bug present in
  libfive's stdlib), revolve_z, elongate.
- **Ring and Mirror nodes**; extended vector font (new glyphs, kerning
  fixes, masking fixes). *(fork era)*

## UI

- **Loading splash**: opening a large graph (hundreds of nodes takes
  a few seconds to evaluate) now shows a splash card naming the file
  instead of a frozen blank window. Shown from startup and forced to
  paint before the blocking load begins, so it's visible for the whole
  wait rather than flashing once at the end.
- **Fixed a crash on File > New (and graph clear) after loading a
  file**: the per-node canvas inspector cards were only deferred-
  deleted when their nodes were freed, so a repaint in the gap
  dereferenced a dead node (both the card's own paint and the
  floating-label sweep read the node). Cards now detach synchronously
  the moment their node goes away - unregistered from the label
  registry and hidden - so no paint path can touch freed memory.
- **Ctrl-C from the launching terminal** runs the normal quit path,
  including the unsaved-changes prompt, instead of killing the
  process (self-pipe SIGINT handler into the Qt event loop; Python
  is initialized without its own signal handlers so it can't reclaim
  SIGINT). A freshly loaded file is marked clean, so quitting right
  after opening never prompts.

- **Analytics overlay**: Render ▸ Analytics overlay integrates the
  model and pins volume, center of mass (with an in-scene crosshair
  marker), and tight bounds to the viewport.
- **Recent files & directory memory**: File ▸ Open Recent, and all
  export dialogs open in the last export directory (kept separate
  from the model directory).
- **Export nodes glow at low zoom**: render/mesh/export nodes get
  amber skeleton cards with halos and role-based floating labels
  ("Mesh 0" instead of "m0") with top priority — the nodes new users
  hunt for announce themselves.
- **Export dialog shows real dimensions**: model size in mm replaces
  the voxel-count readout; a live "min feature ≈" line (mm/µm) shows
  what the chosen resolution can resolve; simplification deviation
  follows at a tenth of a voxel (empirically tuned on real prints).
- **Saner export defaults**: mesh dialogs start at 7 voxels/unit,
  vector/raster at 60 (the old content-scaled default suggested 1 -
  a useless mesh), and the dialog now shows the minimum feature size
  the chosen resolution can resolve.
- **Render ▸ Export image**: render the current view to a PNG from
  the viewport window - orientation-matched, optional section cut,
  size and transparent-background options. Same engine as the
  headless `--render`.
- **Parallel viewport rendering**: the raymarched render fans out
  across physical cores (~3× on heavy models), in the viewport,
  `--render`, and heightmap export. Output is byte-identical to
  single-threaded.
- **No more zoom clipping**: the render's depth range now always
  covers the whole model (voxel budget capped instead), so extreme
  close-ups stop guillotining geometry. Deep zoom is now a quality
  target, not a wall.
- **Viewport chrome**: the section slider is permanent (left edge,
  View ▸ Section slider), the key-light gizmo has a tooltip and a
  View ▸ Light control toggle, and a corner eye button shows/hides
  bounding boxes (lit = visible; the old "Hide UI" menu item is now
  an honest "Bounding boxes" checkbox).
- **Cross-section view**: press X in a viewport for a slider that
  pulls a screen-parallel cut plane through the model — rotate to aim
  the cut, slide to expose walls, voids, and clearances. Also
  available headlessly (`--render --section 0.5`).
- **Warm dark theme**: new palette (amber accent, warm charcoal ramp),
  Fusion style + app-wide stylesheet so menus/dialogs/scrollbars match
  the canvas. Fixed `Colors::adjust()` so hover highlights work on dark
  fills.
- **Canvas facelift**: gradient wires with hover glow and shadows,
  circular 14px ports (up from 10px squares) with hover rings, node
  drop shadows + gradient fills + amber selection, zoom-faded two-level
  dot grid, tinted rubber-band.
- **Type-tinted node headers**: title bands (and zoomed-out cards) tint
  by output datum type — shape green, float yellow, etc. — and
  shape-consuming shape-producers (CSG/deforms/transforms) tint violet,
  so operators read differently from sources at a glance.
- **Floating node labels**: far zoom-out fades in web-map style
  labels — fixed screen size, overlap-avoiding, custom-named nodes
  prioritized (amber-edged) over auto-named ones, leader lines when
  nudged. Big graphs read like a labeled map.
- **Zoom LOD**: below ~32% zoom nodes collapse to solid type-tinted name
  cards and wires simplify to cosmetic strokes; graphs read as a labeled
  map instead of gray fuzz. Zoom-out unclamped deep (0.02).
- **Fuzzy add-node palette**: double-click canvas or press Tab for a
  type-to-filter popup over all nodes (prefix > word > substring >
  subsequence); right-click keeps the browsable menu.
- **Sticky wire connections**: release a wire away from a port and it
  stays live following the cursor (canvas pans/selects underneath);
  left-click a valid port to complete, right-click/Escape to cancel.
  Classic drag-to-connect still works.
- **Enhanced render mode** (new default, menu-selectable alongside
  Shaded/Height-map): hemispheric ambient, key light, depth-buffer AO,
  fresnel rim, gamma-correct output. Key-light direction is draggable
  via a trackball gizmo in the viewport corner (persists in config).
- **Viewport chrome**: near-black gradient backdrop, labeled X/Y/Z
  axes, monospace coordinate readout.
- **Eased hover/selection glows** (~120ms) on nodes and wires.
- **Canvas opens centered on the loaded graph** (zooming out to fit),
  instead of at the origin.
- Script editor uses the system fixed-width font (was hardcoded
  Courier), configurable via the config file.

## Fixes

- **Stale-UI-header builds fixed at the root**: AUTOUIC under Ninja
  ran one build behind on .ui edits (mismatched dialog layouts,
  startup segfaults, ritual clean rebuilds). The .ui files are now
  explicitly wrapped (qt6_wrap_ui); a plain `ninja` is always
  sufficient.
- **Undo-after-open crash**: the undo stack survived file loads,
  so undoing after File ▸ Open replayed commands into the freed
  previous graph. The stack now clears on every load.

## App & infrastructure

- **App-level test suite**: CTest drives the real binary over every
  example (validate, render, deterministic resave round-trip) - CI
  now asserts behavior, not just compilation.
- **Headless render**: `--render out.png [--view iso|top|front]
  [--size N]` draws the model's shapes shaded on a transparent
  background with no display or GL — gallery scripts, CI previews,
  and file-manager thumbnailers can see models directly.
- **Protocol 7 + deterministic saves**: output datums no longer
  serialize their computed-geometry reprs (scripts regenerate them on
  load; only the uid-anchoring stub remains). Files shrink 8-35%,
  git diffs show real changes only, and load→save round-trips are
  byte-identical. `--resave FILE` batch-migrates old files.
- **Field analytics**: `--analyze` integrates the model and reports
  volume (or area for 2D), center of mass, and tight bounds as JSON
  — guess-and-check counterweighting ends here. `--node NAME` scopes
  to a single node.
- **Node library as JSON**: `--describe-nodes` dumps every node's
  name, category, typed inputs (with defaults), and outputs — the
  vocabulary for agents and the seed for generated wiki docs.
- **Per-node renders**: `--render --node NAME` draws one node's
  geometry alone, construction or not — visual bisection when the
  final shape is wrong.
- **Headless CLI verbs**: `--validate` (script/datum errors to stderr,
  exit code) and `--export FILE [--resolution R] [--detect-features]`
  run a model's export node with no display or dialogs — CI, batch
  jobs, and AI agents can drive exports directly.
- **Live reload**: the open .sb reloads automatically when it changes
  on disk (only while the session has no unsaved edits), so external
  tools and agents can edit a model while it's open.

- **Qt6 port** (Qt 6.11 verified; per-monitor HiDPI, modern platform
  support). C++17 across targets.
- **Autosave** with menu toggle. *(fork era — the founding feature)*
- **Config system**: INI at `~/.config/Stibium/Stibium.ini`
  (platform-appropriate elsewhere); autosave interval, editor font,
  last-used directory, and enhanced-mode key-light direction.
- **Test suites wired into CTest** (SbGraphTest + SbFabTest, 237
  assertions).
- **Removed the update checker** (polled upstream's frozen releases;
  dropped the QtNetwork dependency with it).
- Renamed to Stibium in all user-facing surfaces (About box, dialogs);
  binary/install-path plumbing rename deferred to first release.
- CMake 4 compatibility; modernized BUILDING.md (Qt6, Arch recipe).

## Bug fixes

- `exp(constant)` parsed as `abs(constant)`: the exp node's
  constant-fold called the wrong math function (upstream-era typo in
  `node_c.c`; only fired when exp received a literal constant).
- SIGFPE crash exporting models with a sub-voxel-thin axis
  (`region.c` division by zero).
- `new`/`free` mismatch in the math-string parser (UB on every parse
  since 2015).
- Mesher hardening: buffer-overflow guard now returns instead of
  printing, degenerate-contour guard, UB-free erase idioms.
- Segfault printing Python errors (worked around fork-era; proper fix
  tracked in TODO).

## Removed

- Update checker (polled upstream's frozen releases; dropped the
  QtNetwork dependency with it).

## Identity & licensing

- Renamed to **Stibium** in all user-facing surfaces, with full
  Antimony/kokopelli lineage in the About box. Binary/install plumbing
  renames deferred to the first packaged release.
- New work licensed GPLv3; inherited MIT preserved in
  `THIRD_PARTY_LICENSES.md`.

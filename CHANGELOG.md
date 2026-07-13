# Changes since Antimony

Everything Stibium has added or fixed relative to upstream Antimony
(`mkeeter/antimony`, `develop`). Unversioned until the first tagged
release; newest work at the top of each section.

## Geometry & export

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

## App & infrastructure

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

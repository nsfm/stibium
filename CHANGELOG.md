# Changes since Antimony

Everything Stibium has added or fixed relative to upstream Antimony
(`mkeeter/antimony`, `develop`). Unversioned until the first tagged
release; newest work at the top of each section.

## Geometry & export

- **`mod` and `floor` opcodes** in the math engine (prefix `M`/`F`,
  infix `mod()`/`floor()`), with float, interval, gradient, and region
  backends and CTest coverage.
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
  dot grid, tinted rubber-band, zoom clamped to [0.08, 4].
- **Enhanced render mode** (new default, menu-selectable alongside
  Shaded/Height-map): hemispheric ambient, key light, depth-buffer AO,
  fresnel rim, gamma-correct output.
- **Viewport chrome**: near-black gradient backdrop, labeled X/Y/Z
  axes, monospace coordinate readout.
- **Canvas opens centered on the loaded graph** (zooming out to fit),
  instead of at the origin.
- Script editor uses the system fixed-width font (was hardcoded
  Courier), configurable via the config file.

## App & infrastructure

- **Qt6 port** (Qt 6.11 verified; per-monitor HiDPI, modern platform
  support). C++17 across targets.
- **Autosave** with menu toggle. *(fork era — the founding feature)*
- **Config system**: INI at `~/.config/Stibium/Stibium.ini`
  (platform-appropriate elsewhere); autosave interval and editor font
  as first settings.
- **Test suites wired into CTest** (SbGraphTest + SbFabTest, 237
  assertions).
- CMake 4 compatibility.

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

# Shape library roadmap

Gap analysis of this fork's node library vs. libfive's stdlib, kokopelli
(Antimony's predecessor), Inigo Quilez's SDF catalog, fogleman/sdf, Curv,
and BOSL2/OpenSCAD expectations. Compiled 2026-07-12. Companion to the
campaign list in TODO.md.

## Ground truth

**Opcodes available** (`lib/fab/inc/fab/tree/node/opcodes.h`):
`+ - * / min max pow abs square sqrt sin cos tan asin acos atan atan2 neg
exp` plus coordinates and constants. Note **atan2 exists** (full
float/interval/gradient backends) but is reachable only through the infix
`=expr;` syntax — no prefix character. `clamp(x,a,b)` composes as
`min(max(x,a),b)` today.

**Genuinely absent:** `mod`, `floor`, `round`, `sign`, `log`, noise.

**Cost of a new opcode:** ~12 mechanical touchpoints, all switch/table
additions — opcodes.h/.c, eval.c, math_f/i/g/r.c, node_c.c, printers.c,
printers_ss.cpp, parser.c, v2syntax.l/.y. Each opcode needs float,
interval, AND gradient implementations (rendering and meshing depend on
all three). `mod`'s interval case needs care at the wrap boundary.

## Phase 1 — pure .node / shapes.py additions (no C changes)

In value-per-effort order:

1. **Primitive gap-fill (one-liners):** half-space/plane (arbitrary
   normal), 3D capsule (round-capped segment between two points), capped
   cone / frustum (two radii), rounded cylinder & cone, hex prism
   (nut-trap core), cut sphere (dome), slab, elongate, exact-SDF box.
2. **2D stroke kit:** line segment (round caps), **polyline / polygon
   from a vertex list** (biggest 2D gap), quadratic bezier segment.
   Unlocks stroke-based drawing for the photolithography workflow.
3. **2D shape set:** star (n points, inner/outer r), pie/sector, arc,
   annulus, trapezoid, rhombus, oriented box, vesica, ellipse
   (scaled-circle approx), crescent, cross.
4. **Functional-parts kit (parametric CSG):** teardrop hole
   (overhang-safe horizontal holes), countersink / counterbore /
   clearance hole with screw-size presets (M2–M12 lookup table),
   hex nut trap + square-nut slot, heat-set insert boss, PCB standoff
   (tube + fillet base), snap-fit cantilever clip + lid lip, dovetail
   (male/female with clearance), hex-grille vent panel (trig lattice).
5. **Threads and gears — possible TODAY, no new opcodes:** a thread is
   `r − profile(cos(nθ − 2πz/pitch))` in cylindrical coordinates via
   atan2 (trig is inherently periodic, so no mod needed). ISO metric
   presets + printer-clearance parameter. Involute spur gear via one
   tooth + N rotated copies (kokopelli-era demos did this), rack as a
   prism array. Caveat: test interval evaluation at the atan2 branch cut
   early; orient the seam away from geometry.
6. **Deform set:** general bend (arc), twirl (localized twist with
   falloff, from libfive), revolve about Z (odd gap — x and y exist),
   twist-extrude combo node, trig displacement (knurling, ripples) —
   underexploited with existing opcodes.
7. **CSG extras:** blend_rough (libfive), loft_between (4-point loft),
   pipe/groove/tongue seam ops (bead along the intersection curve).
8. **Infix function node** exposing the `=expr;` syntax so atan2 (and
   user formulas) are reachable without writing prefix strings.

## Phase 2 — the mod/floor opcode PR (highest-leverage change)

Add `mod` + `floor` together (`round`/`fract` derive; take `sign` and
native `clamp` opportunistically in the same PR). Add prefix chars +
infix tokens. This single change unlocks:

- **Domain repetition**: infinite repeat x/y/z/xyz, exact finite repeat,
  radial/polar repeat, mirror-alternating repeat — all **O(1) field cost
  regardless of copy count** (the existing Array nodes union N copies at
  O(N), and have a known crash FIXME at shapes.py:532 that mod-repetition
  sidesteps entirely).
- Per-cell indexed variation (floor(p/s) as cell ID driving
  size/rotation/hash) — recursive/fractal tilings for litho masks.
- Voxelize/pixelate, brick/parquet patterns, stairs-blend CSG.
- fract-sin hash — poor-man's noise for speckle/grain before Phase 3.

Interval-arithmetic caveat: when the input interval spans a period
boundary, mod must return the full period range.

## Phase 3 — heavy items

1. **Noise opcode** (value or simplex; interval = ±amplitude bounds,
   gradient analytic) → woodgrain/marble/organic displacement as nodes,
   with seed inputs (reproducible-unique).
2. Image → field import (litho masks, engraving, lithophanes).
3. FreeType text (stroke font exists; outline font is the upgrade).
4. Sweep along arbitrary path (partially composable from polyline+loft).
5. PCB footprint library — only if kokopelli-style board import (outline
   keep-out + mounting-hole pattern + rounded-trace polylines) proves
   insufficient.

## Kokopelli heritage note

Antimony's predecessor had a PCB library (pads, footprints, traces as
rounded segment chains, board outline + drills) that never made the
migration. The useful extract for enclosure work: board-outline import,
mounting-hole pattern node, and the 2D polyline op (Phase 1.2).

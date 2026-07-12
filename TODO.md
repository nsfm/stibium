# TODO — potential upgrades

Running list of improvements considered for this fork. Larger items first.

## Meshing / export

- **Indexed mesh storage in the mesher.** `Mesher` holds the whole mesh as
  `std::list<Triangle>`, and with feature detection on, `remove_dupes()` /
  `prune_flags()` build a `std::map` + `std::set` over every vertex and
  triangle on top of it (~275 bytes/triangle measured, ~2.2x the
  detect-off path; 43M triangles ≈ 12 GB). Moving to flat vectors with an
  indexed vertex buffer would roughly halve memory, remove the OOM crash
  on large exports, and speed up meshing considerably (pointer-chasing a
  list dominates cache behavior). `4x5_ground_glass_holder` takes ~2 min
  to mesh today.
- **Parallel meshing.** `triangulate_region()` recursively subdivides into
  independent sub-regions — a natural fork/join. The mesher is currently
  single-threaded; the machine has cores to spare during export.
- **Adaptive meshing (libfive-style dual contouring).** The uniform grid
  is why raw exports are enormous. libfive (same author, MPL-2.0 core) does
  manifold dual contouring over an adaptive octree — fewer, feature-aligned
  triangles at the source instead of decimating after the fact. Big
  engineering; the meshoptimizer post-pass covers most of the value.
- **Report simplification stats after export** ("2.1M → 48k triangles,
  max deviation 0.009") — the numbers exist, there's just no UI channel.
- **Batch shrink tool for existing STLs** (weld + error-bounded simplify,
  mirrors the export pipeline). Working prototype exists; needs a decision
  on where it lives (antimony-models/scripts/?) and overwrite-vs-suffix.

## Nodes / shape library

- Port more of the libfive stdlib (MPL-2.0, license-compatible):
  gyroid / TPMS infill, polar + linear arrays beyond stock, rounded
  primitives (capsule, rounded box), elongate, bend, twirl.
- 2D-specific ops for the vector-graphics workflow (pure-2D fillet/chamfer
  variants; text-on-path; stroke/outline of a 2D field).

## UI

- **Canvas annotations**: comment boxes + named highlighted zones behind
  nodes (QGraphicsObject at low z, persisted as canvas metadata in the .sb
  JSON alongside node positions). Both sticky notes and named zones.
- Enable "Detect features" by default in the resolution dialog — DONE
  (2026-07-12), label de-experimentalized.

## Tooling

- **Headless CLI renderer** (`antimony-render model.sb -o front.png`):
  render core (`render16`/`shaded8`) is pure C, no Qt/GL needed. Unlocks
  a generated gallery for antimony-models and wiki illustrations.
  Reference harness for driving lib/fab headless already exists (parse
  math string → Region → triangulate/render).
- **Getting-started wiki** (GitHub Pages). Nothing comparable exists
  anywhere; node reference could be partially generated from the .node
  files (they're self-describing Python).

## Known small bugs / cleanups

- `v2parse()` leaks scanner/parser/locals on the parse-failure early
  return (v2parser.cpp) — pre-existing, low priority.
- `RenderTask::render()` (viewport/render/task.cpp) missing-field warnings
  for Region initializers; harmless but noisy.
- STL writer (`lib/fab/src/formats/stl.c`) writes the triangle count via
  `sizeof(float)` — works because sizeof(int)==sizeof(float), but it's
  fragile-looking; trivial cleanup.

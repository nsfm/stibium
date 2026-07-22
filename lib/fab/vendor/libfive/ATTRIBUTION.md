# Vendored libfive (Dual Contouring mesher subset)

This directory contains a **subset** of the libfive CAD kernel, vendored so
Stibium can run libfive's flagship Dual Contouring mesher head-to-head against
Stibium's own `stibnite` adaptive-Delaunay mesher.

- **Upstream:** <https://github.com/libfive/libfive>
- **Commit:** `c9e97343e0af998cd1696e85583eccba95532b96` (2025-11-12)
- **License:** Mozilla Public License, version 2.0 (the `libfive/` core
  layer). Full text in [`LICENSE.MPL2`](LICENSE.MPL2). Every source and
  header file retains its original per-file MPL-2.0 notice.
- **Copyright:** (c) 2015-2025 Matt Keeter and libfive contributors.

MPL-2.0 is compatible with Stibium's AGPL-3.0-or-later: MPL-covered files stay
under the MPL (file-level copyleft), and they are redistributed here unmodified
except where noted below.

## What is vendored

Only the pieces needed to go from a `libfive::Tree` to a Dual-Contouring mesh:

- `include/` -- the libfive public + internal headers (whole tree, for the
  templated brep/eval machinery).
- `src/tree/`, `src/oracle/`, `src/eval/` -- the expression tree and the
  interval / array / derivative / feature evaluators.
- `src/render/brep/{edge_tables,manifold_tables,neighbor_tables,progress}.cpp`
  and `src/render/brep/dc/*`, `src/render/brep/vol/*` -- the DC octree,
  manifold tables, dual walk, and volume acceleration structure.

## NOT vendored (deliberately)

- The `ISO_SIMPLEX` and `HYBRID` meshers (`src/render/brep/simplex/`,
  `.../hybrid/`). Stibium ports libfive's **default** (`DUAL_CONTOURING`).
- Studio (GPL), the language bindings, the stdlib, the solver, contours (2D).

## Modifications

- `src/render/brep/mesh_dc.cpp` is libfive's `src/render/brep/mesh.cpp` with
  the `ISO_SIMPLEX` and `HYBRID` branches (and their includes) removed, so the
  simplex/hybrid translation units need not be compiled. The `DUAL_CONTOURING`
  path is byte-for-byte upstream. The change is marked in-file with a
  `NOTE (Stibium vendor trim)` comment.

The Stibium<->libfive bridge itself (`libfive_bridge.cpp/.h`) is Stibium's own
work (AGPL-3.0-or-later), not part of libfive.

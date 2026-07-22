# Attributions

Stibium builds on ideas and prior art. This file credits that lineage.

For the licenses of code Stibium actually redistributes (inherited source,
vendored libraries, bundled fonts), see
[`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md).

---

## Direct lineage

Stibium is a continuation of **[Antimony](https://github.com/mkeeter/antimony)**
by **[Matt Keeter](https://mattkeeter.com)**, itself a spiritual successor to
**[kokopelli](https://github.com/mkeeter/kokopelli)** by way of
[fabserver](http://kokompe.cba.mit.edu). The functional-representation model,
the node-graph editor, and the "no meshes until export" philosophy all descend
from that work. Inherited Antimony/kokopelli code is redistributed under its
original MIT license (see `THIRD_PARTY_LICENSES.md`).

---

## Algorithmic concepts

Stibium's geometry kernel (`lib/fab`) owes a real intellectual debt to two
later projects by Matt Keeter.

### libfive

- **Project:** <https://libfive.com> · <https://github.com/libfive/libfive>
- **Author:** Matt Keeter and contributors
- **What we learned from it:** the modern f-rep kernel design - evaluating
  signed-distance/interval expression trees over a spatial hierarchy, pruning
  empty regions by interval arithmetic, and the feature-preserving
  dual-contouring family of meshing built on that evaluation. These ideas
  inform Stibium's evaluator and meshing paths.

### fidget

- **Project:** <https://www.mattkeeter.com/projects/fidget/> ·
  <https://github.com/mkeeter/fidget>
- **Author:** Matt Keeter
- **What we learned from it:** the "tape" view of expression evaluation (a
  flattened instruction stream for an implicit surface) plus interval-based
  tape shortening, register scheduling, and the GPU tape pipeline. Stibium's
  export path reimplements these concepts against its own VM.

---

_Miscredited or missing a debt? Open an issue._

# Third-party licenses

Stibium's own new work is licensed under the **GNU Affero General Public
License v3.0 or later** (see [`LICENSE`](LICENSE)). This file accounts for
every third-party work that Stibium redistributes together with the license
each is carried under. All of them are compatible with AGPL-3.0-or-later.

---

## Inherited source

### Antimony

- **Role:** Stibium is a continuation of Antimony; inherited code is spread
  throughout the tree.
- **Upstream:** <https://github.com/mkeeter/antimony>
- **License:** MIT
- **Copyright:** © 2013–2022 Matthew Keeter and other contributors
- **Full text:** reproduced below.

Antimony in turn includes code from
[kokopelli](https://github.com/mkeeter/kokopelli):

- **Copyright:** © 2012–2013 Massachusetts Institute of Technology;
  © 2013 Matthew Keeter
- **License:** MIT (same text below).

```
Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
```

---

## Vendored libraries

### meshoptimizer

- **Role:** mesh simplification / optimization for export.
- **Location:** [`app/vendor/meshoptimizer/`](app/vendor/meshoptimizer/)
- **Upstream:** <https://github.com/zeux/meshoptimizer>
- **Vendored version:** 1.0 (`MESHOPTIMIZER_VERSION 1000`)
- **License:** MIT
- **Copyright:** © 2016–2025 Arseny Kapoulkine
- **Full text:** [`app/vendor/meshoptimizer/LICENSE.md`](app/vendor/meshoptimizer/LICENSE.md)

### Eigen

- **Role:** linear algebra for the `fab` geometry kernel.
- **Location:** [`lib/fab/vendor/Eigen/`](lib/fab/vendor/Eigen/)
- **Upstream:** <https://eigen.tuxfamily.org>
- **Vendored version:** 3.2.4
- **License:** MPL-2.0 (primary). A few files carry BSD or LGPL-2.1
  licenses; see the upstream note in
  [`lib/fab/vendor/Eigen/COPYING.README`](lib/fab/vendor/Eigen/COPYING.README).
  All are AGPL-compatible. (Defining `EIGEN_MPL2_ONLY` would restrict the
  build to the MPL-2.0/BSD subset if a pure-MPL2 guarantee is ever wanted.)
- **Full text:** [`lib/fab/vendor/Eigen/COPYING.MPL2`](lib/fab/vendor/Eigen/COPYING.MPL2)
  (and the sibling `COPYING.*` files).

### libfive (Dual Contouring mesher subset)

- **Role:** vendored subset of the libfive CAD kernel, used to run libfive's
  flagship Dual Contouring mesher head-to-head against Stibium's stibnite
  mesher (experimental `STIBIUM_EXPORT_LIBFIVE=1` export path). Built as an
  isolated shared library; not part of the default mesher.
- **Location:** [`lib/fab/vendor/libfive/`](lib/fab/vendor/libfive/)
- **Upstream:** <https://github.com/libfive/libfive>
- **Vendored commit:** `c9e97343e0af998cd1696e85583eccba95532b96` (2025-11-12)
- **License:** MPL-2.0 (the `libfive/` core layer; the GPL Studio is **not**
  vendored). Full text:
  [`lib/fab/vendor/libfive/LICENSE.MPL2`](lib/fab/vendor/libfive/LICENSE.MPL2).
  Every vendored file retains its per-file MPL-2.0 notice.
- **Copyright:** © 2015–2025 Matt Keeter and libfive contributors.
- **Modifications & provenance:**
  [`lib/fab/vendor/libfive/ATTRIBUTION.md`](lib/fab/vendor/libfive/ATTRIBUTION.md)
  (one file trimmed to DC-only; the Stibium↔libfive bridge is Stibium's own
  AGPL work). MPL-2.0 is file-level copyleft and AGPL-3.0-compatible.

### stb_truetype

- **Role:** glyph outline extraction for the Text (font) node.
- **Location:** [`lib/fab/vendor/stb/`](lib/fab/vendor/stb/)
- **Upstream:** <https://github.com/nothings/stb>
- **Vendored version:** stb_truetype v1.26
- **License:** dual-licensed — **MIT or public domain (Unlicense)**, your choice.
- **Copyright:** © Sean Barrett
- **Full text:** in the header of
  [`lib/fab/vendor/stb/stb_truetype.h`](lib/fab/vendor/stb/stb_truetype.h).

---

## Bundled fonts

Five fonts ship for the Text (font) node, each under its own open license
(OFL 1.1, Apache-2.0, or the Bitstream Vera license). They are redistributed
pristine and unmodified; a baked SDF grid is rendered output, not a modified
version of the font software. The per-font breakdown, Reserved-Font-Name
analysis, and full license texts live with the fonts:

- **Index & rationale:** [`py/fonts/README.md`](py/fonts/README.md)
- **License texts:** [`py/fonts/licenses/`](py/fonts/licenses/)

| Font              | License                                  |
| ----------------- | ---------------------------------------- |
| Overpass          | OFL 1.1 + Apache-2.0                     |
| Archivo Black     | OFL 1.1                                  |
| B612 Mono         | OFL 1.1                                  |
| DejaVu Sans Mono  | Bitstream Vera + public-domain additions |
| Saira Stencil One | OFL 1.1                                  |

---

## System dependencies (not redistributed)

These are linked from the system at build time and are **not** bundled in
this repository, so their license texts are not reproduced here; they are
listed for completeness. Notably **CGAL** (GPLv3, optional, header-only)
backs the constrained Delaunay triangulation path when present.

See [`BUILDING.md`](BUILDING.md) for the full build-dependency list.

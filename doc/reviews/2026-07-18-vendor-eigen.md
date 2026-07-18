# Vendored Eigen review — keep 3.2.4

**Date:** 2026-07-18
**Scope:** `lib/fab/vendor/Eigen` (linear algebra, used by the mesher's feature-point QEF solves)
**Verdict:** **KEEP the vendored 3.2.4 copy as-is.** No functional need, no measurable perf upside, and a switch would almost certainly break every bit-exact golden for zero benefit.
**Mode:** read-only investigation. Nothing was modified.

---

## 1. What is vendored

| Property | Value |
|---|---|
| Version | **3.2.4** (2015) — `EIGEN_WORLD/MAJOR/MINOR = 3/2/4`, from `vendor/Eigen/src/Core/util/Macros.h` |
| Footprint | **330 files, 3.7 MB** |
| Provenance | committed whole in `af863f37 "Move Eigen into fab library"` (single import, never patched since) |
| License | MPL2 primary, with BSD/LGPL/GPL/MINPACK component notices (`COPYING.*`) — compatible with the GPLv3 repo |
| System alternative | `eigen 5.0.1-2` installed at `/usr/include/eigen3` (Arch package) |

### Build wiring
`lib/fab/CMakeLists.txt`:

```cmake
target_include_directories(SbFab SYSTEM PRIVATE vendor)      # line 131-132
target_include_directories(SbFabTest SYSTEM PRIVATE ../../vendor ...)  # line 159-160
```

- `vendor/` is added as an **`-isystem`** dir (`SYSTEM PRIVATE`), so `#include <Eigen/Dense>` resolves to `vendor/Eigen/Dense`, and Eigen's own template warnings are suppressed. Also means any *deprecation warnings* from a newer Eigen would be silently swallowed too.
- System Eigen lives at `/usr/include/eigen3`, which is **not** on the default compiler include path, so there is zero risk of accidentally picking up the system copy. The vendored one wins unambiguously.
- Relevant global flags (`CMakeLists.txt` lines 5-7): `-O3 -DRELEASE`, **no `-march`, no `-mavx`, no `-mfma`, no `-ffast-math`.** The target is baseline x86-64 (SSE2). `SbFab` additionally gets `-fno-trapping-math` (line 90). No `EIGEN_*` tuning macros are defined anywhere in `lib/fab/src` or `inc`.

---

## 2. Actual usage map

Eigen is included in exactly **two translation units**, both in the mesher's triangulate module. Everything else that "uses `Eigen::`" in the grep is vendored Eigen internals.

**Includes** (`grep 'include <Eigen'`, real code only):
- `lib/fab/src/tree/triangulate/delaunay.cpp:24-25` → `<Eigen/Dense>`, `<Eigen/SVD>`

**Type alias:**
- `lib/fab/inc/fab/tree/triangulate/triangle.h:8` → `typedef Eigen::Vector3d Vec3f;` (used throughout the mesher; note the name lies — it's `double`)

### 2a. `delaunay.cpp` — the hot QEF path (float / `Vector3f`)

The performance- and numerics-critical user. Adaptive-Delaunay stage-D feature solve, ~line 1429-1456:

- `Eigen::MatrixXf A(fc.n, 3)`, `Eigen::VectorXf b(fc.n)` — `fc.n` is a `uint8_t` (≤255 rows, typically a small fan), 3 columns.
- `Eigen::JacobiSVD<Eigen::MatrixXf> svd(A, ComputeThinU | ComputeThinV)`
- `svd.singularValues()`, `svd.matrixU()`, `svd.matrixV()`, `.transpose()`, `.asDiagonal()`
- **Manual clamped pseudo-inverse** (does *not* use `svd.solve`): `cutoff = sv(0) * 0.1f`; `inv(q) = 1/sv(q)` only where `sv(q) > cutoff`; then `x = centroid + V * inv.asDiagonal() * (Uᵀ b)`. The `0.1` relative singular-value cutoff is the QEF regularizer whose exact float behavior downstream referees depend on.
- Junction detector, ~line 1818-1823: `Eigen::Matrix3f m; ... Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> es(m); es.eigenvalues()` (3×3 covariance, ascending eigenvalues, ratio test `ev[1]/ev[2] > 0.2f`).
- `Eigen::Vector3f`, `::Zero()`, `.dot()`, `.row()`, scalar `/`.

### 2b. `mesher.cpp` — the DMC feature solve (double / `Vector3d`)

`lib/fab/src/tree/triangulate/mesher.cpp:414-437`:
- `Eigen::MatrixX3d A`, `Eigen::VectorXd B` (N rows, 3 cols; N = contour size, small).
- `Eigen::JacobiSVD<Eigen::MatrixX3d> svd(A, ComputeFullU | ComputeFullV)`
- For edge features: `svd.setThreshold(singular.minCoeff()/singular.maxCoeff() * 1.01)`
- `new_pt = svd.solve(B) + center` — here the SVD **threshold path is live** (unlike delaunay's manual clamp).

**Total real API surface:** `Matrix{Xf,X3d,3f}`, `Vector{Xf,Xd,3f,3d}`, `JacobiSVD` (thin+full, runtime option flags), `SelfAdjointEigenSolver`, plus `singularValues / matrixU / matrixV / solve / setThreshold / asDiagonal / dot / row / transpose / Zero`. Small and stable.

---

## 3. What an upgrade to 3.4.x / 5.x would actually change

### API compatibility — essentially a non-issue (verified against the installed 5.0.1 headers)
The entire surface above is stable public API from 3.2 through 5.x. The one theoretical break is the **runtime-options `JacobiSVD` constructor** (`JacobiSVD<Mat> svd(A, ComputeThinU|ComputeThinV)`), which 3.4+ deprecated in favor of a template `Options` parameter. Checked the actual system header:

- `/usr/include/eigen3/Eigen/src/SVD/JacobiSVD.h:593` — the runtime-options constructor **still exists** in 5.0.1.
- Line 591: its `EIGEN_DEPRECATED` marker is **commented out** — *"TODO: re-enable after fixing a few 3p libraries that error on deprecation warnings."*

So the current call sites would compile against 5.0.1 **without even a deprecation warning** (and `-isystem` would suppress one anyway). `SelfAdjointEigenSolver`, `setThreshold`, `solve`, etc. are unchanged. **API is not the blocker.**

### Numerical behavior — this is the whole problem
- `JacobiSVD` was refactored substantially between 3.2 and 5.x (SVD base reworked, `BDCSVD` added, `compute_impl` paths, helper functions like the 2×2 real preconditioner and `numext::hypot` changed). The *algorithm* (two-sided Jacobi) is the same and equally stable, but **bit-for-bit output is not guaranteed to match**, even at the same SSE2 baseline, because operation ordering and inlined helpers differ.
- No documented change *loosens or tightens* the default JacobiSVD threshold in a way that would alter results at the tolerances this code uses — `delaunay.cpp` bypasses the built-in threshold entirely (manual `0.1` clamp), and `mesher.cpp` sets its own. So the *logic* is stable; only the *low bits* are at risk.
- **FMA is the classic ULP-breaker** (`a*b+c` fused ≠ separate mul+add), but the current build has **no `-mfma`/`-march=native`**, so a newer Eigen at these flags would *not* newly enable FMA. That lowers, but does not eliminate, the drift risk — plain SSE2 kernel reordering can still move a bit.

### Performance — no meaningful upside for this workload
- Matrices are **N×3 with tiny N** (≤255, usually a handful). `JacobiSVD` on such matrices is dominated by scalar Jacobi rotations, not by vectorized GEMM.
- The AVX2/FMA/blocking improvements added in 3.3+ target **large dense products**. They do essentially nothing for 3-column solves, which fall through to scalar/near-scalar paths.
- Realistic expected speedup from upgrading Eigen here: **unmeasurable** against total mesher time. This is not where cycles go.

### Compile-time / footprint
- Included in only 2 TUs, so header bloat barely matters. Modern Eigen is a larger tree (~5-6 MB) — a repo-size cost, not a build-time one at this include count.

---

## 4. System Eigen switch — feasibility and risk

- **Installed:** `eigen 5.0.1-2` at `/usr/include/eigen3`. A `find_package(Eigen3)` switch is *mechanically trivial* — drop the `vendor` include for Eigen, add `find_package(Eigen3 REQUIRED)` + `Eigen3::Eigen`, and the call sites compile unchanged (per §3).
- **Why it's the wrong move for this campaign:**
  1. **Bit-exact goldens will trip.** `lib/fab/tests/mesher.cpp` dumps meshes via `canonical_dump` using `snprintf(..., "%a %a %a ...")` — **hexadecimal float, i.e. exact bit representation** (line 142). The dump is canonicalized (sorted), so topology/ordering churn is absorbed, but **any single-ULP change to a feature-point coordinate changes the `%a` string and fails the diff.** With a refactored SVD feeding those coordinates, a bit-identical result is unlikely. The `[golden]` cube/sphere dumps and the `dmesh` goldens are the exact tripwire.
  2. **The 627,666-assertion suite** is mostly tolerance/logic-based and would likely survive low-bit drift — but the goldens above are the sensitive gate, and they are maximally sensitive by design.
  3. **A system dep floats the numerics.** 5.0.1 today, 5.1/6.0 after the next `pacman -Syu` — each upgrade can silently re-perturb the SVD and re-break goldens on a machine that didn't change a line of Stibium. That is the *opposite* of what a bit-exact regression campaign wants. Vendoring **pins** the numerical baseline; the system package un-pins it.
  4. **No offsetting benefit:** no perf win (§3), no security surface (Eigen is compile-time header math over trusted in-process data, not a parser), no API need.

This lands exactly where the owner already sits on the sibling vendored lib: *not eager to introduce a system dependency for a marginal upgrade* — except here it isn't even marginal-upgrade-for-a-cost, it's cost-for-no-upgrade.

---

## 5. Recommendation

**Keep the vendored Eigen 3.2.4 as-is.** Ranked:

1. **KEEP vendored 3.2.4 (recommended).** It compiles clean under the current C++17 toolchain (build is green), the API surface is tiny and stable, there is no perf or security forcing function, and it pins the SVD numerics the golden-mesh provenance chain depends on. Age alone is not a defect for a header-only pinned math kernel.
2. **Upgrade vendored copy in-place to 3.4.x (only if a real forcing function appears)** — e.g. a future compiler rejects 3.2.4, or a concrete Eigen bug bites. If so: swap in a pinned 3.4.x tarball (not system), treat it as a deliberate **numerical-baseline-shift commit** that **regenerates every golden in the same commit** under fixed flags, reviewed as such. Keeps numerics pinned while modernizing.
3. **Switch to system Eigen 5.0.1 — not recommended.** Trivial to wire, but floats the numerics against the OS package, breaks bit-exact goldens for no benefit, and adds a system dependency the owner is explicitly reluctant to take on.

Do **not** do a casual "dependency bump." Any Eigen change here is a numerical event, not a version-hygiene chore.

---

## 6. Verification checklist (exact commands)

Confirm the findings, and — *if* an upgrade is ever attempted — measure the actual golden drift before committing.

```bash
cd /home/nate/code/stibium

# --- Confirm vendored version & footprint ---
grep -E 'EIGEN_(WORLD|MAJOR|MINOR)_VERSION' lib/fab/vendor/Eigen/src/Core/util/Macros.h
find lib/fab/vendor/Eigen -type f | wc -l          # -> 330
du -sh lib/fab/vendor/Eigen                          # -> 3.7M

# --- Confirm the real usage is only 2 TUs ---
grep -rn 'include <Eigen' lib/fab/src lib/fab/inc
grep -rn 'Eigen::' lib/fab/src/tree/triangulate/delaunay.cpp \
                   lib/fab/src/tree/triangulate/mesher.cpp \
                   lib/fab/inc/fab/tree/triangulate/triangle.h

# --- Confirm build wiring (-isystem vendor) ---
grep -n 'SYSTEM PRIVATE\|vendor' lib/fab/CMakeLists.txt

# --- Confirm system Eigen availability & version ---
pacman -Q eigen
grep -E 'EIGEN_(WORLD|MAJOR|MINOR)_VERSION' /usr/include/eigen3/Eigen/src/Core/util/Macros.h

# --- Confirm the goldens are bit-exact hex floats ---
grep -n '%a' lib/fab/tests/mesher.cpp

# ============================================================
# ONLY if trialing an upgrade — measure golden drift safely:
# ============================================================
# 1. On the current tree, regenerate the reference goldens and stash them:
#    (build first per doc; plain ninja is safe)
ninja -C build SbFabTest
(cd build && ./lib/fab/SbFabTest "[golden]")        # writes golden_*.txt to CWD
mkdir -p /tmp/eigen-baseline && cp build/golden_*.txt /tmp/eigen-baseline/

# 2. Point the two TUs at the candidate Eigen (in-place vendored swap or
#    a temporary -isystem to /usr/include/eigen3), rebuild, regenerate:
ninja -C build SbFabTest
(cd build && ./lib/fab/SbFabTest "[golden]")

# 3. Diff — ANY output means the SVD moved bits; goldens must be
#    consciously re-baselined, not silently accepted:
for f in build/golden_*.txt; do diff -u /tmp/eigen-baseline/$(basename $f) $f; done

# 4. Full gate (627,666 assertions) must stay green regardless:
(cd build && ./lib/fab/SbFabTest)
```

**Interpretation:** empty diffs in step 3 = bit-identical, upgrade is numerically free (unlikely but possible). Any diff = a real numerical baseline shift that has to be reviewed and the goldens regenerated in the same commit — never a quiet bump.

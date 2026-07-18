# Vendored Catch review — 2026-07-18

Read-only investigation. Scope: the vendored Catch test framework, its usage
across the suite, and whether/how to modernize it. A separate worker should be
able to execute the migration from this doc.

## TL;DR

- Vendored **CATCH v1.1 build 3, generated 2015-05-21** (single header,
  `vendor/catch/catch.hpp`, 9,460 lines / 330 KB). Added to the repo 2015-07-05.
  ~11 years old.
- It works and produces no warnings (it's marked `pragma GCC system_header`), so
  it is not *broken*. The cost is **compile time**: the framework re-parses on
  every test TU and the implementation recompiles in each `main.cpp`.
- **Recommendation: upgrade to Catch2 v3** (Arch system package `catch2`,
  currently 3.15.2). It preserves the Catch tag semantics **exactly** —
  including the `[.hidden]` convention and the comma-OR filter that the whole
  regression harness (`SbFabTest "[.dmesh],[.dtrace],[.dchain]"`) depends on.
- **doctest is NOT drop-in for this codebase** despite its reputation. It does
  not implement Catch's `[.tag]`-hides-by-default convention or `[tag]`-in-name
  filtering, and its `WARN` macro has different semantics from the 37 `WARN(...)`
  message-stream call sites here. Migrating to doctest would mean rewriting the
  hidden-tag convention *and* every documented invocation string. Rejected on the
  "preserve exact invocation patterns" constraint.

## 1. What's vendored

| Fact | Value |
|---|---|
| Path | `vendor/catch/catch.hpp` (+ `vendor/catch/LICENSE_1_0.txt`) |
| Version | `CATCH v1.1 build 3 (master branch)`, generated `2015-05-21` |
| Size | 330 KB, 9,460 lines, single header |
| Added | commit `119a55ab` "Add Catch unit testing library", 2015-07-05 |
| License | Boost Software License 1.0 (compatible with the GPLv3 fork) |
| Style | header-only; implementation gated behind `CATCH_CONFIG_MAIN`/`CATCH_CONFIG_RUNNER` |

Include path is wired via `target_include_directories(... SYSTEM PRIVATE
../../vendor)` in both test targets, so tests include it as
`<catch/catch.hpp>` (the `catch/` subdir under `vendor/`).

## 2. Usage map

**Test targets (2), both custom-main:**

- `SbFabTest` — `lib/fab/CMakeLists.txt:136`. Sources: `analytics, contour,
  grid, main, mesh, mesher, parser, aa, dmesh, shape, tape` (+ `gpu` when
  `OpenGL_EGL_FOUND`). Links `SbFab`. `add_test(NAME SbFabTest ...)`. C++17.
- `SbGraphTest` — `lib/graph/CMakeLists.txt:29`. Sources: `datum, graph, link,
  main, node, script, subgraph`. Links `SbGraph`. C++17.

Binaries: `build/lib/fab/SbFabTest` (118 MB), `build/lib/graph/SbGraphTest`
(23 MB). No PCH anywhere in the project (checked — no `target_precompile_headers`).

**Both `main.cpp` files** use the runner form:

```cpp
#define CATCH_CONFIG_RUNNER
#include <catch/catch.hpp>
...
int main(int argc, char** argv) {
    /* project preInit + Py_Initialize */
    return Catch::Session().run(argc, argv);
}
```

**Macro / API inventory** (19 test TUs, counted):

| API | Count | v3 status |
|---|---|---|
| `REQUIRE` | 522 | identical |
| `TEST_CASE` | 137 | identical |
| `CAPTURE` | 128 | identical |
| `SECTION` | 46 | identical |
| `Approx` | 38 | **moved to `Catch::Approx`** (must qualify) |
| `WARN` | 37 | identical (message-stream form) |
| `CHECK` | 32 | identical |
| `CATCH_CONFIG_RUNNER` | 2 | identical (with new include) |
| `REQUIRE_FALSE` | 1 | identical |

No version-1-specific internals are used (no reporter subclassing, no
`Catch::Detail`, no custom matchers, no `ConsoleReporter` tweaks). `WARN` is
used as `WARN("msg " << expr)` — the Catch message-stream form, which Catch2 v3
keeps verbatim. `Approx` is used **unqualified** (`Approx(x).epsilon(...)`).

**Include forms in the wild:** mix of `<catch/catch.hpp>` (16 files) and
`"catch/catch.hpp"` (`gpu.cpp`, `dmesh.cpp`, `aa.cpp`). Both resolve via the
`vendor` include dir.

**Hidden-tag convention (the load-bearing part).** Tests use Catch's
leading-dot "hidden" tags so they don't run by default and are opt-in by name:
`[.aa] [.aafuzz] [.bench] [.dchain] [.dmesh] [.dmeshBC] [.dmeshSTL] [.dmeshVS]
[.dtrace] [.gpu] [.gpubench] [.gpubig] [.gpudbg] [.gpudiv] [.gpudump] [.gpulib]
[.gpu2] [.gpu2bench]` plus visible tags `[grid] [mesh] [tape] [fuzzer] [golden]`.
The harness runs e.g. `SbFabTest "[.dmesh],[.dtrace],[.dchain]"` (comma = OR).
This convention is referenced throughout `doc/*.md`, `CHANGELOG.md`, `TODO.md`,
and the project memory as the regression gate.

## 3. Options assessed

### (a) Keep as-is
- **Not broken.** No compiler warnings (`system_header` suppresses them), no
  C++17 incompatibility observed — it compiles clean with GCC under `-std=c++17`.
- **The only real cost is compile time** (measured, this machine, `-O0`):
  - Framework-only TU with one trivial test (what each of the ~17 non-main test
    files pays on every rebuild that touches them): **~0.64 s** of pure Catch
    parse overhead per TU.
  - `main.cpp` with `CATCH_CONFIG_RUNNER` (compiles the whole framework
    implementation): **~3.25 s**. This recompiles whenever `main.cpp` or its
    deps change, and at the campaign's real `-O2`+flags it is materially worse.
  - Preprocessed size of a bare `#include "catch.hpp"`: **89,781 lines**.
  - Net: order **~15 s** of Catch-attributable parse/compile per clean test
    build, and a fixed ~0.64 s tax on every incremental single-file test rebuild.
- Verdict: acceptable if untouched, but this is the binary Nate rebuilds
  constantly during mesher work, so the fixed per-TU tax is paid often.

### (b) Upgrade to Catch2 v3 — RECOMMENDED
- System package present and current: `catch2 3.15.2-1` in `extra`
  (BSL-1.0, header + compiled lib). Not yet installed; one `pacman -S catch2`.
- **Preserves tag semantics exactly** — `[.hidden]`, `[tag]`-in-name filtering,
  comma-OR command line. The invocation strings in every doc and the harness do
  **not** change.
- **Compile-time win:** v3 is no longer a monolithic single header. Each test TU
  includes only `<catch2/catch_test_macros.hpp>` (far smaller than the 330 KB
  v1 header), and the framework implementation is a **prebuilt system library**
  (`libCatch2.a`/`Catch2::Catch2`) — the ~3.25 s impl compile in each `main.cpp`
  disappears entirely (linked, not compiled).
- **Breakage is small and mechanical** (see §4): include lines, `Approx` →
  `Catch::Approx`, the two `main.cpp` runner includes, and CMake link lines.
  `WARN`, `TEST_CASE`, `SECTION`, `REQUIRE`, `CHECK`, `CAPTURE` are unchanged.

### (c) Migrate to doctest — REJECTED for this repo
- `doctest 2.5.2-1` is available in `extra` (MIT, single header) and compiles
  ~10× faster than even Catch2 v3, which is attractive for a constantly-rebuilt
  binary.
- **But it breaks the two things the constraints forbid changing:**
  1. **Hidden-tag convention.** doctest has no `[.tag]`-hides-by-default rule and
     does not parse `[tag]` substrings inside the test name for filtering (it
     uses `--test-case=<name-glob>`, `--subcase=`, and `* doctest::skip()` /
     `* doctest::test_suite("x")` decorators). Every `[.dmesh]`-style tag and
     every `SbFabTest "[...]"` invocation string in the docs/harness would have
     to be redesigned and rewritten.
  2. **`WARN` semantics.** In Catch, `WARN(msg << expr)` streams a message. In
     doctest, `WARN(expr)` is a *soft assertion* level. The 37 call sites here
     are message streams and would silently change meaning; they'd need to
     become `MESSAGE(...)`.
  - Also: `SECTION` → `SUBCASE` (46 sites), `Approx` → `doctest::Approx`, and a
    different custom-main API (`DOCTEST_CONFIG_IMPLEMENT` + `doctest::Context`).
- Net: doctest is a bigger, riskier migration *here specifically*, and it
  violates "preserve the hidden-tag invocation patterns." Revisit only if the
  team decides to redesign the tag convention wholesale.

### (d) Other options
- **FetchContent Catch2 v3** instead of the system package — pins an exact
  version, survives a machine without `catch2` installed, but adds a
  configure-time download and its own build. Reasonable fallback if reproducible
  builds across machines matter more than using the distro package. Same source
  changes as (b).
- **Bump the vendored single header to Catch2 v2's amalgamated
  `catch_amalgamated`** — keeps the vendored/offline model and the exact tag
  semantics, smaller source delta than v3 (v2 still ships a single header and
  keeps `Approx` unqualified). Downside: v2 is EOL upstream and you keep paying
  the in-tree impl-compile cost. Only pick this if staying vendored/offline is a
  hard requirement.

## 4. Migration plan — Catch2 v3 (recommended)

Scope: **21 files** touch Catch (19 test TUs + 2 `main.cpp`, which overlap the
list) plus **2 CMakeLists**. All changes are mechanical; ~5 minutes of hand work
concentrated in the two `main.cpp` and the `Approx` qualification.

**Step 0 — install / provide the dep**
```
sudo pacman -S catch2      # provides Catch2::Catch2 + Catch2::Catch2WithMain
```

**Step 1 — includes (sed-able across all test files)**
Replace the framework include. Non-main files:
```
# in every lib/{fab,graph}/tests/*.cpp EXCEPT main.cpp
#include <catch/catch.hpp>   ->   #include <catch2/catch_test_macros.hpp>
"catch/catch.hpp"            ->   <catch2/catch_test_macros.hpp>
```
Any TU that uses `Approx` also needs `#include <catch2/catch_approx.hpp>`
(the `Approx` files are `contour.cpp` and any other of the 38-site set — grep
`Approx` to get the exact list; add the include only where used).

**Step 2 — `Approx` qualification (38 sites)**
v3 moved `Approx` into the `Catch` namespace and it must be qualified. Two clean
options:
- Per-file: add `using Catch::Approx;` after the includes in each TU that uses
  it (smallest diff, no call-site churn), **or**
- Sed `Approx(` → `Catch::Approx(` at the 38 call sites.
Prefer the `using` — it's one line per file and leaves the assertions untouched.

**Step 3 — the two `main.cpp` runners (hand edit)**
```cpp
#define CATCH_CONFIG_RUNNER          // DELETE this line
#include <catch/catch.hpp>           // ->
#include <catch2/catch_session.hpp>
...
return Catch::Session().run(argc, argv);   // UNCHANGED
```
`Catch::Session().run(argc,argv)` is identical in v3. Only the include and the
now-unneeded `CATCH_CONFIG_RUNNER` define change. (Because we keep our own
`main` with `preInit`/`Py_Initialize`, link `Catch2::Catch2`, **not**
`Catch2WithMain`.)

**Step 4 — CMake (2 files)**
`lib/fab/CMakeLists.txt` and `lib/graph/CMakeLists.txt`:
```cmake
find_package(Catch2 3 REQUIRED)                 # near the top of each, or in root
...
target_link_libraries(SbFabTest   SbFab   Catch2::Catch2)
target_link_libraries(SbGraphTest SbGraph Catch2::Catch2)
```
Remove `vendor` from the test targets' include dirs **only if** nothing else in
those targets needs it (`SbFabTest` also lists `../../vendor` and Boost/Python —
keep the dir, just stop relying on it for catch). Leaving the `vendor/catch`
header in place is harmless; deleting `vendor/catch/` is optional cleanup once
the build is green.

**Tag / invocation mapping — NO CHANGE.** `[.dmesh]`, `[.dtrace]`, `[.dchain]`,
`[.gpu]`, `[golden]`, `[tape]`, `[fuzzer]`, comma-OR, leading-dot-hides: all
identical in Catch2 v3. `SbFabTest "[.dmesh],[.dtrace],[.dchain]"` runs exactly
as before. This is the whole reason v3 wins over doctest.

**What needs human eyes (not sed):**
- The two `main.cpp` runner edits (Step 3).
- Confirming which TUs need `catch_approx.hpp` (grep `Approx`, add include there).
- One full compile to surface any incidental v1→v3 include-transitivity gaps
  (v1's single header dragged in `<algorithm>`/`<vector>`/etc. that a TU may
  have been leaning on implicitly; v3's granular headers don't). Fix = add the
  missing standard include to the offending TU. Expect 0–3 of these.

## 5. Risks

- **Implicit-include regressions** (medium likelihood, low effort): a test TU
  that relied on the fat v1 header pulling in a std header may fail to compile
  under granular v3 headers. Fix per-TU. Caught by the first build.
- **`Approx` epsilon/margin defaults** (low): v1 and v3 share the same relative
  epsilon default; the 38 sites all set `.epsilon(...)` explicitly anyway, so
  numeric tolerance behavior is preserved. Verify by running the suite and
  confirming assertion count is unchanged.
- **Assertion-count drift** (low but must be checked): the gate is *exactly*
  627,666 assertions / 55 test cases. A macro that changed evaluation (it
  shouldn't — none of ours did) would move the count. Verify explicitly (§6).
- **Reporter/output format**: v3's console reporter formatting differs slightly
  from v1. If any script greps SbFabTest stdout for a specific format, check it.
  (Golden dumps are written to files, not parsed from stdout, so `[golden]` is
  unaffected.)
- **System-package availability**: relying on `pacman catch2` means a fresh
  machine needs the package. If that's a problem, use option (d) FetchContent.

## 6. Verification checklist (exact commands)

Run from repo root. Anchor against the current baseline **before** touching
anything (capture the golden dumps and assertion counts on HEAD first).

```sh
# 0. Baseline on the CURRENT (v1) build, for comparison
cd /home/nate/code/stibium
ninja -C build SbFabTest SbGraphTest
build/lib/fab/SbFabTest  | tail -3      # note "N assertions in M test cases"
build/lib/graph/SbGraphTest | tail -3

# --- apply the migration ---

# 1. Configure + build must be clean
sudo pacman -S --needed catch2
cmake --build build --target SbFabTest SbGraphTest   # or: ninja -C build
#   (if using a fresh configure) cmake -B build -G Ninja && ninja -C build

# 2. Default (non-hidden) suite: assertion + case counts must MATCH baseline
build/lib/fab/SbFabTest        # expect the SAME "assertions in N test cases"
build/lib/graph/SbGraphTest    # as the v1 baseline captured in step 0

# 3. Regression gate — hidden-tag invocations must run and pass unchanged
build/lib/fab/SbFabTest "[.dmesh],[.dtrace],[.dchain]"
build/lib/fab/SbFabTest "[golden]"        # then cmp dumps vs canonical
build/lib/fab/SbFabTest "[tape]"
build/lib/fab/SbFabTest "[.aa]"
# (GPU tags only if built with EGL: build/lib/fab/SbFabTest "[.gpu]")

# 4. Golden determinism: dumps must be byte-identical to pre-migration
#    (run [golden] in a scratch dir before and after, then:)
cmp golden_cube_plain.txt   <baseline copy>
cmp golden_sphere_detect.txt <baseline copy>

# 5. ctest wiring still discovers both targets
ctest --test-dir build -N | grep -E 'SbFabTest|SbGraphTest'

# 6. Compile-time sanity (optional, confirms the win): time a clean rebuild
#    of SbFabTest before/after — main.cpp should no longer compile the framework.
```

**Pass criteria:** identical assertion/test-case counts in step 2, all hidden-tag
runs green in step 3, byte-identical goldens in step 4. Any count delta is a
stop-and-investigate — the migration must be assertion-count-neutral.

## 7. Other decade-old vendored code in the tree

(meshoptimizer excluded per instructions — it's `app/vendor/meshoptimizer`,
v1.0, and actively used/recently touched 2026-07-12.)

- **`lib/fab/vendor/Eigen` — Eigen 3.2.4 (2015), 3.7 MB.** ~11 years old, same
  vintage as Catch. System Arch package is **eigen 5.0.1**. This is the most
  significant stale vendored dep: it's a core numerical dependency of `SbFab`
  (not just tests), so an upgrade is higher-value but also higher-risk (API and
  numerical-behavior surface is much larger than a test framework). Worth a
  dedicated review — do **not** bundle it with the Catch change. Note Vec3f is
  `Eigen::Vector3d` (project memory); any Eigen bump must preserve exact
  numerical/meshing determinism (the golden dumps are the referee).
- No other vendored trees found. `vendor/` holds only `catch/`; `app/vendor/`
  only `meshoptimizer/`; `lib/fab/vendor/` only `Eigen/`. (`app/undo/undo_catcher.h`
  matched a "catch" grep but is unrelated project code — an undo-stack catcher.)

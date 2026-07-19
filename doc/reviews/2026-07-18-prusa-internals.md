# PrusaSlicer internals: how `--info` computes mesh health, and what it means for our pinch splits

**Date:** 2026-07-18
**Source:** PrusaSlicer `version_2.9.4` (installed CLI: `PrusaSlicer-2.9.4+arch8`), cloned to
`/tmp/PrusaSlicer`. All code citations are file:line in that tree.
**Referee:** `prusa-slicer --info` on freshly-exported stibnite meshes (screws + bino, r1).
**Rule of the house:** measured, not predicted. Every number is an export + a referee pass.

---

## 0. Executive answer

There are **two completely different stat engines** in PrusaSlicer, chosen by file format:

- **STL** → `admesh` exact + nearby **position weld**, then an orientation flood-fill.
  Coincident duplicate vertices are re-fused by *bit-exact float32 position*, so our
  index-space pinch split is **100% invisible on STL** (measured: split STL ≡ plain STL,
  identical stats). `facets_reversed` is a flood-fill artifact; `open_edges` is inflated by
  greedy 4-fan pairing.
- **3MF** → index-preserving `its` builder, **no weld, no flood**. Our splits survive.
  `facets_reversed` is **structurally 0** (the pass that mints it never runs).

**The artifact theory is confirmed** for `facets_reversed` (an admesh flood-fill artifact:
150 → 0 the instant you switch to 3MF) and **partially** for `open_edges` (292 → 73 on bino
for the same geometry). But the residual is **not** all artifact: our indexed mesh genuinely
carries 274 non-manifold edges, and *any* honest slicer topology test trips on them. Naive
pinch-split does not fix this — it trades 274 non-manifold edges for **2202 open edges**
(measured), because it un-shares the seams without re-closing the sheets.

---

## 1. The STL path (admesh) — exact-position weld + flood-fill

`load_stl` → `TriangleMesh::ReadSTLFile(path, repair=true)` (default repair)
→ `trianglemesh_repair_on_import(stl)`.
Files: `src/libslic3r/Format/STL.cpp:29`, `src/libslic3r/TriangleMesh.cpp:218`,
`TriangleMesh.cpp:98`.

### 1a. Exact weld — the key mechanism (`bundled_deps/admesh/admesh/connect.cpp`)

`stl_check_facets_exact` (`connect.cpp:432`) builds edge connectivity by hashing each
half-edge on the **raw float32 bits of its two endpoint positions**:

- `HashEdge::load_exact` (`connect.cpp:54`): `memcpy` the 6 floats (two `stl_vertex`) into the
  hash key; the only normalization is **negative-zero → positive-zero** (`connect.cpp:72-84`).
  Match is `memcmp` of all 24 bytes (`connect.cpp:43`). **This is a bit-exact position match,
  not index and not tolerance.** Two vertices at the same position but different indices hash
  identically and pair.
- Pairing is **greedy and order-dependent** (`insert_edge`, `connect.cpp:177`): first
  half-edge with a given key is stored; the next matching one is paired and the stored entry is
  *deleted*. So a 4-incident edge (4 half-edges, identical key) pairs as (f1,f2) then (f3,f4) —
  in insertion order, i.e. **geometrically arbitrary**. Degenerate facets (two equal vertices)
  are removed first (`connect.cpp:444-455`).

**Consequence for us:** our pinch split mints coincident duplicate vertices (one per sheet) at
*bit-identical* positions. `load_exact` re-fuses them by position. The split cannot survive an
STL round-trip. **Measured:** `bino_split.stl` and `bino_plain.stl` produce *byte-identical*
prusa stats (292 / 150), and `screws_split.stl` ≡ `screws_plain.stl` (manifold=yes / 0).

### 1b. Nearby weld — the tolerance ladder (`connect.cpp:481`)

Runs only if `connected_facets_3_edge < number_of_facets` after exact (i.e. some edge is still
unpaired), from `trianglemesh_repair_on_import` (`TriangleMesh.cpp:119-138`):

- `tolerance = shortest_edge` (the **smallest Chebyshev edge length in the mesh**, tracked in
  `load_exact`, `connect.cpp:58-59`); `increment = bounding_diameter / 10000`; **2 iterations**
  (`TriangleMesh.cpp:119-134`).
- `HashEdge::load_nearby` (`connect.cpp:87`) snaps each endpoint to a grid cell of size
  `tolerance` and matches by cell; on a match it **moves** vertices together
  (`match_neighbors_nearby` → `change_vertices`, `connect.cpp:277-426`) — destructive.

For a diagonal ~100 mm, `increment ≈ 0.01 mm`, so tolerance can climb to ~`shortest_edge +
0.02 mm` by iteration 2. This is the vice for any epsilon strategy (§4.3).

### 1c. Orientation flood-fill — where `facets_reversed` comes from (`normals.cpp`)

`stl_fix_normal_directions` (`normals.cpp:116`) floods orientation from **facet 0** across the
(position-welded, greedy-paired) neighbor graph. Any neighbor whose `which_vertex_not > 2`
(wound opposite across the shared edge) is flipped by `reverse_facet`, which does
`++stl->stats.facets_reversed` (`normals.cpp:36`). On a contradiction (a Möbius loop from
mis-pairing) it reverts the whole part (`normals.cpp:162-166`). `stl_calculate_volume`
(`util.cpp:290`) additionally flips *all* facets if total volume < 0 (not counted).

**So `facets_reversed` counts flips the flood-fill needed to make the greedy-paired graph
consistent — it is downstream of the arbitrary 4-fan pairing, not a property of our oriented
mesh.** Proof: it is **0** on the index path (§2), which does no flood-fill.

### 1d. `backwards_edges` and the `open_edges` formula

`stl_verify_neighbors` (`util.cpp`, near `backwards_edges`): for each connected edge, if the
neighbor's `which_vertex_not >= 3` (paired but reversed winding), `++backwards_edges`.
Then (`TriangleMesh.cpp:232-235`):

```
facets_w_1_bad_edge = connected_facets_2_edge - connected_facets_3_edge
facets_w_2_bad_edge = connected_facets_1_edge - connected_facets_2_edge
facets_w_3_bad_edge = number_of_facets       - connected_facets_1_edge
open_edges = backwards_edges + w1 + 2*w2 + 3*w3
```

`manifold()` ⇔ `open_edges == 0` (`TriangleMesh.hpp:130`). Because even-incidence fans fully
pair (no "bad edges"), the whole STL `open_edges` on our meshes is essentially `backwards_edges`
— i.e. mis-oriented greedy pairings. This is why bino reads `open == backwards` in the autopsy.

### 1e. 4+-incident edges, summarized

admesh has **no concept of a >2-fan**. It greedily pairs half-edges by exact position in
arbitrary order, then the flood-fill and `verify_neighbors` report the fallout as
`facets_reversed` / `backwards_edges` / `open_edges`. There is no "non-manifold edge" counter;
the non-manifoldness is laundered into those three numbers.

---

## 2. The 3MF path (`its`) — index-preserving, no weld

3MF import builds an `indexed_triangle_set` straight from the file's `<vertex>`/`<triangle>`
arrays and constructs `TriangleMesh(std::move(its), mesh_stats)`
(`src/libslic3r/Format/3mf.cpp:2569-2605`). That constructor (`TriangleMesh.cpp:90`) calls
`fill_initial_stats` (`TriangleMesh.cpp:64`), which computes:

```
number_of_parts = its_number_of_patches(its, face_neighbors)   // index-edge connectivity
open_edges      = its_num_open_edges(face_neighbors)            // count of neighbor slots == -1
```

`repaired_errors` (which holds `facets_reversed`, `backwards_edges`) is taken **verbatim from
the 3MF metadata** and is **0** for a fresh stibnite export. **No admesh, no weld, no
flood-fill runs.** There is no `merge_vertices` / `require_manifold` / repair call anywhere in
the 3MF import path (grepped).

### 2a. The neighbor index is keyed on vertex INDICES (`MeshSplitImpl.hpp:221`)

`create_face_neighbors_index`: for each face edge `(a,b)` (vertex **indices**), it scans other
faces incident to index `a`, matches one whose third vertex gives edge `(b,a)` (opposite
winding), assigns it, and **`break`s** — at most one neighbor per edge
(`MeshSplitImpl.hpp:246-259`). `its_num_open_edges` then counts unpaired slots
(`TriangleMesh.cpp:1633`). `its_number_of_patches` flood-fills across index-shared edges
(`MeshSplitImpl.hpp:205`; header note at `TriangleMesh.hpp:275`: "shared edge defined with 2
shared vertex indices").

**Consequence:** coincident-but-distinct-index vertices are **distinct** here. Our pinch split
survives. But a 4-fan that still shares the *same two indices* is left with 2 unpaired
half-edges (only one pair made) → counted as open, and the fan acts as a **patch barrier**
(inflating `number_of_parts`). The comment at `TriangleMesh.hpp:256` says it outright: *"This
function will happily create non-manifolds if more than two faces share the same vertex position
or edge."*

---

## 3. The empirical matrix (measured, r1)

`prusa-slicer --info`. `reversed` blank = not printed = 0. Both bino 3MFs are geometrically
identical (122 858 distinct vertex positions); the split only adds duplicate-index vertices.

| model / file            | manifold | open_edges | facets_reversed | parts | idx nonmanifold edges | idx boundary edges |
|-------------------------|:--------:|:----------:|:---------------:|:-----:|:---------------------:|:------------------:|
| **bino** plain **STL**  | no       | 292        | **150**         | 0     | — (welded away)       | —                  |
| bino split **STL**      | no       | 292        | 150             | 0     | — (welded away)       | —                  |
| **bino** plain **3MF**  | no       | **73**     | **0**           | 227   | 274 (272×4-way, 2×6)  | 0                  |
| bino split **3MF**      | no       | **2202**   | 0               | 746   | **0**                 | **2202**           |
| **screws** plain **STL**| **yes**  | 0          | 0               | 1     | — (welded away)       | —                  |
| screws split **STL**    | yes      | 0          | 0               | 1     | — (welded away)       | —                  |
| screws plain **3MF**    | no       | 8          | 0               | 2     | 1                     | 0                  |
| screws split **3MF**    | no       | 8          | 0               | 2     | 0                     | 8                  |

Reads straight off the grid:

1. **`facets_reversed` is a pure STL/admesh artifact.** 150 on STL → 0 on 3MF for identical
   geometry. It is the orientation flood-fill (§1c) reconciling arbitrary 4-fan pairings.
   Confirms and closes the MESH-WAR open question ("the 14→150 reversed delta"): it's neither
   the seal nor decimation *per se* — it is admesh's flood-fill on a non-manifold input, and it
   vanishes the moment prusa is not asked to weld.
2. **Our split is invisible on STL.** split STL ≡ plain STL, bit-for-bit, on both models.
   Exact-position weld (§1a) re-fuses the coincident copies. This is the mechanism behind the
   TIER-A note "coincident copies make prusa WORSE (it re-welds them)".
3. **3MF cuts the artifact but exposes real non-manifoldness.** bino open 292 → 73 and reversed
   150 → 0, but the 274 genuine index-nonmanifold edges remain (as 73 open + 227 parts).
4. **Naive pinch-split is worse in 3MF, and we now know exactly why.** It drives idx-nonmanifold
   274 → 0 but mints **2202 boundary (open) edges** (the "~2198 geometric edges" of the TIER-A
   note, reproduced here as an exact `open==boundary` count). The split un-shares each seam edge
   per sheet but the sheets **do not re-close** along the seam → every split half-edge whose
   partner went to the other sheet copy becomes a boundary. `incidence_hist` proves it: plain
   `{2: 192284, 4: 272, 6: 2}` → split `{1: 2202, 2: 191733}`. Non-manifold became boundary,
   one-for-one.
5. **STL can *accidentally* win (screws).** One 4-fan that greedy-pairs *consistently* yields
   manifold=yes on STL, while 3MF's at-most-one-neighbor leaves it 8-open. So 3MF is not a
   universal improvement — it removes the reversed artifact and helps tangled models, but it
   never opportunistically 2-manifolds a same-index fan the way a lucky greedy weld can.

---

## 4. Verdict and the cheapest path to `manifold=yes / 0 / 0`

### 4.1 Is the artifact theory correct?

**Yes for `facets_reversed`** (100% STL flood-fill artifact — 0 on 3MF). **Yes for the bulk of
STL `open_edges`** (292 → 73 is weld/greedy-pairing inflation). **No for the residual:** the
plain-3MF 73 open / 274 index-nonmanifold edges are *real* — our indexed mesh is genuinely
non-manifold at the self-contact seams (field oracle in the bino autopsy already confirmed those
seams are real geometry, not tessellation phantoms). prusa's index test is an honest 2-manifold
test and it is right to fail our mesh. The correctness guard (weld fan-claim) is defensible; the
"prusa is lying" framing is only half true — it lies about *reversed*, it tells the truth about
*non-manifold*.

### 4.2 Why `manifold=yes/0/0` is a topology problem, not a format/epsilon problem

prusa's 3MF manifold test only requires **every index-edge to have exactly one opposite-wound
partner**. That is reachable *without moving any vertex* if — and only if — each 4-fan is split
into two edges whose two faces are a consistently-wound, **locally-closed** pair, walked
consistently along the whole seam chain. The current PINCHSPLIT does the duplication but not the
consistent walk, so chains leave dangling boundaries (2202). This is precisely the
**"seam-consistent sheet labeling"** MESH-WAR already named as the missing design.

The autopsy's split of the seam population matters here: **149 doubled-sheet** edges (two near-
parallel walls kissing) resolve to two closed manifolds cleanly; **167 genuine junctions**
(true 4-way self-contact) can be given a consistent 2-face pairing too (prusa's test is
local/topological and does not check self-intersection), *provided the chain walk closes*. So
0-open is reachable in principle for both classes via a correct split — the blocker is the
consistent-labeling walk, not geometry.

### 4.3 What defeats the STL weld, and why epsilon is the wrong tool

- **Exact weld** is bit-exact float32. Any separation larger than a float32 ULP at the
  coordinate magnitude (~6e-6 mm at 50 mm) defeats it.
- **But nearby weld** (§1b) can re-fuse: its grid tolerance climbs to ≈ `shortest_edge +
  bounding_diameter/10000` ≈ `shortest_edge + 0.02 mm` over its 2 iterations. Our real walls are
  0.012–0.5 mm (autopsy §2). **An epsilon large enough to reliably beat nearby (>~0.02 mm) is
  large enough to fuse real 0.012 mm walls.** Epsilon-separation on STL is caught in a vice.
- **And epsilon inherits the same open-seam problem as the naive split** (separating the sheets
  opens the seam boundary unless the sheets independently close), *plus* it perturbs geometry off
  the exact isosurface. It is strictly dominated by an index-space split exported to 3MF, which
  preserves the split with **zero** geometric perturbation and **no** weld to fight.

### 4.4 Ranked recommendation

1. **Ship 3MF, not STL, as the slicer-facing format — free, do it now.** For identical geometry
   it eliminates the `facets_reversed` artifact (150 → 0) and roughly quarters `open_edges`
   (292 → 73 on bino) with zero geometry change and no weld to undo our work. It converts the
   scary headline ("no / 292 open / 150 reversed") into an honest small residual
   ("no / 73 open / 0 reversed / 227 parts"). Caveat (screws): keep an STL fallback for meshes
   whose only defect is a single consistently-weldable fan, where STL can read manifold=yes.
2. **Make PINCHSPLIT seam-closure-consistent, then export via 3MF — the path to 0-open.** The
   current split is topologically incomplete (274 nm → 2202 open). Fix the seam-chain walk so
   each sheet re-closes (consistent A/B labeling along each of the 251 chains; the 226 isolated
   single edges are trivial local edits, the 31-edge chain needs the walk). Do it on the `its`
   representation — 3MF preserves it; STL would weld it away. Target: idx-nonmanifold 0 **and**
   idx-boundary 0 → prusa manifold=yes / 0 / 0. This is autopsy "Cure 2" with the closure bug
   named.
3. **Decouple `facets_reversed` worry from the roadmap.** Since it is 0 on 3MF, the decimation
   winding-veto (autopsy Cure 1) is only needed if we keep shipping STL. Under recommendation 1
   it becomes moot for the slicer referee (still worth it for internal razor hygiene, separate
   metric).
4. **Do NOT pursue epsilon separation as the primary route** (§4.3): vice'd against nearby-weld,
   perturbs geometry, and still opens seams. Only consider a *local* epsilon (< 0.1·sp, gap <
   0.1 mm) for the ≤28 sub-feature seam edges as a fallback where a clean split walk is awkward
   — matching the autopsy's local-weld carve-out.
5. **Re-evaluate the weld fan-claim once export is 3MF + correct split.** It is a legitimate
   correctness guard, but its cost was denominated in STL manifold-ness, which stops mattering
   under rec. 1. Check whether it is still load-bearing after the split walk lands.

### 4.5 One caveat to carry forward

prusa's 3MF manifold test is *topological only* — it does not test self-intersection. A
consistent split that makes genuine junctions read manifold=yes still leaves two surfaces
touching at a curve (a real self-contact). That is fine for FDM slicing (the slice contours are
still well-defined) and is what the guarantees-ladder RUNG 2 (tri-tri sweep) would separately
certify. Do not conflate "prusa manifold=yes" with "no self-contact"; they are different claims.

---

## 5. Reproduction

```
# exports (r1); PINCHSPLIT=1 enables our split pass; .3mf suffix picks 3MF writer
STIBIUM_EXPORT_DMESH=1 [STIBIUM_DMESH_PINCHSPLIT=1] \
  build/app/stibium --export /tmp/X.{stl,3mf} --resolution 1 examples/torture/zeiss_id02_bino.sb
prusa-slicer --info /tmp/X.{stl,3mf} | grep -E 'manifold|open_edges|facets_reversed|number_of_parts'
```
Index-level analysis: unzip the `.3mf`, parse `3D/3dmodel.model`, histogram undirected
index-edge incidence (`>2` = non-manifold, `==1` = boundary/open). Script inline in this
session's scratchpad.

Code tree: `/tmp/PrusaSlicer` @ `version_2.9.4`. Key files:
`bundled_deps/admesh/admesh/{connect,normals,util}.cpp`,
`src/libslic3r/{TriangleMesh.cpp,MeshSplitImpl.hpp,Format/3mf.cpp,Format/STL.cpp,Model.cpp}`.

# Bino slicer-hygiene autopsy: the residual 336 open / 150 reversed

**Date:** 2026-07-18
**Referee:** `prusa-slicer --info` (formal judge) + `tools/openedges_exact.py`,
`tools/overlaps.py`, and `/tmp/autopsy.py` (bit-exact-weld edge classifier).
**Model:** `examples/torture/zeiss_id02_bino.sb` @ resolution 1 (r1), sp ≈ 0.33 mm
(from the seal comment's own calibration: `3e-4 sp` = "a tenth of a micron at
product scale" → sp ≈ 0.33 mm; the 0.1 mm feature bar ≈ 0.30 sp).
**Rule of the house:** measured, not predicted. Every number below is from an
export + a referee pass, not a theory.

---

## 0. The headline: there are ZERO geometric holes

`openedges_exact.py` (bit-exact float weld, signed directed-edge count) on the
shipped baseline STL:

```
128310 tris, 64252 welded verts, 0 open edges, 0 unbalanced-multi
```

The mesh is a **closed, orientation-coherent, bit-exact-watertight** surface.
prusa nonetheless reports:

```
manifold = no   open_edges = 336   backwards_edges = 336   facets_reversed = 150
```

`open_edges == backwards_edges == 336` is the tell. prusa's admesh topology
builder assumes **exactly two facets per edge**. It is not counting holes; it is
counting edges it could not pair into a 2-manifold. The whole residual is
**non-manifold seams + orientation folds**, and both live entirely in the
close-perimeter pocket region (`x∈[-14,11], y∈[37,62], z∈[57,74]`) — the
"CLOSE-PERIMETER can of worms" named in MESH-WAR.

---

## 1. Attribution matrix (the A/B grid)

Every row is a full r1 export + prusa + autopsy. `nm` = fat edges (>2 incident
half-edges, on bit-exact weld); `fold180` = 2-incident anti-parallel neighbor
pairs (thin fins); `pr_open`/`pr_rev` = prusa `open_edges`/`facets_reversed`.

| config (env)                | tris    | nm  | fold180 | pr_open | pr_rev |
|-----------------------------|---------|-----|---------|---------|--------|
| **baseline** (all on)       | 128310  | 318 | 168     | **336** | **150**|
| `WINDING=0`                 | 128310  | 318 | 168     | 336     | 150    |
| `SEAL=0`                    | 128800  | 502 | 693     | 492     | 162    |
| `SEAL=0 WINDING=0` (raw)    | 128800  | 502 | 693     | 492     | 162    |
| `DECIMATE=0`                | 245258  | 571 | 482     | 440     | **6**  |

Two exact identities fall straight out of the grid:

- **`WINDING=0` ≡ baseline, bit-for-bit** (336/150, 318 nm, 168 fold). The
  winding pass flips **zero** facets on bino (`WINDING: … flipped` never
  prints). It is a **complete no-op here** — neither helps nor hurts.
- **`SEAL=0` ≡ `SEAL=0 WINDING=0`** (492/162). Same fact from the other side:
  with seal off, winding still changes nothing.

The seal debug line on the shipped path:
```
SEAL: 117 near-coincident vertices merged, 488 twin facets annihilated
```

---

## 2. Attribution per pass

### SEAL — the workhorse (owns the twin/coincident class)
Turning it off: open 336→492, nm 318→502, fold180 168→693, reversed 150→162,
and 243 exact coincident-triple twin pairs reappear (`coincident triangle sets:
243` under `SEAL=0`, `0` with seal on). Seal annihilates 244 twin pairs (488
facets) and merges 117 vertices. **It cures every seam whose two sheets are
bit-close** (< eps = 3e-4·sp = 1e-4 mm).

**Why 318 survive it:** the residual seams are not numerical twins. Measured
inter-sheet apex gap on the 316 four-incident edges:

```
min 0.01196 mm   median 0.24229 mm   max 0.49782 mm
apex-gap ≤ 1e-2 mm: 0   ≤ 3e-2: 28   ≤ 1e-1: 85   ≤ 3e-1: 198   ≤ 1e0: 316
```

**Not one is within 100× of the weld eps.** The smallest gap (0.012 mm = 0.036
sp) is 120× the eps; the median (0.73 sp) is a feature-scale wall separation.
Seal *cannot* reach these without a weld radius that would fuse real 0.1 mm
walls. This is a **topology** problem, not a weld-epsilon tuning problem.

### WINDING — a no-op on bino (owns nothing here)
Zero flips. The pass floods orientation across clean 2-facet edges and votes one
sign per component, but only on **decisive** components (`comp.size() >= 8` and
`|vote| > 0.5·mag`). The 150 misoriented facets (see below) sit in **fragments
bounded by the 318 nm edges** — nm/pinch edges are flood barriers, so each fold
is its own <8-tri island that never meets the gate. The pass is architecturally
blind to exactly the defect that remains.

### DECIMATE — mints the reversed facets (owns the 150 class)
This is the surprise and it resolves the MESH-WAR open question ("the 14→150
reversed delta may predate the seal"). It does **not** predate the seal — it is
**minted by flat-decimation**:

- `DECIMATE=0`: **facets_reversed = 6** (essentially orientation-coherent).
- baseline (decimate on): **facets_reversed = 150**.

Decimation is a genuine tradeoff on bino: it **helps** the seam class (nm
571→318, open 440→336, fold 482→168, tris 245K→128K) by collapsing close-
perimeter slivers, but its flat-collapse **mints 144 orientation-incoherent
facets** that prusa then reverses. The RAZOR independent-set fix (2026-07-18)
killed the *strict-interior fights*, but prusa's `facets_reversed` — winding
incoherence vs. neighbors — is a different metric and still reads 150.

---

## 3. Anatomy classes (counts + coordinates)

### Class A — non-manifold pinch seams  → prusa open 336 / backwards 336
318 fat edges, all **net winding 0**:
- incidence histogram: `{4: 316, 6: 2}` (316 four-way edges, 2 six-way triple junctions)
- connectivity: **251 chains** — **226 isolated single edges** (scattered pinch
  points), one **31-edge chain** (a real seam loop), plus a 7, a 6, and a
  handful of 2–3-edge stubs.
- geometry split (by pairing the 4 incident normals): **149 doubled-sheet**
  (two near-parallel pairs separated by the wall gap) + **167 genuine junction**
  (4-way fan / true self-contact).
- mechanism: prusa pairs only 2 of the 4 half-edges at each edge; the leftover
  pair is reported as **both** an open edge and a backwards edge → open==back==336.

Representative sites (edge midpoints, from `/tmp/autopsy.py`):
```
junction:  (4.624, 54.739, 67.406)  (11.487, 55.540, 69.118)  (-11.533, 58.323, 71.630)
doubled:   (11.242, 58.084, 68.614) (-11.284, 57.306, 71.631) (5.981, 58.589, 69.618)
```
Densest 2 mm bins: `(-12,56,70)` 18, `(-6,50,70)` 13, `(-6,56,70)` 13, `(12,60,70)` 11.

**Field oracle confirms these are real geometry** (`STIBIUM_DMESH_FPROBE`, ends
the is-it-real argument). Every seam midpoint sits on the isosurface
(`|f|/|g| = 0.02–0.06`, i.e. ≈0 distance). Stepping ±0.2 mm across the junction
at `(4.62,54.74,67.41)`: the field is solidly **interior on both sides**
(`f=-0.54`, `|f|/|g|≈0.045 mm` inside) while the seam line itself is a near-zero
pinch ridge (`f=+0.12`). That is a genuine thin-feature self-contact, not a
tessellation phantom — matching the sub-sp apex gaps of §2.

### Class B — decimation-minted orientation folds → prusa facets_reversed 150
150 facets whose winding disagrees with their neighbors, created by
flat-decimation collapses, invisible pre-decimation (6). Bounded into sub-8-tri
fragments by Class-A edges → winding pass can't vote them. Co-located with
Class A (same pocket region, `z∈[52,74]`).

### Class C — thin fins (fold180) → contribute to both
168 anti-parallel 2-incident pairs; min triangle area 1.5e-6, median 1.5e-3 mm².
These are near-zero-area flaps at the seam. `DECIMATE=0` has 482 of them (coherently
wound → only 6 reversed), so fins ≠ reversed; they are the sliver substrate that
decimation both thins and mis-collapses.

---

## 4. Ranked cures with referee plans

### Cure 1 (biggest lever, source fix) — orientation-guard the decimation collapse
**Target: facets_reversed 150 → ~6.** The 150 are minted by `decimate_flats`.
Extend the independent-set collapse rule with a **winding veto**: reject any
collapse that would leave an incident facet wound against the field gradient (or
against its surviving neighbors) at its new centroid — the same gradient-vs-normal
oracle `repair_winding` already computes, applied *per candidate collapse* instead
of post-hoc. Pre-decimation coherence (6) proves the target is reachable without
touching geometry.
- **Expected:** pr_rev 150→≤10; nm/open unchanged (Class A untouched); tris held.
- **Referee:** `prusa-slicer --info` (rev<10?), `tools/tilt.py` on screws (no
  regression — screws are already manifold=yes/0/76-class), depth referee
  (`_CHIP_DEBUG` worst-chip must hold 0.097 sp).

### Cure 2 (owns the open count) — non-manifold edge split at the pipeline tail
**Target: open 336 → ~0, nm 318 → ~0.** Add a tail pass after seal: for every
edge with 4 incident half-edges, **split the shared edge per sheet** (duplicate
the two edge vertices, reassign each 2-facet sheet to its own copy). 
- **149 doubled-sheet edges** resolve cleanly — each sheet keeps its opposite-wound
  pair → 2-manifold, still bit-watertight (the sheets were already closed
  elsewhere; splitting only un-shares the contact edge).
- **167 genuine junctions**: same duplication yields two touching-but-distinct
  manifold copies (prusa manifold=yes; the geometry genuinely self-contacts).
- **28 sub-feature edges** (gap ≤ 0.03 mm = 0.09 sp, below the 0.1 mm bar): safe to
  *weld* the sheets instead — a raised local seal only where gap < 0.1·sp.
- **Cost:** ~+600 verts, geometry byte-identical. The 226 isolated pinch points are
  independent local edits; the one 31-edge chain needs a consistent split walk.
- **Referee:** `openedges_exact.py` (still 0), `/tmp/autopsy.py` (nm→0),
  `prusa-slicer --info` (open→0, manifold=yes).

### Cure 3 (best-of-both, larger) — decimation redesign
Decimation *helps* Class A (571→318) and *hurts* Class B (6→150). A collapse
scheme that thins close-perimeter slivers **without** minting folds captures the
nm reduction AND the pre-dec reversed floor. Larger than Cure 1; do Cure 1 first
and only escalate if Cure 2's split leaves a stubborn nm floor that decimation
would have thinned.

### What will NOT work
- **Raising `STIBIUM_DMESH_SEAL`.** Gaps are 0.012–0.5 mm; a weld that reaches them
  fuses real walls. Measured, dead on arrival.
- **Loosening the WINDING gate.** The gate exists because ungated fragment flips
  minted backwards edges (632 flips → +12 opens, per the pass comment). The folds
  are in nm-bounded islands the flood can't safely vote. Fix the source (Cure 1),
  don't unlock the vote.
- **Higher resolution (r2).** These defects are coincident/topological
  (resolution-immune per MESH-WAR); r2 halves *absolute* defect size but the nm
  seams and orientation folds are counts, not sizes. Must be topological.

### Suggested order
1. **Cure 1** (decimation winding-veto) — one source change, kills 150 reversed,
   referee-cheap, screws already clean so low regression risk.
2. **Cure 2** (nm edge split) — kills the 336 open; the harder but well-scoped one.
3. Re-prusa. If manifold=yes / 0 / <10, the bino class joins screws as DONE.
   Then and only then consider Cure 3 / an r2 validation dish.

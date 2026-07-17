# Mesh revision ledger

One row per reviewed revision.  STLs live in build/zeiss_dmesh/
(KEEP - Nate reviews rev-by-rev); renders in renders/.  Stats are
the dmesh export line: tris / open / non-manifold / constrained /
worst chip depth (sp).  "Verdict" is Nate's eyes - the formal
referee.  Update this table with EVERY rev handed over for review.

## Full zeiss (examples/torture/zeiss_id02_merged.sb, m20, r1)

| rev | commit | change vs previous | tris | open | nm | cstr | worst | verdict |
|-----|--------|--------------------|------|------|----|------|-------|---------|
| d0 | pre-892e2190 | no density, manual era | ~400K | 0 | 20 | ~12K | 0.774 | baseline |
| d1 | - | manual global density knob | 407K | 0 | 50 | ~12K | 0.609 | old champ |
| autod14 | ade2b437 | stage-D auto-density + retreat loop | 726K | 0 | 168 | ~12K | 0.573 | collars at ungated ideal |
| autod16 | 983d4187 | phantom oracle + projection leak fix | ~1M | 0 | 187 | ~13K | 0.776* | *honest instruments read higher |
| autod17 | 6afdbd9f | repair keep-out retired | 1.04M | 0 | 124 | ~13K | 0.212 | reference for a day |
| autod19 | 5c31a6a5 | ring-eater dedupe fix + strips default on | 2.6M | 0 | 640 | ~13K | 0.324 | "cleanest yet, ~97%" (2.5x tris sanctioned) |
| autod20 | b152a810 | shallow channel + 2nd ring-eater fix | 2.79M | 4 | 965 | 0! | 0.576 | trust-gate cliff: ALL law torched |
| autod21 | 112b11a1 | trust gate fallback-only | 3.37M | 12 | 963 | 16.1K | 0.576 | shipped worst attempt (retreat vandalism) |
| autod22 | (autopsy) | graduated rollback + ship-best | 3.04M | 4† | 967 | 16.1K | 0.576 | †phantom: 0 geometric |
| autod23 | 8730c4e8 | holes counted geometrically | 3.04M | 0 | 967 | 16.1K | 0.576 | **"cleanest ever"**: sharp shallow engravings, clean steps, sharp flush screws; printable grade |
| autod24 | (A/B) | STRIPS=0 | 979K | 0 | 128 | 15.8K | 0.604 | ratty joint bands - strips load-bearing |
| autod25 | fa454dc0 | strip gap bar 0.5 sp | 2.83M | 0 | 985 | 15.7K | 0.638* | *renders identical; dial kept, default 0 |
| autod26 | (A/B) | strip gap bar 0.25 sp | 2.78M | 0 | 991 | 15.7K | 0.576 | -10% size but fuzz at fine features; autod23 still best |
| autod27 | da10e554 | kink-law referee fix (grid-aligned creases keep law) | 1.14M | 0 | 191 | 16.0K | 0.211 | best numbers ever BUT spikier rims, ultrafine detail loss - churn was subsidizing detail density (densdiff: the 2M delta sits ON detail areas) |
| autod28 | c0f0aa4f | DC bias suppressed in crease band (air-chord fix) | 1.13M | 0 | 175 | 16.0K | 0.294 | splits 2906 -> 1041; plinth referee: "solved the air chords... perfect" (Nate); full-model slicer verdict pending |

## Bino referee (examples/torture/zeiss_id02_bino.sb, fast loop)

| rev | change | tris | open | nm | cstr | worst | verdict |
|-----|--------|------|------|----|------|-------|---------|
| v15 | dedupe + strip cores | - | 0 | 588 | - | - | "by FAR the cleanest" - crowned, strips became default |
| v16 | committed HEAD baseline (pre-shallow) | 716K | 4 | 588 | 2836 | 0.334 | secretly 4-open |
| v17 | dup-guard point-to-segment 0.1 sp | 711K | 0 | 693 | 3042 | 0.310 | +53 chains, watertight restored by coverage |
| v18 | + shallow seed channel | 783K | 0 | 684 | 3210 | 0.287 | shallow engravings sharp (Nate, on autod23) |
| v19 | per-chain trust conviction (REFUTED) | 742K | 14 | 657 | 3028 | 0.347 | strips oracle-verified law - reverted |
| v20 | graduated rollback (no-op on bino) | 783K | 0 | 684 | 3210 | 0.287 | = v18 |
| v21 | strip gap bar 0.5 | 761K | 0 | 690 | 3244 | 0.351 | marginal refund; default stayed 0 |
| v22 | kink-law referee fix | 275K | 0 | 199 | 3224 | 0.167 | best numbers ever BUT rims spikier (churn subsidy withdrawn) |
| v23 | DC crease-band suppression | 267K | 0 | 183 | 3224 | 0.238 | splits HALVED (1511 -> 653); worst up - phantom roofs were hiding real divots |

## Open questions the next rev must answer

- Selective detail density to replace the churn subsidy: per-leaf
  CHAIN-SEGMENT CROWDING as a stage-D driver (fine detail = dense
  polylines; the densdiff hotspots are all chain-crowded regions).
- Quilting on cones (lamp model): chainless-curvature density
  signal (Nate's design seed, 2026-07-18).
- Thin cone caps (knob model): thin-geometry tangle class.

## Plinth verdict (2026-07-18, c0f0aa4f) - AIR CHORDS SOLVED

(An earlier version of this entry declared the teeth "real
geometry" - REFUTED by Nate's model reading + the z=9.5 wall
probe.  The record stands corrected:)

The air-chords were the DC all-surface-cell bias (built the same
morning the tracer landed, for UNCONSTRAINED corners) roofing
over the junction's AIR WEDGE beside a now-constrained crease -
volume-adding facets wall-to-face across the crease fence,
labeled "visually benign" in its own comment.  Convicted by
elimination: repair oracle (0 drops), snap (nosnap identical),
density (hidden + level 3: nothing), real geometry (FPROBE: all
vertices on-surface AND no ledge exists - Nate's "simple
cylinder on cube" reading correct).  DC=0 test: 88 teeth -> 0.
FIX: bias suppressed within STIBIUM_DMESH_DC_BAND (default 0.75
cells) of constrained segments; centroid-only there.  NATE'S
EYES on the plinth: "solved the air chords... perfect."
Along the way: kink-law referee fix (model-wide), law-blind
crowding dial (refuted, parked), repair post-projection oracle
(armor), FPROBE point instrument (the closer).
Still open from this arc: detail-density for the autod23/27
tradeoff (chain crowding), cone quilting (chainless curvature),
thin cone caps (knob model).

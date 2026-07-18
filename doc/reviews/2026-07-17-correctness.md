# Opus correctness review - 2026-07-17 (read-only; findings NOT yet acted on except #1's decimation half, which corroborated Nate's gash sighting and is being fixed in the live loop)

## Findings (ranked)
1. **Decimation runs AFTER the quality metrics and has no link-condition
   check** (moderate; certain/likely).  open/nm counts describe the
   pre-decimation mesh; ledger columns lie about the shipped STL.  The
   naive nearest-neighbour collapse can also mint non-manifold edges or
   doubled tris that the in-plane orientation guard cannot see.
   Prescribed test: screws decimate on/off + independent edge sweep on
   the FINAL STL.  [CORROBORATES Nate's razor-gash sighting on autod31.]
2. **tri_key packs vertex indices into 21-bit fields** (moderate;
   certain overflow, likely reachable): silent collision above 2,097,151
   vertices in the manifold pass - wrong-triangle sheet assignment,
   possible nm/hole minting.  Multi-million-tri exports are already
   close.  Cheap guard: assert verts/3 < 2^21 at the manifold block.
3. CCDT near-coincidence weld overwrites a sign witness's info to
   surface unconditionally (minor/speculative) - siblings guard this;
   this one path doesn't.  Instrument before fixing.
4. delaunay_chains symmetry check is dead code ("|| true") - fallback
   path only; delete and confirm csg counts.
5. Kink-referee zero-grad bar is ABSOLUTE (1e-6 model units) unlike
   every neighbouring spacing-relative tolerance - scale-dependence
   risk; test by scaling the plinth 100x both ways.
6. Unguarded integer shifts on extreme knob values (operator foot-guns
   only; defaults safe).

## Clean areas verified sound
Two-pass drill-down carryover (no OTHER pass-2 empty-state readers);
retreat-loop state channels; prov pointer lifetime (no dt.remove
anywhere, Compact_container stability); snap-tent manifoldness and
winding; the sheet[k] operator[] trap's reconciliation; weld stencil
indexing; local_step band lookup boundary behaviour; NaN/Inf
discipline throughout the eval consumers.

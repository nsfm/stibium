# Opus performance review - 2026-07-17 (read-only; not yet acted on)

## The 20 GB, explained
DT vertex ~500-600 B per inserted point at 3D Delaunay's ~6.5 cells/vertex
(more on the CCDT path).  20 GB implies TENS OF MILLIONS of inserted
points.  "Point count is the entire memory story, and it is also the
675 s insert story.  Every finding that reduces inserted-point count
pays twice."

## Findings, ranked by impact/risk
1. **SAMPLE THINNING (the big lever, time AND memory)**: the DT holds
   every signed sample - culled-box far-field corners plus the FULL
   64-corner lattice of every near-surface leaf - though refinement only
   cares about vertices bounding inside<->outside edges.  Drop interior
   witnesses not lattice-adjacent to a sign change; keep a one-ring
   shell around the surface band, coarsen by octree level away from it.
   "The decimation lesson applied UPSTREAM."  Gate on the open-edge
   referee.  Plausibly a large fraction of both 675 s and 20 GB.
2. **FREE PASS-1 STATE BEFORE PASS 2** (the cheap first win): c.soup
   .samples/.surface, c.seen_samples, c.edge_index, c.edges all stay
   resident while pass 2 samples (only crease_leaves is moved).  Clear +
   shrink_to_fit after the move - roughly halves sampling-phase RSS.
3. **CCDT CONCURRENT INSERTION IS NOT FEASIBLE** - the conforming TDS
   bases are not Parallel_tag bases and constraint bookkeeping is not
   thread-safe.  The only true-parallel route: bulk samples into a plain
   concurrent DT + CCDT only for the crease neighborhood - large
   refactor, research item, not a quick fix.  Near-term insert-time win
   is finding 1.
4. box_dense_factor: LINEAR scan over all dense boxes per point in hot
   paths (corridor live by default) - grid/BVH it, purely mechanical.
5. Strip detection O(P^2) over chain vertices (~1e8 pair tests at
   attempt 0) with an inner dense_boxes scan - grid at strip_r -> O(P).
6. Dup-trace guard + chain mediation are O(n^2)-ish in the tracer -
   spatial index collapses both.  (Segment referee already well-batched,
   no finding.)
7. Refinement loop full edge sweep x 48 rounds + unreserved seen set -
   cheap: reserve(); structural: incremental re-detection (only edges
   near fresh inserts), gate on open-edge referee parity.
8. Transient copies: pts double-copy held through insert; vgrid map
   never reserved (multi-million rehash); ship-best deep-copies per
   attempt (swap instead).  Peak-shaving.
9. Cascade/rollback paths rebuild the ENTIRE DT per retry (up to 8 x 4
   nested) - the tail that turns 25 min into an hour.  Longer-term:
   incremental constraint re-insertion.

## Reviewer's sequencing
2 first (low-risk memory), then 1 (big lever, referee-gated), then
4/5/6 (mechanical indexing), 3 + 7 as research items, 8-9 peak-shaving.

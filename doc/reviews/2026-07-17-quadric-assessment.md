# Opus assessment: Fast-Quadric-Mesh-Simplification as dependency (2026-07-17)

## Bottom line
DON'T vendor sp4cerat - app/vendor/meshoptimizer (MIT, already in
tree, already called at export_mesh.cpp:208-244) is the same
algorithm class (QEM edge collapse) with the knobs sp4cerat lacks:
- meshopt_SimplifyLockBorder - watertightness protection
- meshopt_simplifyWithAttributes + vertex_lock - hard-lock the
  constrained-crease set (the sharp-edge law CANNOT erode)
- deviation-bound targeting (model units via simplifyScale) vs
  sp4cerat's count-only target (blind to the 0.1mm bar)

sp4cerat: sound classic Garland-Heckbert, but NO crease
preservation, file-scope global state (not thread-safe),
count-target API, "capable hobby header" quality.  License would
be fine (MIT->GPLv3, notice retention - same as meshoptimizer);
redundancy is the blocker, not law.

## What Nate's experiment proved
zeiss_autod31_decimate.stl (727,058 tris) had ONLY the in-house
conservative flat pass - the meshopt QEM pass wasn't running on
the DMesh path (deviation unset).  The sp4cerat result (363,528 =
50.002%, its default reduce_fraction) with "no visual defects" is
a GREEN LIGHT for QEM on this geometry - deliverable from
already-vendored code.

## Integration plan (5-9 hrs total)
1. Enable meshopt on DMesh exports (env deviation, e.g.
   STIBIUM_DMESH_SIMPLIFY=<mm>) - 1-2h
2. Crease-lock: simplifyWithAttributes + vertex_lock built from
   the mesher's constrained vertices + LockBorder - 3-5h
   (plumb the crease set through DMesh)
3. Referee: [.dmesh*] suite, open_edges==0 post-simplify,
   --resave determinism, extract renders, Nate's eyes - 1-2h

## Risks (any QEM)
Crease erosion at aggressive ratios (cured by vertex_lock);
border handling (LockBorder); determinism (meshopt already inside
the --resave envelope).  Re-referee watertightness, never assume.

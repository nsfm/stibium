# Opus architecture review - 2026-07-17 (read-only; not yet acted on)

## What the design gets right
- The public seam is genuinely good: five entry points, DSoup/DChains/DMesh as contracts, soup-in/mesh-out as a real testing seam. Diagnostics-in-struct discipline is what makes the referee ledger possible.
- Damage-is-measured-not-predicted is the right architecture, honestly implemented (retreat loop + ship-best).
- Provenance discipline is real (geometric vs index-space open counting).
- Knobs are measured into their defaults with the measurement in the comment - the file is its own lab notebook.
- Debt is in PRESENTATION, not algorithm: mesh_impl is ~2,980 lines in one function template; tracer ~870.

## DO SOON (low-risk, high-leverage)
1. Extract the instrument/debug layer behind one Debug struct (10 debug envs, 30+ getenv sites). 0.5-1 day, near-zero risk.
2. Consolidate ~26 env knobs into DMeshConfig resolved once per export. 1-1.5 days. Tiers:
   - permanent product dials: AUTODENSE(_MAX), CCDT, MANIFOLD, SNAP(_SURF), DECIMATE, STRIPS, TANGLE_DOT, SEP, LIVE, SHALLOW, DC(_BAND), REPAIR, ROUNDS, TRACE(_LOCAL), SEEDGATE, DENSE
   - frozen experiment switches: WELD, PAIR_FALLBACK, SDB, KEEPOUT, CROWD_LAW, STRIP_GAP, SLIDE
   - candidate-dead: HIDDEN (pending a graze-vs-feature oracle that doesn't exist)
   Note: static const char* read-once caching is thread-UNSAFE for the planned eval threading.
3. Extract a geom utility header: the 7-tap gradient probe is copied at SIXTEEN sites; point-to-segment clamp at 8. A gradient-robustness fix currently has to be made 16 times. 0.5-1 day.

## DO WHEN STABLE (wait for the quality bar to stop moving)
4. Split mesh_impl into named phase functions over a shared MeshState (3-5 days, medium-high risk - the phases share mutable CGAL state through enclosing scope; extraction would DISCOVER the implicit contracts).
5. Formalize the two-pass Collector handoff as a SurveyState sub-struct moved as one unit - "any future per-leaf state pass 2 needs will silently read empty unless someone remembers to move it" (the shallow-channel bug's exact shape). 1-2 days.
6. Regression tests for the referee-only bug classes - HIGHEST-VALUE GAP. Each mapped to a documented bug:
   - two-pass handoff liveness (pass-2 tseeds > 0 on a shallow-crease model) - cheapest, do first
   - retreat ships watertight AND ships-best (the autod21 regression is one '<' away)
   - geometric-vs-index open counting (guards weld_ids semantics)
   - grid-aligned spacing tie-flip (the one-ulp 16K->45K tri blowup)
   - trust-gate cliff (constraints survive local chord-sag rejects)
   2-4 days; mostly authoring small deterministic CI models (zeiss too slow for CI).
7. Optional: promote the referee ledger to a JSON-emitting harness - deferred deliberately; automating the numbers must not create pressure to trust them over the eyeball referee (the house rule that saved the campaign twice).

## Reviewer's caveat
Could not build/run (machine busy) - config-consolidation risk rests on grep evidence that no knob mutates at runtime; wants a before/after [.dmesh*] run when the machine is free.

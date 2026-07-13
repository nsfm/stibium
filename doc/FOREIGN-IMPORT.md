# Importing libfive Studio / Fidget projects into Stibium

*Recon, 2026-07-13. Companion to doc/LIBFIVE-RECON.md. Grounded in
shallow clones of both repos, Stibium's parser/opcodes, and a working
converter prototype measured against every model both projects ship.*

## The formats

**libfive Studio `.io`** is raw Guile Scheme source, nothing else —
save is literally `out << editor->getScript()`
(`studio/src/window.cpp:569`). No metadata; bounds/resolution are
Scheme calls inside the script. The files are arbitrary programs by
design: of the 8 shipped examples, effectively none are pure
straight-line (macros in `pawn.io`, recursion + `map`/`apply` in
`menger.io`, unicode infix operators, rational literals). Studio can
also run the pane as Python — same story.

**Fidget `.vm`** is a flat SSA tape — the post-flatten DAG, one op per
line, named operands, last line is the output:

```
# A quarter-circle in the lower-left quadrant
y var-y
x var-x
mxy max x y
y2 square y
x2 square x
r2 add x2 y2
f const 0.5
circle sub r2 f
out max mxy circle
```

Reader at `fidget-core/src/context/mod.rs:861`. Opcodes: leaves
`var-x var-y var-z const`; unary `abs neg sqrt square floor ceil round
sin cos tan asin acos atan ln not exp`; binary `add mul min max div
atan2 sub compare mod and or`. Fidget's `.rhai` scripts are the
parametric analogue of `.io` (straight-line vs programmatic split,
`draw()` output, `.remap()` methods).

## Op-vocabulary fit

Everything maps directly onto our opcode set except four exotics with
algebraic rewrites and three genuinely unmappable ops:

- `ceil(x)` → `-floor(-x)`; `round(x)` → `floor(x+0.5)`;
  `recip(x)` → `1/x`; `nth_root(x,n)` → `pow(x, 1/n)`.
- `atan2` exists as an opcode but has no V1 prefix char — emit V2
  `atan2(a,b)` embedded via `=…;`.
- **No equivalent**: `compare` (sign), boolean `and/or/not`, libfive's
  `NANFILL`, general oracles (mesh oracles ≈ our `grid`).

**Empirical result: across all six shipped `.vm` models, zero
occurrences of the unmappable ops.** Frozen import has 100% op
coverage on the real corpus.

## Blow-up and parser-stack measurements (prototype, real models)

`.vm` is a shared DAG; our math strings are trees, so shared
subexpressions inline textually. Measured expansion:

| model | DAG nodes | tree nodes | blow-up | max depth | V1 bytes |
|---|---|---|---|---|---|
| quarter | 9 | 11 | 1.2× | 5 | 14 B |
| tanglecube | 20 | 28 | 1.4× | 9 | 33 B |
| hi | 56 | 97 | 1.7× | 10 | 184 B |
| colonnade | 808 | 2,782 | 3.4× | 51 | 4.9 KB |
| bear | 657 | 17,063 | 26× | 73 | 48 KB |
| prospero | 7,866 | 26,018 | 3.3× | 677 | 74 KB |

Blow-up is polynomial on this corpus (worst 26×); biggest string is
74 KB — trivial for a JSON field. Parser stack scales with left-spine
depth, not length: worst observed 677 against lemon's 4096. **No
practical size wall**; a converter should still count depth and warn
on pathological inputs.

## The libfive escape hatch: don't parse Scheme

libfive serializes its own trees: binary `libfive_tree_save/load`
(`libfive.h:285`, DAG-deduped ids — no blow-up) and text
`libfive_tree_print` (S-expression, expands sharing). So the strong
`.io` path is: run the file through libfive's own interpreter (Guile
or Python binding), dump the evaluated tree, and translate the flat
tree exactly like a `.vm`. Macros, recursion, unicode operators —
all sidestepped, because their own evaluator did the work. Reject
`COMPARE`/`NANFILL`/`ORACLE` trees with a clear error.

## Recommendation (build in this order)

1. **[S] `--import-vm file.vm -o out.sb`** — Fidget tape → one frozen
   script node (`output('shape', Shape('<v1 string>', bounds))`).
   ~200 lines of Python; ships prospero/bear/menger day one. `.vm`
   carries no bounds: default box or interval-probe pass.
2. **[S-M] libfive frozen import** via the tree escape hatch — a thin
   adapter invoking libfive to evaluate `.io`, then the same
   translator as (1). Effort is dependency plumbing, not logic.
3. **[M] Straight-line Rhai/Scheme → script-node transpiler**,
   contingent on our stdlib replication (the libfive-raid quick wins),
   for the minority of files that are declarative — the only path
   yielding *parametric* graphs. Detect macros/recursion and fall
   back to (2).

**Don't build**: a Scheme macro expander or Rhai interpreter inside
Stibium. Freeze the programmatic files instead.

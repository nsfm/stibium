# Agent Surface — design plan

How AI agents (and scripts, and CI) become first-class Stibium users.
Written 2026-07-13, the night the headless verbs landed, by the agent
who will be using this surface — biases disclosed accordingly.

## Where we are

The read/verify loop shipped:

| verb | role | status |
|---|---|---|
| `--validate` | conscience — script/datum errors, exit code | ✅ |
| `--render` | eyes — shaded PNG, iso/top/front | ✅ |
| `--export` | hands — 3MF/STL/SVG/DXF/PNG | ✅ |
| live reload | the user watches the agent work | ✅ |

What's missing is **write**: an agent today edits `.sb` JSON directly,
which works (the format is honest JSON and connections are readable
Python expressions like `s0.shape`) but is fragile at the edges.

## The layers, in build order

### 1. Deterministic serialization (foundation, small)

Stable node/datum/key ordering on save, so `.sb` files diff cleanly in
git and an agent's read-modify-write cycle doesn't shuffle unrelated
lines. Also: normalize float formatting.

One deliberate change beyond ordering: **output-datum `expr` strings
are derived state** (serialized `fab.types.Shape(...)` reprs that get
recomputed on load). The writer should either omit them or write a
placeholder — they're the noisiest, least-human, least-agent part of
the format, and removing them makes diffs mean something. Needs a
protocol bump with a loader that tolerates their absence.

### 2. Machine-readable node reference (small)

`antimony --describe-nodes` → JSON: every `.node` file's name,
category, title, inputs (name, type, default), outputs, and doc line.
The `.node` files are already self-describing scripts; this is a
parser and a dump. Same generator can feed human docs later.

This is the agent's *vocabulary*. Without it, every session starts by
re-reading `py/nodes/`; with it, one tool call loads the whole menu.

### 3. Python authoring library (the workhorse)

A small pure-Python package (`stibium.author` or similar) that builds
graphs and writes `.sb`:

```python
g = Graph()
prof = g.node("gear_profile", teeth=9, module=1.2)
gear = g.node("extrude_z", shape=prof.shape, z0=0, z1=0.8)
g.node("export_mesh", shape=gear.out)
g.save("gear.sb")
```

Why Python instead of CLI verbs like `--add-node`/`--link`: agents are
*better* at writing programs than at issuing imperative command
sequences — a program shows the whole intended structure at once, can
be reviewed before it runs, and re-running it is idempotent. CLI graph
surgery (`--add-node`, one flag-soup invocation per edit) is the worst
of both worlds: stateful, order-dependent, unreviewable. Skip it.

The library needs no engine: it emits `.sb` JSON (layer 1 makes that
stable) with node scripts taken from the node library (layer 2 tells
it what exists). Verification is the already-shipped loop:
`--validate` then `--render` then iterate.

Node *placement* (inspector positions) should be auto-laid-out
left-to-right by dependency depth so generated files open readably.

### 4. MCP server on the live session (the endgame)

A socket surface over a running GUI session exposing roughly:

- `describe_nodes()` — layer 2, live
- `get_graph()` — nodes, datums, links, errors, as JSON
- `add_node(script) / set_datum(node, name, expr) / link / delete`
- `get_errors()` — structured `--validate`
- `render(view, size)` — the viewport the user is looking at
- `undo() / redo()` — the same undo stack the user has

The user watches the model change live; the agent sees what the user
sees. The file-watcher reload already covers the disk-based version of
this loop; the MCP server removes the save/reload round-trip and makes
errors/undo shared state. Build it as a thin adapter over the same
operations layer 3 uses — not a second implementation.

Transport: local socket, JSON-RPC or MCP proper. Auth: filesystem
permissions (it's a local session tool).

## Patterns that help the agent reason (asked, and answering honestly)

- **Names are load-bearing.** `a0`/`t3` tells me nothing;
  `hub`, `tooth_profile`, `mount_holes` is a model I can navigate a
  week later. The authoring library should push named nodes hard, and
  the canvas-annotations TODO item (sticky notes / zones) becomes
  agent-readable documentation *inside* the file.
- **Sub-renders over descriptions.** `--render --node NAME` (render
  one node's output alone) would let an agent bisect a wrong model
  visually — the same way a human clicks through nodes when the
  final shape looks off. Cheap to add to the existing render verb;
  the single highest-value follow-up on the read side.
- **Structured errors.** `--validate --json` (node, datum, line,
  message) once anything programmatic consumes it.
- **Determinism everywhere.** Same input → byte-same file. Agents
  (and git) trust diffs; every nondeterministic byte is noise in the
  loop. This is layer 1 and it's why it's first.

## Order and cost

1 and 2 are each an evening. 3 is a weekend, and its API should be
designed by writing the *test models first* (the gear above, a
bracket, a mask) and making them pleasant. 4 is a project; do it only
after 3's operations feel right, and reuse them.

Not rushing this is correct — but 1 and 2 have no design risk and
unblock everything, so they can land anytime.

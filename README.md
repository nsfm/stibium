# Stibium

*Stibium* is a computer-aided design (CAD) tool from a parallel universe
in which CAD software evolved from Lisp machines rather than drafting
tables.

Models are **functional representations**: pure math describing, for
every point in space, how far inside or outside the shape it is. You
build them in a **node graph**: boxes wired to boxes, every parameter
live, every node a small Python script you can open and edit. There are
no meshes until the moment you export one, so nothing ever gets
approximated, degenerate, or non-manifold while you work. Zoom in
forever; booleans always succeed; a thread is an equation.

It is a continuation of [Antimony](https://github.com/mkeeter/antimony)
by [Matt Keeter](https://mattkeeter.com/projects/antimony), itself a
spiritual successor to [kokopelli](https://github.com/mkeeter/kokopelli)
by way of [fabserver](http://kokompe.cba.mit.edu).

## What Stibium adds

Antimony went into maintenance mode in 2016. Stibium picks the project
back up and keeps going. Since the fork:

**Geometry in, geometry out**
- **Mesh import**: drop an STL (scan, vendor part, old export) into
  the graph and it becomes a solid distance field like everything
  else: subtract it, blend it, lattice it, re-export it. Sign comes
  from winding numbers, so imperfect scans still work.
- **Fidget import**: `stibium model.vm --import-vm out.sb` translates
  math tapes from Keeter's successor project
  [Fidget](https://github.com/mkeeter/fidget) into Stibium projects.
- **Modern export**: 3MF (typically 6× smaller than STL), SVG and DXF
  vector output with sharp-corner recovery for laser cutting and
  photolithography, and error-bounded mesh simplification. Meshing
  and rendering run on all cores.

**A deeper modeling vocabulary**
- O(1) **domain repetition**: infinite arrays, polar rings, mirrored
  tilings, and self-similar scale recursion at constant cost (the
  classic Array nodes union N copies and weep).
- **Parts library**: provably-mating ISO threads, involute gears,
  printable clearance holes with M2–M12 tables, dovetails, snap clips,
  hex vents, heat-set bosses.
- Smooth (log-sum-exp) and quadratic blends, chamfers and fillets,
  exact-distance primitives, gyroid lattices, bends and twirls.

**Headless and scriptable**
- The `.sb` format is deterministic, diffable JSON; saves are
  byte-stable.
- Headless verbs close the loop without a display: `--validate`
  (exit code), `--render` (see the model), `--analyze` (volume, center
  of mass, tight bounds as JSON), `--export`, `--resave`,
  `--describe-nodes` (the node library as JSON).
- Headless verbs mean CI for your CAD models.

**Quality of life**
- Cross-section preview, analytics overlay, recent files, export
  dialogs that suggest sane resolutions and show real dimensions, a
  canvas that stays smooth on 300-node graphs, and a seven-year
  backlog of crash fixes cleared.

The full ledger is in [CHANGELOG.md](CHANGELOG.md).

## Try it

[Build from source](BUILDING.md) (Linux and Mac).

Open an example to get oriented: `examples/showcase_gear.sb` is a
lightened involute gear built the way a user would build it, and
`examples/import_bead.sb` shows a scanned mesh being drilled by an
f-rep cylinder.

## Documentation

- [USAGE.md](doc/USAGE.md) - the tour
- [SCRIPTING.md](doc/SCRIPTING.md) - writing your own nodes

## Support

If you hit an issue or have a question,
[open an issue](https://github.com/nsfm/stibium/issues).

For Antimony history and design rationale, see
[Matt Keeter's writeup](https://mattkeeter.com/projects/antimony).

## License

Stibium's new work is released under the [GPLv3](LICENSE).

Code inherited from Antimony remains under its original MIT License;
see [THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md).

Copyright (c) 2013-2015 Matthew Keeter and other contributors (Antimony)

Antimony includes code from
[kokopelli](https://github.com/mkeeter/kokopelli), which is
© 2012-2013 Massachusetts Institute of Technology
© 2013 Matthew Keeter

#!/usr/bin/env python3
"""cyl_referee.py <cyl.stl> <chips.txt>

Acceptance referee for examples/detail_spacing_test_cyl.sb: bins
chip/air-chord defects by ring rim (each rim circle at r=10) so
defect-vs-gap reads as a curve.  Rings (z-extents): r0 20.5-22,
r1 15-16.5, r2 11.5-13, r3 9-10.5, r4 7-8.5; gaps above each:
8/4/2/1/0.5."""
import sys, collections

RIMS = [  # (label, z of rim circle)
    ('r0.top(gap8)', 22.0), ('r0.bot(gap4)', 20.5),
    ('r1.top(gap4)', 16.5), ('r1.bot(gap2)', 15.0),
    ('r2.top(gap2)', 13.0), ('r2.bot(gap1)', 11.5),
    ('r3.top(gap1)', 10.5), ('r3.bot(gap0.5)', 9.0),
    ('r4.top(gap0.5)', 8.5), ('r4.bot(base)', 7.0),
]

rows = [tuple(map(float, l.split())) for l in open(sys.argv[2])]
byrim = collections.defaultdict(lambda: [0, 0, 0.0])
other = 0
for x, y, z, dep, ln in rows:
    best, bd = None, 0.6
    for name, rz in RIMS:
        if abs(z - rz) < bd:
            best, bd = name, abs(z - rz)
    if best is None:
        other += 1
        continue
    e = byrim[best]
    e[0] += dep < 0
    e[1] += dep > 0
    e[2] = max(e[2], abs(dep))
print('rim              chips  airs  worst(sp)')
for name, rz in RIMS:
    e = byrim.get(name, [0, 0, 0.0])
    print('  %-15s %5d %5d  %.3f' % (name, e[0], e[1], e[2]))
print('  unbinned defects:', other)

#!/usr/bin/env python3
"""quiltdens.py <census.txt> <cx> <cy> <cz> [R]

DENS-map autopsy around a coordinate: prints every DENS leaf row
within R (default 6) of the point, plus GATE live/sep and nearby
hot CELL residuals - the offline referee for density-shadow
(quilted-patch) hunts.  Leaves NOT printed are base-level by
omission: a quilted patch with no DENS row at all means no trigger
ever considered the leaf (not even as a smooth pocket)."""
import sys

path, cx, cy, cz = sys.argv[1], *map(float, sys.argv[2:5])
R = float(sys.argv[5]) if len(sys.argv) > 5 else 6.0

dens, gate, cells = [], [], []
for ln in open(path):
    p = ln.split()
    if not p:
        continue
    if p[0] in ('DENS', 'GATE', 'CELL') and len(p) >= 5:
        try:
            x, y, z = float(p[1]), float(p[2]), float(p[3])
        except ValueError:
            continue
        if (x-cx)**2 + (y-cy)**2 + (z-cz)**2 > R*R:
            continue
        if p[0] == 'DENS':
            dens.append((x, y, z, int(p[4]), p[5]))
        elif p[0] == 'GATE':
            gate.append((x, y, z, p[4], p[6]))
        else:
            cells.append((x, y, z, p[5]))

print('DENS leaves within %.1f: %d' % (R, len(dens)))
for x, y, z, lv, tag in sorted(dens, key=lambda d: d[3]):
    print('  DENS (%8.3f %8.3f %8.3f) level %d %s' % (x, y, z, lv, tag))
print('GATE leaves within %.1f: %d' % (R, len(gate)))
for x, y, z, live, sep in gate:
    print('  GATE (%8.3f %8.3f %8.3f) %s %s' % (x, y, z, live, sep))
print('CELL rows within %.1f: %d (nr>=0.03: %d)'
      % (R, len(cells), sum(1 for c in cells
                            if float(c[3].split('=')[1]) >= 0.03)))

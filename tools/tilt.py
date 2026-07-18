#!/usr/bin/env python3
"""tilt.py <a.stl> [b.stl ...] - screws error-by-definition benchmark.

The screws model has ONLY axis-aligned surfaces (vertical walls,
horizontal flats): every triangle whose normal is not within TOL of
+-x/+-y/+-z is error by definition.  Reports tilted count, tilted
area fraction, and the area-weighted tilt angle (degrees off the
nearest axis) - lower is better on every column.
"""
import sys, struct
import numpy as np

TOL_DEG = 1.0

def read_stl(path):
    with open(path, 'rb') as f:
        head = f.read(80)
        if head[:5] == b'solid':
            # try ascii
            f.seek(0)
            txt = f.read()
            if b'facet' in txt[:400]:
                vs = []
                for line in txt.decode(errors='ignore').splitlines():
                    t = line.split()
                    if t[:1] == ['vertex']:
                        vs.append([float(t[1]), float(t[2]), float(t[3])])
                v = np.array(vs).reshape(-1, 3, 3)
                return v
            f.seek(80)
        n = struct.unpack('<I', f.read(4))[0]
        data = np.frombuffer(f.read(n * 50), dtype=np.uint8).reshape(n, 50)
        tri = data[:, 12:48].copy().view('<f4').reshape(n, 3, 3)
        return tri

for path in sys.argv[1:]:
    v = read_stl(path)
    e1 = v[:, 1] - v[:, 0]
    e2 = v[:, 2] - v[:, 0]
    nrm = np.cross(e1, e2)
    area2 = np.linalg.norm(nrm, axis=1)
    ok = area2 > 0
    nrm = nrm[ok] / area2[ok][:, None]
    area = area2[ok] / 2
    # angle off nearest axis: max |component| of the unit normal
    off = np.degrees(np.arccos(np.clip(np.abs(nrm).max(axis=1), -1, 1)))
    tilted = off > TOL_DEG
    total = area.sum()
    ta = area[tilted].sum()
    aw = (off * area).sum() / total if total > 0 else 0
    print(f"{path}: {len(area)} tris, {tilted.sum()} tilted "
          f"({100 * ta / total:.2f}% of area), "
          f"area-weighted tilt {aw:.3f} deg")

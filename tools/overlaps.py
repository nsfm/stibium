#!/usr/bin/env python3
"""overlaps.py <a.stl> [...] - z-fighting census.

Counts (on bit-exact welded ids):
  dup-tri:   triangle vertex-triples appearing more than once
  fat-edge:  edges with >2 incident triangles
  copl-pair: adjacent tri pairs whose dihedral is ~180 deg folded
             (normals anti-parallel across a shared edge = doubled
             sheet) or ~0 with overlapping span (coincident).
"""
import sys, struct
from collections import defaultdict
import numpy as np

def read_stl(path):
    with open(path, 'rb') as f:
        f.seek(80)
        n = struct.unpack('<I', f.read(4))[0]
        d = np.frombuffer(f.read(n * 50), dtype=np.uint8).reshape(n, 50)
        return d[:, 12:48].copy().view('<f4').reshape(n, 3, 3)

for path in sys.argv[1:]:
    v = read_stl(path)
    flat = v.reshape(-1, 3)
    uniq, inv = np.unique(flat, axis=0, return_inverse=True)
    tri = inv.reshape(-1, 3)

    tkey = defaultdict(int)
    for t in tri:
        tkey[tuple(sorted(t))] += 1
    dup = sum(k - 1 for k in tkey.values() if k > 1)

    e1 = v[:, 1] - v[:, 0]
    e2 = v[:, 2] - v[:, 0]
    nrm = np.cross(e1, e2)
    ln = np.linalg.norm(nrm, axis=1)
    ok = ln > 0
    nrm[ok] /= ln[ok][:, None]

    einc = defaultdict(list)
    for ti, t in enumerate(tri):
        for a, b in ((0, 1), (1, 2), (2, 0)):
            p, q = t[a], t[b]
            einc[(min(p, q), max(p, q))].append(ti)
    fat = sum(1 for lst in einc.values() if len(lst) > 2)

    fold = 0
    for (p, q), lst in einc.items():
        if len(lst) != 2:
            continue
        d = float(np.dot(nrm[lst[0]], nrm[lst[1]]))
        if d < -0.9999:
            fold += 1
    print(f"{path}: dup-tri {dup}, fat-edge {fat}, "
          f"fold180 {fold}, tris {len(tri)}")

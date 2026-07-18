#!/usr/bin/env python3
"""sheets.py <a.stl> [...] - near-coincident sheet census.

For each triangle, find non-adjacent triangles whose plane is
nearly parallel (|dot| > 0.999), whose centroid-to-plane distance
is tiny (< eps), and whose centroid projects within the other's
bounding sphere.  These pairs z-fight in a renderer at certain
angles while every adjacency census reads clean.
"""
import sys, struct
from collections import defaultdict
import numpy as np

EPS = 5e-3     # plane separation bar (model units)

def read_stl(path):
    with open(path, 'rb') as f:
        f.seek(80)
        n = struct.unpack('<I', f.read(4))[0]
        d = np.frombuffer(f.read(n * 50), dtype=np.uint8).reshape(n, 50)
        return d[:, 12:48].copy().view('<f4').reshape(n, 3, 3)

for path in sys.argv[1:]:
    v = read_stl(path)
    n_t = len(v)
    cen = v.mean(axis=1)
    e1 = v[:, 1] - v[:, 0]
    e2 = v[:, 2] - v[:, 0]
    nrm = np.cross(e1, e2)
    ln = np.linalg.norm(nrm, axis=1)
    ok = ln > 0
    nrm[ok] /= ln[ok][:, None]
    rad = np.linalg.norm(v - cen[:, None, :], axis=2).max(axis=1)

    # shared-vertex adjacency to exclude
    flat = v.reshape(-1, 3)
    uniq, inv = np.unique(flat, axis=0, return_inverse=True)
    tri = inv.reshape(-1, 3)
    vert2tri = defaultdict(set)
    for ti, t in enumerate(tri):
        for x in t:
            vert2tri[x].add(ti)

    cell = max(rad.mean() * 2, 1e-6)
    grid = defaultdict(list)
    for ti in range(n_t):
        key = tuple((cen[ti] // cell).astype(np.int64))
        grid[key].append(ti)

    pairs = 0
    seen = set()
    for ti in range(n_t):
        adj = set()
        for x in tri[ti]:
            adj |= vert2tri[x]
        kx, ky, kz = (cen[ti] // cell).astype(np.int64)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for tj in grid.get((kx+dx, ky+dy, kz+dz), ()):
                        if tj <= ti or tj in adj:
                            continue
                        if np.dot(nrm[ti], nrm[tj]) < 0.999:
                            continue
                        d = abs(np.dot(cen[tj] - cen[ti], nrm[ti]))
                        if d > EPS:
                            continue
                        if np.linalg.norm(cen[tj] - cen[ti]) > \
                           rad[ti] + rad[tj]:
                            continue
                        pairs += 1
                        seen.add(ti)
                        seen.add(tj)
    print(f"{path}: near-sheet pairs {pairs}, tris involved "
          f"{len(seen)} of {n_t}")

#!/usr/bin/env python3
"""openedges_exact.py <a.stl> [...] - independent watertight check.

Welds vertices bit-exactly (float triples), then counts geometric
open edges: directed edges without a matching opposite.  0 = law.
"""
import sys, struct
from collections import defaultdict
import numpy as np

def read_stl(path):
    with open(path, 'rb') as f:
        f.seek(80)
        n = struct.unpack('<I', f.read(4))[0]
        data = np.frombuffer(f.read(n * 50), dtype=np.uint8).reshape(n, 50)
        return data[:, 12:48].copy().view('<f4').reshape(n, 3, 3)

for path in sys.argv[1:]:
    v = read_stl(path)
    flat = v.reshape(-1, 3)
    uniq, inv = np.unique(flat, axis=0, return_inverse=True)
    tri = inv.reshape(-1, 3)
    count = defaultdict(int)
    for a, b, c in tri:
        for p, q in ((a, b), (b, c), (c, a)):
            count[(min(p, q), max(p, q))] += (1 if p < q else -1)
    open_edges = [e for e, k in count.items() if k != 0]
    nm = sum(1 for e, k in count.items() if abs(k) > 1)
    print(f"{path}: {len(tri)} tris, {len(uniq)} welded verts, "
          f"{len(open_edges)} open edges, {nm} unbalanced-multi")
    for e in open_edges[:8]:
        print(f"   open at {uniq[e[0]]} - {uniq[e[1]]}")

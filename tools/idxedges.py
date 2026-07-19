#!/usr/bin/env python3
"""Index-edge incidence histogram for a 3MF: >2 = non-manifold, ==1 = boundary/open.
Also reports distinct vertex positions (weld-invariance check)."""
import sys, zipfile, re
from collections import Counter

path = sys.argv[1]
with zipfile.ZipFile(path) as z:
    xml = z.read('3D/3dmodel.model').decode()

verts = re.findall(r'<vertex x="([^"]+)" y="([^"]+)" z="([^"]+)"', xml)
tris = re.findall(r'<triangle v1="(\d+)" v2="(\d+)" v3="(\d+)"', xml)

ec = Counter()
for a, b, c in tris:
    a, b, c = int(a), int(b), int(c)
    for u, v in ((a, b), (b, c), (c, a)):
        ec[(min(u, v), max(u, v))] += 1

hist = Counter(ec.values())
nm = sum(v for k, v in hist.items() if k > 2)
bd = hist.get(1, 0)
npos = len(set(verts))
print(f"{path}: verts={len(verts)} distinct_pos={npos} tris={len(tris)}")
print(f"  incidence_hist={dict(sorted(hist.items()))}")
print(f"  nonmanifold_edges={nm} boundary_edges={bd}")

#!/usr/bin/env python3
"""fightpix.py <axis> <a.stl> [...] - z-fighting pixel census.

Orthographic projection along +x/-x/+y/-y/+z/-z.  A FIGHTING pixel
has >= 2 front-facing triangles within EPS depth of the visible
surface.  Reports count + the worst scar clusters in world coords.
"""
import sys, struct
from collections import defaultdict
import numpy as np

RES = 900
EPS = 2e-3

AX = {'+x': (0, 1, 2, +1), '-x': (0, 1, 2, -1),
      '+y': (1, 0, 2, +1), '-y': (1, 0, 2, -1),
      '+z': (2, 0, 1, +1), '-z': (2, 0, 1, -1)}

def read_stl(path):
    with open(path, 'rb') as f:
        f.seek(80)
        n = struct.unpack('<I', f.read(4))[0]
        d = np.frombuffer(f.read(n * 50), dtype=np.uint8).reshape(n, 50)
        return d[:, 12:48].copy().view('<f4').reshape(n, 3, 3)

axis = sys.argv[1]
d_, u_, v_, sgn = AX[axis]

for path in sys.argv[2:]:
    V = read_stl(path)
    nrm = np.cross(V[:, 1] - V[:, 0], V[:, 2] - V[:, 0])
    ln = np.linalg.norm(nrm, axis=1)
    ok = ln > 0
    front = np.zeros(len(V), bool)
    front[ok] = (nrm[ok, d_] * sgn) > 1e-12 * ln[ok]

    lo = V.reshape(-1, 3).min(axis=0)
    hi = V.reshape(-1, 3).max(axis=0)
    span = max(hi[u_] - lo[u_], hi[v_] - lo[v_])
    px = span / RES

    depth = np.full((RES + 1, RES + 1), -np.inf)
    count = np.zeros((RES + 1, RES + 1), np.int32)
    dstore = defaultdict(list)

    tris = np.nonzero(front)[0]
    for t in tris:
        P = V[t]
        us = (P[:, u_] - lo[u_]) / px
        vs = (P[:, v_] - lo[v_]) / px
        ws = P[:, d_] * sgn
        u0, u1 = int(np.floor(us.min())), int(np.ceil(us.max()))
        v0, v1 = int(np.floor(vs.min())), int(np.ceil(vs.max()))
        u0 = max(u0, 0); v0 = max(v0, 0)
        u1 = min(u1, RES); v1 = min(v1, RES)
        if u1 < u0 or v1 < v0:
            continue
        gu, gv = np.meshgrid(np.arange(u0, u1 + 1),
                             np.arange(v0, v1 + 1), indexing='ij')
        # barycentric
        x1, y1 = us[0], vs[0]
        x2, y2 = us[1], vs[1]
        x3, y3 = us[2], vs[2]
        det = (y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3)
        if abs(det) < 1e-12:
            continue
        l1 = ((y2 - y3) * (gu - x3) + (x3 - x2) * (gv - y3)) / det
        l2 = ((y3 - y1) * (gu - x3) + (x1 - x3) * (gv - y3)) / det
        l3 = 1 - l1 - l2
        inside = (l1 > 1e-4) & (l2 > 1e-4) & (l3 > 1e-4)
        w = l1 * ws[0] + l2 * ws[1] + l3 * ws[2]
        ii, jj = gu[inside], gv[inside]
        for i2, j2, w2 in zip(ii.ravel(), jj.ravel(),
                              w[inside].ravel()):
            dstore[(i2, j2)].append(w2)

    fight = 0
    spots = []
    for (i2, j2), ws2 in dstore.items():
        if len(ws2) < 2:
            continue
        ws2 = sorted(ws2, reverse=True)
        if ws2[0] - ws2[1] < EPS:
            fight += 1
            spots.append((i2, j2, ws2[0]))
    print(f"{path} [{axis}]: fighting pixels {fight} "
          f"of {len(dstore)} covered")
    # cluster spots into bins of 20px, report top 4 in world coords
    bins = defaultdict(int)
    rep = {}
    for i2, j2, w2 in spots:
        k = (i2 // 20, j2 // 20)
        bins[k] += 1
        rep[k] = (i2, j2, w2)
    top = sorted(bins.items(), key=lambda kv: -kv[1])[:4]
    for k, n2 in top:
        i2, j2, w2 = rep[k]
        wc = [0., 0., 0.]
        wc[u_] = lo[u_] + i2 * px
        wc[v_] = lo[v_] + j2 * px
        wc[d_] = w2 * sgn
        print(f"   cluster {n2} px near ({wc[0]:.3f}, {wc[1]:.3f}, "
              f"{wc[2]:.3f})")

#!/usr/bin/env python3
"""areamap.py <axis> <a.stl> <out.png> [xmin xmax ymin ymax zmin zmax]

Triangle-AREA heatmap: orthographic z-buffer along +x/-x/+y/-y/+z/-z,
each pixel colored by log-area of the front-most triangle.  Low-density
(decimated / quilted) patches glow bright; dense regions stay dark.
Optional world-box crop.  The instrument for the cylinder-quilting
splotch hunt (quilted areas = far fewer tris, polygonal patches).
"""
import sys, struct
import numpy as np
from PIL import Image

RES = 1100
AX = {'+x': (0, 1, 2, +1), '-x': (0, 1, 2, -1),
      '+y': (1, 0, 2, +1), '-y': (1, 0, 2, -1),
      '+z': (2, 0, 1, +1), '-z': (2, 0, 1, -1)}


def read_stl(path):
    with open(path, 'rb') as f:
        hdr = f.read(80)
        if hdr[:5] == b'solid' and b'\n' in hdr:
            f.seek(0)
            pts = []
            for ln in f.read().decode(errors='replace').split('\n'):
                p = ln.split()
                if p[:1] == ['vertex']:
                    pts.append([float(v) for v in p[1:4]])
            return np.array(pts).reshape(-1, 3, 3)
        n = struct.unpack('<I', f.read(4))[0]
        d = np.frombuffer(f.read(n * 50), dtype=np.uint8)
        d = d.reshape(n, 50)[:, 12:48].copy().view('<f4')
        return d.reshape(n, 3, 3).astype(np.float64)


axis, path, out = sys.argv[1], sys.argv[2], sys.argv[3]
tri = read_stl(path)
if len(sys.argv) >= 10:
    box = np.array(sys.argv[4:10], dtype=float).reshape(3, 2)
    c = tri.mean(axis=1)
    m = np.all((c >= box[:, 0]) & (c <= box[:, 1]), axis=1)
    tri = tri[m]
d, u, v, sgn = AX[axis]
areas = 0.5 * np.linalg.norm(
    np.cross(tri[:, 1] - tri[:, 0], tri[:, 2] - tri[:, 0]), axis=1)

lo = tri[:, :, [u, v]].reshape(-1, 2).min(axis=0)
hi = tri[:, :, [u, v]].reshape(-1, 2).max(axis=0)
span = (hi - lo).max()
scale = (RES - 4) / span

zbuf = np.full((RES, RES), -1e30)
abuf = np.zeros((RES, RES))
for t in range(len(tri)):
    P = tri[t]
    uu = (P[:, u] - lo[0]) * scale + 2
    vv = (P[:, v] - lo[1]) * scale + 2
    ww = P[:, d] * sgn
    iu0, iu1 = int(uu.min()), int(np.ceil(uu.max()))
    iv0, iv1 = int(vv.min()), int(np.ceil(vv.max()))
    if iu1 < 0 or iv1 < 0 or iu0 >= RES or iv0 >= RES:
        continue
    iu0, iv0 = max(iu0, 0), max(iv0, 0)
    iu1, iv1 = min(iu1, RES - 1), min(iv1, RES - 1)
    gx, gy = np.meshgrid(np.arange(iu0, iu1 + 1),
                         np.arange(iv0, iv1 + 1), indexing='ij')
    det = (vv[1]-vv[2])*(uu[0]-uu[2]) + (uu[2]-uu[1])*(vv[0]-vv[2])
    if abs(det) < 1e-12:
        continue
    l0 = ((vv[1]-vv[2])*(gx-uu[2]) + (uu[2]-uu[1])*(gy-vv[2])) / det
    l1 = ((vv[2]-vv[0])*(gx-uu[2]) + (uu[0]-uu[2])*(gy-vv[2])) / det
    l2 = 1 - l0 - l1
    inside = (l0 >= -1e-9) & (l1 >= -1e-9) & (l2 >= -1e-9)
    if not inside.any():
        continue
    depth = l0*ww[0] + l1*ww[1] + l2*ww[2]
    zs = zbuf[iu0:iu1+1, iv0:iv1+1]
    as_ = abuf[iu0:iu1+1, iv0:iv1+1]
    won = inside & (depth > zs)
    zs[won] = depth[won]
    as_[won] = areas[t]

covered = abuf > 0
la = np.zeros_like(abuf)
la[covered] = np.log10(abuf[covered])
if covered.any():
    p2, p98 = np.percentile(la[covered], [2, 98])
    la = np.clip((la - p2) / max(p98 - p2, 1e-9), 0, 1)
img = np.zeros((RES, RES, 3), dtype=np.uint8)
img[..., 0] = (la * 255).astype(np.uint8)
img[..., 1] = ((1 - la) * 180 * covered).astype(np.uint8)
img[..., 2] = (covered * 40).astype(np.uint8)
Image.fromarray(np.flipud(img.swapaxes(0, 1))).save(out)
n = len(tri)
print(f"{path} [{axis}]: {n} tris, area p50/p90/p99 = "
      f"{np.percentile(areas,50):.5f}/{np.percentile(areas,90):.5f}/"
      f"{np.percentile(areas,99):.5f} -> {out}")

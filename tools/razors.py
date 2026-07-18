#!/usr/bin/env python3
"""razors.py <a.stl> [...] - razor census: triangles whose min
altitude is below a bar (default 0.02 * sp estimated from median
edge length... no - absolute flag).  Razor = z-fighting sliver."""
import sys, struct
import numpy as np

def read_stl(path):
    with open(path, 'rb') as f:
        f.seek(80)
        n = struct.unpack('<I', f.read(4))[0]
        d = np.frombuffer(f.read(n * 50), dtype=np.uint8).reshape(n, 50)
        return d[:, 12:48].copy().view('<f4').reshape(n, 3, 3)

for path in sys.argv[1:]:
    v = read_stl(path)
    e0 = np.linalg.norm(v[:, 1] - v[:, 0], axis=1)
    e1 = np.linalg.norm(v[:, 2] - v[:, 1], axis=1)
    e2 = np.linalg.norm(v[:, 0] - v[:, 2], axis=1)
    area = 0.5 * np.linalg.norm(
        np.cross(v[:, 1] - v[:, 0], v[:, 2] - v[:, 0]), axis=1)
    longest = np.maximum(e0, np.maximum(e1, e2))
    h = 2 * area / np.maximum(longest, 1e-30)
    med = np.median(longest)
    for bar in (0.005, 0.01, 0.02):
        print(f"{path}: h<{bar}: {(h < bar).sum()}", end="  ")
    print(f"(tris {len(v)}, median long edge {med:.3f})")

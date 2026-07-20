#!/usr/bin/env python3
"""ladder_referee.py <flat.stl> [chips.txt]

One-command acceptance referee for examples/detail_spacing_test.sb:
corner-vertex exactness for every groove (the junction law test),
tilted-tri census by groove, and chip totals from an optional
CHIP_DUMP file.  PASS = every corner EXACT, no tilted tris off
the noise floor."""
import sys, struct, math, collections

GROOVES = {  # name -> (left wall x, right wall x)
    'c1(8/8)':   (-12.0, -10.5),
    'c2(8/4)':   (-2.5,  -1.0),
    'c3(4/2)':   (3.0,   4.5),
    'c4(2/1)':   (6.5,   8.0),
    'c5(1/0.5)': (9.0,   10.5),
    'c6(0.5/-)': (11.0,  12.5),
}
RIM_Z, FLOOR_Z, YE = 2.5, 1.0, 10.0

def load(path):
    with open(path, 'rb') as f:
        f.read(80)
        n = struct.unpack('<I', f.read(4))[0]
        return [struct.unpack('<12f', f.read(50)[:48])[3:12]
                for _ in range(n)]

tris = load(sys.argv[1])
verts = set()
for t in tris:
    for v in (t[0:3], t[3:6], t[6:9]):
        verts.add((round(v[0], 4), round(v[1], 4), round(v[2], 4)))

fails = 0
tiltx = collections.Counter()
for t in tris:
    p = [t[0:3], t[3:6], t[6:9]]
    u = [p[1][i]-p[0][i] for i in range(3)]
    v = [p[2][i]-p[0][i] for i in range(3)]
    n = (u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2],
         u[0]*v[1]-u[1]*v[0])
    l = math.sqrt(sum(q*q for q in n))
    if l < 1e-12:
        continue
    if max(abs(q)/l for q in n) < math.cos(math.radians(5)):
        tiltx[round(sum(q[0] for q in p)/3)] += 1

print('== corner-vertex exactness (all four corners per groove):')
for name, (xl, xr) in GROOVES.items():
    worst = 0.0
    for x in (xl, xr):
        for y in (-YE, YE):
            for z in (RIM_Z, FLOOR_Z):
                c = (x, y, z)
                d = min((math.dist(c, v) for v in verts
                         if abs(v[0]-x) < 1 and abs(v[1]-y) < 1
                         and abs(v[2]-z) < 1), default=9)
                worst = max(worst, d)
    ok = worst < 1e-3
    fails += not ok
    print('  %-10s worst corner miss %.4f  %s'
          % (name, worst, 'EXACT' if ok else '<-- FAIL'))
print('== tilted tris by x:', dict(sorted(tiltx.items())) or 'NONE')
fails += sum(tiltx.values()) > 0
if len(sys.argv) > 2:
    rows = [tuple(map(float, l.split()))
            for l in open(sys.argv[2])]
    ch = [r for r in rows if r[3] < 0]
    ai = [r for r in rows if r[3] > 0]
    print('== chips %d (worst %.3f), air-chords %d (worst %.3f)'
          % (len(ch), min((r[3] for r in ch), default=0),
             len(ai), max((r[3] for r in ai), default=0)))
print('VERDICT:', 'PASS' if fails == 0 else 'FAIL (%d)' % fails)

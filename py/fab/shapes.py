import math
import functools
import operator

from fab.types import Shape, Transform

def preserve_color(f):
    """ Function decorator that preserves the color of the first input
        (which is expected to be a fab.types.Shape)
    """
    def p(s, *args, **kwargs):
        color = s._r, s._g, s._b
        return set_color(f(s, *args, **kwargs), s._r, s._g, s._b)
    return p

def union(a, b):
    return a | b

def intersection(a, b):
    return a & b

def difference(a, b):
    return a & ~b

@preserve_color
def offset(a, o):
    """ Assumes a linear distance field for bounds calculations!
    """
    if o < 0:
        return Shape('-%sf%g' % (a.math, o),
                a.bounds.xmin, a.bounds.ymin, a.bounds.zmin,
                a.bounds.xmax, a.bounds.ymax, a.bounds.zmax)
    else:
        return Shape('-%sf%g' % (a.math, o),
                a.bounds.xmin - o, a.bounds.ymin - o, a.bounds.zmin - o,
                a.bounds.xmax + o, a.bounds.ymax + o, a.bounds.zmax + o)

def clearance(a, b, o):
    return b | (a & ~offset(b, o))

@preserve_color
def shell(a, o):
    return a & ~offset(a, -o)

@preserve_color
def buffer(a):
    return a

def set_color(a, r, g, b):
    """ Applies a given color to an input shape a and returns it.
    """
    q = Shape(a.math, a.bounds)
    q._r, q._g, q._b = r, g, b
    return q

@preserve_color
def invert(a):
    """ Inverts a shape within its existing bounds.
    """
    if a.bounds.is_bounded_xyz():
        return cube(a.bounds.xmin, a.bounds.xmax,
                    a.bounds.ymin, a.bounds.zmax,
                    a.bounds.zmin, a.bounds.zmax) & (~a)
    elif a.bounds.is_bounded_xy():
        return rectangle(a.bounds.xmin, a.bounds.xmax,
                         a.bounds.ymin, a.bounds.ymax) & (~a)
    else:
        return Shape()

################################################################################

def circle(x, y, r):
    """ Defines a circle from a center point and radius.
    """
    # sqrt((X-x)**2 + (Y-y)**2) - r
    r = abs(r)
    return Shape(
            '-r+q%sq%sf%g' % (('-Xf%g' % x) if x else 'X',
                              ('-Yf%g' % y) if y else 'Y', r),
            x - r, y - r, x + r, y + r)

def circle_edge(x0, y0, x1, y1):
    """ Defines a circle from two points on its radius.
    """
    xmid = (x0+x1)/2.0
    ymid = (y0+y1)/2.0
    r = math.sqrt((xmid-x0)**2 +(ymid-y0)**2)
    return circle(xmid, ymid, r)

def polygon_radius(x, y, r, N):
    """ Makes a polygon with a center-to-vertex distance r
        The polygon is oriented so that the bottom is always flat.
    """
    # Find the center-to-edge distance
    r_ = -r * math.cos(math.pi / N)
    # Make an offset half-region shape
    half = Shape('-f%gY' % r_)
    # Take the union of a bunch of rotated half-region shapes
    p = functools.reduce(operator.and_,
            [rotate(half, 360./N * i) for i in range(N)])
    # Apply appropriate bounds and return
    return Shape(p.math, -r, -r, r, r)

################################################################################

def triangle(x0, y0, x1, y1, x2, y2):
    """ Defines a triangle from three points.
    """
    # Find the angles of the points about the center
    xm = (x0 + x1 + x2) / 3.
    ym = (y0 + y1 + y2) / 3.
    angles = [math.atan2(y - ym, x - xm) for x, y in [(x0,y0), (x1,y1), (x2,y2)]]

    # Sort the angles so that the smallest one is first
    if angles[1] < angles[0] and angles[1] < angles[2]:
        angles = [angles[1], angles[2], angles[0]]
    elif angles[2] < angles[0] and angles[2] < angles[1]:
        angles = [angles[2], angles[0], angles[1]]

    # Enforce that points must be in clockwise order by swapping if necessary
    if angles[2] > angles[1]:
        x0, y0, x1, y1 = x1, y1, x0, y0

    def edge(x, y, dx, dy):
        # dy*(X-x)-dx*(Y-y)
        return '-*f%(dy)g-Xf%(x)g*f%(dx)g-Yf%(y)g' % locals()

    e0 = edge(x0, y0, x1-x0, y1-y0)
    e1 = edge(x1, y1, x2-x1, y2-y1)
    e2 = edge(x2, y2, x0-x2, y0-y2)

    # -min(e0, min(e1, e2))
    return Shape(
            'ni%(e0)si%(e1)s%(e2)s' % locals(),
            min(x0, x1, x2), min(y0, y1, y2),
            max(x0, x1, x2), max(y0, y1, y2))

def right_triangle(x, y, w, h):
   # max(max(x-X,y-Y),X-(x*(Y-y)+(x+w)*(y+h-Y))/h)
   ws = math.copysign(1,w)
   hs = math.copysign(1,h)
   return Shape(
      'aa*f%(ws)g-f%(x)gX*f%(hs)g-f%(y)gY*f%(ws)g-X/+*f%(x)g-Yf%(y)g*+f%(x)gf%(w)g-+f%(y)gf%(h)gYf%(h)g' % locals(),
       x, y, x + w, y + h)

################################################################################

def rectangle(xmin, xmax, ymin, ymax):
    # max(max(xmin - X, X - xmax), max(ymin - Y, Y - ymax)
    return Shape(
            'aa-f%(xmin)gX-Xf%(xmax)ga-f%(ymin)gY-Yf%(ymax)g' % locals(),
            xmin, ymin, xmax, ymax)

def rounded_rectangle(xmin, xmax, ymin, ymax, r):
    """ Returns a rectangle with rounded corners.
        r is a roundedness fraction between 0 (not rounded)
        and 1 (completely rounded)
    """
    r *= min(xmax - xmin, ymax - ymin)/2
    return (
        rectangle(xmin, xmax, ymin+r, ymax-r) |
        rectangle(xmin+r, xmax-r, ymin, ymax) |
        circle(xmin+r, ymin+r, r) |
        circle(xmin+r, ymax-r, r) |
        circle(xmax-r, ymin+r, r) |
        circle(xmax-r, ymax-r, r)
    )

################################################################################

def tab(x, y, width, height, angle=0, chamfer=0.2):
    tab = rectangle(-width/2, width/2, 0, height)
    cutout = triangle(width/2 - chamfer*height, height,
                      width/2, height,
                      width/2, height - chamfer*height)
    tab &= ~(cutout | reflect_x(cutout))

    return move(rotate(tab, angle), x, y)

################################################################################

def slot(x, y, width, height, angle=0, chamfer=0.2):
    slot = rectangle(-width/2, width/2, -height, 0)
    inset = triangle(width/2, 0,
                     width/2 + height * chamfer, 0,
                     width/2, -chamfer*height)
    slot |= inset | reflect_x(inset)

    return move(rotate(slot, angle), x, y)

################################################################################

@preserve_color
def move(part, dx, dy, dz=0):
    return part.map(Transform(
        '-Xf%g' % dx, '-Yf%g' % dy, '-Zf%g' % dz,
        '+Xf%g' % dx, '+Yf%g' % dy, '+Zf%g' % dz))

translate = move

@preserve_color
def origin_xy(a, x0, y0, x1, y1):
    return move(a, x1 - x0, y1 - y0)

@preserve_color
def origin_xyz(a, x0, y0, z0, x1, y1, z1):
    return move(a, x1 - x0, y1 - y0, z1 - z0)

################################################################################

@preserve_color
def recenter(part, x, y, z):
    if not math.isinf(part.bounds.xmax) and not math.isinf(part.bounds.xmin):
        dx = x - (part.bounds.xmax + part.bounds.xmin) / 2
    else:
        dx = 0
    if not math.isinf(part.bounds.ymax) and not math.isinf(part.bounds.ymin):
        dy = y - (part.bounds.ymax + part.bounds.ymin) / 2
    else:
        dy = 0
    if not math.isinf(part.bounds.zmax) and not math.isinf(part.bounds.zmin):
        dz = z - (part.bounds.zmax + part.bounds.zmin) / 2
    else:
        dz = 0
    return move(part, dx, dy, dz)

################################################################################

@preserve_color
def rotate(part, angle, x0=0, y0=0):
    p = move(part, -x0, -y0, 0)
    angle *= math.pi/180
    ca, sa = math.cos(angle), math.sin(angle)
    nca, nsa = math.cos(-angle), math.sin(-angle)

    return move(p.map(Transform(
        '+*f%(ca)gX*f%(sa)gY'    % locals(),
        '+*f%(nsa)gX*f%(ca)gY'   % locals(),
        '+*f%(nca)gX*f%(nsa)gY'  % locals(),
        '+*f%(sa)gX*f%(nca)gY'   % locals())),
        x0, y0, 0)

################################################################################

@preserve_color
def reflect_x(part, x0=0):
    # X' = 2*x0-X
    # X  = 2*x0-X'
    return part.map(Transform(
        '-*f2f%gX' % x0, '',
        '-*f2f%gX' % x0, ''))

@preserve_color
def reflect_y(part, y0=0):
    # Y' = 2*y0-Y
    return part.map(Transform(
        '', '-*f2f%gY' % y0,
        '', '-*f2f%gY' % y0))

@preserve_color
def reflect_z(part, z0=0):
    # Z' = 2*z0-Z
    return part.map(Transform(
        '', '', '-*f2f%gZ' % z0,
        '', '', '-*f2f%gZ' % z0))

@preserve_color
def reflect_xy(part, x0=0, y0=0):
   # X' = x0 + (Y-y0)
   # Y' = y0 + (X-x0)
   # X = x0 + (Y'-y0)
   # Y = y0 + (X'-x0)
   return part.map(Transform(
      '+f%(x0)g-Yf%(y0)g' % locals(),
      '+f%(y0)g-Xf%(x0)g' % locals(),
      '+f%(x0)g-Yf%(y0)g' % locals(),
      '+f%(y0)g-Xf%(x0)g' % locals()))

@preserve_color
def reflect_yz(part, y0=0, z0=0):
   # Y' = y0 + (Z-z0)
   # Z' = z0 + (Y-y0)
   # Y = y0 + (Z'-z0)
   # Z = z0 + (Y'-y0)
   return part.map(Transform(
      'X',
      '+f%(y0)g-Zf%(z0)g' % locals(),
      '+f%(z0)g-Yf%(y0)g' % locals(),
      'X',
      '+f%(y0)g-Zf%(z0)g' % locals(),
      '+f%(z0)g-Yf%(y0)g' % locals()))

@preserve_color
def reflect_xz(part, x0=0, z0=0):
   # X' = x0 + (Z-z0)
   # Z' = z0 + (X-x0)
   # X = x0 + (Z'-z0)
   # Z = z0 + (X'-x0)
   return part.map(Transform(
      '+f%(x0)g-Zf%(z0)g' % locals(),
      'Y',
      '+f%(z0)g-Xf%(x0)g' % locals(),
      '+f%(x0)g-Zf%(z0)g' % locals(),
      'Y',
      '+f%(z0)g-Xf%(x0)g' % locals()))

################################################################################

@preserve_color
def scale_x(part, x0, sx):
    # X' = x0 + (X-x0)/sx
    return part.map(Transform(
        '+f%(x0)g/-Xf%(x0)gf%(sx)g' % locals()
                if x0 else '/Xf%g' % sx,
        'Y',
        '+f%(x0)g*f%(sx)g-Xf%(x0)g' % locals()
                if x0 else '*Xf%g' % sx,
        'Y'))

@preserve_color
def scale_y(part, y0, sy):
    # Y' = y0 + (Y-y0)/sy
    return part.map(Transform(
        'X',
        '+f%(y0)g/-Yf%(y0)gf%(sy)g' % locals()
                if y0 else '/Yf%g' % sy,
        'X',
        '+f%(y0)g*f%(sy)g-Yf%(y0)g' % locals()
                if y0 else '*Yf%g' % sy))

@preserve_color
def scale_z(part, z0, sz):
    # Z' = z0 + (Y-y0)/sz
    # Z  = (Z'-z0)*sz + z0
    return part.map(Transform(
        'X', 'Y',
        '+f%(z0)g/-Zf%(z0)gf%(sz)g' % locals()
                if z0 else '/Zf%g' % sz,
        'X', 'Y',
        '+f%(z0)g*f%(sz)g-Zf%(z0)g' % locals()
                if z0 else '*Zf%g' % sz))

@preserve_color
def scale_xy(part, x0, y0, sx, sy=None):
    # X' = x0 + (X-x0)/sx
    # Y' = y0 + (Y-y0)/sy
    # X  = (X'-x0)*sx + x0
    # Y  = (Y'-y0)*sy + y0
    if sy is None:
        sy = sx
    return part.map(Transform(
        '+f%(x0)g/-Xf%(x0)gf%(sx)g' % locals(),
        '+f%(y0)g/-Yf%(y0)gf%(sy)g' % locals(),
        '+f%(x0)g*f%(sx)g-Xf%(x0)g' % locals(),
        '+f%(y0)g*f%(sy)g-Yf%(y0)g' % locals()))

@preserve_color
def scale_xyz(part, x0, y0, z0, sx, sy, sz):
   # X' = x0 + (X-x0)/sx
   # Y' = y0 + (Y-y0)/sy
   # Z' = z0 + (Z-z0)/sz
   # X = x0 + (X'-x0)*sx
   # Y = y0 + (Y'-y0)*sy
   # Z = z0 + (Z'-z0)*sz
   return part.map(Transform(
      '+f%(x0)g/-Xf%(x0)gf%(sx)g' % locals(),
      '+f%(y0)g/-Yf%(y0)gf%(sy)g' % locals(),
      '+f%(z0)g/-Zf%(z0)gf%(sz)g' % locals(),
      '+f%(x0)g*-Xf%(x0)gf%(sx)g' % locals(),
      '+f%(y0)g*-Yf%(y0)gf%(sy)g' % locals(),
      '+f%(z0)g*-Zf%(z0)gf%(sz)g' % locals()))

@preserve_color
def scale_cos_xy_z(part, x0, y0, z0, z1, amp, off, t0, t1):
   dz = z1-z0
   t0 = math.radians(t0)
   t1 = math.radians(t1)
   # X' = x0 + (X-x0)/(off+amp*math.cos(theta0+(theta1-theta0)*(Z-z0)/dz))
   # X = x0 + (X'-x0)*(off+amp*math.cos(theta0+(theta1-theta0)*(Z-z0)/dz))
   # Y' = y0 + (Y-y0)/(off+amp*math.cos(theta0+(theta1-theta0)*(Z-z0)/dz))
   # Y = y0 + (Y'-y0)*(off+amp*math.cos(theta0+(theta1-theta0)*(Z-z0)/dz))
   return part.map(Transform(
      '/+f%(x0)g-Xf%(x0)g+f%(off)g*f%(amp)gc+f%(t0)g/*-f%(t1)gf%(t0)g-Zf%(z0)gf%(dz)g' % locals(),
      '/+f%(y0)g-Yf%(y0)g+f%(off)g*f%(amp)gc+f%(t0)g/*-f%(t1)gf%(t0)g-Zf%(z0)gf%(dz)g' % locals(),
      'Z',
      '*+f%(x0)g-Xf%(x0)g+f%(off)g*f%(amp)gc+f%(t0)g/*-f%(t1)gf%(t0)g-Zf%(z0)gf%(dz)g' % locals(),
      '*+f%(y0)g-Yf%(y0)g+f%(off)g*f%(amp)gc+f%(t0)g/*-f%(t1)gf%(t0)g-Zf%(z0)gf%(dz)g' % locals(),
      'Z'))

@preserve_color
def scale_cos_x_y(part, x0, y0, y1, amp, off, t0, t1):
   dy = y1 - y0
   t0 = math.radians(t0)
   t1 = math.radians(t1)
   # X' = x0 + (X-x0)/(off+amp*math.cos(theta0+(theta1-theta0)*(Y-y0)/dy))
   # X = x0 + (X'-x0)*(off+amp*math.cos(theta0+(theta1-theta0)*(Y-y0)/dy))
   return part.map(Transform(
      '/+f%(x0)g-Xf%(x0)g+f%(off)g*f%(amp)gc+f%(t0)g/*-f%(t1)gf%(t0)g-Yf%(y0)gf%(dy)g' % locals(),
      'Y',
      '*+f%(x0)g-Xf%(x0)g+f%(off)g*f%(amp)gc+f%(t0)g/*-f%(t1)gf%(t0)g-Yf%(y0)gf%(dy)g' % locals(),
      'Y'))

def scale_z_r(part, x0, y0, z0, r0, s0, r1, s1):
   dr = r1 - r0
   # Z' = z0 + (Z-z0)*dr/((s1-s0)*sqrt((X-x0)^2+(Y-y0)^2)-s1*r0+s0*r1)
   # Z = z0 + (Z'-z0)*((s1-s0)*sqrt((X-x0)^2+(Y-y0)^2)-s1*r0+s0*r1)/dr
   return part.map(Transform(
      'X', 'Y',
      '+f%(z0)g/*-Zf%(z0)gf%(dr)g+-*-f%(s1)gf%(s0)gr+q-Xf%(x0)gq-Yf%(y0)g*f%(s1)gf%(r0)g*f%(s0)gf%(r1)g' % locals(),
      'X', 'Y',
      '+f%(z0)g/*-Zf%(z0)g+-*-f%(s1)gf%(s0)gr+q-Xf%(x0)gq-Yf%(y0)g*f%(s1)gf%(r0)g*f%(s0)gf%(r1)gf%(dr)g' % locals()))

################################################################################

@preserve_color
def extrude_z(part, zmin, zmax):
    # max(part, max(zmin-Z, Z-zmax))
    return Shape(
            'am__f1%sa-f%gZ-Zf%g' % (part.math, zmin, zmax),
            part.bounds.xmin, part.bounds.ymin, zmin,
            part.bounds.xmax, part.bounds.ymax, zmax)

################################################################################

def loft_xy_z(a, b, zmin, zmax):
    """ Creates a blended loft between two shapes.

        Input shapes should be 2D (in the XY plane).
        The resulting loft will be shape a at zmin and b at zmax.
    """
    # ((z-zmin)/(zmax-zmin))*b + ((zmax-z)/(zmax-zmin))*a
    # In the prefix string below, we add caps at zmin and zmax then
    # factor out the division by (zmax - zmin)
    dz = zmax - zmin
    a_, b_ = a.math, b.math
    return Shape(('aa-Zf%(zmax)g-f%(zmin)g' +
                  'Z/+*-Zf%(zmin)g%(b_)s' +
                     '*-f%(zmax)gZ%(a_)sf%(dz)g') % locals(),
                min(a.bounds.xmin, b.bounds.xmin),
                min(a.bounds.ymin, b.bounds.ymin), zmin,
                max(a.bounds.xmax, b.bounds.xmax),
                max(a.bounds.ymax, b.bounds.ymax), zmax)

################################################################################

@preserve_color
def shear_x_y(part, ymin, ymax, dx0, dx1):
    dx = dx1 - dx0
    dy = ymax - ymin

    # X' = X-dx0-dx*(Y-ymin)/dy
    # X  = X'+dx0+(dx)*(Y-ymin)/dy
    return part.map(Transform(
            '--Xf%(dx0)g/*f%(dx)g-Yf%(ymin)gf%(dy)g' % locals(),
            'Y',
            '++Xf%(dx0)g/*f%(dx)g-Yf%(ymin)gf%(dy)g' % locals(),
            'Y'))

@preserve_color
def shear_xy_z(part, zmin, zmax, dx0, dy0, dx1, dy1):
    dx = dx1 - dx0
    dy = dy1 - dy0
    dz = zmax - zmin

    # X' = X-dx0-dx*(Z-zmin)/dz
    # Y' = Y-dy0-dy*(Z-zmin)/dz
    # X  = X'+dx0+(dx)*(Y-ymin)/dy
    # Y  = Y'+dy0+(dy)*(Y-ymin)/dy
    return part.map(Transform(
            '--Xf%(dx0)g/*f%(dx)g-Zf%(zmin)gf%(dz)g' % locals(),
            '--Yf%(dy0)g/*f%(dy)g-Zf%(zmin)gf%(dz)g' % locals(),
            '++Xf%(dx0)g/*f%(dx)g-Zf%(zmin)gf%(dz)g' % locals(),
            '++Yf%(dy0)g/*f%(dy)g-Zf%(zmin)gf%(dz)g' % locals()))

@preserve_color
def shear_cos_xy_z(part, z0, z1, ampx, offx, ampy, offy, t0, t1):
   dz = z1-z0
   t0 = math.radians(t0)
   t1 = math.radians(t1)
   # X' = X-(offx+ampx*math.cos(theta0+(theta1-theta0)*(Z-z0)/dz))
   # X = X'+(offx+ampx*math.cos(theta0+(theta1-theta0)*(Z-z0)/dz))
   # Y' = Y-(offy+ampy*math.cos(theta0+(theta1-theta0)*(Z-z0)/dz))
   # Y = Y'+(offy+ampy*math.cos(theta0+(theta1-theta0)*(Z-z0)/dz))
   return part.map(Transform(
      '-X+f%(offx)g*f%(ampx)gc+f%(t0)g/*-f%(t1)gf%(t0)g-Zf%(z0)gf%(dz)g' % locals(),
      '-Y+f%(offy)g*f%(ampy)gc+f%(t0)g/*-f%(t1)gf%(t0)g-Zf%(z0)gf%(dz)g' % locals(),
      'Z',
      '+X+f%(offx)g*f%(ampx)gc+f%(t0)g/*-f%(t1)gf%(t0)g-Zf%(z0)gf%(dz)g' % locals(),
      '+Y+f%(offy)g*f%(ampy)gc+f%(t0)g/*-f%(t1)gf%(t0)g-Zf%(z0)gf%(dz)g' % locals(),
      'Z'))

@preserve_color
def shear_cos_x_y(part,y0,y1,amp,off,t0,t1):
   dy = y1-y0
   t0 = math.radians(t0)
   t1 = math.radians(t1)
   # X' = X-(off+amp*math.cos(theta0+(theta1-theta0)*(Y-y0)/dy))
   # X = X'+(off+amp*math.cos(theta0+(theta1-theta0)*(Y-y0)/dy))
   return part.map(Transform(
      '-X+f%(off)g*f%(amp)gc+f%(t0)g/*-f%(t1)gf%(t0)g-Yf%(y0)gf%(dy)g' % locals(),
      'Y',
      '+X+f%(off)g*f%(amp)gc+f%(t0)g/*-f%(t1)gf%(t0)g-Yf%(y0)gf%(dy)g' % locals(),
      'Y'))
################################################################################

@preserve_color
def taper_x_y(part, x0, y0, y1, s0, s1):
    dy = y1 - y0
    ds = s1 - s0
    s0y1 = s0 * y1
    s1y0 = s1 * y0

    #   X'=x0+(X-x0)*(y1-y0)/(Y*(s1-s0)+s0*y1-s1*y0))
    #   X=(X'-x0)*(Y*(s1-s0)+s0*y1-s1*y0)/(y1-y0)+x0
    return part.map(Transform(
        '+f%(x0)g/*-Xf%(x0)gf%(dy)g-+*Yf%(ds)gf%(s0y1)gf%(s1y0)g' % locals(),
        'Y',
        '+f%(x0)g*-Xf%(x0)g/-+*Yf%(ds)gf%(s0y1)gf%(s1y0)gf%(dy)g' % locals(),
        'Y'))

@preserve_color
def iterate2d(part, i, j, dx, dy):
    return iterate3d(part, i, j, 1, dx, dy, 1)

@preserve_color
def iterate3d(part, i, j, k, dx, dy, dz):
    """ Tiles a part in the X, Y, and Z directions.
    """
    if i < 1 or j < 1 or k < 1:
        raise ValueError("Invalid value for iteration")

    #FIXME: if you try to do this all at once, segfault
    xiterate = functools.reduce(operator.or_, [move(part, a*dx, 0, 0) for a in range(i)])
    yiterate = functools.reduce(operator.or_, [move(xiterate, 0, b*dy, 0) for b in range(j)])
    ziterate = functools.reduce(operator.or_, [move(yiterate, 0, 0, c*dz) for c in range(k)])
    return ziterate

@preserve_color
def iterate_polar(part, x, y, n):
    """ Tiles a part by rotating it n times about x,y
    """

    if n < 1:
        raise ValueError("Invalid count for iteration")

    return functools.reduce(operator.or_,
            [rotate(part, 360./n * i, x, y)
             for i in range(n)])

################################################################################

def blend(a, b, amount):
    joint = a | b

    # sqrt(abs(a)) + sqrt(abs(b)) - amount
    fillet = Shape('-+rb%srb%sf%g' % (a.math, b.math, amount),
                   joint.bounds)
    return joint | fillet

def chamfer_union(a, b, r):
    """ Union with a 45-degree chamfer of width r along the seam.
        Assumes distance-like fields near the surfaces.
    """
    if r <= 0:
        return a | b
    # min(min(a, b), (a + b - r) * sqrt(1/2))
    s = 'ii%s%s*f0.7071067811865476+%s-%sf%g' % (a.math, b.math,
                                                  a.math, b.math, r)
    return Shape(s, (a | b).bounds)

def chamfer_intersection(a, b, r):
    """ Intersection with a 45-degree chamfer of width r along the seam.
    """
    if r <= 0:
        return a & b
    # max(max(a, b), (a + b + r) * sqrt(1/2))
    s = 'aa%s%s*f0.7071067811865476++%s%sf%g' % (a.math, b.math,
                                                  a.math, b.math, r)
    return Shape(s, (a & b).bounds)

def chamfer_difference(a, b, r):
    """ Difference with a 45-degree chamfer of width r along the cut edge.
    """
    if r <= 0:
        return a & ~b
    # max(max(a, -b), (a - b + r) * sqrt(1/2))
    s = 'aa%sn%s*f0.7071067811865476++%sn%sf%g' % (a.math, b.math,
                                                    a.math, b.math, r)
    return Shape(s, (a & ~b).bounds)

def fillet_union(a, b, r):
    """ Union with a rounded fillet of radius r along the seam.
        Assumes distance-like fields near the surfaces.
    """
    if r <= 0:
        return a | b
    # min(a, b) - max(r - |a - b|, 0)^2 / (4*r)
    s = '-i%s%s/qaf0-f%gb-%s%sf%g' % (a.math, b.math, r,
                                      a.math, b.math, 4*r)
    return Shape(s, (a | b).bounds)

def fillet_intersection(a, b, r):
    """ Intersection with a rounded fillet of radius r along the seam.
    """
    if r <= 0:
        return a & b
    # max(a, b) + max(r - |a - b|, 0)^2 / (4*r)
    s = '+a%s%s/qaf0-f%gb-%s%sf%g' % (a.math, b.math, r,
                                      a.math, b.math, 4*r)
    return Shape(s, (a & b).bounds)

def fillet_difference(a, b, r):
    """ Difference with a rounded fillet of radius r along the cut edge.
    """
    if r <= 0:
        return a & ~b
    # max(a, -b) + max(r - |a + b|, 0)^2 / (4*r)
    s = '+a%sn%s/qaf0-f%gb+%s%sf%g' % (a.math, b.math, r,
                                       a.math, b.math, 4*r)
    return Shape(s, (a & ~b).bounds)

def morph(a, b, weight):
    """ Morphs between two shapes.
    """
    # shape = weight*a+(1-weight)*b
    s = "+*f%g%s*f%g%s" % (weight, a.math, 1-weight, b.math)
    return Shape(s, (a | b).bounds)

################################################################################

def cylinder(x, y, zmin, zmax, r):
    return extrude_z(circle(x, y, r), zmin, zmax)

def cylinder_x(xmin, xmax, y, z,r):
   from fab.types import Shape, Transform
   # max(sqrt((Y-y)^2+(Z-z)^2)-r,max(xmin-X,X-xmax))
   return Shape(
      'a-r+q-Yf%(y)gq-Zf%(z)gf%(r)ga-f%(xmin)gX-Xf%(xmax)g' % locals(),
      xmin, y-r, z-r, xmax, y+r, z+r)

def cylinder_y(x, ymin, ymax, z, r):
   from fab.types import Shape, Transform
   # max(sqrt((X-x)^2+(Z-z)^2)-r,max(ymin-Y,Y-ymax))
   return Shape(
      'a-r+q-Xf%(x)gq-Zf%(z)gf%(r)ga-f%(ymin)gY-Yf%(ymax)g' % locals(),
      x-r, ymin, z-r, x+r, ymax,z+r)

def sphere(x, y, z, r):
    return Shape('-r++q%sq%sq%sf%g' % (('-Xf%g' % x) if x else 'X',
                                  ('-Yf%g' % y) if y else 'Y',
                                  ('-Zf%g' % z) if z else 'Z',
                                  r),
            x - r, y - r, z - r, x + r, y + r, z + r)

def cube(xmin, xmax, ymin, ymax, zmin, zmax):
    return extrude_z(rectangle(xmin, xmax, ymin, ymax), zmin, zmax)

def rounded_cube(xmin, xmax, ymin, ymax, zmin, zmax, r):
    """ Returns a cube with rounded corners and edges
        r is a roundedness fraction between 0 (not rounded)
        and 1 (completely rounded)
    """
    r *= min([xmax - xmin, ymax - ymin, zmax - zmin])/2
    s = (
        extrude_z(rectangle(xmin + r, xmax - r, ymin + r, ymax - r),
                  zmin, zmax) |
        extrude_z(rectangle(xmin, xmax, ymin + r, ymax - r) |
                  rectangle(xmin + r, xmax - r, ymin, ymax),
                  zmin + r, zmax - r)
    )
    for i in range(8):
        s |= sphere((xmin + r) if (i & 1) else (xmax - r),
                    (ymin + r) if (i & 2) else (ymax - r),
                    (zmin + r) if (i & 4) else (zmax - r), r)
    for i in range(4):
        s |= cylinder(
                (xmin + r) if (i & 1) else (xmax - r),
                (ymin + r) if (i & 2) else (ymax - r),
                zmin + r, zmax - r, r)
        s |= cylinder_x(
                xmin + r, xmax - r,
                (ymin + r) if (i & 1) else (ymax - r),
                (zmin + r) if (i & 2) else (zmax - r), r)
        s |= cylinder_y(
                (xmin + r) if (i & 1) else (xmax - r),
                ymin + r, ymax - r,
                (zmin + r) if (i & 2) else (zmax - r), r)
    return s

def cone(x, y, z0, z1, r):
    flipped = z1 < z0
    if flipped:
        z1 = 2*z0 - z1
    cyl = cylinder(x, y, z0, z1, r)
    out = taper_xy_z(cyl, x, y, z0, z1, 1.0, 0.0)
    return reflect_z(out, z0) if flipped else out

def pyramid(xmin, xmax, ymin, ymax, z0, z1):
    flipped = z1 < z0
    if flipped:
        z1 = 2*z0 - z1
    c = cube(xmin, xmax, ymin, ymax, z0, z1)
    out = taper_xy_z(c, (xmin+xmax)/2, (ymin+ymax)/2, z0, z1, 1, 0)
    return reflect_z(out, z0) if flipped else out

def torus_x(x, y, z, R, r):
   # sqrt((R - sqrt((Y-y)^2+(Z-z)^2))^2 + (X-x)^2)-r
   return move(Shape(
      '-r+q-f%(R)gr+qYqZqXf%(r)g' % locals(),
       -r, -(R + r), -(R + r), r, R + r, R + r), x, y, z)

def torus_y(x, y, z, R, r):
   # sqrt((R - sqrt((X-x)^2+(Z-z)^2))^2 + (Y-y)^2)-r
   return move(Shape(
      '-r+q-f%(R)gr+qXqZqYf%(r)g' % locals(),
       -(R+r), -r, -(R+r), R + r, r, R + r), x, y, z)

def torus_z(x, y, z, R, r):
   return move(Shape(
      '-r+q-f%(R)gr+qXqYqZf%(r)g' % locals(),
       -(R+r), -(R+r), -r, R + r, R + r, r), x, y, z)
################################################################################

# 3D shapes and operations

@preserve_color
def rotate_x(part, angle, y0=0, z0=0):
    p = move(part, 0, -y0, -z0)
    angle *= math.pi/180
    ca, sa = math.cos(angle), math.sin(angle)
    nca, nsa = math.cos(-angle), math.sin(-angle)

    return move(p.map(Transform(
        '', '+*f%(ca)gY*f%(sa)gZ'  % locals(),
            '+*f%(nsa)gY*f%(ca)gZ' % locals(),

        'X', '+*f%(nca)gY*f%(nsa)gZ' % locals(),
             '+*f%(sa)gY*f%(nca)gZ' % locals())),
        0, y0, z0)

@preserve_color
def rotate_y(part, angle, x0=0, z0=0):

    p = move(part, -x0, 0, -z0)
    angle *= math.pi/180
    ca, sa = math.cos(angle), math.sin(angle)
    nca, nsa = math.cos(-angle), math.sin(-angle)

    return move(p.map(Transform(
            '+*f%(ca)gX*f%(sa)gZ'  % locals(), 'Y',
            '+*f%(nsa)gX*f%(ca)gZ' % locals(),
            '+*f%(nca)gX*f%(nsa)gZ' % locals(), 'Y',
            '+*f%(sa)gX*f%(nca)gZ' % locals())),
            x0, 0, z0)

rotate_z = rotate

################################################################################

@preserve_color
def shear_x_z(part, z0, z1, dx0, dx1):
    #   X' = X-dx0-(dx1-dx0)*(Z-z0)/(z1-z0)
    #   X = X'+dx0+(dx1-dx0)*(Z-z0)/(z1-z0)
    return part.map(Transform(
        '--Xf%(dx0)g/*f%(dx)g-Zf%(z0)gf%(dz)g' % locals(), '', '',
        '++Xf%(dx0))g/*f%(dx)g-Zf%(z0)gf%(dz)g' % locals(), '', ''))

################################################################################

@preserve_color
def taper_xy_z(part, x0, y0, z0, z1, s0, s1):

    dz = z1 - z0

    # X' =  x0 +(X-x0)*dz/(s1*(Z-z0) + s0*(z1-Z))
    # Y' =  y0 +(Y-y0)*dz/(s1*(Z-z0) + s0*(z1-Z))
    # X  = (X' - x0)*(s1*(Z-z0) + s0*(z1-Z))/dz + x0
    # Y  = (Y' - y0)*(s1*(Z-z0) + s0*(z1-Z))/dz + y0
    return part.map(Transform(
        '+f%(x0)g/*-Xf%(x0)gf%(dz)g+*f%(s1)g-Zf%(z0)g*f%(s0)g-f%(z1)gZ'
            % locals(),
        '+f%(y0)g/*-Yf%(y0)gf%(dz)g+*f%(s1)g-Zf%(z0)g*f%(s0)g-f%(z1)gZ'
            % locals(),
        '',
        '+/*-Xf%(x0)g+*f%(s1)g-Zf%(z0)g*f%(s0)g-f%(z1)gZf%(dz)gf%(x0)g'
            % locals(),
        '+/*-Yf%(y0)g+*f%(s1)g-Zf%(z0)g*f%(s0)g-f%(z1)gZf%(dz)gf%(y0)g'
            % locals(),
        ''))

################################################################################

@preserve_color
def revolve_y(a):
    ''' Revolve a part in the XY plane about the Y axis. '''
    #   X' = +/- sqrt(X**2 + Z**2)
    pos = a.map(Transform('r+qXqZ', '', '', '', '', ''))
    neg = a.map(Transform('nr+qXqZ', '', '', '', '', ''))
    m = max(abs(a.bounds.xmin), abs(a.bounds.xmax))
    return Shape((pos | neg).math, -m, a.bounds.ymin, -m,
                                    m, a.bounds.ymax,  m)


@preserve_color
def revolve_x(a):
    ''' Revolve a part in the XY plane about the X axis. '''
    #   Y' = +/- sqrt(Y**2 + Z**2)
    pos = a.map(Transform('', 'r+qYqZ', '', '', '', ''))
    neg = a.map(Transform('', 'nr+qYqZ', '', '', '', ''))
    m = max(abs(a.bounds.ymin), abs(a.bounds.ymax))
    return Shape((pos | neg).math, a.bounds.xmin, -m, -m,
                                   a.bounds.xmax,  m,  m)

@preserve_color
def revolve_xy_x(a, y):
    """ Revolves the given shape about the x-axis
        (offset by the given y value)
    """
    return move(revolve_x(move(a, 0, -y)), 0, y)

@preserve_color
def revolve_xy_y(a, x):
    """ Revolves the given shape about the y-axis
        (offset by the given x value)
    """
    return move(revolve_y(move(a, -x, 0)), x, 0)

################################################################################

@preserve_color
def attract(part, x, y, z, r):

    # Shift the part so that it is centered
    part = move(part, -x, -y, -z)

    # exponential fallout value
    # x*(1 - exp(-sqrt(x**2 + y**2 + z**2) / r))
    d = '+f1xn/r++qXqYqZf%g' % r
    p = part.map(Transform(
        '*X'+d, '*Y'+d, '*Z'+d, '', '', ''))

    b = r/math.e
    return move(Shape(
        p.math,
        part.bounds.xmin - b, part.bounds.ymin - b, part.bounds.zmin - b,
        part.bounds.xmax + b, part.bounds.ymax + b, part.bounds.zmax + b),
        x, y, z)

@preserve_color
def repel(part, x, y, z, r):
    # Shift the part so that it is centered
    part = move(part, -x, -y, -z)

    # exponential fallout value
    # x*(1 - exp(-sqrt(x**2 + y**2 + z**2) / r))
    d = '-f1xn/r++qXqYqZf%g' % r
    p = part.map(Transform('*X'+d, '*Y'+d, '*Z'+d, '', '', ''))

    b = r/math.e
    return move(Shape(
        p.math,
        part.bounds.xmin - b, part.bounds.ymin - b, part.bounds.zmin - b,
        part.bounds.xmax + b, part.bounds.ymax + b, part.bounds.zmax + b),
        x, y, z)

################################################################################

@preserve_color
def cylinder_y_wrap(part, radius):
    tx = "(X / %(radius)f)" % locals()
    tz = "(Z / %(radius)f)" % locals()
    dist = "(sqrt( (%(tx)s)**2 + (%(tz)s)**2 ))" % locals()
    angle = "(atan2( %(tx)s, %(tz)s ))" % locals()

    Xfn = "=%(angle)s * %(radius)f;" % locals()
    Yfn = "_"
    Zfn = "=(%(dist)s - 1) * %(radius)f;" % locals()

    angle = "(X / %(radius)f)" % locals()
    r = "(%(radius)f + Z)" % locals()
    sina = "(sin(%(angle)s))" % locals()
    cosa = "(cos(%(angle)s))" % locals()

    Xfn_inv = "=%(r)s * %(sina)s;" % locals()
    Yfn_inv = ""
    Zfn_inv = "=%(r)s * %(cosa)s;" % locals()

    return part.map(Transform(  Xfn, Yfn, Zfn,
                                Xfn_inv, Yfn_inv, Zfn_inv))

@preserve_color
def twist_xy_z(part, x, y, z0, z1, t0, t1):
    # First, we'll move and scale so that the relevant part of the model
    # is at x=y=0 and scaled so that z is between 0 and 1.
    p1 = scale_z(move(part, -x, -y, -z0), 0, 1.0/(z1 - z0))

    t0 = math.pi * t0 / 180.0
    t1 = math.pi * t1 / 180.0

    # X' =  X*cos(t1*z + t0*(1-z)) + Y*sin(t1*z + t0*(1-z))
    # Y' = -X*sin(t1*z + t0*(1-z)) + Y*cos(t1*z + t0*(1-z))
    # X =  X*cos(t1*z + t0*(1-z)) - Y*sin(t1*z + t0*(1-z))
    # Y =  X*sin(t1*z + t0*(1-z)) + Y*cos(t1*z + t0*(1-z))
    p2 = p1.map(Transform(
        '+*Xc+*f%(t1)gZ*f%(t0)g-f1Z*Ys+*f%(t1)gZ*f%(t0)g-f1Z' % locals(),
        '+n*Xs+*f%(t1)gZ*f%(t0)g-f1Z*Yc+*f%(t1)gZ*f%(t0)g-f1Z' % locals(),
        '-*Xc+*f%(t1)gZ*f%(t0)g-f1Z*Ys+*f%(t1)gZ*f%(t0)g-f1Z' % locals(),
        '+*Xs+*f%(t1)gZ*f%(t0)g-f1Z*Yc+*f%(t1)gZ*f%(t0)g-f1Z' % locals()))

    return move(scale_z(p2, 0, z1 - z0), x, y, z0)

################################################################################

def function_prefix_xy(fn, xmin, xmax, ymin, ymax):
    """ Takes an arbitrary prefix math-string and makes it a function.
        Returns the function intersected with the given bounding rectangle.
    """
    return Shape(fn) & rectangle(xmin, xmax, ymin, ymax)

def function_prefix_xyz(fn, xmin, xmax, ymin, ymax, zmin, zmax):
    """ Takes an arbitrary prefix math-string and makes it a function.
        Returns the function intersected with the given bounding cube.
    """
    return Shape(fn) & cube(xmin, xmax, ymin, ymax, zmin, zmax)

################################################################################

def text(text, x, y, height=1, align='LB'):
    if text == '':
        return Shape()
    dx, dy = 0, -1
    text_shape = None

    for line in text.split('\n'):
        line_shape = None

        for c in line:
            if not c in _glyphs.keys():
                print('Warning:  Unknown character "%s"' % c)
            else:
                chr_math = move(_glyphs[c], dx, dy)
                if line_shape is None:  line_shape  = chr_math
                else:                   line_shape |= chr_math
                dx += _widths[c] + 0.1
        dx -= 0.1

        if line_shape is not None:
            if align[0] == 'L':
                pass
            elif align[0] == 'C':
                line_shape = move(line_shape, -dx / 2, 0)
            elif align[0] == 'R':
                line_shape = move(line_shape, -dx, 0)

            if text_shape is None:  text_shape  = line_shape
            else:                   text_shape |= line_shape

        dy -= 1.55
        dx = 0
    dy += 1.55
    if text_shape is None:  return None

    if align[1] == 'T':
        pass
    elif align[1] == 'B':
        text_shape = move(text_shape, 0, -dy,)
    elif align[1] == 'C':
        text_shape = move(text_shape, 0, -dy/2)

    if height != 1:
        text_shape = scale_xy(text_shape, 0, 0, height)

    return move(text_shape, x, y)


_glyphs = {}
_widths = {}

shape = triangle(0, 0, 0.35, 1, 0.1, 0)
shape |= triangle(0.1, 0, 0.35, 1, 0.45, 1)
shape |= triangle(0.35, 1, 0.45, 1, 0.8, 0)
shape |= triangle(0.7, 0, 0.35, 1, 0.8, 0)
shape |= rectangle(0.2, 0.6, 0.3, 0.4)
_widths['A'] = 0.8
_glyphs['A'] = shape


shape = circle(0.25, 0.275, 0.275)
shape &= ~circle(0.25, 0.275, 0.175)
shape = shear_x_y(shape, 0, 0.35, 0, 0.1)
shape |= rectangle(0.51, 0.61, 0, 0.35)
shape = move(shape, -0.05, 0)
_widths['a'] = 0.58
_glyphs['a'] = shape


shape = circle(0.3, 0.725, 0.275)
shape &= ~circle(0.3, 0.725, 0.175)
shape |= circle(0.3, 0.275, 0.275)
shape &= ~circle(0.3, 0.275, 0.175)
shape &= rectangle(0.3, 1, 0, 1)
shape |= rectangle(0, 0.1, 0, 1)
shape |= rectangle(0.1, 0.3, 0, 0.1)
shape |= rectangle(0.1, 0.3, 0.45, 0.55)
shape |= rectangle(0.1, 0.3, 0.9, 1)
_widths['B'] = 0.575
_glyphs['B'] = shape


shape = circle(0.25, 0.275, 0.275)
shape &= ~circle(0.25, 0.275, 0.175)
shape &= rectangle(0.25, 1, 0, 0.275) | rectangle(0, 1, 0.275, 1)
shape |= rectangle(0, 0.1, 0, 1)
shape |= rectangle(0.1, 0.25, 0, 0.1)
_widths['b'] = 0.525
_glyphs['b'] = shape


shape = circle(0.3, 0.7, 0.3) & ~circle(0.3, 0.7, 0.2)
shape |= circle(0.3, 0.3, 0.3) & ~circle(0.3, 0.3, 0.2)
shape &= ~rectangle(0, 0.6, 0.3, 0.7)
shape &= ~triangle(0.3, 0.5, 1.5, 1.5, 1.5, -0.5)
shape &= ~rectangle(0.3, 0.6, 0.2, 0.8)
shape |= rectangle(0, 0.1, 0.3, 0.7)
_widths['C'] = 0.57
_glyphs['C'] = shape


shape = circle(0.275, 0.275, 0.275)
shape &= ~circle(0.275, 0.275, 0.175)
shape &= ~triangle(0.275, 0.275, 1, 1, 1, -0.55)
_widths['c'] = 0.48
_glyphs['c'] = shape


shape = circle(0.1, 0.5, 0.5) & ~circle(0.1, 0.5, 0.4)
shape &= rectangle(0, 1, 0, 1)
shape |= rectangle(0, 0.1, 0, 1)
_widths['D'] = 0.6
_glyphs['D'] = shape


shape = reflect_x(_glyphs['b'], _widths['b']/2)
_widths['d'] = _widths['b']
_glyphs['d'] = shape


shape = rectangle(0, 0.1, 0, 1)
shape |= rectangle(0.1, 0.6, 0.9, 1)
shape |= rectangle(0.1, 0.6, 0, 0.1)
shape |= rectangle(0.1, 0.5, 0.45, 0.55)
_widths['E'] = 0.6
_glyphs['E'] = shape


shape = circle(0.275, 0.275, 0.275)
shape &= ~circle(0.275, 0.275, 0.175)
shape &= ~triangle(0.1, 0.275, 0.75, 0.275, 0.6, 0)
shape |= rectangle(0.05, 0.55, 0.225, 0.315)
shape &=  circle(0.275, 0.275, 0.275)
_widths['e'] = 0.55
_glyphs['e'] = shape


shape = rectangle(0, 0.1, 0, 1)
shape |= rectangle(0.1, 0.6, 0.9, 1)
shape |= rectangle(0.1, 0.5, 0.45, 0.55)
_widths['F'] = 0.6
_glyphs['F'] = shape


shape = circle(0.4, 0.75, 0.25) & ~circle(0.4, 0.75, 0.15)
shape &= rectangle(0, 0.4, 0.75, 1)
shape |= rectangle(0, 0.4, 0.45, 0.55)
shape |= rectangle(0.15, 0.25, 0, 0.75)
_widths['f'] = 0.4
_glyphs['f'] = shape


shape = circle(0.275, -0.1, 0.275)
shape &= ~circle(0.275, -0.1, 0.175)
shape &= rectangle(0, 0.55, -0.375, -0.1)
shape |= circle(0.275, 0.275, 0.275) & ~circle(0.275, 0.275, 0.175)
shape |= rectangle(0.45, 0.55, -0.1, 0.55)
_widths['g'] = 0.55
_glyphs['g'] = shape


shape = circle(0.3, 0.7, 0.3) & ~circle(0.3, 0.7, 0.2)
shape |= circle(0.3, 0.3, 0.3) & ~circle(0.3, 0.3, 0.2)
shape &= ~rectangle(0, 0.6, 0.3, 0.7)
shape |= rectangle(0, 0.1, 0.3, 0.7)
shape |= rectangle(0.5, 0.6, 0.3, 0.4)
shape |= rectangle(0.3, 0.6, 0.4, 0.5)
_widths['G'] = 0.6
_glyphs['G'] = shape


shape = rectangle(0, 0.1, 0, 1)
shape |= rectangle(0.5, 0.6, 0, 1)
shape |= rectangle(0.1, 0.5, 0.45, 0.55)
_widths['H'] = 0.6
_glyphs['H'] = shape


shape = circle(0.275, 0.275, 0.275)
shape &= ~circle(0.275, 0.275, 0.175)
shape &= rectangle(0, 0.55, 0.275, 0.55)
shape |= rectangle(0, 0.1, 0, 1)
shape |= rectangle(0.45, 0.55, 0, 0.275)
_widths['h'] = 0.55
_glyphs['h'] = shape


shape = rectangle(0, 0.5, 0, 0.1)
shape |= rectangle(0, 0.5, 0.9, 1)
shape |= rectangle(0.2, 0.3, 0.1, 0.9)
_widths['I'] = 0.5
_glyphs['I'] = shape


shape = rectangle(0.025, 0.125, 0, 0.55)
shape |= circle(0.075, 0.7, 0.075)
_widths['i'] = 0.15
_glyphs['i'] = shape


shape = circle(0.275, 0.275, 0.275)
shape &= ~circle(0.275, 0.275, 0.175)
shape &= rectangle(0, 0.55, 0, 0.275)
shape |= rectangle(0.45, 0.55, 0.275, 1)
_widths['J'] = 0.55
_glyphs['J'] = shape


shape = circle(0.0, -0.1, 0.275)
shape &= ~circle(0.0, -0.1, 0.175)
shape &= rectangle(0, 0.55, -0.375, -0.1)
shape |= rectangle(0.175, 0.275, -0.1, 0.55)
shape |= circle(0.225, 0.7, 0.075)
_widths['j'] = 0.3
_glyphs['j'] = shape


shape = rectangle(0, 0.6, 0, 1)
shape &= ~triangle(0.1, 1, 0.5, 1, 0.1, 0.6)
shape &= ~triangle(0.5, 0, 0.1, 0, 0.1, 0.4)
shape &= ~triangle(0.6, 0.95, 0.6, 0.05, 0.18, 0.5)
_widths['K'] = 0.6
_glyphs['K'] = shape


shape = rectangle(0, 0.5, 0, 1)
shape &= ~triangle(0.1, 1, 0.5, 1, 0.1, 0.45)
shape &= ~triangle(0.36, 0, 0.1, 0, 0.1, 0.25)
shape &= ~triangle(0.6, 1, 0.5, 0.0, 0.18, 0.35)
shape &= ~triangle(0.1, 1, 0.6, 1, 0.6, 0.5)
_widths['k'] = 0.5
_glyphs['k'] = shape


shape = rectangle(0, 0.6, 0, 0.1)
shape |= rectangle(0, 0.1, 0, 1)
_widths['L'] = 0.6
_glyphs['L'] = shape


shape = rectangle(0.025, 0.125, 0, 1)
_widths['l'] = 0.15
_glyphs['l'] = shape


shape = rectangle(0, 0.1, 0, 1)
shape |= rectangle(0.7, 0.8, 0, 1)
shape |= triangle(0, 1, 0.1, 1, 0.45, 0)
shape |= triangle(0.45, 0, 0.35, 0, 0, 1)
shape |= triangle(0.7, 1, 0.8, 1, 0.35, 0)
shape |= triangle(0.35, 0, 0.8, 1, 0.45, 0)
_widths['M'] = 0.8
_glyphs['M'] = shape


shape = circle(0.175, 0.35, 0.175) & ~circle(0.175, 0.35, 0.075)
shape |= circle(0.425, 0.35, 0.175) & ~circle(0.425, 0.35, 0.075)
shape &= rectangle(0, 0.65, 0.35, 0.65)
shape |= rectangle(0, 0.1, 0, 0.525)
shape |= rectangle(0.25, 0.35, 0, 0.35)
shape |= rectangle(0.5, 0.6, 0, 0.35)
_widths['m'] = 0.6
_glyphs['m'] = shape


shape = rectangle(0, 0.1, 0, 1)
shape |= rectangle(0.5, 0.6, 0, 1)
shape |= triangle(0, 1, 0.1, 1, 0.6, 0)
shape |= triangle(0.6, 0, 0.5, 0, 0, 1)
_widths['N'] = 0.6
_glyphs['N'] = shape


shape = circle(0.275, 0.275, 0.275)
shape &= ~circle(0.275, 0.275, 0.175)
shape &= rectangle(0, 0.55, 0.325, 0.55)
shape |= rectangle(0, 0.1, 0, 0.55)
shape |= rectangle(0.45, 0.55, 0, 0.325)
_widths['n'] = 0.55
_glyphs['n'] = shape


shape = circle(0.3, 0.7, 0.3) & ~circle(0.3, 0.7, 0.2)
shape |= circle(0.3, 0.3, 0.3) & ~circle(0.3, 0.3, 0.2)
shape &= ~rectangle(0, 0.6, 0.3, 0.7)
shape |= rectangle(0, 0.1, 0.3, 0.7)
shape |= rectangle(0.5, 0.6, 0.3, 0.7)
_widths['O'] = 0.6
_glyphs['O'] = shape


shape = circle(0.275, 0.275, 0.275)
shape &= ~circle(0.275, 0.275, 0.175)
_widths['o'] = 0.55
_glyphs['o'] = shape


shape = circle(0.3, 0.725, 0.275)
shape &= ~circle(0.3, 0.725, 0.175)
shape &= rectangle(0.3, 1, 0, 1)
shape |= rectangle(0, 0.1, 0, 1)
shape |= rectangle(0.1, 0.3, 0.45, 0.55)
shape |= rectangle(0.1, 0.3, 0.9, 1)
_widths['P'] = 0.575
_glyphs['P'] = shape


shape = circle(0.275, 0.275, 0.275)
shape &= ~circle(0.275, 0.275, 0.175)
shape |= rectangle(0, 0.1, -0.375, 0.55)
_widths['p'] = 0.55
_glyphs['p'] = shape


shape = circle(0.3, 0.7, 0.3) & ~circle(0.3, 0.7, 0.2)
shape |= circle(0.3, 0.3, 0.3) & ~circle(0.3, 0.3, 0.2)
shape &= ~rectangle(0, 0.6, 0.3, 0.7)
shape |= rectangle(0, 0.1, 0.3, 0.7)
shape |= rectangle(0.5, 0.6, 0.3, 0.7)
shape |= triangle(0.5, 0.1, 0.6, 0.1, 0.6, 0)
shape |= triangle(0.5, 0.1, 0.5, 0.3, 0.6, 0.1)
_widths['Q'] = 0.6
_glyphs['Q'] = shape


shape = circle(0.275, 0.275, 0.275) & ~circle(0.275, 0.275, 0.175)
shape |= rectangle(0.45, 0.55, -0.375, 0.55)
_widths['q'] = 0.55
_glyphs['q'] = shape


shape = circle(0.3, 0.725, 0.275)
shape &= ~circle(0.3, 0.725, 0.175)
shape &= rectangle(0.3, 1, 0, 1)
shape |= rectangle(0, 0.1, 0, 1)
shape |= rectangle(0.1, 0.3, 0.45, 0.55)
shape |= rectangle(0.1, 0.3, 0.9, 1)
shape |= triangle(0.3, 0.5, 0.4, 0.5, 0.575, 0)
shape |= triangle(0.475, 0.0, 0.3, 0.5, 0.575, 0)
_widths['R'] = 0.575
_glyphs['R'] = shape


shape = circle(0.55, 0, 0.55) & ~scale_x(circle(0.55, 0, 0.45), 0.55, 0.8)
shape &= rectangle(0, 0.55, 0, 0.55)
shape = scale_x(shape, 0, 0.7)
shape |= rectangle(0, 0.1, 0, 0.55)
_widths['r'] = 0.385
_glyphs['r'] = shape


shape = circle(0.275, 0.725, 0.275)
shape &= ~circle(0.275, 0.725, 0.175)
shape &= ~rectangle(0.275, 0.55, 0.45, 0.725)
shape |= reflect_x(reflect_y(shape, 0.5), .275)
_widths['S'] = 0.55
_glyphs['S'] = shape


shape = circle(0.1625, 0.1625, 0.1625)
shape &= ~scale_x(circle(0.165, 0.165, 0.0625), 0.165, 1.5)
shape &= ~rectangle(0, 0.1625, 0.1625, 0.325)
shape |= reflect_x(reflect_y(shape, 0.275), 0.1625)
shape = scale_x(shape, 0, 1.5)
_widths['s'] = 0.4875
_glyphs['s'] = shape


shape = rectangle(0, 0.6, 0.9, 1) | rectangle(0.25, 0.35, 0, 0.9)
_widths['T'] = 0.6
_glyphs['T'] = shape


shape = circle(0.4, 0.25, 0.25) & ~circle(0.4, 0.25, 0.15)
shape &= rectangle(0, 0.4, 0, 0.25)
shape |= rectangle(0, 0.4, 0.55, 0.65)
shape |= rectangle(0.15, 0.25, 0.25, 1)
_widths['t'] = 0.4
_glyphs['t'] = shape

# Fix kerning for lowercase t
shape = circle(0.4, 0.25, 0.25) & ~circle(0.4, 0.25, 0.15)
shape &= rectangle(0, 0.4, 0, 0.25)
shape |= rectangle(0, 0.4, 0.55, 0.65)
shape |= rectangle(0.15, 0.25, 0.25, 1)
shape = translate(shape, -.125, 0, 0)
_widths['𝑡'] = 0.275
_glyphs['𝑡'] = shape

shape = circle(0.3, 0.3, 0.3) & ~circle(0.3, 0.3, 0.2)
shape &= rectangle(0, 0.6, 0, 0.3)
shape |= rectangle(0, 0.1, 0.3, 1)
shape |= rectangle(0.5, 0.6, 0.3, 1)
_widths['U'] = 0.6
_glyphs['U'] = shape


shape = circle(0.275, 0.275, 0.275) & ~circle(0.275, 0.275, 0.175)
shape &= rectangle(0, 0.55, 0, 0.275)
shape |= rectangle(0, 0.1, 0.275, 0.55)
shape |= rectangle(0.45, 0.55, 0, 0.55)
_widths['u'] = 0.55
_glyphs['u'] = shape


shape = triangle(0, 1, 0.1, 1, 0.35, 0)
shape |= triangle(0.35, 0, 0.25, 0, 0, 1)
shape |= reflect_x(shape, 0.3)
_widths['V'] = 0.6
_glyphs['V'] = shape


shape = triangle(0, 0.55, 0.1, 0.55, 0.35, 0)
shape |= triangle(0.35, 0, 0.25, 0, 0, 0.55)
shape |= reflect_x(shape, 0.3)
_widths['v'] = 0.6
_glyphs['v'] = shape


shape = triangle(0, 1, 0.1, 1, 0.25, 0)
shape |= triangle(0.25, 0, 0.15, 0, 0, 1)
shape |= triangle(0.15, 0, 0.35, 1, 0.45, 1)
shape |= triangle(0.45, 1, 0.25, 0, 0.15, 0)
shape |= reflect_x(shape, 0.4)
_widths['W'] = 0.8
_glyphs['W'] = shape


shape = triangle(0, 0.55, 0.1, 0.55, 0.25, 0)
shape |= triangle(0.25, 0, 0.15, 0, 0, 0.55)
shape |= triangle(0.15, 0, 0.35, 0.5, 0.45, 0.5)
shape |= triangle(0.45, 0.5, 0.25, 0, 0.15, 0)
shape |= reflect_x(shape, 0.4)
_widths['w'] = 0.8
_glyphs['w'] = shape


shape = triangle(0, 1, 0.125, 1, 0.8, 0)
shape |= triangle(0.8, 0, 0.675, 0, 0, 1)
shape |= reflect_x(shape, 0.4)
_widths['X'] = 0.8
_glyphs['X'] = shape


shape = triangle(0, 0.55, 0.125, 0.55, 0.55, 0)
shape |= triangle(0.55, 0, 0.425, 0, 0, 0.55)
shape |= reflect_x(shape, 0.275)
_widths['x'] = 0.55
_glyphs['x'] = shape


shape = triangle(0, 1, 0.1, 1, 0.45, 0.5)
shape |= triangle(0.45, 0.5, 0.35, 0.5, 0, 1)
shape |= reflect_x(shape, 0.4)
shape |= rectangle(0.35, 0.45, 0, 0.5)
_widths['Y'] = 0.8
_glyphs['Y'] = shape


shape = triangle(0, 0.55, 0.1, 0.55, 0.325, 0)
shape |= triangle(0.325, 0, 0.225, 0, 0, 0.55)
shape |= reflect_x(shape, 0.275) | move(reflect_x(shape, 0.275), -0.225, -0.55)
shape &= rectangle(0, 0.55, -0.375, 0.55)
_widths['y'] = 0.55
_glyphs['y'] = shape


shape = rectangle(0, 0.6, 0, 1)
shape &= ~triangle(0, 0.1, 0, 0.9, 0.45, 0.9)
shape &= ~triangle(0.6, 0.1, 0.15, 0.1, 0.6, 0.9)
_widths['Z'] = 0.6
_glyphs['Z'] = shape


shape = rectangle(0, 0.6, 0, 0.55)
shape &= ~triangle(0, 0.1, 0, 0.45, 0.45, 0.45)
shape &= ~triangle(0.6, 0.1, 0.15, 0.1, 0.6, 0.45)
_widths['z'] = 0.6
_glyphs['z'] = shape


shape = Shape("f1.0", 0, 0.55, 0, 1)
_widths[' '] = 0.55
_glyphs[' '] = shape


shape = circle(0.075, 0.075, 0.075)
shape = scale_y(shape, 0.075, 3)
shape &= rectangle(0.0, 0.15, -0.15, 0.075)
shape &= ~triangle(0.075, 0.075, 0.0, -0.15, -0.5, 0.075)
shape |= circle(0.1, 0.075, 0.075)
_widths[','] = 0.175
_glyphs[','] = shape


shape = circle(0.075, 0.075, 0.075)
_widths['.'] = 0.15
_glyphs['.'] = shape


shape = rectangle(0, 0.1, 0.55, 0.8)
_widths["'"] = 0.1
_glyphs["'"] = shape

shape = rectangle(0, 0.1, 0.55, 0.8) | rectangle(0.2, 0.3, 0.55, 0.8)
_widths['"'] = 0.3
_glyphs['"'] = shape


shape = circle(0.075, 0.15, 0.075) | circle(0.075, 0.45, 0.075)
_widths[':'] = 0.15
_glyphs[':'] = shape


shape = circle(0.075, 0.15, 0.075)
shape = scale_y(shape, 0.15, 3)
shape &= rectangle(0.0, 0.15, -0.075, 0.15)
shape &= ~triangle(0.075, 0.15, 0.0, -0.075, -0.5, 0.15)
shape |= circle(0.075, 0.45, 0.075)
shape |= circle(0.1, 0.15, 0.075)
_widths[';'] = 0.15
_glyphs[';'] = shape


shape = rectangle(0.025, 0.125, 0.3, 1)
shape |= circle(0.075, 0.075, 0.075)
_widths['!'] = 0.1
_glyphs['!'] = shape


shape = rectangle(0.05, 0.4, 0.35, 0.45)
_widths['-'] = 0.45
_glyphs['-'] = shape


shape = circle(0, 0.4, 0.6) & ~scale_x(circle(0, 0.4, 0.5), 0, 0.7)
shape &= rectangle(0, 0.6, -0.2, 1)
shape = scale_x(shape, 0, 1/2.)
_widths[')'] = 0.3
_glyphs[')'] = shape


shape = circle(0.6, 0.4, 0.6) & ~scale_x(circle(0.6, 0.4, 0.5), 0.6, 0.7)
shape &= rectangle(0, 0.6, -0.2, 1)
shape = scale_x(shape, 0, 1/2.)
_widths['('] = 0.3
_glyphs['('] = shape


shape = rectangle(0, 0.3, 0, 1)
shape &= ~circle(0, 1, 0.2)
shape &= ~rectangle(0, 0.2, 0, 0.7)
_widths['1'] = 0.3
_glyphs['1'] = shape


shape = circle(0.275, .725, .275)
shape &= ~circle(0.275, 0.725, 0.175)
shape &= ~rectangle(0, 0.55, 0, 0.725)
shape |= rectangle(0, 0.55, 0, 0.1)
shape |= triangle(0, 0.1, 0.45, 0.775, 0.55, 0.725)
shape |= triangle(0, 0.1, 0.55, 0.725, 0.125, 0.1)
_widths['2'] = 0.55
_glyphs['2'] = shape


shape = circle(0.3, 0.725, 0.275)
shape &= ~circle(0.3, 0.725, 0.175)
shape |= circle(0.3, 0.275, 0.275)
shape &= ~circle(0.3, 0.275, 0.175)
shape &= ~rectangle(0, 0.275, 0.275, 0.725)
_widths['3'] = 0.55
_glyphs['3'] = shape


shape = triangle(-0.10, 0.45, 0.4, 1, 0.4, 0.45)
shape |= rectangle(0.4, 0.5, 0, 1)
shape &= ~triangle(0.4, 0.85, 0.4, 0.55, 0.1, 0.55)
shape &= rectangle(0, 0.5, 0, 1)
_widths['4'] = 0.5
_glyphs['4'] = shape


shape = circle(0.325, 0.325, 0.325) & ~circle(0.325, 0.325, 0.225)
shape &= ~rectangle(0, 0.325, 0.325, 0.65)
shape |= rectangle(0, 0.325, 0.55, 0.65)
shape |= rectangle(0, 0.1, 0.55, 1)
shape |= rectangle(0.1, 0.65, 0.9, 1)
_widths['5'] = 0.65
_glyphs['5'] = shape


shape = circle(0.275, 0.725, 0.275) & ~scale_y(circle(0.275, 0.725, 0.175), .725, 1.2)
shape &= rectangle(0, 0.55, 0.725, 1)
shape &= ~triangle(0.275, 0.925, 0.55, 0.9, 0.55, 0.725)
shape = scale_y(shape, 1, 2)
shape = scale_x(shape, 0, 1.1)
shape &= ~rectangle(0.275, 0.65, 0., 0.7)
shape |= rectangle(0, 0.1, 0.275, 0.45)
shape |= circle(0.275, 0.275, 0.275) & ~circle(0.275, 0.275, 0.175)
_widths['6'] = 0.55
_glyphs['6'] = shape


shape = rectangle(0, 0.6, 0.9, 1)
shape |= triangle(0, 0, 0.475, 0.9, 0.6, 0.9)
shape |= triangle(0, 0, 0.6, 0.9, 0.125, 0)
_widths['7'] = 0.6
_glyphs['7'] = shape


shape = circle(0.3, 0.725, 0.275)
shape &= ~circle(0.3, 0.725, 0.175)
shape |= circle(0.3, 0.275, 0.275)
shape &= ~circle(0.3, 0.275, 0.175)
_widths['8'] = 0.55
_glyphs['8'] = shape


shape = reflect_x(reflect_y(_glyphs['6'], 0.5), _widths['6']/2)
_widths['9'] = _widths['6']
_glyphs['9'] = shape


shape = circle(0.5, 0.5, 0.5) & ~scale_x(circle(0.5, 0.5, 0.4), 0.5, 0.7**0.5)
shape = scale_x(shape, 0, 0.7)
_widths['0'] = 0.7
_glyphs['0'] = shape


shape = rectangle(0., 0.5, 0.45, 0.55)
shape |= rectangle(0.2, 0.3, 0.25, 0.75)
_widths['+'] = 0.55
_glyphs['+'] = shape


shape = triangle(0, 0, 0.425, 1, 0.55, 1)
shape |= triangle(0, 0, 0.55, 1, 0.125, 0)
_widths['/'] = 0.55
_glyphs['/'] = shape


shape = circle(0.275, 0.725, 0.275) & ~circle(0.275, 0.725, 0.175)
shape &= ~rectangle(0, 0.275, 0.45, 0.725)
shape |= rectangle(0.225, 0.325, 0.3, 0.55)
shape |= circle(0.275, 0.075, 0.075)
_widths['?'] = 0.55
_glyphs['?'] = shape


shape = rectangle(0.5, 0.6, 0, 0.8)
shape |= rectangle(0, 0.1, 0, 0.8)
shape |= rectangle(0.1, 0.6, 0, 0.1)
shape |= rectangle(0.1, 0.6, 0.4, 0.5)
shape = shear_x_y(shape, 0, 0.8, 0, 0.1)
_widths['#'] = 0.75
_glyphs['#'] = shape


shapeA = circle(0.25, 0.275, 0.275)
shapeA &= ~circle(0.25, 0.275, 0.175)
shapeA = shear_x_y(shapeA, 0, 0.35, 0, 0.1)
shapeA |= rectangle(0.51, 0.61, 0, 0.35)
shapeA = move(shapeA, -0.1, 0)
shapeA = scale_x(shapeA, 0, 0.8)
shapeE = circle(0.275, 0.275, 0.275)
shapeE &= ~circle(0.275, 0.275, 0.175)
shapeE &= ~triangle(0.1, 0.275, 0.75, 0.275, 0.6, 0)
shapeE |= rectangle(0.05, 0.55, 0.225, 0.315)
shapeE &=  circle(0.275, 0.275, 0.275)
shapeE = scale_x(shapeE, 0, 0.9)
shapeE = move(shapeE, 0.32, 0)
_widths['æ'] = 1
_glyphs['æ'] = shapeA | shapeE


del shapeA
del shapeE
del shape

################################################################################
# 3D primitives and deformations (phase 1 node campaign)
################################################################################

################################################################################

# already defined above.
#
# Math strings are prefix notation.  Opcodes:
#   + - * /   binary
#   i=min a=max p=pow
#   b=abs q=square r=sqrt n=neg x=exp s=sin c=cos t=tan
#   X Y Z coords, f<number> constant, m<xf><yf><zf><body> coordinate remap
################################################################################

################################################################################
# 3D primitives
################################################################################

def half_space(px, py, pz, nx, ny, nz):
    """ Half-space bounded by a plane through (px,py,pz) with the given
        outward normal (nx,ny,nz).  Solid on the side the normal points
        AWAY from (field = signed distance, positive along +normal).
        The normal is normalised internally, so the field is a true
        distance and the region is unbounded.
    """
    n = math.sqrt(nx*nx + ny*ny + nz*nz)
    if n == 0:
        raise ValueError("half_space normal must be non-zero")
    nx, ny, nz = nx/n, ny/n, nz/n
    # nx*(X-px) + ny*(Y-py) + nz*(Z-pz)
    return Shape(
        '++*f%g-Xf%g*f%g-Yf%g*f%g-Zf%g' % (nx, px, ny, py, nz, pz))


def capsule(x0, y0, z0, x1, y1, z1, r):
    """ Round-capped segment (exact SDF) between two 3D points, radius r.
        IQ sdCapsule: h = clamp(dot(pa,ba)/dot(ba,ba), 0, 1);
        length(pa - ba*h) - r.
    """
    r = abs(r)
    bx, by, bz = x1 - x0, y1 - y0, z1 - z0
    bb = bx*bx + by*by + bz*bz
    if bb == 0:                       # degenerate segment -> sphere
        return sphere(x0, y0, z0, r)
    pax, pay, paz = '-Xf%g' % x0, '-Yf%g' % y0, '-Zf%g' % z0
    # D = dot(pa, ba)
    D = '++*f%g%s*f%g%s*f%g%s' % (bx, pax, by, pay, bz, paz)
    # H = clamp(D/bb, 0, 1) = min(max(D/bb, 0), 1)
    H = 'ia/%sf%gf0f1' % (D, bb)
    cx = '-%s*f%g%s' % (pax, bx, H)
    cy = '-%s*f%g%s' % (pay, by, H)
    cz = '-%s*f%g%s' % (paz, bz, H)
    # sqrt(cx^2 + cy^2 + cz^2) - r
    return Shape(
        '-r++q%sq%sq%sf%g' % (cx, cy, cz, r),
        min(x0, x1) - r, min(y0, y1) - r, min(z0, z1) - r,
        max(x0, x1) + r, max(y0, y1) + r, max(z0, z1) + r)


def capped_cone_z(z0, z1, r0, r1):
    """ Z-axis frustum: radius r0 at z0, r1 at z1, centred on the Z axis.
        Field is CSG-style (not a normalised distance): the lateral wall
        is sqrt(X^2+Y^2) - radius(Z) intersected with the end caps.
    """
    if z1 == z0:
        raise ValueError("capped_cone_z needs z0 != z1")
    k = (r1 - r0) / (z1 - z0)
    rad = '+f%g*f%g-Zf%g' % (r0, k, z0)         # r0 + k*(Z-z0)
    lateral = '-r+qXqY%s' % rad                 # sqrt(X^2+Y^2) - rad
    caps = 'a-f%gZ-Zf%g' % (z0, z1)             # max(z0-Z, Z-z1)
    rm = max(abs(r0), abs(r1))
    return Shape('a%s%s' % (lateral, caps),
                 -rm, -rm, min(z0, z1), rm, rm, max(z0, z1))


def rounded_cylinder_z(x, y, z0, z1, r, rr):
    """ Cylinder along Z centred at (x,y) from z0 to z1, radius r, with the
        top and bottom rim rounded by radius rr.  Exact distance to the
        revolved rounded-rectangle profile.
    """
    r = abs(r)
    H = (z1 - z0) / 2.0
    zc = (z0 + z1) / 2.0
    rr = min(abs(rr), r, abs(H))
    rho = 'r+q-Xf%gq-Yf%g' % (x, y)             # sqrt((X-x)^2+(Y-y)^2)
    qx = '-%sf%g' % (rho, r - rr)               # rho - (r-rr)
    qy = '-b-Zf%gf%g' % (zc, H - rr)            # |Z-zc| - (H-rr)
    inner = 'ia%s%sf0' % (qx, qy)               # min(max(qx,qy), 0)
    outer = 'r+qa%sf0qa%sf0' % (qx, qy)         # len(max(qx,0), max(qy,0))
    # inner + outer - rr
    return Shape('-+%s%sf%g' % (inner, outer, rr),
                 x - r, y - r, min(z0, z1), x + r, y + r, max(z0, z1))


def hex_prism_z(x, y, z0, z1, r):
    """ Hexagonal prism along Z, centred at (x,y), flat top/bottom edges
        (nut-trap orientation).  r is the centre-to-vertex (circum) radius.
    """
    return extrude_z(move(polygon_radius(0, 0, r, 6), x, y),
                     min(z0, z1), max(z0, z1))


def tri_prism_z(x, y, z0, z1, r):
    """ Triangular prism along Z, centred at (x,y), flat bottom edge.
        r is the centre-to-vertex (circum) radius.
    """
    return extrude_z(move(polygon_radius(0, 0, r, 3), x, y),
                     min(z0, z1), max(z0, z1))


def cut_sphere(x, y, z, r, h):
    """ Dome: a sphere of radius r centred at (x,y,z) cut by the plane
        Z = z + h, keeping the lower portion.  h in [-r, r]; h>=r keeps the
        whole sphere, h<=-r is empty.  Field is CSG-style (sign-correct).
    """
    r = abs(r)
    # max(sphere, Z-(z+h))
    sph = '-r++q-Xf%gq-Yf%gq-Zf%gf%g' % (x, y, z, r)
    plane = '-Zf%g' % (z + h)
    zt = z + max(-r, min(r, h))
    return Shape('a%s%s' % (sph, plane),
                 x - r, y - r, z - r, x + r, y + r, zt)


def slab(xmin=None, xmax=None, ymin=None, ymax=None, zmin=None, zmax=None):
    """ Axis-aligned region bounded by any subset of the six planes.
        Pass None (or +/-inf) for an axis end to leave it open.  With all
        six finite this is identical to `cube`; with open ends it is an
        unbounded clipping region (a slab / half-space stack).
    """
    def fin(v):
        return v is not None and not math.isinf(v)

    cons = []
    if fin(xmin): cons.append('-f%gX' % xmin)   # xmin - X
    if fin(xmax): cons.append('-Xf%g' % xmax)   # X - xmax
    if fin(ymin): cons.append('-f%gY' % ymin)
    if fin(ymax): cons.append('-Yf%g' % ymax)
    if fin(zmin): cons.append('-f%gZ' % zmin)
    if fin(zmax): cons.append('-Zf%g' % zmax)

    if not cons:
        return Shape()                           # all of space

    m = cons[0]
    for c in cons[1:]:
        m = 'a' + m + c                          # fold with max

    inf = float('inf')
    return Shape(m,
                 xmin if fin(xmin) else -inf,
                 ymin if fin(ymin) else -inf,
                 zmin if fin(zmin) else -inf,
                 xmax if fin(xmax) else inf,
                 ymax if fin(ymax) else inf,
                 zmax if fin(zmax) else inf)


def box_exact(xmin, xmax, ymin, ymax, zmin, zmax):
    """ Axis-aligned box with a single-formula EXACT signed distance field
        (IQ sdBox): length(max(|p-c|-b, 0)) + min(max component, 0).
        A drop-in exact-field alternative to `cube`.
    """
    cx, cy, cz = (xmin + xmax)/2.0, (ymin + ymax)/2.0, (zmin + zmax)/2.0
    bx, by, bz = (xmax - xmin)/2.0, (ymax - ymin)/2.0, (zmax - zmin)/2.0
    qx = '-b-Xf%gf%g' % (cx, bx)                 # |X-cx| - bx
    qy = '-b-Yf%gf%g' % (cy, by)
    qz = '-b-Zf%gf%g' % (cz, bz)
    outside = 'r++qa%sf0qa%sf0qa%sf0' % (qx, qy, qz)
    inside = 'ia%sa%s%sf0' % (qx, qy, qz)        # min(max(qx,max(qy,qz)),0)
    return Shape('+%s%s' % (outside, inside),
                 xmin, ymin, zmin, xmax, ymax, zmax)


################################################################################
# Deforms
################################################################################

@preserve_color
def bend_x(part, radius, x0=0):
    """ Arc bend driven by the X coordinate (IQ cheap-bend), curving the
        shape in the XY plane about the line X=x0, Y=0.  `radius` is the
        bend radius (angle at position X is (X-x0)/radius).
    """
    if radius == 0:
        return part
    k = 1.0 / radius
    A = '*f%g-Xf%g' % (k, x0)                    # angle = k*(X-x0)
    U = '-Xf%g' % x0
    # X' = x0 + (cos(A)*U - sin(A)*Y);  Y' = sin(A)*U + cos(A)*Y
    Xf = '+f%g-*c%s%s*s%sY' % (x0, A, U, A)
    Yf = '+*s%s%s*c%sY' % (A, U, A)
    p = part.map(Transform(Xf, Yf, '', '', '', ''))

    b = part.bounds
    Rmax = max(math.hypot(xx - x0, yy)
               for xx in (b.xmin, b.xmax) for yy in (b.ymin, b.ymax))
    return Shape(p.math,
                 x0 - Rmax, -Rmax, b.zmin,
                 x0 + Rmax,  Rmax, b.zmax)


@preserve_color
def bend_y(part, radius, y0=0):
    """ Arc bend driven by the Y coordinate, curving the shape in the XY
        plane about the line Y=y0, X=0.  `radius` is the bend radius.
    """
    if radius == 0:
        return part
    k = 1.0 / radius
    A = '*f%g-Yf%g' % (k, y0)                    # angle = k*(Y-y0)
    V = '-Yf%g' % y0
    # X' = cos(A)*X - sin(A)*V;  Y' = y0 + sin(A)*X + cos(A)*V
    Xf = '-*c%sX*s%s%s' % (A, A, V)
    Yf = '+f%g+*s%sX*c%s%s' % (y0, A, A, V)
    p = part.map(Transform(Xf, Yf, '', '', '', ''))

    b = part.bounds
    Rmax = max(math.hypot(xx, yy - y0)
               for xx in (b.xmin, b.xmax) for yy in (b.ymin, b.ymax))
    return Shape(p.math,
                 -Rmax, y0 - Rmax, b.zmin,
                  Rmax, y0 + Rmax, b.zmax)


def _twirl_angle(amount, radius):
    """ Falloff angle expression: radians(amount) * exp(-norm/radius),
        where norm = sqrt(X^2+Y^2+Z^2) in the centred frame. """
    a = math.radians(amount)
    norm = 'r++qXqYqZ'
    return '*f%gxn/%sf%g' % (a, norm, radius)


def _twirl_bounds(part, plane_axes):
    """ Conservative bounds for a norm-preserving twirl.  plane_axes names
        the two axes that rotate; their extent collapses to the enclosing
        disc, the third axis is left unchanged. """
    b = part.bounds
    lo = {'x': b.xmin, 'y': b.ymin, 'z': b.zmin}
    hi = {'x': b.xmax, 'y': b.ymax, 'z': b.zmax}
    u, v = plane_axes
    R = max(math.hypot(uu, vv)
            for uu in (lo[u], hi[u]) for vv in (lo[v], hi[v]))
    lo2, hi2 = dict(lo), dict(hi)
    lo2[u], hi2[u], lo2[v], hi2[v] = -R, R, -R, R
    return (lo2['x'], lo2['y'], lo2['z'], hi2['x'], hi2['y'], hi2['z'])


@preserve_color
def twirl_x(part, amount, radius, x=0, y=0, z=0):
    """ Localised twist about the X axis with exponential falloff
        (libfive centered_twirl).  `amount` in degrees, `radius` sets the
        falloff distance from the centre point (x,y,z).
    """
    part = move(part, -x, -y, -z)
    A = _twirl_angle(amount, radius)
    # X unchanged; Y' = cos(A)*Y + sin(A)*Z;  Z' = cos(A)*Z - sin(A)*Y
    Yf = '+*c%sY*s%sZ' % (A, A)
    Zf = '-*c%sZ*s%sY' % (A, A)
    p = part.map(Transform('', Yf, Zf, '', '', ''))
    return move(Shape(p.math, *_twirl_bounds(part, ('y', 'z'))), x, y, z)


@preserve_color
def twirl_y(part, amount, radius, x=0, y=0, z=0):
    """ Localised twist about the Y axis with exponential falloff. """
    part = move(part, -x, -y, -z)
    A = _twirl_angle(amount, radius)
    # Y unchanged; X' = cos(A)*X - sin(A)*Z;  Z' = cos(A)*Z + sin(A)*X
    Xf = '-*c%sX*s%sZ' % (A, A)
    Zf = '+*c%sZ*s%sX' % (A, A)
    p = part.map(Transform(Xf, '', Zf, '', '', ''))
    return move(Shape(p.math, *_twirl_bounds(part, ('x', 'z'))), x, y, z)


@preserve_color
def twirl_z(part, amount, radius, x=0, y=0, z=0):
    """ Localised twist about the Z axis with exponential falloff. """
    part = move(part, -x, -y, -z)
    A = _twirl_angle(amount, radius)
    # Z unchanged; X' = cos(A)*X + sin(A)*Y;  Y' = cos(A)*Y - sin(A)*X
    Xf = '+*c%sX*s%sY' % (A, A)
    Yf = '-*c%sY*s%sX' % (A, A)
    p = part.map(Transform(Xf, Yf, '', '', '', ''))
    return move(Shape(p.math, *_twirl_bounds(part, ('x', 'y'))), x, y, z)


@preserve_color
def revolve_z(a):
    """ Revolve a 2D shape in the XZ plane (X = radius, Z = height) about
        the Z axis.  Mirrors revolve_x / revolve_y.
    """
    #   X' = +/- sqrt(X**2 + Y**2)
    pos = a.map(Transform('r+qXqY', '', '', '', '', ''))
    neg = a.map(Transform('nr+qXqY', '', '', '', '', ''))
    m = max(abs(a.bounds.xmin), abs(a.bounds.xmax))
    return Shape((pos | neg).math, -m, -m, a.bounds.zmin,
                                    m,  m, a.bounds.zmax)


@preserve_color
def elongate(part, hx, hy, hz):
    """ Stretch the centre of a shape by inserting a slab of half-width
        (hx,hy,hz) per axis: q = p - clamp(p, -h, h), evaluate at q.
    """
    def cl(coord, h):
        if h == 0:
            return coord
        # coord - clamp(coord, -h, h) = coord - min(max(coord,-h), h)
        return '-%sia%sf%gf%g' % (coord, coord, -h, h)

    Xf, Yf, Zf = cl('X', hx), cl('Y', hy), cl('Z', hz)
    p = part.map(Transform(Xf, Yf, Zf, '', '', ''))
    b = part.bounds
    return Shape(p.math,
                 b.xmin - hx, b.ymin - hy, b.zmin - hz,
                 b.xmax + hx, b.ymax + hy, b.zmax + hz)

################################################################################
# 2D stroke and shape kit (phase 1 node campaign)
################################################################################

################################################################################

# Prefix math conventions (see existing shapes): operators precede operands.
#   + - * /  binary          i=min a=max p=pow
#   b=abs q=square r=sqrt n=neg   x=exp s/c/t=sin/cos/tan  S/C/T=asin/acos/atan
#   X Y Z coords   f<number> constant   (negative consts emit as f-1.5, ok)
# 2D shapes: field is independent of Z; z bounds are ±inf (Shape 4-arg form,
# exactly like circle/rectangle).
################################################################################

def _sub(coord, v):
    """coord if v == 0 else (coord - v), as a prefix subexpression."""
    return coord if v == 0 else '-%sf%g' % (coord, v)


def _min_chain(exprs):
    """Prefix string for min(exprs[0], exprs[1], ...)."""
    acc = exprs[-1]
    for e in reversed(exprs[:-1]):
        acc = 'i' + e + acc
    return acc


def _max_chain(exprs):
    """Prefix string for max(exprs[0], exprs[1], ...)."""
    acc = exprs[-1]
    for e in reversed(exprs[:-1]):
        acc = 'a' + e + acc
    return acc


################################################################################
# Stroke kit
################################################################################

def segment(x0, y0, x1, y1, w):
    """ Round-capped line segment between (x0,y0) and (x1,y1), full width w.
        Exact signed-distance field (IQ sdSegment):
            d = length(pa - ba*clamp(dot(pa,ba)/dot(ba,ba),0,1)) - w/2
    """
    bax, bay = x1 - x0, y1 - y0
    dd = bax * bax + bay * bay
    w2 = abs(w) / 2.0
    if dd == 0:                       # degenerate: a disk
        return circle(x0, y0, w2)

    Xx, Yy = _sub('X', x0), _sub('Y', y0)
    # h = clamp( (bax*(X-x0)+bay*(Y-y0)) / dd , 0, 1)
    num = '+*f%g%s*f%g%s' % (bax, Xx, bay, Yy)
    h = 'ia/%sf%gf0f1' % (num, dd)               # min(max(u,0),1)
    px = '-%s*f%g%s' % (Xx, bax, h)
    py = '-%s*f%g%s' % (Yy, bay, h)
    d = '-r+q%sq%sf%g' % (px, py, w2)
    return Shape(d,
                 min(x0, x1) - w2, min(y0, y1) - w2,
                 max(x0, x1) + w2, max(y0, y1) + w2)


def polyline(points, w):
    """ Open polyline through a list of [x,y] points, stroked with round joins
        and caps at full width w (union of segments; python-side unroll).
    """
    pts = [list(p) for p in points]
    if len(pts) == 0:
        return Shape()
    if len(pts) == 1:
        return circle(pts[0][0], pts[0][1], abs(w) / 2.0)
    segs = [segment(pts[i][0], pts[i][1], pts[i + 1][0], pts[i + 1][1], w)
            for i in range(len(pts) - 1)]
    return functools.reduce(operator.or_, segs)


def polygon(points):
    """ Filled polygon from an ordered vertex list (list of [x,y]).
        Generalizes `triangle`: the field is -min over the signed edge
        half-planes, so the result is the intersection of the edge half-spaces
        -- i.e. exact for CONVEX polygons.  Non-convex inputs are clipped to
        the intersection of their edge lines.  Winding is normalized to
        clockwise automatically (like triangle).  Field is non-euclidean
        (un-normalized edge distances), matching triangle.
    """
    # Drop consecutive duplicate vertices (zero-length edges break the field).
    pts = []
    for p in points:
        p = [float(p[0]), float(p[1])]
        if not pts or pts[-1] != p:
            pts.append(p)
    if len(pts) > 1 and pts[0] == pts[-1]:
        pts.pop()
    if len(pts) < 3:
        raise ValueError("polygon needs at least 3 distinct vertices")

    # Signed area > 0  =>  counter-clockwise; flip to clockwise.
    area = 0.0
    n = len(pts)
    for i in range(n):
        x0, y0 = pts[i]
        x1, y1 = pts[(i + 1) % n]
        area += x0 * y1 - x1 * y0
    if area > 0:
        pts = pts[::-1]

    def edge(x, y, dx, dy):
        # dy*(X-x) - dx*(Y-y)
        return '-*f%g%s*f%g%s' % (dy, _sub('X', x), dx, _sub('Y', y))

    edges = []
    for i in range(n):
        x0, y0 = pts[i]
        x1, y1 = pts[(i + 1) % n]
        edges.append(edge(x0, y0, x1 - x0, y1 - y0))

    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]
    return Shape('n' + _min_chain(edges),
                 min(xs), min(ys), max(xs), max(ys))


def _bezier_quad_flatten(p0, p1, p2, tol, out):
    """ Recursive adaptive subdivision of a quadratic bezier.
        Appends interior/end points (not p0) to `out`.
    """
    # Perpendicular distance from control point p1 to chord p0->p2.
    ax, ay = p2[0] - p0[0], p2[1] - p0[1]
    L = math.hypot(ax, ay)
    if L == 0:
        dev = math.hypot(p1[0] - p0[0], p1[1] - p0[1])
    else:
        dev = abs(ax * (p0[1] - p1[1]) - ay * (p0[0] - p1[0])) / L
    # For a quadratic, max curve deviation from the chord is <= dev/2.
    if dev <= 2 * tol or L < 1e-12:
        out.append(list(p2))
        return
    # de Casteljau split at t=0.5
    m01 = [(p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2]
    m12 = [(p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2]
    mid = [(m01[0] + m12[0]) / 2, (m01[1] + m12[1]) / 2]
    _bezier_quad_flatten(p0, m01, mid, tol, out)
    _bezier_quad_flatten(mid, m12, p2, tol, out)


def bezier_quad_points(x0, y0, x1, y1, x2, y2, tol=0.02):
    """ Returns the flattened polyline point list for a quadratic bezier
        with control points p0,p1,p2 and chordal tolerance `tol`.
    """
    p0, p1, p2 = [x0, y0], [x1, y1], [x2, y2]
    pts = [list(p0)]
    _bezier_quad_flatten(p0, p1, p2, tol, pts)
    return pts


def bezier_quad(x0, y0, x1, y1, x2, y2, w, tol=0.02):
    """ Stroked quadratic bezier (control points p0,p1,p2), full width w.
        APPROXIMATION: the curve is adaptively flattened to a polyline
        (chordal tolerance `tol`) and stroked as a union of round-capped
        segments.  Choice: polyline approximation (not the exact closed form).
    """
    return polyline(bezier_quad_points(x0, y0, x1, y1, x2, y2, tol), w)


################################################################################
# Shape set
################################################################################

def star(x, y, n, r_inner, r_outer, angle=90):
    """ n-pointed star centered at (x,y), alternating between r_outer (points)
        and r_inner (valleys).  `angle` (degrees) orients the first point.
        Built as a fan of triangles from the center (star is star-shaped
        about its center), unrolled python-side.
    """
    n = int(n)
    a0 = math.radians(angle)
    verts = []
    for k in range(2 * n):
        r = r_outer if (k % 2 == 0) else r_inner
        a = a0 + k * math.pi / n
        verts.append((x + r * math.cos(a), y + r * math.sin(a)))
    tris = [triangle(x, y, verts[k][0], verts[k][1],
                     verts[(k + 1) % (2 * n)][0], verts[(k + 1) % (2 * n)][1])
            for k in range(2 * n)]
    return functools.reduce(operator.or_, tris)


def _wedge(Xx, Yy, a0, a1):
    """ Angular-wedge field (negative inside) for the sector from a0 to a1
        (radians, CCW), given prefix (X-x) and (Y-y) subexpressions.
    """
    span = (a1 - a0) % (2 * math.pi)
    # left ray a0 : inside is cross(d0,p) >= 0  -> field = sin a0*Xx - cos a0*Yy
    hp0 = '-*f%g%s*f%g%s' % (math.sin(a0), Xx, math.cos(a0), Yy)
    # right ray a1: inside is cross(d1,p) <= 0 -> field = cos a1*Yy - sin a1*Xx
    hp1 = '-*f%g%s*f%g%s' % (math.cos(a1), Yy, math.sin(a1), Xx)
    # span <= pi : intersection (max); span > pi : union (min)
    return ('a' if span <= math.pi + 1e-12 else 'i') + hp0 + hp1


def pie(x, y, r, a0, a1):
    """ Pie / circular sector: disk of radius r clipped to the angular range
        [a0, a1] (degrees, CCW).  Constructive (circle intersected with an
        angular wedge); sign-correct, field is not an exact SDF.
    """
    ra0, ra1 = math.radians(a0), math.radians(a1)
    span = (ra1 - ra0) % (2 * math.pi)
    disk = circle(x, y, r)
    if span == 0:                     # full turn
        return disk
    Xx, Yy = _sub('X', x), _sub('Y', y)
    field = 'a' + disk.math + _wedge(Xx, Yy, ra0, ra1)
    return Shape(field, x - r, y - r, x + r, y + r)


sector = pie


def arc(x, y, r, a0, a1, w):
    """ Stroked circular arc of radius r spanning [a0,a1] (degrees, CCW),
        full stroke width w, with round caps.  Constructive: annulus of
        half-width w/2 clipped to the wedge, unioned with round end caps.
        Sign-correct; field is not an exact SDF.
    """
    w2 = abs(w) / 2.0
    ra0, ra1 = math.radians(a0), math.radians(a1)
    span = (ra1 - ra0) % (2 * math.pi)
    Xx, Yy = _sub('X', x), _sub('Y', y)
    # ring: |dist - r| - w/2
    ring = '-b-r+q%sq%sf%gf%g' % (Xx, Yy, r, w2)
    if span == 0:                     # full ring, no caps
        return Shape(ring,
                     x - r - w2, y - r - w2, x + r + w2, y + r + w2)
    body = 'a%s%s' % (ring, _wedge(Xx, Yy, ra0, ra1))
    cap0 = circle(x + r * math.cos(ra0), y + r * math.sin(ra0), w2).math
    cap1 = circle(x + r * math.cos(ra1), y + r * math.sin(ra1), w2).math
    field = _min_chain([body, cap0, cap1])
    return Shape(field,
                 x - r - w2, y - r - w2, x + r + w2, y + r + w2)


def annulus(x, y, r_outer, r_inner):
    """ Annulus / 2D ring between r_inner and r_outer.  Exact SDF:
            d = |length(p) - rmid| - halfwidth
    """
    rmid = (r_outer + r_inner) / 2.0
    half = (r_outer - r_inner) / 2.0
    Xx, Yy = _sub('X', x), _sub('Y', y)
    field = '-b-r+q%sq%sf%gf%g' % (Xx, Yy, rmid, half)
    return Shape(field, x - r_outer, y - r_outer, x + r_outer, y + r_outer)


ring2d = annulus


def trapezoid(width_bottom, width_top, height, x=0, y=0):
    """ Isosceles trapezoid centered at (x,y): bottom edge width_bottom at
        y-height/2, top edge width_top at y+height/2.  Convex -> built via
        `polygon` (non-euclidean field, sign-correct).
    """
    hb, ht, h2 = width_bottom / 2.0, width_top / 2.0, height / 2.0
    verts = [[x - hb, y - h2], [x + hb, y - h2],
             [x + ht, y + h2], [x - ht, y + h2]]
    return polygon(verts)


def rhombus(rx, ry, x=0, y=0):
    """ Axis-aligned rhombus (diamond) with half-diagonals rx (horizontal) and
        ry (vertical), centered at (x,y).  Convex -> built via `polygon`.
    """
    verts = [[x + rx, y], [x, y + ry], [x - rx, y], [x, y - ry]]
    return polygon(verts)


def oriented_box(x0, y0, x1, y1, thickness):
    """ Sharp rectangle whose long axis runs from (x0,y0) to (x1,y1), with the
        given thickness (across the axis).  Exact SDF (IQ sdOrientedBox).
    """
    l = math.hypot(x1 - x0, y1 - y0)
    th2 = abs(thickness) / 2.0
    if l == 0:
        return rectangle(x0 - th2, x0 + th2, y0 - th2, y0 + th2)
    dx, dy = (x1 - x0) / l, (y1 - y0) / l
    cx, cy = (x0 + x1) / 2.0, (y0 + y1) / 2.0
    l2 = l / 2.0
    Cx, Cy = _sub('X', cx), _sub('Y', cy)
    qx = '+*f%g%s*f%g%s' % (dx, Cx, dy, Cy)        #  dx*Cx + dy*Cy
    qy = '+*f%g%s*f%g%s' % (-dy, Cx, dx, Cy)       # -dy*Cx + dx*Cy
    ax = '-b%sf%g' % (qx, l2)                       # |qx| - l/2
    ay = '-b%sf%g' % (qy, th2)                      # |qy| - th/2
    term1 = 'r+qa%sf0qa%sf0' % (ax, ay)            # length(max(a,0))
    term2 = 'ia%s%sf0' % (ax, ay)                  # min(max(ax,ay),0)
    field = '+%s%s' % (term1, term2)
    # bounds from the four corners
    px, py = -dy, dx
    corners = [(cx + sx * l2 * dx + sy * th2 * px,
                cy + sx * l2 * dy + sy * th2 * py)
               for sx in (-1, 1) for sy in (-1, 1)]
    xs = [c[0] for c in corners]
    ys = [c[1] for c in corners]
    return Shape(field, min(xs), min(ys), max(xs), max(ys))


def vesica(x, y, r, d):
    """ Vesica / lens: intersection of two disks of radius r whose centers are
        at (x-d, y) and (x+d, y).  Requires d < r.  Sign-correct; the max of
        two circle SDFs is exact outside and on the boundary.
    """
    return circle(x - d, y, r) & circle(x + d, y, r)


def crescent(x, y, r_outer, r_inner, offset):
    """ Crescent / moon: disk of radius r_outer at (x,y) with a disk of radius
        r_inner at (x+offset, y) subtracted.  Constructive (difference).
    """
    return circle(x, y, r_outer) & ~circle(x + offset, y, r_inner)


moon = crescent


def cross(arm_length, thickness, rounding=0, x=0, y=0):
    """ Plus / cross centered at (x,y): a horizontal and a vertical bar, each
        of half-length arm_length and half-thickness `thickness`.  `rounding`
        > 0 rounds the interior (reflex) corners via a fillet of that radius
        (sdCross semantics); the four outer tips stay square.
    """
    L, T = arm_length, thickness
    harm = rectangle(x - L, x + L, y - T, y + T)
    varm = rectangle(x - T, x + T, y - L, y + L)
    if rounding > 0:
        return fillet_union(harm, varm, rounding)
    return harm | varm


def ellipse_approx(x, y, a, b):
    """ Axis-aligned ellipse with semi-axes a (x) and b (y), centered (x,y).
        APPROXIMATION: this is a unit circle scaled by (a,b), so the field is
        NON-EUCLIDEAN -- sqrt((X/a)^2+(Y/b)^2)-1 is not true distance to the
        ellipse (error grows with eccentricity).  Zero-set is exact; use for
        boolean geometry, not for uniform offsets.
    """
    Xx, Yy = _sub('X', x), _sub('Y', y)
    field = '-r+q/%sf%gq/%sf%gf1' % (Xx, a, Yy, b)
    return Shape(field, x - a, y - b, x + a, y + b)

################################################################################
# Functional parts kit (phase 1 node campaign)
################################################################################

# ============================================================================
# Design notes
# ------------
# * Positioning/orientation is baked directly into the prefix math string
#   (via _remap / coordinate substitution and by building primitives at their
#   final coordinates) rather than through move()/rotate() coordinate maps.
#   This keeps the fields exactly composable and lets the compound shapes be
#   sign-sampled by a plain prefix evaluator.
# * NEGATIVE shapes: the screw-hole / nut-trap / groove features model the
#   material to be REMOVED.  Wire them into a Difference node's `b` input.
#   Each such docstring says "NEGATIVE".  The bosses/standoffs/clips/tongues
#   are POSITIVE solids (wire into Union).
# ============================================================================


# ---------------------------------------------------------------------------
# tiny string helpers (no coordinate maps -> fully evaluable / bakeable)
# ---------------------------------------------------------------------------

def _off(coord, val):
    """ Prefix for (coord - val), collapsing to bare coord when val == 0. """
    return coord if val == 0 else '-%sf%g' % (coord, val)


def _remap(math_str, ux, uy, uz):
    """ Substitute the coordinate opcodes X, Y, Z in a prefix math string
        with the prefix sub-expressions ux, uy, uz (done simultaneously).

        Only the single-letter coordinate opcodes X/Y/Z are touched; constants
        (f...) and every other opcode - including the trig opcodes S/C/T - are
        left untouched, so this is a safe static coordinate remap that produces
        a valid prefix string with no `m` (map) node.
    """
    tmp = math_str.replace('X', '\x00').replace('Y', '\x01').replace('Z', '\x02')
    return tmp.replace('\x00', ux).replace('\x01', uy).replace('\x02', uz)


def _cone_z(x, y, zb, zt, r_bottom, r_top):
    """ Truncated cone (frustum) aligned with +Z, evaluable (no map).
        radius(Z) = r_bottom + (r_top-r_bottom)*(Z-zb)/(zt-zb).
        Returns a solid: negative inside, positive outside.
    """
    k = (r_top - r_bottom) / (zt - zb)
    radius = '+f%g*f%g%s' % (r_bottom, k, _off('Z', zb))
    radial = '-r+q%sq%s%s' % (_off('X', x), _off('Y', y), radius)
    caps = 'a-f%gZ-Zf%g' % (zb, zt)
    rmax = max(abs(r_bottom), abs(r_top))
    return Shape('a%s%s' % (radial, caps),
                 x - rmax, y - rmax, zb, x + rmax, y + rmax, zt)


def _extrude_prism(profile2d, amin, amax):
    """ Extrude a 2D (X,Y-only) profile along +Z between amin and amax WITHOUT
        emitting a coordinate map (unlike extrude_z, which pins the body's Z to
        1 via `m__f1`).  Map-free so the result can be safely coordinate-swapped
        by _remap to orient the prism along an arbitrary axis.
    """
    caps = 'a-f%gZ-Zf%g' % (amin, amax)
    return Shape('a%s%s' % (profile2d.math, caps),
                 profile2d.bounds.xmin, profile2d.bounds.ymin, amin,
                 profile2d.bounds.xmax, profile2d.bounds.ymax, amax)


_EPS = 1e-3   # small overshoot so subtracted negatives break the surface cleanly


# ---------------------------------------------------------------------------
# 1. Teardrop hole  (overhang-safe horizontal hole)
# ---------------------------------------------------------------------------

def teardrop_profile(r):
    """ 2D teardrop cross-section in the XY plane, circle centred at origin,
        peak toward +Y.  circle(r) unioned with a 90-degree (45-deg walls)
        triangular peak whose flanks are tangent at +/-45 degrees, so it
        prints without support when the peak points up.
    """
    r = abs(r)
    t = r * math.sqrt(0.5)          # tangent point offset (r*sin45 == r*cos45)
    apex = r * math.sqrt(2.0)       # peak height above centre
    peak = triangle(-t, t, t, t, 0.0, apex)
    return circle(0, 0, r) | peak


def teardrop_hole(x, y, z, r, length, axis='X'):
    """ Overhang-safe horizontal through-hole: a teardrop-section prism.

        NEGATIVE shape - wire into a Difference `b` input to cut a supportless
        horizontal bolt/rod hole in an FDM print.

        Parameters
        ----------
        x, y, z : centre of the circular cross-section (mm).
        r       : bore radius (mm).
        length  : hole length along `axis`, centred on (x,y,z).
        axis    : 'X' or 'Y' (horizontal); 'Z' also accepted.  The teardrop
                  peak always points toward +Z (up) so the 45-degree flanks
                  self-support during printing.
    """
    prof = teardrop_profile(r)
    tube = _extrude_prism(prof, -length / 2.0, length / 2.0)
    m = tube.math
    apex = r * math.sqrt(2.0)
    if axis == 'X':
        # profile-X <- worldY, profile-Y(peak) <- worldZ, length <- worldX
        mm = _remap(m, _off('Y', y), _off('Z', z), _off('X', x))
        bounds = (x - length / 2, y - r, z - r,
                  x + length / 2, y + r, z + apex)
    elif axis == 'Y':
        mm = _remap(m, _off('X', x), _off('Z', z), _off('Y', y))
        bounds = (x - r, y - length / 2, z - r,
                  x + r, y + length / 2, z + apex)
    elif axis == 'Z':
        mm = _remap(m, _off('X', x), _off('Y', y), _off('Z', z))
        bounds = (x - r, y - r, z - length / 2,
                  x + apex, y + r, z + length / 2)
    else:
        raise ValueError("axis must be 'X', 'Y' or 'Z'")
    return Shape(mm, *bounds)


# ---------------------------------------------------------------------------
# 2. Screw-hole family  (all NEGATIVE - subtract from your part)
# ---------------------------------------------------------------------------
#
# Close-fit ("close") clearance-hole diameters, mm.  ISO 273 "close" column /
# common maker-space practice (also matches ASME/DIN normal-to-close).  HIGH
# confidence for M2-M12.
CLEARANCE_CLOSE_D = {
    'M2': 2.4, 'M2.5': 2.9, 'M3': 3.4, 'M4': 4.5, 'M5': 5.5,
    'M6': 6.6, 'M8': 9.0, 'M10': 11.0, 'M12': 13.5,
}

# Countersunk flat-head THEORETICAL head diameter (dk, sharp-cornered) for a
# 90-degree head - DIN 965 / ISO 7046.  The theoretical (sharp) diameter is the
# correct value for cutting the cone recess (it is where the cone meets the top
# face).  Actual machined dk runs ~0.3-0.6 mm smaller.  MEDIUM-HIGH confidence.
COUNTERSINK_HEAD_D = {
    'M2': 4.0, 'M2.5': 5.0, 'M3': 6.0, 'M4': 8.0, 'M5': 10.0,
    'M6': 12.0, 'M8': 16.0, 'M10': 20.0, 'M12': 24.0,
}
COUNTERSINK_ANGLE = 90.0   # included angle, DIN 965

# Socket-head cap-screw head diameter dk (max) and head height k (max):
# DIN 912 / ISO 4762.  HIGH confidence.
SOCKET_HEAD_D = {
    'M2': 3.8, 'M2.5': 4.5, 'M3': 5.5, 'M4': 7.0, 'M5': 8.5,
    'M6': 10.0, 'M8': 13.0, 'M10': 16.0, 'M12': 18.0,
}
SOCKET_HEAD_H = {
    'M2': 2.0, 'M2.5': 2.5, 'M3': 3.0, 'M4': 4.0, 'M5': 5.0,
    'M6': 6.0, 'M8': 8.0, 'M10': 10.0, 'M12': 12.0,
}


def _clear_d(size, clearance_d):
    if clearance_d is not None:
        return clearance_d
    if size not in CLEARANCE_CLOSE_D:
        raise ValueError('unknown screw size %r (have %s)'
                         % (size, sorted(CLEARANCE_CLOSE_D)))
    return CLEARANCE_CLOSE_D[size]


def clearance_hole(x, y, z_top, length, size='M3', clearance_d=None):
    """ Plain cylindrical clearance shaft, drilled DOWN from the top face at
        z_top by `length` (bore axis == Z).

        NEGATIVE shape.  Diameter from the ISO 273 "close" preset for `size`
        (override with clearance_d).  Default configuration is M3.
    """
    d = _clear_d(size, clearance_d)
    return cylinder(x, y, z_top - length, z_top + _EPS, d / 2.0)


def counterbore_hole(x, y, z_top, length, size='M3', bore_depth=None,
                     clearance_d=None, head_d=None):
    """ Clearance shaft + cylindrical recess for a socket-head cap screw
        (DIN 912 / ISO 4762).  Recess is at the top face (z_top), sunk by
        bore_depth (default = DIN 912 head height for `size`).

        NEGATIVE shape.  Default configuration is M3.
    """
    d = _clear_d(size, clearance_d)
    hd = head_d if head_d is not None else SOCKET_HEAD_D[size]
    depth = bore_depth if bore_depth is not None else SOCKET_HEAD_H[size]
    shaft = cylinder(x, y, z_top - length, z_top + _EPS, d / 2.0)
    recess = cylinder(x, y, z_top - depth, z_top + _EPS, hd / 2.0)
    return shaft | recess


def countersink_hole(x, y, z_top, length, size='M3', angle=None,
                     clearance_d=None, head_d=None):
    """ Clearance shaft + conical recess for a flat-head screw
        (DIN 965 / ISO 7046, 90-degree head by default).

        NEGATIVE shape.  The cone runs from the theoretical head diameter at
        the top face (z_top) down to the clearance diameter; cone depth is set
        by the included `angle`.  Default configuration is M3.
    """
    d = _clear_d(size, clearance_d)
    hd = head_d if head_d is not None else COUNTERSINK_HEAD_D[size]
    ang = COUNTERSINK_ANGLE if angle is None else angle
    half = math.radians(ang / 2.0)
    cone_depth = (hd / 2.0 - d / 2.0) / math.tan(half)
    shaft = cylinder(x, y, z_top - length, z_top + _EPS, d / 2.0)
    cone = _cone_z(x, y, z_top - cone_depth, z_top + _EPS,
                   d / 2.0, hd / 2.0 + _EPS)
    return shaft | cone


# ---------------------------------------------------------------------------
# hex-prism primitive  (nut-trap / vent core)   -- local helper, merge dedupes
# ---------------------------------------------------------------------------

def _hex_prism_af(x, y, zmin, zmax, across_flats, rotation=0):
    """ Regular hexagonal prism aligned with +Z, specified by its
        width-across-flats (wrench size).  rotation is in degrees; at
        rotation 0 two flats are horizontal (normals along +/-Y).

        Built as the intersection of three |n.p| <= apothem slabs, so the
        field is baked (no coordinate map).
    """
    a = across_flats / 2.0            # apothem == half of across-flats
    terms = []
    for th in (30.0, 90.0, 150.0):
        ang = math.radians(th + rotation)
        nx, ny = math.cos(ang), math.sin(ang)
        lin = '+*f%g%s*f%g%s' % (nx, _off('X', x), ny, _off('Y', y))
        terms.append('-b%sf%g' % (lin, a))
    hex2d = 'aa%s%s%s' % (terms[0], terms[1], terms[2])
    corner = a / math.cos(math.radians(30.0))   # centre-to-vertex
    prof = Shape(hex2d, x - corner, y - corner, x + corner, y + corner)
    return extrude_z(prof, zmin, zmax)


# ---------------------------------------------------------------------------
# 3. Hex nut trap (+ side-entry slot)
# ---------------------------------------------------------------------------
#
# Hex-nut width across flats, mm - ISO 4032 style nut (== spanner size).
# HIGH confidence.
NUT_WIDTH_AF = {
    'M2': 4.0, 'M2.5': 5.0, 'M3': 5.5, 'M4': 7.0, 'M5': 8.0,
    'M6': 10.0, 'M8': 13.0, 'M10': 16.0, 'M12': 18.0,
}
# Standard hex-nut thickness m (max), mm - ISO 4032.  Used for the default
# pocket depth.  MEDIUM-HIGH confidence.
NUT_THICKNESS = {
    'M2': 1.6, 'M2.5': 2.0, 'M3': 2.4, 'M4': 3.2, 'M5': 4.7,
    'M6': 5.2, 'M8': 6.8, 'M10': 8.4, 'M12': 10.8,
}


def nut_trap_hex(x, y, z0, depth=None, size='M3', clearance=0.2,
                 slot_length=0.0, slot_dir='X'):
    """ Hexagonal pocket sized for a captive nut, optionally with a
        rectangular side-entry slide channel so the nut can be pushed in
        from the edge.

        NEGATIVE shape.  Pocket runs from z0 upward by `depth` (default =
        ISO 4032 nut thickness for `size`).  `clearance` is added to the
        across-flats (printer fit, default 0.2 mm).  slot_length > 0 adds a
        channel of the same width running in +slot_dir ('X' or 'Y').
        Default configuration is M3.
    """
    af = NUT_WIDTH_AF[size] + clearance
    d = depth if depth is not None else NUT_THICKNESS[size]
    pocket = _hex_prism_af(x, y, z0, z0 + d, af)
    if slot_length and slot_length > 0:
        if slot_dir == 'X':
            chan = rectangle(x, x + slot_length, y - af / 2.0, y + af / 2.0)
        elif slot_dir == 'Y':
            chan = rectangle(x - af / 2.0, x + af / 2.0, y, y + slot_length)
        else:
            raise ValueError("slot_dir must be 'X' or 'Y'")
        pocket = pocket | extrude_z(chan, z0, z0 + d)
    return pocket


# ---------------------------------------------------------------------------
# 4. Heat-set insert boss  &  PCB standoff  (POSITIVE solids)
# ---------------------------------------------------------------------------
#
# Heat-set threaded-insert bore diameters, mm (brass tapered inserts, e.g.
# CNC Kitchen / McMaster style).  These are the recommended boss BORE (hole)
# diameters.  MEDIUM confidence - varies a little by brand/knurl.
HEATSET_BORE_D = {
    'M2': 3.2, 'M2.5': 3.6, 'M3': 4.0, 'M4': 5.6, 'M5': 6.4,
}
# Typical standard-length insert lengths, mm (bore depth = length + 1).
# MEDIUM/LOW confidence - brand dependent; override `bore_depth` for exact.
HEATSET_LEN = {
    'M2': 4.0, 'M2.5': 4.0, 'M3': 5.7, 'M4': 5.8, 'M5': 6.4,
}


def heatset_boss(x, y, z0, height, size='M3', wall=1.6,
                 bore_depth=None, fillet=0.0):
    """ Cylindrical boss with a blind bore sized for a heat-set insert.

        POSITIVE solid - wire into a Union with your body.

        Boss radius = insert-bore radius + `wall`.  Boss runs from z0 up by
        `height`; the bore is drilled from the top down by
        (insert length + 1) unless bore_depth is given.  If fillet > 0 a
        rounded foot of that radius is blended at the base (fillet_union-ready
        for joining to a wall/floor).  Default configuration is M3.
    """
    bore_r = HEATSET_BORE_D[size] / 2.0
    boss_r = bore_r + wall
    depth = bore_depth if bore_depth is not None else HEATSET_LEN[size] + 1.0
    z_top = z0 + height
    boss = cylinder(x, y, z0, z_top, boss_r)
    if fillet and fillet > 0:
        foot = cylinder(x, y, z0, z0 + fillet, boss_r + fillet)
        boss = fillet_union(boss, foot, fillet)
    bore = cylinder(x, y, z_top - depth, z_top + _EPS, bore_r)
    return boss & ~bore


# Pilot-hole diameters for thread-forming screws driven straight into plastic,
# mm (~0.8 x nominal - a common self-tapping-into-3D-print rule of thumb).
# LOW-MEDIUM confidence; tune to your material.
PILOT_D = {
    'M2': 1.6, 'M2.5': 2.0, 'M3': 2.5, 'M4': 3.3, 'M5': 4.2,
    'M6': 5.0, 'M8': 6.8, 'M10': 8.5, 'M12': 10.2,
}


def standoff_pcb(x, y, z0, height, size='M3', mode='clearance',
                 wall=1.5, hole_d=None, fillet=0.0):
    """ PCB standoff: a boss with a hole running all the way through.

        POSITIVE solid - wire into a Union.

        mode='clearance' -> ISO 273 close clearance bore (bolt passes through).
        mode='pilot'     -> smaller self-tapping pilot bore.
        Override with hole_d.  fillet > 0 adds a rounded foot.
        Default configuration is M3, clearance.
    """
    if hole_d is not None:
        d = hole_d
    elif mode == 'pilot':
        d = PILOT_D[size]
    elif mode == 'clearance':
        d = CLEARANCE_CLOSE_D[size]
    else:
        raise ValueError("mode must be 'clearance' or 'pilot'")
    boss_r = d / 2.0 + wall
    z_top = z0 + height
    boss = cylinder(x, y, z0, z_top, boss_r)
    if fillet and fillet > 0:
        foot = cylinder(x, y, z0, z0 + fillet, boss_r + fillet)
        boss = fillet_union(boss, foot, fillet)
    hole = cylinder(x, y, z0 - _EPS, z_top + _EPS, d / 2.0)
    return boss & ~hole


# ---------------------------------------------------------------------------
# 5. Snap-fit cantilever clip  &  lid tongue/groove lip
# ---------------------------------------------------------------------------

def snap_clip(x, y, z, arm_length=10.0, arm_thickness=2.0, arm_width=4.0,
              barb_depth=1.2, barb_angle=45.0):
    """ Cantilever snap-fit clip: a rectangular arm with a hooked barb at the
        free end.  POSITIVE solid.

        The arm root is at (x, y, z); the arm runs along +X (length), its
        thickness is along +Y, its width along +Z.  The barb protrudes +Y at
        the tip with a lead-in ramp of `barb_angle` (from the arm face) and a
        square retention face.
    """
    arm = rectangle(x, x + arm_length, y, y + arm_thickness)
    tip = x + arm_length
    ramp = barb_depth / math.tan(math.radians(barb_angle))
    barb = triangle(tip, y + arm_thickness,
                    tip, y + arm_thickness + barb_depth,
                    tip - ramp, y + arm_thickness)
    return extrude_z(arm | barb, z, z + arm_width)


def lid_lip(x0, x1, y0, y1, width=1.5, height=2.0, kind='tongue',
            clearance=0.2, z0=0.0):
    """ Perimeter tongue/groove ring for box+lid registration around the
        rectangular footprint (x0,x1,y0,y1) - the ring is centred on that
        border line.

        kind='tongue' -> POSITIVE raised ring (Union onto the lid).
        kind='groove' -> NEGATIVE ring, widened by `clearance` (Difference
                         from the box rim).
        Ring wall thickness = `width`; it rises `height` from z0.
    """
    w = width if kind == 'tongue' else (width + clearance)
    h = w / 2.0
    outer = rectangle(x0 - h, x1 + h, y0 - h, y1 + h)
    inner = rectangle(x0 + h, x1 - h, y0 + h, y1 - h)
    ring = outer & ~inner
    if kind == 'tongue':
        return extrude_z(ring, z0, z0 + height)
    elif kind == 'groove':
        return extrude_z(ring, z0 - _EPS, z0 + height)
    else:
        raise ValueError("kind must be 'tongue' or 'groove'")


# ---------------------------------------------------------------------------
# 6. Sliding dovetail (male pin / female socket)
# ---------------------------------------------------------------------------

def dovetail(x, y, z, width=10.0, height=6.0, depth=20.0, angle=15.0,
             kind='male', clearance=0.2):
    """ Sliding dovetail prism: a trapezoid cross-section (wider at the top,
        i.e. undercut) extruded along +Z (the slide direction).

        kind='male'   -> POSITIVE pin.
        kind='female' -> NEGATIVE socket, enlarged by `clearance` on the
                         cross-section (Difference from a block).

        (x, y, z) is the centre of the bottom face.  `width` is the wide
        (top) face width, `height` the pin height, `angle` the flank angle
        from vertical (degrees).
    """
    w = width if kind == 'male' else (width + clearance)
    h = height if kind == 'male' else (height + clearance)
    nb = w / 2.0 - h * math.tan(math.radians(angle))   # narrow (bottom) half
    base = rectangle(x - w / 2.0, x + w / 2.0, y, y + h)
    left = triangle(x - w / 2.0, y, x - nb, y, x - w / 2.0, y + h)
    right = triangle(x + nb, y, x + w / 2.0, y, x + w / 2.0, y + h)
    prof = base & ~(left | right)
    if kind == 'male':
        return extrude_z(prof, z, z + depth)
    elif kind == 'female':
        return extrude_z(prof, z - _EPS, z + depth + _EPS)
    else:
        raise ValueError("kind must be 'male' or 'female'")


# ---------------------------------------------------------------------------
# 7. Hex-pattern vent grille  (trig quasi-lattice)
# ---------------------------------------------------------------------------

def vent_grille_hex(x0, x1, y0, y1, z0, z1, cell=5.0, strut=1.0):
    """ Ventilation panel perforated with a hexagonal quasi-lattice, built
        from three 60-degree-rotated cosines thresholded (no `mod` opcode
        needed - trig is inherently periodic; cf. the gyroid node).

        POSITIVE solid (panel with holes).  The strut network is the region
        where T(x,y) < threshold, with
            T = cos(k.d0.p) + cos(k.d1.p) + cos(k.d2.p),  k = 2*pi/cell
        for the three lattice directions d0,d1,d2 at 0/60/120 degrees.  The
        struts are then extruded through z0..z1 and clipped to the panel rect.

        `cell` is the lattice period; `strut` is the strut width, mapped
        heuristically to the threshold (APPROXIMATE - the field is not a true
        distance field, so realised strut width is nominal, not exact).
    """
    k = 2.0 * math.pi / cell
    c1 = 'c*f%gX' % k
    c2 = 'c*f%g+*f0.5X*f0.8660254Y' % k
    c3 = 'c*f%g+*nf0.5X*f0.8660254Y' % k
    T = '+%s+%s%s' % (c1, c2, c3)                       # T in [-1.5, 3]
    frac = min(max(strut / cell, 0.0), 1.0)
    thr = -1.5 + 4.5 * frac                             # struts where T < thr
    struts2d = Shape('-%sf%g' % (T, thr), x0, y0, x1, y1)
    panel = cube(x0, x1, y0, y1, z0, z1)
    return extrude_z(struts2d, z0, z1) & panel

################################################################################
# ISO threads, involute gears, racks (phase 1 node campaign)
################################################################################


# --------------------------------------------------------------------------
# ISO metric coarse-pitch table (nominal diameter mm -> pitch mm)
# --------------------------------------------------------------------------
ISO_COARSE_PITCH = {
    2: 0.40, 2.5: 0.45, 3: 0.50, 4: 0.70, 5: 0.80,
    6: 1.00, 8: 1.25, 10: 1.50, 12: 1.75,
}


def iso_pitch(nominal_d):
    """Coarse pitch for a nominal metric diameter, or a sane fallback."""
    if nominal_d in ISO_COARSE_PITCH:
        return ISO_COARSE_PITCH[nominal_d]
    # nearest tabulated key
    key = min(ISO_COARSE_PITCH, key=lambda k: abs(k - nominal_d))
    return ISO_COARSE_PITCH[key]


# --------------------------------------------------------------------------
# THREAD FIELD
# --------------------------------------------------------------------------
def _iso_H(pitch):
    """Height of the ISO fundamental (sharp) triangle: H = (sqrt3/2) * P."""
    return math.sqrt(3.0) / 2.0 * pitch


def _thread_radii(d, pitch, clearance):
    """Return (r_crest, r_root, r_mid, amp) for the clamped-triangle profile."""
    H = _iso_H(pitch)
    rmaj = d / 2.0
    r_crest = rmaj + clearance            # crest flat (external major radius)
    r_root = rmaj - 5.0 * H / 8.0 + clearance   # root flat (5H/8 depth)
    r_mid = rmaj - 3.0 * H / 8.0 + clearance    # midline of the sharp V
    amp = H / 2.0                         # sharp-V half amplitude
    return r_crest, r_root, r_mid, amp


def _thread_field_infix(d, pitch, n_starts, clearance):
    """The infix ``=...;`` field string for the helical thread surface.

    Field is < 0 inside the solid, uses atan2 for theta and the seam-free
    asin(sin(...)) triangle wave for the profile.
    """
    r_crest, r_root, r_mid, amp = _thread_radii(d, pitch, clearance)
    K = 2.0 * math.pi / pitch             # axial angular frequency
    two_over_pi = 2.0 / math.pi

    def c(v):
        # parenthesised full-precision constant (handles negatives cleanly)
        return '(%r)' % v

    tri = 'asin(sin(%s*atan2(Y,X) - %s*Z))*%s' % (c(n_starts), c(K),
                                                   c(two_over_pi))
    sharp = '%s + %s*(%s)' % (c(r_mid), c(amp), tri)
    R = 'min(max(%s, %s), %s)' % (sharp, c(r_root), c(r_crest))
    field = '=sqrt(X*X + Y*Y) - (%s);' % R
    return field, r_crest, r_root


def _cap(zmin, zmax):
    """Prefix expression max(zmin - Z, Z - zmax) -- the axial end caps."""
    return 'a-f%rZ-Zf%r' % (zmin, zmax)


def iso_thread_external(d, pitch=None, length=10.0, n_starts=1,
                        clearance=0.0, zmin=0.0):
    """ISO metric *external* thread (a threaded rod).

    d          nominal (major) diameter [mm]
    pitch      thread pitch [mm]; None -> ISO coarse lookup
    length     threaded length [mm]
    n_starts   number of thread starts (integer; keeps the seam continuous)
    clearance  radial offset [mm]; 0 for the rod itself
    """
    if pitch is None:
        pitch = iso_pitch(d)
    zmax = zmin + length
    field, r_crest, r_root = _thread_field_infix(d, pitch, n_starts, clearance)
    math_str = 'a' + field + _cap(zmin, zmax)   # max(field, caps)
    return Shape(math_str,
                 -r_crest, -r_crest, zmin,
                  r_crest,  r_crest, zmax)


def thread_rod_iso(nominal_d, length=10.0, n_starts=1, clearance=0.0):
    """Convenience: ISO coarse threaded rod by nominal diameter (M-size)."""
    return iso_thread_external(nominal_d, iso_pitch(nominal_d), length,
                               n_starts, clearance)


def iso_thread_internal_negative(d, pitch=None, length=10.0, n_starts=1,
                                 clearance=0.2, zmin=0.0):
    """ISO metric *internal* thread as a NEGATIVE solid.

    Subtract this from a body to cut a tapped hole that mates with the
    matching external rod.  It is the external-thread field grown radially
    by ``clearance`` (default 0.2 mm on radius) so the printed nut is a
    touch larger than the bolt everywhere.

    Geometry: a core shaft at (minor radius + clearance) with helical
    thread ridges out to (major radius + clearance).  Same handedness and
    phase as :func:`iso_thread_external`, so the two thread together.
    """
    if pitch is None:
        pitch = iso_pitch(d)
    return iso_thread_external(d, pitch, length, n_starts, clearance, zmin)


def tapped_hole_negative(nominal_d, length=10.0, n_starts=1, clearance=0.2):
    """Convenience negative for an ISO coarse tapped hole (subtract me)."""
    return iso_thread_internal_negative(nominal_d, iso_pitch(nominal_d),
                                        length, n_starts, clearance)


# --------------------------------------------------------------------------
# INVOLUTE SPUR GEAR
# --------------------------------------------------------------------------
def _inv(a):
    """Involute function inv(a) = tan(a) - a."""
    return math.tan(a) - a


def gear_geometry(module, teeth, pressure_angle=20.0):
    """Return the standard radii for an involute spur gear."""
    m = float(module)
    z = int(teeth)
    a = math.radians(pressure_angle)
    r_pitch = m * z / 2.0
    r_base = r_pitch * math.cos(a)
    r_add = r_pitch + m           # addendum = m
    r_ded = r_pitch - 1.25 * m    # dedendum = 1.25 m
    return dict(m=m, z=z, alpha=a, r_pitch=r_pitch, r_base=r_base,
                r_add=r_add, r_ded=r_ded, r_root=max(r_ded, 0.05 * m))


def _gear_flank_points(g, samples=12):
    """Sampled involute points for the +theta flank, root -> tip.

    Tooth centred on the +X axis; the -theta flank is the mirror.
    """
    r_base, r_add = g['r_base'], g['r_add']
    r_root, alpha = g['r_root'], g['alpha']
    z = g['z']
    ht = math.pi / (2.0 * z)          # half tooth angle at the pitch circle
    inv_p = _inv(alpha)

    def flank_angle(r):
        # clamp for r just under r_base (numerical) -> acos arg <= 1
        ar = math.acos(min(1.0, r_base / r))
        return ht + inv_p - _inv(ar)

    r_start = max(r_base, r_root)
    pts = []
    # radial extension below the base circle (documented undercut approx)
    if r_base > r_root:
        a_edge = ht + inv_p           # flank_angle(r_base)
        pts.append((r_root * math.cos(a_edge), r_root * math.sin(a_edge)))
    for i in range(samples):
        r = r_start + (r_add - r_start) * i / (samples - 1)
        th = flank_angle(r)
        pts.append((r * math.cos(th), r * math.sin(th)))
    return pts


def spur_gear_2d(module=2.0, teeth=20, pressure_angle=20.0, samples=12):
    """2D involute spur-gear cross section (in the XY plane)."""
    g = gear_geometry(module, teeth, pressure_angle)
    right = _gear_flank_points(g, samples)        # +theta, root->tip
    left = [(x, -y) for (x, y) in right]          # -theta, root->tip

    # boundary vertex loop, angle strictly increasing: -a_edge .. +a_edge
    verts = left + right[::-1]

    def collinear(p0, p1, p2):
        return abs((p1[0] - p0[0]) * (p2[1] - p0[1]) -
                   (p1[1] - p0[1]) * (p2[0] - p0[0])) < 1e-12

    tooth = None
    O = (0.0, 0.0)
    for i in range(len(verts) - 1):
        p1, p2 = verts[i], verts[i + 1]
        if collinear(O, p1, p2):
            continue                              # zero-area radial sliver
        tri = triangle(0.0, 0.0, p1[0], p1[1], p2[0], p2[1])
        tooth = tri if tooth is None else (tooth | tri)

    teeth_all = iterate_polar(tooth, 0.0, 0.0, g['z'])
    gear = teeth_all | circle(0.0, 0.0, g['r_root'])
    gear = gear & circle(0.0, 0.0, g['r_add'])

    # tighten bounds
    ra = g['r_add']
    return Shape(gear.math, -ra, -ra, ra, ra)


def spur_gear(module=2.0, teeth=20, pressure_angle=20.0, thickness=6.0,
              samples=12, zmin=0.0):
    """3D involute spur gear (2D section extruded in Z)."""
    g2 = spur_gear_2d(module, teeth, pressure_angle, samples)
    return extrude_z(g2, zmin, zmin + thickness)


# --------------------------------------------------------------------------
# RACK (involute limit case: straight trapezoidal teeth)
# --------------------------------------------------------------------------
def rack_2d(module=2.0, n_teeth=8, pressure_angle=20.0, base_height=None):
    """2D involute rack cross section.  Pitch line at Y=0, teeth point +Y."""
    m = float(module)
    n = int(n_teeth)
    a = math.radians(pressure_angle)
    p = math.pi * m                      # circular pitch
    add = m
    ded = 1.25 * m
    if base_height is None:
        base_height = 1.5 * m

    half_top = p / 4.0 - add * math.tan(a)   # tooth half-width at addendum
    half_bot = p / 4.0 + ded * math.tan(a)   # tooth half-width at dedendum

    # one tooth centred at X=0 (two triangles form the trapezoid)
    bl = (-half_bot, -ded)
    br = (half_bot, -ded)
    tr = (half_top, add)
    tl = (-half_top, add)
    tooth = (triangle(bl[0], bl[1], br[0], br[1], tr[0], tr[1]) |
             triangle(bl[0], bl[1], tr[0], tr[1], tl[0], tl[1]))

    total = n * p
    x0 = -total / 2.0 + p / 2.0          # centre the row
    teeth = functools.reduce(operator.or_,
                             [move(tooth, x0 + i * p, 0.0)
                              for i in range(n)])

    body = rectangle(-total / 2.0, total / 2.0,
                            -ded - base_height, 0.0)
    rack = teeth | body
    return Shape(rack.math,
                 -total / 2.0, -ded - base_height,
                  total / 2.0, add)


def rack(module=2.0, n_teeth=8, pressure_angle=20.0, thickness=6.0,
         base_height=None, zmin=0.0):
    """3D involute rack (2D section extruded in Z)."""
    r2 = rack_2d(module, n_teeth, pressure_angle, base_height)
    return extrude_z(r2, zmin, zmin + thickness)


################################################################################
# Domain repetition (O(1) field cost regardless of copy count)
################################################################################

def _rep_expr(coord, spacing, center):
    """ mod-space fold: coord' = mod(coord - off, spacing) + off,
        where off centers the cell on `center`. """
    off = center - spacing / 2.
    return '+M-%sf%g' % (coord, off) + 'f%g' % spacing + 'f%g' % off

def _rep_mirror_expr(coord, spacing, center):
    """ mirrored fold: neighboring cells are reflections.
        coord' = center + |mod(coord - center + s/2, 2s) - s| - s/2 """
    s = spacing
    t = 'M+-%sf%g' % (coord, center) + 'f%g' % (s / 2.) + 'f%g' % (2 * s)
    return '+f%g' % center + '-b-%s' % t + 'f%g' % s + 'f%g' % (s / 2.)

def _rep_finite_expr(coord, spacing, count, start):
    """ finite fold: coord' = coord - s * clamp(floor((coord-start)/s),
        0, count-1); cells [start, start + count*s) collapse onto the
        first cell. """
    k = 'F/-%sf%g' % (coord, start) + 'f%g' % spacing
    k = 'ia%sf0f%g' % (k, count - 1)
    return '-%s*f%g%s' % (coord, spacing, k)

def _remap_shape(a, xexpr='_', yexpr='_', zexpr='_',
                 xmin=None, xmax=None, ymin=None, ymax=None,
                 zmin=None, zmax=None):
    b = a.bounds
    pick = lambda v, d: d if v is None else v
    return Shape('m%s%s%s%s' % (xexpr, yexpr, zexpr, a.math),
                 pick(xmin, b.xmin), pick(ymin, b.ymin),
                 pick(zmin, b.zmin), pick(xmax, b.xmax),
                 pick(ymax, b.ymax), pick(zmax, b.zmax))

def repeat_x(a, spacing, x0=0):
    """ Repeats a shape infinitely along X. The shape should live within
        one cell of the given spacing, centered on x0.
    """
    if spacing <= 0:
        return a
    return _remap_shape(a, xexpr=_rep_expr('X', spacing, x0),
                        xmin=float('-inf'), xmax=float('inf'))

def repeat_y(a, spacing, y0=0):
    if spacing <= 0:
        return a
    return _remap_shape(a, yexpr=_rep_expr('Y', spacing, y0),
                        ymin=float('-inf'), ymax=float('inf'))

def repeat_z(a, spacing, z0=0):
    if spacing <= 0:
        return a
    return _remap_shape(a, zexpr=_rep_expr('Z', spacing, z0),
                        zmin=float('-inf'), zmax=float('inf'))

def repeat_xy(a, spacing_x, spacing_y, x0=0, y0=0):
    """ Infinite 2D grid repetition. """
    if spacing_x <= 0 or spacing_y <= 0:
        return a
    return _remap_shape(a, xexpr=_rep_expr('X', spacing_x, x0),
                        yexpr=_rep_expr('Y', spacing_y, y0),
                        xmin=float('-inf'), xmax=float('inf'),
                        ymin=float('-inf'), ymax=float('inf'))

def repeat_mirror_x(a, spacing, x0=0):
    """ Infinite repetition along X where neighboring cells are mirror
        images (seamless for asymmetric shapes). """
    if spacing <= 0:
        return a
    return _remap_shape(a, xexpr=_rep_mirror_expr('X', spacing, x0),
                        xmin=float('-inf'), xmax=float('inf'))

def repeat_mirror_y(a, spacing, y0=0):
    if spacing <= 0:
        return a
    return _remap_shape(a, yexpr=_rep_mirror_expr('Y', spacing, y0),
                        ymin=float('-inf'), ymax=float('inf'))

def repeat_x_finite(a, spacing, count, x0=0):
    """ Repeats a shape `count` times along X at the given spacing,
        starting from the cell centered on x0. O(1) field cost (unlike
        the Array nodes, which union N copies).
    """
    count = int(count)
    if spacing <= 0 or count < 2:
        return a
    b = a.bounds
    return _remap_shape(a,
            xexpr=_rep_finite_expr('X', spacing, count, x0 - spacing/2.),
            xmin=b.xmin, xmax=b.xmax + spacing * (count - 1))

def repeat_y_finite(a, spacing, count, y0=0):
    count = int(count)
    if spacing <= 0 or count < 2:
        return a
    b = a.bounds
    return _remap_shape(a,
            yexpr=_rep_finite_expr('Y', spacing, count, y0 - spacing/2.),
            ymin=b.ymin, ymax=b.ymax + spacing * (count - 1))

def repeat_polar(a, n, x=0, y=0):
    """ Repeats a shape n times around the point (x, y). The source
        shape should live in the sector straddling the +X direction.
        O(1) field cost (unlike iterate_polar, which unions N copies).
    """
    n = int(n)
    if n < 2:
        return a
    import math as _m
    sector = 2 * _m.pi / n
    theta = ('mod(atan2(Y-%g,X-%g)+%g,%g)-%g'
             % (y, x, sector / 2., sector, sector / 2.))
    r = 'sqrt((X-%g)*(X-%g)+(Y-%g)*(Y-%g))' % (x, x, y, y)
    xexpr = '=%g+%s*cos(%s);' % (x, r, theta)
    yexpr = '=%g+%s*sin(%s);' % (y, r, theta)
    b = a.bounds
    rad = max(abs(b.xmin - x), abs(b.xmax - x),
              abs(b.ymin - y), abs(b.ymax - y))
    return _remap_shape(a, xexpr=xexpr, yexpr=yexpr,
                        xmin=x - rad, xmax=x + rad,
                        ymin=y - rad, ymax=y + rad)

def repeat_scale_xy(a, factor, x=0, y=0, r0=1.0):
    """ Infinite self-similar (scale) repetition about the point (x, y):
        the pattern is copied at every scale factor**k, k in Z. The
        source shape should live in the annulus r0 <= r < r0*factor.
        This is log-space domain repetition - recursion depth is
        unlimited at O(1) field cost. (Note: the remapped field is
        non-euclidean away from the source band; sign is exact.)
    """
    if factor <= 1 or r0 <= 0:
        return a
    import math as _m
    lnk = _m.log(factor)
    r = 'sqrt((X-%g)*(X-%g)+(Y-%g)*(Y-%g))' % (x, x, y, y)
    q = 'exp(mod(log(%s/%g),%g))*%g/%s' % (r, r0, lnk, r0, r)
    xexpr = '=%g+(X-%g)*%s;' % (x, x, q)
    yexpr = '=%g+(Y-%g)*%s;' % (y, y, q)
    rmax = r0 * factor * 1.0
    b = a.bounds
    return Shape('m%s%s_%s' % (xexpr, yexpr, a.math),
                 x - rmax, y - rmax, b.zmin, x + rmax, y + rmax, b.zmax)

def iterate_scaled(part, n, factor, dx=0, dy=0, x0=0, y0=0):
    """ Finite copy -> translate -> scale chain: n copies, each scaled
        by `factor` about (x0, y0) relative to the previous and offset
        by (dx, dy) scaled to its generation. Automates the manual
        recursive-feature loop. O(n) field cost - use repeat_scale_xy
        for unlimited depth about a fixed point.
    """
    n = int(n)
    out = part
    s = 1.0
    px, py = 0.0, 0.0
    for i in range(1, n):
        px += dx * s
        py += dy * s
        s *= factor
        copy = scale_xy(part, x0, y0, s)
        out |= move(copy, px, py, 0)
    return out

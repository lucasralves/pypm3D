import typing as tp
import math
import numpy as np
import matplotlib.pyplot as plt


class Vec3D(tp.NamedTuple):
    x: float
    y: float
    z: float

def division(a: float, b: float) -> float:
    if math.fabs(b) < 1e-12:
        return a * b / (b * b + 1e-12)
    else:
        return a / b


# Source
def quad_source_potential(x: float, y: float, z: float,
                          x1: float, y1: float,
                          x2: float, y2: float,
                          x3: float, y3: float,
                          x4: float, y4: float,
                          r1: float, r2: float, r3: float, r4: float,
                          e1: float, e2: float, e3: float, e4: float,
                          h1: float, h2: float, h3: float, h4: float,
                          d12: float, d23: float, d34: float, d41: float,
                          m12: float, m23: float, m34: float, m41: float) -> float:
    
    phi = (0.25 / math.pi) * ( (((x - x1) * (y2 - y1) - (y - y1) * (x2 - x1)) / d12) * math.log((r1 + r2 + d12) / (r1 + r2 - d12)) + \
                               (((x - x2) * (y3 - y2) - (y - y2) * (x3 - x2)) / d23) * math.log((r2 + r3 + d23) / (r2 + r3 - d23)) + \
                               (((x - x3) * (y4 - y3) - (y - y3) * (x4 - x3)) / d34) * math.log((r3 + r4 + d34) / (r3 + r4 - d34)) + \
                               (((x - x4) * (y1 - y4) - (y - y4) * (x1 - x4)) / d41) * math.log((r4 + r1 + d41) / (r4 + r1 - d41)) - \
                               math.fabs(z) * ( math.atan((m12 * e1 - h1) / z * r1) - math.atan((m12 * e2 - h2) / z * r2) + \
                                                math.atan((m23 * e2 - h2) / z * r2) - math.atan((m23 * e3 - h3) / z * r3) + \
                                                math.atan((m34 * e3 - h3) / z * r3) - math.atan((m34 * e4 - h4) / z * r4) + \
                                                math.atan((m41 * e4 - h4) / z * r4) - math.atan((m41 * e1 - h1) / z * r1) ) )

    return phi

def tri_source_potential(x: float, y: float, z: float,
                         x1: float, y1: float,
                         x2: float, y2: float,
                         x3: float, y3: float,
                         r1: float, r2: float, r3: float,
                         e1: float, e2: float, e3: float,
                         h1: float, h2: float, h3: float,
                         d12: float, d23: float, d31: float,
                         m12: float, m23: float, m31: float,) -> float:
    
    phi = (0.25 / math.pi) * ( (((x - x1) * (y2 - y1) - (y - y1) * (x2 - x1)) / d12) * math.log((r1 + r2 + d12) / (r1 + r2 - d12)) + \
                               (((x - x2) * (y3 - y2) - (y - y2) * (x3 - x2)) / d23) * math.log((r2 + r3 + d23) / (r2 + r3 - d23)) + \
                               (((x - x3) * (y1 - y3) - (y - y3) * (x1 - x3)) / d31) * math.log((r3 + r1 + d31) / (r3 + r1 - d31)) - \
                               math.fabs(z) * ( math.atan((m12 * e1 - h1) / z * r1) - math.atan((m12 * e2 - h2) / z * r2) + \
                                                math.atan((m23 * e2 - h2) / z * r2) - math.atan((m23 * e3 - h3) / z * r3) + \
                                                math.atan((m31 * e3 - h3) / z * r3) - math.atan((m31 * e1 - h1) / z * r1) ) )

    return phi


# Panel
def panel(n: int, p: Vec3D, p1: Vec3D, p2: Vec3D, p3:Vec3D, p4: Vec3D) -> tp.List:
    
    r1 = math.sqrt((p.x - p1.x) * (p.x - p1.x) + (p.y - p1.y) * (p.y - p1.y) + p.z * p.z)
    r2 = math.sqrt((p.x - p2.x) * (p.x - p2.x) + (p.y - p2.y) * (p.y - p2.y) + p.z * p.z)
    r3 = math.sqrt((p.x - p3.x) * (p.x - p3.x) + (p.y - p3.y) * (p.y - p3.y) + p.z * p.z)
    d12 = math.sqrt((p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y))
    d23 = math.sqrt((p3.x - p2.x) * (p3.x - p2.x) + (p3.y - p2.y) * (p3.y - p2.y))
    m12 = division((p2.y - p1.y), (p2.x - p1.x))
    m23 = division((p3.y - p2.y), (p3.x - p2.x))
    e1 = (p.x - p1.x) * (p.x - p1.x) + p.z * p.z
    e2 = (p.x - p2.x) * (p.x - p2.x) + p.z * p.z
    e3 = (p.x - p3.x) * (p.x - p3.x) + p.z * p.z
    h1 = (p.x - p1.x) * (p.y - p1.y)
    h2 = (p.x - p2.x) * (p.y - p2.y)
    h3 = (p.x - p3.x) * (p.y - p3.y)

    if n == 4:

        r4 = math.sqrt((p.x - p4.x) * (p.x - p4.x) + (p.y - p4.y) * (p.y - p4.y) + p.z * p.z)
        d34 = math.sqrt((p4.x - p3.x) * (p4.x - p3.x) + (p4.y - p3.y) * (p4.y - p3.y))
        d41 = math.sqrt((p1.x - p4.x) * (p1.x - p4.x) + (p1.y - p4.y) * (p1.y - p4.y))
        m34 = division((p4.y - p3.y), (p4.x - p3.x))
        m41 = division((p1.y - p4.y), (p1.x - p4.x))
        e4 = (p.x - p4.x) * (p.x - p4.x) + p.z * p.z
        h4 = (p.x - p4.x) * (p.y - p4.y)

        source_pot = quad_source_potential(p.x, p.y, p.z, p1.x, p1.y, p2.x, p2.y, p3.x, p3.y, p4.x, p4.y, r1, r2, r3, r4, e1, e2, e3, e4, h1, h2, h3, h4, d12, d23, d34, d41, m12, m23, m34, m41)

    else:
        
        d31 = math.sqrt((p1.x - p3.x) * (p1.x - p3.x) + (p1.y - p3.y) * (p1.y - p3.y))
        m31 = division((p1.y - p3.y), (p1.x - p3.x))

        source_pot = tri_source_potential(p.x, p.y, p.z, p1.x, p1.y, p2.x, p2.y, p3.x, p3.y, r1, r2, r3, e1, e2, e3, h1, h2, h3, d12, d23, d31, m12, m23, m31)
       

    return source_pot


def plot_xy_yz(n_sides: int, p1: Vec3D, p2: Vec3D, p3: Vec3D, p4: Vec3D, z0: float, x0: float):

    # Grid
    n = 100
    l = 5
    x = np.linspace(-l, l, num=n)
    y = np.linspace(-l, l, num=n)
    z = np.linspace(-l, l, num=n)

    XY, YX = np.meshgrid(x, y)
    YZ, ZY = np.meshgrid(y, z)

    # Parameters
    source_pot_xy = np.empty((n, n), dtype=np.double)
    source_pot_yz = np.empty((n, n), dtype=np.double)

    # Fill parameters
    for i in range(n):
        for j in range(n):

            # xy
            p = Vec3D(XY[i, j], YX[i, j], z0)
            pot_xy = panel(n_sides, p, p1, p2, p3, p4)

            source_pot_xy[i, j] = pot_xy

            # yz
            p = Vec3D(x0, YZ[i, j], ZY[i, j])
            pot_yz = panel(n_sides, p, p1, p2, p3, p4)

            source_pot_yz[i, j] = pot_yz
    
    # xy
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.suptitle('xy')
    
    c1 = ax1.contourf(XY, YX, source_pot_xy)
    ax1.axis('equal')
    fig.colorbar(c1)

    c2 = ax2.contourf(YZ, ZY, source_pot_yz)
    ax2.axis('equal')
    fig.colorbar(c2)

    plt.show()

    # plt.figure()
    # plt.title('Potential - source')
    # plt.contourf(X, Y, source_phi)

    # if n_sides == 4:
    #     plt.plot([p1.x, p2.x, p3.x, p4.x, p1.x], [p1.y, p2.y, p3.y, p4.y, p4.x], 'k')
    # else:
    #     plt.plot([p1.x, p2.x, p3.x, p1.x], [p1.y, p2.y, p3.y, p1.y], 'k')
    
    # plt.colorbar()
    # plt.axis('equal')

    # plt.figure()
    # plt.title('Velocity - source')
    # plt.contourf(X, Y, source_vel[:, :, 2])
    # plt.colorbar()
    # plt.quiver(X, Y, source_vel[:, :, 0], source_vel[:, :, 1])

    # if n_sides == 4:
    #     plt.plot([p1.x, p2.x, p3.x, p4.x, p1.x], [p1.y, p2.y, p3.y, p4.y, p4.x], 'k')
    # else:
    #     plt.plot([p1.x, p2.x, p3.x, p1.x], [p1.y, p2.y, p3.y, p1.y], 'k')
    
    # plt.axis('equal')

    # plt.show()

    return

if __name__ == '__main__':

    p1 = Vec3D(+1, +1, 0)
    p2 = Vec3D(-1, +1, 0)
    p3 = Vec3D(-1, -1, 0)
    p4 = Vec3D(+1, -1, 0)

    plot_xy_yz(4, p1, p2, p3, p4, -0.01, 0.0)
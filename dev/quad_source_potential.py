import math

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

def quad_source_velocity(x: float, y: float, z: float,
                         x1: float, y1: float,
                         x2: float, y2: float,
                         x3: float, y3: float,
                         x4: float, y4: float,
                         r1: float, r2: float, r3: float, r4: float,
                         e1: float, e2: float, e3: float, e4: float,
                         h1: float, h2: float, h3: float, h4: float,
                         d12: float, d23: float, d34: float, d41: float,
                         m12: float, m23: float, m34: float, m41: float) -> float:
    
    u = (0.25 / math.pi) * ( ((y2 - y1) / d12) * math.log((r1 + r2 + d12) / (r1 + r2 - d12)) + \
                             ((y3 - y2) / d23) * math.log((r2 + r3 + d23) / (r2 + r3 - d23)) + \
                             ((y4 - y3) / d34) * math.log((r3 + r4 + d34) / (r3 + r4 - d34)) + \
                             ((y1 - y4) / d41) * math.log((r4 + r1 + d41) / (r4 + r1 - d41)) )
    
    v = (0.25 / math.pi) * ( ((x1 - x2) / d12) * math.log((r1 + r2 - d12) / (r1 + r2 + d12)) + \
                             ((x2 - x3) / d23) * math.log((r2 + r3 - d23) / (r2 + r3 + d23)) + \
                             ((x3 - x4) / d34) * math.log((r3 + r4 - d34) / (r3 + r4 + d34)) + \
                             ((x4 - x1) / d41) * math.log((r4 + r1 - d41) / (r4 + r1 + d41)) )
    
    w = (0.25 / math.pi) * ( math.atan( (m12 * e1 - h1) / (z * r1) ) - math.atan( (m12 * e2 - h2) / (z * r2) ) + \
                             math.atan( (m23 * e2 - h2) / (z * r2) ) - math.atan( (m23 * e3 - h3) / (z * r3) ) + \
                             math.atan( (m34 * e3 - h3) / (z * r3) ) - math.atan( (m34 * e4 - h4) / (z * r4) ) + \
                             math.atan( (m41 * e4 - h4) / (z * r4) ) - math.atan( (m41 * e1 - h1) / (z * r1) ) )

    vel = Vec3D(u, v, w)

    return vel

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

def tri_source_velocity(x: float, y: float, z: float,
                        x1: float, y1: float,
                        x2: float, y2: float,
                        x3: float, y3: float,
                        r1: float, r2: float, r3: float,
                        e1: float, e2: float, e3: float,
                        h1: float, h2: float, h3: float,
                        d12: float, d23: float, d31: float,
                        m12: float, m23: float, m31: float,) -> float:
    
    u = (0.25 / math.pi) * ( ((y2 - y1) / d12) * math.log((r1 + r2 + d12) / (r1 + r2 - d12)) + \
                             ((y3 - y2) / d23) * math.log((r2 + r3 + d23) / (r2 + r3 - d23)) + \
                             ((y1 - y3) / d31) * math.log((r3 + r1 + d31) / (r3 + r1 - d31)) )
    
    v = (0.25 / math.pi) * ( ((x1 - x2) / d12) * math.log((r1 + r2 - d12) / (r1 + r2 + d12)) + \
                             ((x2 - x3) / d23) * math.log((r2 + r3 - d23) / (r2 + r3 + d23)) + \
                             ((x3 - x1) / d31) * math.log((r3 + r1 - d31) / (r3 + r1 + d31)) )
    
    w = (0.25 / math.pi) * ( math.atan( (m12 * e1 - h1) / (z * r1) ) - math.atan( (m12 * e2 - h2) / (z * r2) ) + \
                             math.atan( (m23 * e2 - h2) / (z * r2) ) - math.atan( (m23 * e3 - h3) / (z * r3) ) + \
                             math.atan( (m31 * e3 - h3) / (z * r3) ) - math.atan( (m31 * e1 - h1) / (z * r1) ) )

    vel = Vec3D(u, v, w)

    return vel
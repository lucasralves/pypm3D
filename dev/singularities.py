from typing import NamedTuple
from math import sqrt, fabs, pi, atan2, atan, log


def divide(a: float, b: float) -> float:
    if fabs(b) < 1e-12:
        return a * b / (b * b + 1e-12)
    else:
        return a / b

class Vec3D(NamedTuple):
    x: float
    y: float
    z: float

def quad_panel(p: Vec3D, p1: Vec3D, p2: Vec3D, p3: Vec3D, p4: Vec3D):

    d12 = sqrt( (p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y) )
    d23 = sqrt( (p3.x - p2.x) * (p3.x - p2.x) + (p3.y - p2.y) * (p3.y - p2.y) )
    d34 = sqrt( (p4.x - p3.x) * (p4.x - p3.x) + (p4.y - p3.y) * (p4.y - p3.y) )
    d41 = sqrt( (p1.x - p4.x) * (p1.x - p4.x) + (p1.y - p4.y) * (p1.y - p4.y) )

    C12 = (p2.x - p1.x) / d12
    C23 = (p3.x - p2.x) / d23
    C34 = (p4.x - p3.x) / d34
    C41 = (p1.x - p4.x) / d41

    S12 = (p2.y - p1.y) / d12
    S23 = (p3.y - p2.y) / d23
    S34 = (p4.y - p3.y) / d34
    S41 = (p1.y - p4.y) / d41

    s12_1 = (p1.x - p.x) * C12 + (p1.y - p.y) * S12
    s12_2 = (p2.x - p.x) * C12 + (p2.y - p.y) * S12

    s23_2 = (p2.x - p.x) * C23 + (p2.y - p.y) * S23
    s23_3 = (p3.x - p.x) * C23 + (p3.y - p.y) * S23

    s34_3 = (p3.x - p.x) * C34 + (p3.y - p.y) * S34
    s34_4 = (p4.x - p.x) * C34 + (p4.y - p.y) * S34

    s41_4 = (p4.x - p.x) * C41 + (p4.y - p.y) * S41
    s41_1 = (p1.x - p.x) * C41 + (p1.y - p.y) * S41

    R12 = (p.x - p1.x) * S12 - (p.y - p1.y) * C12
    R23 = (p.x - p2.x) * S23 - (p.y - p2.y) * C23
    R34 = (p.x - p3.x) * S34 - (p.y - p3.y) * C34
    R41 = (p.x - p4.x) * S41 - (p.y - p4.y) * C41

    r1 = sqrt((p.x - p1.x) * (p.x - p1.x) + (p.y - p1.y) * (p.y - p1.y) + p.z * p.z)
    r2 = sqrt((p.x - p2.x) * (p.x - p2.x) + (p.y - p2.y) * (p.y - p2.y) + p.z * p.z)
    r3 = sqrt((p.x - p3.x) * (p.x - p3.x) + (p.y - p3.y) * (p.y - p3.y) + p.z * p.z)
    r4 = sqrt((p.x - p4.x) * (p.x - p4.x) + (p.y - p4.y) * (p.y - p4.y) + p.z * p.z)

    Q12 = log((r1 + r2 + d12) / (r1 + r2 - d12))
    Q23 = log((r2 + r3 + d23) / (r2 + r3 - d23))
    Q34 = log((r3 + r4 + d34) / (r3 + r4 - d34))
    Q41 = log((r4 + r1 + d41) / (r4 + r1 - d41))

    F12 = R12 * fabs(p.z) * (r1 * s12_2 - r2 * s12_1)
    F23 = R23 * fabs(p.z) * (r2 * s23_3 - r3 * s23_2)
    F34 = R34 * fabs(p.z) * (r3 * s34_4 - r4 * s34_3)
    F41 = R41 * fabs(p.z) * (r4 * s41_1 - r1 * s41_4)

    G12 = r1 * r2 * R12 * R12 + p.z * p.z * s12_1 * s12_2
    G23 = r2 * r3 * R23 * R23 + p.z * p.z * s23_2 * s23_3
    G34 = r3 * r4 * R34 * R34 + p.z * p.z * s34_3 * s34_4
    G41 = r4 * r1 * R41 * R41 + p.z * p.z * s41_4 * s41_1

    J12 = atan2(F12, G12)
    J23 = atan2(F23, G23)
    J34 = atan2(F34, G34)
    J41 = atan2(F41, G41)

    R12_x = S12; R12_y = - C12
    R23_x = S23; R23_y = - C23
    R34_x = S34; R34_y = - C34
    R41_x = S41; R41_y = - C41

    r1_x = (p.x - p1.x) / r1; r1_y = (p.y - p1.y) / r1; r1_z = p.z / r1
    r2_x = (p.x - p2.x) / r2; r2_y = (p.y - p2.y) / r2; r2_z = p.z / r2
    r3_x = (p.x - p3.x) / r3; r3_y = (p.y - p3.y) / r3; r3_z = p.z / r3
    r4_x = (p.x - p4.x) / r4; r4_y = (p.y - p4.y) / r4 ;r4_z = p.z / r4

    s12_1_x = s12_2_x = - C12; s12_1_y = s12_2_y = - S12
    s23_2_x = s23_3_x = - C23; s23_2_y = s23_3_y = - S23
    s34_3_x = s34_4_x = - C34; s34_3_y = s34_4_y = - S34
    s41_4_x = s41_1_x = - C41; s41_4_y = s41_1_y = - S41

    F12_x = fabs(p.z) * (R12_x * (r1 * s12_2 - r2 * s12_1) + R12 * (r1_x * s12_2 + r1 * s12_2_x - (r2_x * s12_1 + r2 * s12_1_x)))
    F23_x = fabs(p.z) * (R23_x * (r2 * s23_3 - r3 * s23_2) + R23 * (r2_x * s23_3 + r2 * s23_3_x - (r3_x * s23_2 + r3 * s23_2_x)))
    F34_x = fabs(p.z) * (R34_x * (r3 * s34_4 - r4 * s34_3) + R34 * (r3_x * s34_4 + r3 * s34_4_x - (r4_x * s34_3 + r4 * s34_3_x)))
    F41_x = fabs(p.z) * (R41_x * (r4 * s41_1 - r1 * s41_4) + R41 * (r4_x * s41_1 + r4 * s41_1_x - (r1_x * s41_4 + r1 * s41_4_x)))

    F12_y = fabs(p.z) * (R12_y * (r1 * s12_2 - r2 * s12_1) + R12 * (r1_y * s12_2 + r1 * s12_2_y - (r2_y * s12_1 + r2 * s12_1_y)))
    F23_y = fabs(p.z) * (R23_y * (r2 * s23_3 - r3 * s23_2) + R23 * (r2_y * s23_3 + r2 * s23_3_y - (r3_y * s23_2 + r3 * s23_2_y)))
    F34_y = fabs(p.z) * (R34_y * (r3 * s34_4 - r4 * s34_3) + R34 * (r3_y * s34_4 + r3 * s34_4_y - (r4_y * s34_3 + r4 * s34_3_y)))
    F41_y = fabs(p.z) * (R41_y * (r4 * s41_1 - r1 * s41_4) + R41 * (r4_y * s41_1 + r4 * s41_1_y - (r1_y * s41_4 + r1 * s41_4_y)))

    F12_z = (p.z / fabs(p.z)) * R12 * (r1 * s12_2 - r2 * s12_1)
    F23_z = (p.z / fabs(p.z)) * R23 * (r2 * s23_3 - r3 * s23_2)
    F34_z = (p.z / fabs(p.z)) * R34 * (r3 * s34_4 - r4 * s34_3)
    F41_z = (p.z / fabs(p.z)) * R41 * (r4 * s41_1 - r1 * s41_4)

    G12_x = R12 * R12 * (r1_x * r2 + r1 * r2_x) + 2 * r1 * r2 * R12_x * R12 + p.z * p.z * (s12_1_x * s12_2 + s12_1 * s12_1_x)
    G23_x = R23 * R23 * (r2_x * r3 + r2 * r3_x) + 2 * r2 * r3 * R23_x * R23 + p.z * p.z * (s23_2_x * s23_3 + s23_2 * s23_2_x)
    G34_x = R34 * R34 * (r3_x * r4 + r3 * r4_x) + 2 * r3 * r4 * R34_x * R34 + p.z * p.z * (s34_3_x * s34_4 + s34_3 * s34_3_x)
    G41_x = R41 * R41 * (r4_x * r1 + r4 * r1_x) + 2 * r4 * r1 * R41_x * R41 + p.z * p.z * (s41_4_x * s41_1 + s41_4 * s41_4_x)

    G12_y = R12 * R12 * (r1_y * r2 + r1 * r2_y) + 2 * r1 * r2 * R12_y * R12 + p.z * p.z * (s12_1_y * s12_2 + s12_1 * s12_1_y)
    G23_y = R23 * R23 * (r2_y * r3 + r2 * r3_y) + 2 * r2 * r3 * R23_y * R23 + p.z * p.z * (s23_2_y * s23_3 + s23_2 * s23_2_y)
    G34_y = R34 * R34 * (r3_y * r4 + r3 * r4_y) + 2 * r3 * r4 * R34_y * R34 + p.z * p.z * (s34_3_y * s34_4 + s34_3 * s34_3_y)
    G41_y = R41 * R41 * (r4_y * r1 + r4 * r1_y) + 2 * r4 * r1 * R41_y * R41 + p.z * p.z * (s41_4_y * s41_1 + s41_4 * s41_4_y)
    
    G12_z = R12 * R12 * (r1_z * r2 + r1 * r1_z) + 2 * p.z * s12_1 * s12_2
    G23_z = R23 * R23 * (r2_z * r3 + r2 * r2_z) + 2 * p.z * s23_2 * s23_3
    G34_z = R34 * R34 * (r3_z * r4 + r3 * r3_z) + 2 * p.z * s34_3 * s34_4
    G41_z = R41 * R41 * (r4_z * r1 + r4 * r4_z) + 2 * p.z * s41_4 * s41_1

    J12_x = ((F12_x * G12 - F12 * G12_x) / (G12 * G12)) / (1 + (F12 / G12) * (F12 / G12))
    J23_x = ((F23_x * G23 - F23 * G23_x) / (G23 * G23)) / (1 + (F23 / G23) * (F23 / G23))
    J34_x = ((F34_x * G34 - F34 * G34_x) / (G34 * G34)) / (1 + (F34 / G34) * (F34 / G34))
    J41_x = ((F41_x * G41 - F41 * G41_x) / (G41 * G41)) / (1 + (F41 / G41) * (F41 / G41))

    J12_y = ((F12_y * G12 - F12 * G12_y) / (G12 * G12)) / (1 + (F12 / G12) * (F12 / G12))
    J23_y = ((F23_y * G23 - F23 * G23_y) / (G23 * G23)) / (1 + (F23 / G23) * (F23 / G23))
    J34_y = ((F34_y * G34 - F34 * G34_y) / (G34 * G34)) / (1 + (F34 / G34) * (F34 / G34))
    J41_y = ((F41_y * G41 - F41 * G41_y) / (G41 * G41)) / (1 + (F41 / G41) * (F41 / G41))

    J12_z = ((F12_z * G12 - F12 * G12_z) / (G12 * G12)) / (1 + (F12 / G12) * (F12 / G12))
    J23_z = ((F23_z * G23 - F23 * G23_z) / (G23 * G23)) / (1 + (F23 / G23) * (F23 / G23))
    J34_z = ((F34_z * G34 - F34 * G34_z) / (G34 * G34)) / (1 + (F34 / G34) * (F34 / G34))
    J41_z = ((F41_z * G41 - F41 * G41_z) / (G41 * G41)) / (1 + (F41 / G41) * (F41 / G41))

    if R12 < 0 and R23 < 0 and R34 < 0 and R41 < 0:
        delta = 2 * pi
    else:
        delta = 0.0
    
    # source
    source_pot = R12 * Q12 + R23 * Q23 + R34 * Q34 + R41 * Q41 + fabs(p.z) * (delta + J12 + J23 + J34 + J41)
    source_vel = Vec3D(
        x=(S12 * Q12 + S23 * Q23 + S34 * Q34 + S41 * Q41),
        y=-(C12 * Q12 + C23 * Q23 + C34 * Q34 + C41 * Q41),
        z=(p.z / fabs(p.z)) * (delta + J12 + J23 + J34 + J41),
    )

    # doublet
    doublet_pot = - source_vel.z
    doublet_vel = Vec3D(
        x=-(p.z / fabs(p.z)) * (J12_x + J23_x + J34_x + J41_x),
        y=-(p.z / fabs(p.z)) * (J12_y + J23_y + J34_y + J41_y),
        z=-(p.z / fabs(p.z)) * (J12_z + J23_z + J34_z + J41_z),
    )


    return [source_pot, source_vel, doublet_pot, doublet_vel]

def tri_panel(p: Vec3D, p1: Vec3D, p2: Vec3D, p3: Vec3D):

    d12 = sqrt( (p2.x - p1.x) * (p2.x - p1.x) + (p2.y - p1.y) * (p2.y - p1.y) )
    d23 = sqrt( (p3.x - p2.x) * (p3.x - p2.x) + (p3.y - p2.y) * (p3.y - p2.y) )
    d31 = sqrt( (p1.x - p3.x) * (p1.x - p3.x) + (p1.y - p3.y) * (p1.y - p3.y) )

    C12 = (p2.x - p1.x) / d12
    C23 = (p3.x - p2.x) / d23
    C31 = (p1.x - p3.x) / d31

    S12 = (p2.y - p1.y) / d12
    S23 = (p3.y - p2.y) / d23
    S31 = (p1.y - p3.y) / d31

    s12_1 = (p1.x - p.x) * C12 + (p1.y - p.y) * S12
    s12_2 = (p2.x - p.x) * C12 + (p2.y - p.y) * S12

    s23_2 = (p2.x - p.x) * C23 + (p2.y - p.y) * S23
    s23_3 = (p3.x - p.x) * C23 + (p3.y - p.y) * S23

    s31_3 = (p3.x - p.x) * C31 + (p3.y - p.y) * S31
    s31_1 = (p1.x - p.x) * C31 + (p1.y - p.y) * S31

    R12 = (p.x - p1.x) * S12 - (p.y - p1.y) * C12
    R23 = (p.x - p2.x) * S23 - (p.y - p2.y) * C23
    R31 = (p.x - p3.x) * S31 - (p.y - p3.y) * C31

    r1 = sqrt((p.x - p1.x) * (p.x - p1.x) + (p.y - p1.y) * (p.y - p1.y) + p.z * p.z)
    r2 = sqrt((p.x - p2.x) * (p.x - p2.x) + (p.y - p2.y) * (p.y - p2.y) + p.z * p.z)
    r3 = sqrt((p.x - p3.x) * (p.x - p3.x) + (p.y - p3.y) * (p.y - p3.y) + p.z * p.z)

    Q12 = log((r1 + r2 + d12) / (r1 + r2 - d12))
    Q23 = log((r2 + r3 + d23) / (r2 + r3 - d23))
    Q31 = log((r3 + r1 + d31) / (r3 + r1 - d31))

    F12 = R12 * fabs(p.z) * (r1 * s12_2 - r2 * s12_1)
    F23 = R23 * fabs(p.z) * (r2 * s23_3 - r3 * s23_2)
    F31 = R31 * fabs(p.z) * (r3 * s31_1 - r1 * s31_3)

    G12 = r1 * r2 * R12 * R12 + p.z * p.z * s12_1 * s12_2
    G23 = r2 * r3 * R23 * R23 + p.z * p.z * s23_2 * s23_3
    G31 = r3 * r1 * R31 * R31 + p.z * p.z * s31_3 * s31_1

    J12 = atan2(F12, G12)
    J23 = atan2(F23, G23)
    J31 = atan2(F31, G31)

    R12_x = S12; R12_y = - C12
    R23_x = S23; R23_y = - C23
    R31_x = S31; R31_y = - C31

    r1_x = (p.x - p1.x) / r1; r1_y = (p.y - p1.y) / r1; r1_z = p.z / r1
    r2_x = (p.x - p2.x) / r2; r2_y = (p.y - p2.y) / r2; r2_z = p.z / r2
    r3_x = (p.x - p3.x) / r3; r3_y = (p.y - p3.y) / r3; r3_z = p.z / r3

    s12_1_x = s12_2_x = - C12; s12_1_y = s12_2_y = - S12
    s23_2_x = s23_3_x = - C23; s23_2_y = s23_3_y = - S23
    s31_3_x = s31_1_x = - C31; s31_3_y = s31_1_y = - S31

    F12_x = fabs(p.z) * (R12_x * (r1 * s12_2 - r2 * s12_1) + R12 * (r1_x * s12_2 + r1 * s12_2_x - (r2_x * s12_1 + r2 * s12_1_x)))
    F23_x = fabs(p.z) * (R23_x * (r2 * s23_3 - r3 * s23_2) + R23 * (r2_x * s23_3 + r2 * s23_3_x - (r3_x * s23_2 + r3 * s23_2_x)))
    F31_x = fabs(p.z) * (R31_x * (r3 * s31_1 - r1 * s31_3) + R31 * (r3_x * s31_1 + r3 * s31_1_x - (r1_x * s31_3 + r1 * s31_3_x)))

    F12_y = fabs(p.z) * (R12_y * (r1 * s12_2 - r2 * s12_1) + R12 * (r1_y * s12_2 + r1 * s12_2_y - (r2_y * s12_1 + r2 * s12_1_y)))
    F23_y = fabs(p.z) * (R23_y * (r2 * s23_3 - r3 * s23_2) + R23 * (r2_y * s23_3 + r2 * s23_3_y - (r3_y * s23_2 + r3 * s23_2_y)))
    F31_y = fabs(p.z) * (R31_y * (r3 * s31_1 - r1 * s31_3) + R31 * (r3_y * s31_1 + r3 * s31_1_y - (r1_y * s31_3 + r1 * s31_3_y)))

    F12_z = (p.z / fabs(p.z)) * R12 * (r1 * s12_2 - r2 * s12_1)
    F23_z = (p.z / fabs(p.z)) * R23 * (r2 * s23_3 - r3 * s23_2)
    F31_z = (p.z / fabs(p.z)) * R31 * (r3 * s31_1 - r1 * s31_3)

    G12_x = R12 * R12 * (r1_x * r2 + r1 * r2_x) + 2 * r1 * r2 * R12_x * R12 + p.z * p.z * (s12_1_x * s12_2 + s12_1 * s12_1_x)
    G23_x = R23 * R23 * (r2_x * r3 + r2 * r3_x) + 2 * r2 * r3 * R23_x * R23 + p.z * p.z * (s23_2_x * s23_3 + s23_2 * s23_2_x)
    G31_x = R31 * R31 * (r3_x * r1 + r3 * r1_x) + 2 * r3 * r1 * R31_x * R31 + p.z * p.z * (s31_3_x * s31_1 + s31_3 * s31_3_x)

    G12_y = R12 * R12 * (r1_y * r2 + r1 * r2_y) + 2 * r1 * r2 * R12_y * R12 + p.z * p.z * (s12_1_y * s12_2 + s12_1 * s12_1_y)
    G23_y = R23 * R23 * (r2_y * r3 + r2 * r3_y) + 2 * r2 * r3 * R23_y * R23 + p.z * p.z * (s23_2_y * s23_3 + s23_2 * s23_2_y)
    G31_y = R31 * R31 * (r3_y * r1 + r3 * r1_y) + 2 * r3 * r1 * R31_y * R31 + p.z * p.z * (s31_3_y * s31_1 + s31_3 * s31_3_y)
 
    G12_z = R12 * R12 * (r1_z * r2 + r1 * r1_z) + 2 * p.z * s12_1 * s12_2
    G23_z = R23 * R23 * (r2_z * r3 + r2 * r2_z) + 2 * p.z * s23_2 * s23_3
    G31_z = R31 * R31 * (r3_z * r1 + r3 * r3_z) + 2 * p.z * s31_3 * s31_1

    J12_x = ((F12_x * G12 - F12 * G12_x) / (G12 * G12)) / (1 + (F12 / G12) * (F12 / G12))
    J23_x = ((F23_x * G23 - F23 * G23_x) / (G23 * G23)) / (1 + (F23 / G23) * (F23 / G23))
    J31_x = ((F31_x * G31 - F31 * G31_x) / (G31 * G31)) / (1 + (F31 / G31) * (F31 / G31))

    J12_y = ((F12_y * G12 - F12 * G12_y) / (G12 * G12)) / (1 + (F12 / G12) * (F12 / G12))
    J23_y = ((F23_y * G23 - F23 * G23_y) / (G23 * G23)) / (1 + (F23 / G23) * (F23 / G23))
    J31_y = ((F31_y * G31 - F31 * G31_y) / (G31 * G31)) / (1 + (F31 / G31) * (F31 / G31))

    J12_z = ((F12_z * G12 - F12 * G12_z) / (G12 * G12)) / (1 + (F12 / G12) * (F12 / G12))
    J23_z = ((F23_z * G23 - F23 * G23_z) / (G23 * G23)) / (1 + (F23 / G23) * (F23 / G23))
    J31_z = ((F31_z * G31 - F31 * G31_z) / (G31 * G31)) / (1 + (F31 / G31) * (F31 / G31))

    if R12 < 0 and R23 < 0 and R31 < 0:
        delta = 2 * pi
    else:
        delta = 0.0
    
    # source
    source_pot = R12 * Q12 + R23 * Q23 + R31 * Q31 + fabs(p.z) * (delta + J12 + J23 + J31)
    source_vel = Vec3D(
        x=(S12 * Q12 + S23 * Q23 + S31 * Q31),
        y=-(C12 * Q12 + C23 * Q23 + C31 * Q31),
        z=(p.z / fabs(p.z)) * (delta + J12 + J23 + J31),
    )

    # doublet
    doublet_pot = - source_vel.z
    doublet_vel = Vec3D(
        x=-(p.z / fabs(p.z)) * (J12_x + J23_x + J31_x),
        y=-(p.z / fabs(p.z)) * (J12_y + J23_y + J31_y),
        z=-(p.z / fabs(p.z)) * (J12_z + J23_z + J31_z),
    )


    return [source_pot, source_vel, doublet_pot, doublet_vel]


def test_quad_source():

    import numpy as np
    import matplotlib.pyplot as plt

    p1 = Vec3D(+1, +1, 0)
    p2 = Vec3D(-1, +1, 0)
    p3 = Vec3D(-1, -1, 0)
    p4 = Vec3D(+1, -1, 0)

    n = 100
    l = 5

    x = np.linspace(-l, l, num=n)
    y = np.linspace(-l, l, num=n)
    z = np.linspace(-l, l, num=n)

    xy, yx = np.meshgrid(x, y)
    z0 = -1e-1

    yz, zy = np.meshgrid(x, y)
    x0 = 0.0

    source_pot_xy = np.empty((n, n), dtype=np.double)
    source_vel_xy = np.empty((n, n, 3), dtype=np.double)

    source_pot_yz = np.empty((n, n), dtype=np.double)
    source_vel_yz = np.empty((n, n, 3), dtype=np.double)

    doublet_pot_xy = np.empty((n, n), dtype=np.double)
    doublet_vel_xy = np.empty((n, n, 3), dtype=np.double)

    doublet_pot_yz = np.empty((n, n), dtype=np.double)
    doublet_vel_yz = np.empty((n, n, 3), dtype=np.double)

    for i in range(n):
        for j in range(n):
            
            p = Vec3D(xy[i, j], yx[i, j], z0)
            pot1, vel1, pot2, vel2 = tri_panel(p, p1, p2, p3) # quad_panel(p, p1, p2, p3, p4)
            
            source_pot_xy[i, j] = pot1
            source_vel_xy[i, j, 0] = vel1.x
            source_vel_xy[i, j, 1] = vel1.y
            source_vel_xy[i, j, 2] = vel1.z

            doublet_pot_xy[i, j] = pot2
            doublet_vel_xy[i, j, 0] = vel2.x
            doublet_vel_xy[i, j, 1] = vel2.y
            doublet_vel_xy[i, j, 2] = vel2.z

            p = Vec3D(x0, yz[i, j], zy[i, j])
            pot1, vel1, pot2, vel2 = tri_panel(p, p1, p2, p3) # quad_panel(p, p1, p2, p3, p4)

            source_pot_yz[i, j] = pot1
            source_vel_yz[i, j, 0] = vel1.x
            source_vel_yz[i, j, 1] = vel1.y
            source_vel_yz[i, j, 2] = vel1.z

            doublet_pot_yz[i, j] = pot2
            doublet_vel_yz[i, j, 0] = vel2.x
            doublet_vel_yz[i, j, 1] = vel2.y
            doublet_vel_yz[i, j, 2] = vel2.z
    

    # source
    plt.figure()

    # xy
    plt.subplot(2, 2, 1)
    plt.contourf(xy, yx, source_pot_xy)
    plt.axis('equal')
    plt.colorbar()

    plt.subplot(2, 2, 2)
    plt.contourf(xy, yx, source_vel_xy[:, :, 2])
    plt.colorbar()
    plt.quiver(xy, yx, source_vel_xy[:, :, 0], source_vel_xy[:, :, 1])
    plt.axis('equal')

    # yz
    plt.subplot(2, 2, 3)
    plt.contourf(yz, zy, source_pot_yz)
    plt.axis('equal')
    plt.colorbar()

    plt.subplot(2, 2, 4)
    plt.contourf(yz, zy, source_vel_yz[:, :, 0])
    plt.colorbar()
    plt.quiver(yz, zy, source_vel_yz[:, :, 1], source_vel_yz[:, :, 2])
    plt.axis('equal')

    # doublet
    plt.figure()

    # xy
    plt.subplot(2, 2, 1)
    plt.contourf(xy, yx, doublet_pot_xy)
    plt.axis('equal')
    plt.colorbar()

    plt.subplot(2, 2, 2)
    plt.contourf(xy, yx, doublet_vel_xy[:, :, 2])
    plt.colorbar()
    plt.quiver(xy, yx, doublet_vel_xy[:, :, 0], doublet_vel_xy[:, :, 1])
    plt.axis('equal')

    # yz
    plt.subplot(2, 2, 3)
    plt.contourf(yz, zy, doublet_pot_yz)
    plt.axis('equal')
    plt.colorbar()

    plt.subplot(2, 2, 4)
    plt.contourf(yz, zy, doublet_vel_yz[:, :, 0])
    plt.colorbar()
    plt.quiver(yz, zy, doublet_vel_yz[:, :, 1], doublet_vel_yz[:, :, 2])
    plt.axis('equal')

    plt.show()

    return

if __name__ == '__main__':
    test_quad_source()
import typing as tp
import ctypes as ct
import numpy as np


class Vec3D(ct.Structure):
    _fields_ = [
        ('x', ct.c_double),
        ('y', ct.c_double),
        ('z', ct.c_double),
    ]

class Vec2D(ct.Structure):
    _fields_ = [
        ('x', ct.c_double),
        ('y', ct.c_double),
    ]

def main(sides: int,
         p1: np.ndarray,
         p2: np.ndarray,
         p3: np.ndarray,
         p4: np.ndarray,
         e1: np.ndarray,
         e2: np.ndarray,
         e3: np.ndarray,
         p_avg: np.ndarray,
         p_ctrl: np.ndarray) -> tp.List[np.ndarray]:

    n = p_ctrl.shape[0]

    p1_vec = Vec2D(p1[0], p1[1])
    p2_vec = Vec2D(p2[0], p2[1])
    p3_vec = Vec2D(p3[0], p3[1])
    p4_vec = Vec2D(p4[0], p4[1])

    e1_vec = Vec3D(e1[0], e1[1], e1[2])
    e2_vec = Vec3D(e2[0], e2[1], e2[2])
    e3_vec = Vec3D(e3[0], e3[1], e3[2])

    p_avg_vec = Vec3D(p_avg[0], p_avg[1], p_avg[2])

    doublet = np.empty(n, dtype=np.double)
    source = np.empty(n, dtype=np.double)

    lib = ct.CDLL('./src/pypm3D/bin/potential/potential.so')

    lib.panel_potential.argtypes = [
        ct.c_int,                         # n
        ct.c_int,                         # sides
        Vec2D,                            # p1
        Vec2D,                            # p2
        Vec2D,                            # p3
        Vec2D,                            # p4
        Vec3D,                            # e1
        Vec3D,                            # e2
        Vec3D,                            # e3
        Vec3D,                            # p_avg
        ct.POINTER(ct.c_double),          # p_ctrl
        ct.POINTER(ct.c_double),          # source
        ct.POINTER(ct.c_double),          # doublet
    ]

    lib.panel_potential.restype = None

    lib.panel_potential(
        n,
        sides,
        p1_vec,
        p2_vec,
        p3_vec,
        p4_vec,
        e1_vec,
        e2_vec,
        e3_vec,
        p_avg_vec,
        np.ctypeslib.as_ctypes(p_ctrl.reshape(p_ctrl.size)),
        np.ctypeslib.as_ctypes(source),
        np.ctypeslib.as_ctypes(doublet),
    )

    return [source, doublet]
import ctypes as ct
import numpy as np


def main(nv_te: int,
         nf: int,
         n_sides: np.ndarray,
         p_avg: np.ndarray,
         e1: np.ndarray, e2: np.ndarray, e3: np.ndarray,
         p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray,
         source: np.ndarray,
         doublet: np.ndarray,
         p_ctrls: np.ndarray) -> np.ndarray:
    
    vel_out = np.zeros_like(p_ctrls)
    vel = np.empty(3, dtype=np.double)
    p_ctrl = np.empty(3, dtype=np.double)

    lib = ct.CDLL('./src/pypm3D/modules/aero/utils/bin/lib.so')

    lib.surface_point_velocity.argtypes = [
        ct.c_int,                         # nf
        ct.POINTER(ct.c_int),             # n_sides
        ct.POINTER(ct.c_double),          # p_avg
        ct.POINTER(ct.c_double),          # e1
        ct.POINTER(ct.c_double),          # e2
        ct.POINTER(ct.c_double),          # e3
        ct.POINTER(ct.c_double),          # p1
        ct.POINTER(ct.c_double),          # p2
        ct.POINTER(ct.c_double),          # p3
        ct.POINTER(ct.c_double),          # p4
        ct.POINTER(ct.c_double),          # source
        ct.POINTER(ct.c_double),          # doublet
        ct.POINTER(ct.c_double),          # p_ctrl
        ct.POINTER(ct.c_double),          # vel
    ]

    lib.surface_point_velocity.restype = None

    for i in range(nv_te):

        p_ctrl[:] = p_ctrls[i, :]

        lib.surface_point_velocity(
            nf,
            np.ctypeslib.as_ctypes(np.asarray([x for x in n_sides])),
            np.ctypeslib.as_ctypes(p_avg.reshape(p_avg.size)),
            np.ctypeslib.as_ctypes(e1.reshape(e1.size)),
            np.ctypeslib.as_ctypes(e2.reshape(e2.size)),
            np.ctypeslib.as_ctypes(e3.reshape(e3.size)),
            np.ctypeslib.as_ctypes(p1.reshape(p1.size)),
            np.ctypeslib.as_ctypes(p2.reshape(p2.size)),
            np.ctypeslib.as_ctypes(p3.reshape(p3.size)),
            np.ctypeslib.as_ctypes(p4.reshape(p4.size)),
            np.ctypeslib.as_ctypes(source),
            np.ctypeslib.as_ctypes(doublet),
            np.ctypeslib.as_ctypes(p_ctrl),
            np.ctypeslib.as_ctypes(vel),
        )

        vel_out[i, :] = vel[:]
    
    return vel_out
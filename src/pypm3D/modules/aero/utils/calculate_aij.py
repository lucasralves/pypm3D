import ctypes as ct
import numpy as np


def main(nf: int,
         n_sides: np.ndarray,
         p_avg: np.ndarray,
         p_ctrl: np.ndarray,
         e1: np.ndarray,
         e2: np.ndarray,
         e3: np.ndarray,
         p1: np.ndarray,
         p2: np.ndarray,
         p3: np.ndarray,
         p4: np.ndarray,
         source: np.ndarray,
         a_j: np.ndarray) -> None:
    
    a_j_x = np.empty((nf,), dtype=np.double)
    a_j_y = np.empty((nf,), dtype=np.double)
    a_j_z = np.empty((nf,), dtype=np.double)
    
    lib = ct.CDLL('./src/pypm3D/modules/aero/utils/bin/lib.so')

    lib.get_a_j_coefs.argtypes = [
        ct.c_int,                         # nf
        ct.POINTER(ct.c_int),             # n_sides
        ct.POINTER(ct.c_double),          # p_avg
        ct.POINTER(ct.c_double),          # p_ctrl
        ct.POINTER(ct.c_double),          # e1
        ct.POINTER(ct.c_double),          # e2
        ct.POINTER(ct.c_double),          # e3
        ct.POINTER(ct.c_double),          # p1
        ct.POINTER(ct.c_double),          # p2
        ct.POINTER(ct.c_double),          # p3
        ct.POINTER(ct.c_double),          # p4
        ct.POINTER(ct.c_double),          # source
        ct.POINTER(ct.c_double),          # a_j_x
        ct.POINTER(ct.c_double),          # a_j_y
        ct.POINTER(ct.c_double),          # a_j_z
    ]

    lib.get_a_j_coefs.restype = None

    lib.get_a_j_coefs(
        nf,
        np.ctypeslib.as_ctypes(np.asarray([x for x in n_sides])),
        np.ctypeslib.as_ctypes(p_avg.reshape(p_avg.size)),
        np.ctypeslib.as_ctypes(p_ctrl.reshape(p_ctrl.size)),
        np.ctypeslib.as_ctypes(e1.reshape(e1.size)),
        np.ctypeslib.as_ctypes(e2.reshape(e2.size)),
        np.ctypeslib.as_ctypes(e3.reshape(e3.size)),
        np.ctypeslib.as_ctypes(p1.reshape(p1.size)),
        np.ctypeslib.as_ctypes(p2.reshape(p2.size)),
        np.ctypeslib.as_ctypes(p3.reshape(p3.size)),
        np.ctypeslib.as_ctypes(p4.reshape(p4.size)),
        np.ctypeslib.as_ctypes(source),
        np.ctypeslib.as_ctypes(a_j_x),
        np.ctypeslib.as_ctypes(a_j_y),
        np.ctypeslib.as_ctypes(a_j_z)
    )

    a_j[:, 0] = a_j_x[:]
    a_j[:, 1] = a_j_y[:]
    a_j[:, 2] = a_j_z[:]

    return
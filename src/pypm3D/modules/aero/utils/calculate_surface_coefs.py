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
         a_ij: np.ndarray,
         b_ij: np.ndarray) -> None:
    
    a_ij_aux = np.empty((nf * nf), dtype=np.double)
    b_ij_aux = np.empty((nf * nf), dtype=np.double)

    lib = ct.CDLL('./src/pypm3D/modules/aero/utils/bin/lib.so')

    lib.set_surface_coefs.argtypes = [
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
        ct.POINTER(ct.c_double),          # a_ij
        ct.POINTER(ct.c_double),          # b_ij
    ]

    lib.set_surface_coefs.restype = None

    lib.set_surface_coefs(
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
        np.ctypeslib.as_ctypes(a_ij_aux),
        np.ctypeslib.as_ctypes(b_ij_aux),
    )

    a_ij[:, :] = a_ij_aux.reshape((nf, nf))[:, :]
    b_ij[:, :] = b_ij_aux.reshape((nf, nf))[:, :]
    
    return
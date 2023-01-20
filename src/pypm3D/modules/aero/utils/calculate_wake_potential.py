import ctypes as ct
import numpy as np


def main(nf: int,
         nte: int,
         nw: int,
         faces: np.ndarray,
         vertices: np.ndarray,
         areas: np.ndarray,
         circulations: np.ndarray,
         p_ctrl: np.ndarray,
         c_ik: np.ndarray,
         d_i: np.ndarray) -> None:
    
    c_ik_aux = np.empty((nf * nte), dtype=np.double)
    d_i_aux = np.empty((nf), dtype=np.double)

    lib = ct.CDLL('./src/pypm3D/modules/aero/utils/bin/lib.so')

    lib.calculate_wake_potential.argtypes = [
        ct.c_int,                         # nf
        ct.c_int,                         # nte
        ct.c_int,                         # nw
        ct.POINTER(ct.c_int),             # faces
        ct.POINTER(ct.c_double),          # vertices
        ct.POINTER(ct.c_double),          # areas
        ct.POINTER(ct.c_double),          # circulations
        ct.POINTER(ct.c_double),          # p_ctrl
        ct.POINTER(ct.c_double),          # c_ik
        ct.POINTER(ct.c_double),          # d_i
    ]

    lib.calculate_wake_potential.restype = None

    lib.calculate_wake_potential(
        nf,
        nte,
        nw,
        np.ctypeslib.as_ctypes(faces.reshape(faces.size)),
        np.ctypeslib.as_ctypes(vertices.reshape(vertices.size)),
        np.ctypeslib.as_ctypes(areas),
        np.ctypeslib.as_ctypes(circulations),
        np.ctypeslib.as_ctypes(p_ctrl.reshape(p_ctrl.size)),
        np.ctypeslib.as_ctypes(c_ik_aux),
        np.ctypeslib.as_ctypes(d_i_aux),
    )

    c_ik[:, :] = c_ik_aux.reshape((nf, nte))[:, :]
    d_i[:] = d_i_aux[:]
    
    return
import ctypes as ct
import numpy as np


def main(nf: int,
         nte: int,
         vertices: np.ndarray,
         trailing_edge: np.ndarray,
         p_ctrl: np.ndarray,
         c_kj: np.ndarray) -> None:
    
    c_kj_x = np.empty(nf * nte, dtype=np.double)
    c_kj_y = np.empty(nf * nte, dtype=np.double)
    c_kj_z = np.empty(nf * nte, dtype=np.double)

    lib = ct.CDLL('./src/pypm3D/modules/aero/utils/bin/lib.so')

    lib.get_c_kj_coefs.argtypes = [
        ct.c_int,                         # nf
        ct.c_int,                         # nte
        ct.POINTER(ct.c_double),          # vertices
        ct.POINTER(ct.c_int),             # trailing_edge
        ct.POINTER(ct.c_double),          # p_ctrl
        ct.POINTER(ct.c_double),          # c_ij_x
        ct.POINTER(ct.c_double),          # c_ij_y
        ct.POINTER(ct.c_double),          # c_ij_z
    ]

    lib.get_c_kj_coefs.restype = None

    lib.get_c_kj_coefs(
        nf,
        nte,
        np.ctypeslib.as_ctypes(vertices.reshape(vertices.size)),
        np.ctypeslib.as_ctypes(trailing_edge.reshape(trailing_edge.size)),
        np.ctypeslib.as_ctypes(p_ctrl.reshape(p_ctrl.size)),
        np.ctypeslib.as_ctypes(c_kj_x),
        np.ctypeslib.as_ctypes(c_kj_y),
        np.ctypeslib.as_ctypes(c_kj_z)
    )

    c_kj_x = c_kj_x.reshape((nf, nte))
    c_kj_y = c_kj_y.reshape((nf, nte))
    c_kj_z = c_kj_z.reshape((nf, nte))

    c_kj[:, :, 0] = c_kj_x[:, :]
    c_kj[:, :, 1] = c_kj_y[:, :]
    c_kj[:, :, 2] = c_kj_z[:, :]

    return
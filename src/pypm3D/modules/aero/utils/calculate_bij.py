import ctypes as ct
import numpy as np


def main(nf: int,
         vertices: np.ndarray,
         faces: np.ndarray,
         p_ctrl: np.ndarray,
         b_ij: np.ndarray) -> None:
    
    b_ij_x = np.empty(nf * nf, dtype=np.double)
    b_ij_y = np.empty(nf * nf, dtype=np.double)
    b_ij_z = np.empty(nf * nf, dtype=np.double)

    lib = ct.CDLL('./src/pypm3D/modules/aero/utils/bin/lib.so')

    lib.get_b_ij_coefs.argtypes = [
        ct.c_int,                         # nf
        ct.POINTER(ct.c_double),          # vertices
        ct.POINTER(ct.c_int),             # faces
        ct.POINTER(ct.c_double),          # p_ctrl
        ct.POINTER(ct.c_double),          # b_ij_x
        ct.POINTER(ct.c_double),          # b_ij_y
        ct.POINTER(ct.c_double),          # b_ij_z
    ]

    lib.get_b_ij_coefs.restype = None

    lib.get_b_ij_coefs(
        nf,
        np.ctypeslib.as_ctypes(vertices.reshape(vertices.size)),
        np.ctypeslib.as_ctypes(faces.reshape(faces.size)),
        np.ctypeslib.as_ctypes(p_ctrl.reshape(p_ctrl.size)),
        np.ctypeslib.as_ctypes(b_ij_x),
        np.ctypeslib.as_ctypes(b_ij_y),
        np.ctypeslib.as_ctypes(b_ij_z)
    )

    b_ij_x = b_ij_x.reshape((nf, nf))
    b_ij_y = b_ij_y.reshape((nf, nf))
    b_ij_z = b_ij_z.reshape((nf, nf))

    b_ij[:, :, 0] = b_ij_x[:, :]
    b_ij[:, :, 1] = b_ij_y[:, :]
    b_ij[:, :, 2] = b_ij_z[:, :]

    return
import ctypes as ct
import numpy as np


def main(nv_te: int,
         nte: int,
         section: int,
         vertices: np.ndarray,
         faces: np.ndarray,
         areas: np.ndarray,
         circulations: np.ndarray,
         p_ctrls: np.ndarray) -> np.ndarray:
    
    vel_out = np.zeros_like(p_ctrls)
    vel = np.empty(3, dtype=np.double)
    p_ctrl = np.empty(3, dtype=np.double)

    lib = ct.CDLL('./src/pypm3D/modules/aero/utils/bin/lib.so')

    lib.wake_point_velocity.argtypes = [
        ct.c_int,                         # nte
        ct.c_int,                         # sections
        ct.POINTER(ct.c_double),          # vertices
        ct.POINTER(ct.c_int),             # faces
        ct.POINTER(ct.c_double),          # areas
        ct.POINTER(ct.c_double),          # circulations
        ct.POINTER(ct.c_double),          # p_ctrl
        ct.POINTER(ct.c_double),          # vel
    ]

    lib.wake_point_velocity.restype = None

    for i in range(nv_te):

        p_ctrl[:] = p_ctrls[i, :]

        lib.wake_point_velocity(
            nte,
            section,
            np.ctypeslib.as_ctypes(vertices.reshape(vertices.size)),
            np.ctypeslib.as_ctypes(faces.reshape(faces.size)),
            np.ctypeslib.as_ctypes(areas),
            np.ctypeslib.as_ctypes(circulations),
            np.ctypeslib.as_ctypes(p_ctrl),
            np.ctypeslib.as_ctypes(vel),
        )

        vel_out[i, :] = vel[:]
    
    return vel_out
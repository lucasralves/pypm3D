import typing as tp
import numpy as np


def main(nf: int,
         a_j: np.ndarray,
         b_ij: np.ndarray,
         c_kj: np.ndarray,
         d_j: np.ndarray,
         freestream: np.ndarray,
         e3: np.ndarray,
         doublet: np.ndarray,
         wake_doublet_sec_0: np.ndarray,
         func: tp.Callable[[np.ndarray], None]) -> None:
    
    vel = np.empty((nf, 3), dtype=np.double)
    cp = np.empty((nf,), dtype=np.double)
    transpiration = np.empty((nf,), dtype=np.double)

    vel[:, 0] = freestream[:, 0] + a_j[:, 0] + np.dot(b_ij[:, :, 0], doublet) + np.dot(c_kj[:, :, 0], wake_doublet_sec_0) + d_j[:, 0]
    vel[:, 1] = freestream[:, 1] + a_j[:, 1] + np.dot(b_ij[:, :, 1], doublet) + np.dot(c_kj[:, :, 1], wake_doublet_sec_0) + d_j[:, 1]
    vel[:, 2] = freestream[:, 2] + a_j[:, 2] + np.dot(b_ij[:, :, 2], doublet) + np.dot(c_kj[:, :, 2], wake_doublet_sec_0) + d_j[:, 2]

    cp[:] = 1.0 - (vel[:, 0] * vel[:, 0] + vel[:, 1] * vel[:, 1] + vel[:, 2] * vel[:, 2]) / (freestream[:, 0] * freestream[:, 0] + freestream[:, 1] * freestream[:, 1] + freestream[:, 2] * freestream[:, 2])

    transpiration[:] = vel[:, 0] * e3[:, 0] + vel[:, 1] * e3[:, 1] + vel[:, 2] * e3[:, 2]

    func(vel, cp, transpiration)
    
    return
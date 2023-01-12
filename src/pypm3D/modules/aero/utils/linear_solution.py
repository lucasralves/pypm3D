import typing as tp
import numpy as np


def main(nf: int,
         a_j: np.ndarray,
         b_ij: np.ndarray,
         c_kj: np.ndarray,
         freestream: np.ndarray,
         e3: np.ndarray,
         trailing_edge_faces: np.ndarray,
         func: tp.Callable[[np.ndarray], None]) -> None:
    
    rhs = np.empty((nf,), dtype=np.double)
    lhs = np.empty((nf, nf), dtype=np.double)

    rhs[:] = - ( (freestream[:, 0] + a_j[:, 0]) * e3[:, 0] + (freestream[:, 1] + a_j[:, 1]) * e3[:, 1] + (freestream[:, 2] + a_j[:, 2]) * e3[:, 2] )

    for face in range(nf):
        # Boundary condition
        lhs[face, :] = b_ij[face, :, 0] * e3[face, 0] + b_ij[face, :, 1] * e3[face, 1] + b_ij[face, :, 2] * e3[face, 2]

        # Kutta condition
        lhs[face, trailing_edge_faces[:, 0]] = lhs[face, trailing_edge_faces[:, 0]] + ( c_kj[face, :, 0] * e3[face, 0] + c_kj[face, :, 1] * e3[face, 1] + c_kj[face, :, 2] * e3[face, 2] )
        lhs[face, trailing_edge_faces[:, 1]] = lhs[face, trailing_edge_faces[:, 1]] - ( c_kj[face, :, 0] * e3[face, 0] + c_kj[face, :, 1] * e3[face, 1] + c_kj[face, :, 2] * e3[face, 2] )
    
    sol = np.linalg.solve(lhs, rhs)

    func(sol)

    return
import typing as tp
import numpy as np


def initial_condition(nf: int,
                      a_ij: np.ndarray,
                      b_ij: np.ndarray,
                      source: np.ndarray,
                      func: tp.Callable[[np.ndarray], None]) -> None:

    rhs = np.empty((nf,), dtype=np.double)
    lhs = np.empty((nf, nf), dtype=np.double)

    rhs[:] = - np.dot(a_ij, source)
    for face in range(nf):
        lhs[face, :] = b_ij[face, :]
    
    sol = np.linalg.solve(lhs, rhs)

    func(sol)

    return

def main(nf: int,
         a_ij: np.ndarray,
         b_ij: np.ndarray,
         c_ik: np.ndarray,
         d_i: np.ndarray,
         source: np.ndarray,
         trailing_edge_faces: np.ndarray,
         func: tp.Callable[[np.ndarray, np.ndarray], None]) -> None:

    rhs = np.empty((nf,), dtype=np.double)
    lhs = np.empty((nf, nf), dtype=np.double)

    rhs[:] = - np.dot(a_ij, source) - d_i[:]
    for face in range(nf):
        lhs[face, :] = b_ij[face, :]
        lhs[face, trailing_edge_faces[:, 0]] = lhs[face, trailing_edge_faces[:, 0]] - c_ik[face, :]
        lhs[face, trailing_edge_faces[:, 1]] = lhs[face, trailing_edge_faces[:, 1]] + c_ik[face, :]
     
    surf_doublet = np.linalg.solve(lhs, rhs)
    wake_doublet = - surf_doublet[trailing_edge_faces[:, 0]] + surf_doublet[trailing_edge_faces[:, 1]]

    func(surf_doublet, wake_doublet)

    return
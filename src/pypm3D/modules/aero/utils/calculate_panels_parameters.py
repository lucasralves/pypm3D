import typing as tp
import numpy as np

from pypm3D.models.vertice_model import VerticeModel

def _scalar(nv: int,
           vertices_connection: tp.List[VerticeModel],
           vals: np.ndarray) -> np.ndarray:

    vals_v = np.empty((nv,), dtype=np.double)

    for i in range(nv):
        vals_v[i] = .0
        for j in range(vertices_connection[i].n):
            vals_v[i] = vals_v[i] + vertices_connection[i].coefs[j] * vals[vertices_connection[i].faces[j]]
    
    return vals_v

def _vec(nv: int,
        vertices_connection: tp.List[VerticeModel],
        vals: np.ndarray) -> np.ndarray:

    vals_v = np.empty((nv, 3), dtype=np.double)

    for i in range(nv):
        vals_v[i, :] = .0
        for j in range(vertices_connection[i].n):
            vals_v[i, :] = vals_v[i, :] + vertices_connection[i].coefs[j] * vals[vertices_connection[i].faces[j], :]
    
    return vals_v

def main(nf: int,
         nv: int,
         a_ij: np.ndarray,
         b_ij: np.ndarray,
         faces: np.ndarray,
         source: np.ndarray, doublet: np.ndarray,
         p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray,
         e1: np.ndarray, e2: np.ndarray, e3: np.ndarray,
         freestream: np.ndarray,
         vertices_connection: tp.List[VerticeModel],
         func: tp.Callable[[np.ndarray, np.ndarray, np.ndarray, np.ndarray], None]) -> None:
    
    # doublet
    doublet_v = _scalar(nv, vertices_connection, doublet)

    # Faces parameters
    vel = np.empty((nf, 3), dtype=np.double)
    cp = np.empty((nf,), dtype=np.double)
    transpiration = np.empty((nf,), dtype=np.double)
    potential = np.empty((nf,), dtype=np.double)

    for face in range(nf):
        p0 = np.array([.0, .0, doublet[face]])

        v1 = np.array([p1[face, 0], p1[face, 1], doublet_v[faces[face, 1]]]) - p0
        v2 = np.array([p2[face, 0], p2[face, 1], doublet_v[faces[face, 2]]]) - p0
        v3 = np.array([p3[face, 0], p3[face, 1], doublet_v[faces[face, 3]]]) - p0

        n1 = np.cross(v1, v2)
        n2 = np.cross(v2, v3)

        if faces[face, 0] == 4:
            v4 = np.array([p4[face, 0], p4[face, 1], doublet_v[faces[face, 4]]]) - p0
            n3 = np.cross(v3, v4)
            n4 = np.cross(v4, v1)
            n = n1 + n2 + n3 + n4
        else:
            n3 = np.cross(v3, v1)
            n = n1 + n2 + n3
        
        vel[face, :] = (n[0] / n[2]) * e1[face, :] + (n[1] / n[2]) * e2[face, :] + freestream[face, :] - e3[face, :] * np.dot(freestream[face, :], e3[face, :])
    
    cp[:] = 1.0 - (vel[:, 0] * vel[:, 0] + vel[:, 1] * vel[:, 1] + vel[:, 2] * vel[:, 2]) / (freestream[:, 0] * freestream[:, 0] + freestream[:, 1] * freestream[:, 1] + freestream[:, 2] * freestream[:, 2])
    transpiration[:] = vel[:, 0] * e3[:, 0] + vel[:, 1] * e3[:, 1] + vel[:, 2] * e3[:, 2]
    potential[:] = np.dot(a_ij, source) + np.dot(b_ij, doublet)
    
    # Vertices parameters
    source_v = _scalar(nv, vertices_connection, source)
    vel_v = _vec(nv, vertices_connection, vel)
    cp_v = _scalar(nv, vertices_connection, cp)
    transpiration_v = _scalar(nv, vertices_connection, transpiration)
    potential_v = _scalar(nv, vertices_connection, potential)

    func(vel, cp, transpiration, potential, vel_v, cp_v, transpiration_v, potential_v, source_v, doublet_v)


    
    return
import typing as tp
import numpy as np

from pypm3D.models.vertice_model import VerticeModel


def main(nv: int,
         vertices_connection: tp.List[VerticeModel],
         vel: np.ndarray,
         cp: np.ndarray,
         transpiration: np.ndarray,
         func: tp.Callable[[np.ndarray], None]):

    vel_v = np.empty((nv, 3), dtype=np.double)
    cp_v = np.empty((nv,), dtype=np.double)
    transpiration_v = np.empty((nv,), dtype=np.double)

    for i in range(nv):

        vel_v[i, :] = .0
        cp_v[i] = .0
        transpiration_v[i] = .0

        for j in range(vertices_connection[i].n):
            
            vel_v[i, :] = vel_v[i, :] + vertices_connection[i].coefs[j] * vel[vertices_connection[i].faces[j], :]
            cp_v[i] = cp_v[i] + vertices_connection[i].coefs[j] * cp[vertices_connection[i].faces[j]]
            transpiration_v[i] = transpiration_v[i] + vertices_connection[i].coefs[j] * transpiration[vertices_connection[i].faces[j]]

    func(vel_v, cp_v, transpiration_v)

    return
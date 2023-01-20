import typing as tp
import numpy as np

from pypm3D.models.mesh_model import MeshModel
from pypm3D.models.solution_model import SolutionModel


def _tri_area(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray) -> float:
    return 0.5 * np.abs(p1[0] * p2[1] + p2[0] * p3[1] + p3[0] * p1[1] - (p3[0] * p2[1] + p2[0] * p1[1] + p1[0] * p3[1]))

def forces_and_moments(mesh: MeshModel,
                       sol: SolutionModel,
                       ref_area: float) -> tp.List[float]:

    force = np.array([.0, .0, .0])

    for face in range(mesh.surface.nf):

        if mesh.surface.faces[face, 0] == 4:
            area = _tri_area(mesh.surface.p1[face, :], mesh.surface.p2[face, :], mesh.surface.p3[face, :]) + _tri_area(mesh.surface.p1[face, :], mesh.surface.p3[face, :], mesh.surface.p4[face, :])
        else:
            area = _tri_area(mesh.surface.p1[face, :], mesh.surface.p2[face, :], mesh.surface.p3[face, :])
        
        force[:] = force[:] + mesh.surface.e3[face, :] * sol.cp[face] * area
    
    force[:] = force[:] / ref_area

    return force
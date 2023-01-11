import numpy as np

from pypm3D import models
from pypm3D.modules.mesh import core

def proc_mesh(vertices: np.ndarray,
              faces: np.ndarray,
              trailing_edge: np.ndarray) -> models.MeshModel:

    mesh = models.MeshModel(vertices, faces, trailing_edge)

    core.find_faces_at_trailing_edge(mesh)
    core.calculate_panels_parameters(mesh)
    core.create_internal_panels(mesh)
    core.calculate_scale_factor(mesh)

    return mesh
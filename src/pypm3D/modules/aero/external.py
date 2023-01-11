import typing as tp
import numpy as np

from pypm3D import models
from pypm3D.modules.aero import core

def solve(mesh: models.MeshModel,
          u_func: tp.Callable[[np.ndarray], np.ndarray],
          length: float,
          time_step: float) -> models.AeroModel:

    # Freestream
    freestream = np.empty((mesh.surface.nf, 3), dtype=np.double)
    for i in range(mesh.surface.nf): freestream[i, :] = u_func(mesh.surface.p_avg[i, :])

    # Connection between vertices
    vertices_connection = core.create_vertices_connection(mesh)
      
    # Indulced velocity coefficients
    a_ij = np.empty((mesh.surface.nf, mesh.surface.nf, 3), dtype=np.double)
    b_kj = np.empty((mesh.surface.nf, mesh.surface.nte, 3), dtype=np.double)
    c_kj = np.empty((mesh.surface.nf, mesh.surface.nte, 3), dtype=np.double)
    d_j = np.empty((mesh.surface.nf, 3), dtype=np.double)

    # Output
    aero = models.AeroModel(mesh.surface.nf, mesh.surface.nv, mesh.surface.nte)

    # Wake
    core.create_wake_mesh(mesh, u_func, length, time_step, aero)

    # Calculate coefficients
    core.set_a_ij(mesh, a_ij)
    core.set_b_kj(mesh, b_kj)
    core.set_c_kj(mesh, c_kj)

    # Initial condition
    core.solve_initial_condition(a_ij, b_kj, c_kj, mesh, freestream, aero)

    # Calculate surface and vertice parameters
    core.calculate_surface_parameters(a_ij, b_kj, c_kj, freestream, mesh, aero)
    core.calculate_vertices_parameters(vertices_connection, mesh, aero)
      
    # for i in range(aero.nw):
    #     core.add_wake_section()

    return aero
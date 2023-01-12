import typing as tp
import numpy as np
import functools as ft

from pypm3D.models.mesh_model import MeshModel
from pypm3D.models.solution_model import SolutionModel
from pypm3D.modules.aero.utils import (calculate_source_strength,
                                       calculate_aij,
                                       calculate_bij,
                                       calculate_ckj,
                                       calculate_ckj_and_dj,
                                       linear_solution,
                                       non_linear_solution,
                                       add_wake_section,
                                       calculate_panels_parameters,
                                       calculate_vertices_parameters,
                                       find_vertices_connection)


def solve(mesh: MeshModel,
          u_func: tp.Callable[[np.ndarray], np.ndarray],
          length: float,
          time_step: float,
          kutta_condition: str = 'linear') -> SolutionModel:
    """
    It finds the source and doublet distribution that satisfy
    the boundary layer and Kutta condition.

    Parameters:
    ----------
    - mesh: Object MeshModel
    - u_func: Callable that calulate the freestream at a given point
    - length: wake length
    - time_step: time step used to compute the wake
    - kutta_condition: must be linear or non-linear. The linear condition
      sets the wake circulation as the difference between the upper and
      lower panel. The non-linear condition sets the pressure on the upper
      and lower panel at the trailing edge to be equal.
    """

    # Wrappers
    def _add_source(sol: SolutionModel, vals: np.ndarray) -> None:
        sol.source[:] = vals[:]
        return
    
    def _add_singularities(trailing_edge_faces: np.ndarray, sol: SolutionModel, vals: np.ndarray) -> None:
        sol.doublet[:] = vals[:]
        sol.wake_doublet[:, 0] = vals[trailing_edge_faces[:, 0]] - vals[trailing_edge_faces[:, 1]]
        return
    
    def _add_faces_parameters(sol: SolutionModel, vel: np.ndarray, cp: np.ndarray, transpiration: np.ndarray) -> None:
        sol.vel[:, :] = vel[:, :]
        sol.cp[:] = cp[:]
        sol.transpiration[:] = transpiration[:]
        return
    
    def _add_vertices_parameters(sol: SolutionModel, vel: np.ndarray, cp: np.ndarray, transpiration: np.ndarray) -> None:
        sol.vel_v[:, :] = vel[:, :]
        sol.cp_v[:] = cp[:]
        sol.transpiration_v[:] = transpiration[:]
        return
    
    # Init wake mesh
    mesh.init_wake(u_func, length, time_step)
    
    # Create solution
    sol = SolutionModel()
    sol.init(mesh.surface.nf, mesh.surface.nv, mesh.wake.nte, mesh.wake.nw)

    # Local parameters
    freestream = np.asarray([u_func(x) for x in mesh.surface.p_avg], dtype=np.double)

    a_j = np.empty((mesh.surface.nf, 3), dtype=np.double)
    b_ij = np.empty((mesh.surface.nf, mesh.surface.nf, 3), dtype=np.double)
    c_kj = np.empty((mesh.surface.nf, mesh.surface.nte, 3), dtype=np.double)
    d_j = np.zeros((mesh.surface.nf, 3), dtype=np.double)

    # Calculate source
    func = ft.partial(_add_source, sol)
    calculate_source_strength.main(freestream, mesh.surface.e3, func)

    # Initialize a_j, b_ij, and c_kj
    calculate_aij.main(mesh.surface.nf, mesh.surface.faces[:, 0], mesh.surface.p_avg, mesh.surface.p_ctrl, mesh.surface.e1, mesh.surface.e2, mesh.surface.e3, mesh.surface.p1, mesh.surface.p2, mesh.surface.p3, mesh.surface.p4, sol.source, a_j)
    calculate_bij.main(mesh.surface.nf, mesh.surface.vertices, mesh.surface.faces, mesh.surface.p_ctrl, b_ij)
    calculate_ckj.main(mesh.surface.nf, mesh.surface.nte, mesh.surface.vertices, mesh.surface.trailing_edge, mesh.surface.p_ctrl, c_kj)

    # Solve using linear Kutta condition
    func = ft.partial(_add_singularities, mesh.surface.trailing_edge_faces, sol)
    linear_solution.main(mesh.surface.nf, a_j, b_ij, c_kj, freestream, mesh.surface.e3, mesh.surface.trailing_edge_faces, func)

    # Vertices connection
    vertices_connection = find_vertices_connection.main(mesh.surface.nv, mesh.surface.vertices, mesh.surface.faces)

    for section in range(mesh.wake.nw):

        # Calculate panels parameters
        func = ft.partial(_add_faces_parameters, sol)
        calculate_panels_parameters.main(mesh.surface.nf, a_j, b_ij, c_kj, d_j, freestream, mesh.surface.e3, sol.doublet, sol.wake_doublet[:, 0], func)

        # Calculate vertices parameters
        func = ft.partial(_add_vertices_parameters, sol)
        calculate_vertices_parameters.main(mesh.surface.nv, vertices_connection, sol.vel, sol.cp, sol.transpiration, func)

        break
        
        # Add wake section
        add_wake_section.main()

        # Calcula c_kj and d_j
        calculate_ckj_and_dj.main()

        # Solve using non linear Kutta condition
        # if kutta_condition == 'linear':
        #     linear_solution.main()
        # elif kutta_condition == 'non-linear':
        #     non_linear_solution.main()
    
    # # Calculate panels parameters
    # calculate_panels_parameters.main()

    # # Calculate vertices parameters
    # calculate_vertices_parameters.main()

    return sol
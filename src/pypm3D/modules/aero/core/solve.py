import typing as tp
import numpy as np
import functools as ft

from pypm3D.models.mesh_model import MeshModel
from pypm3D.models.solution_model import SolutionModel
from pypm3D.modules.aero.utils import warnings, calculate_surface_coefs, calculate_source_strength, linear_solution, find_vertices_connection, calculate_panels_parameters, add_wake_section, calculate_wake_potential, non_linear_solution


def solve(mesh: MeshModel,
          u_func: tp.Callable[[np.ndarray], np.ndarray],
          length: float,
          time_step: float,
          kutta_condition: str = 'linear',
          verbose: bool = False) -> SolutionModel:
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
    def add_surface_and_wake_doublet(sol: SolutionModel, mesh: MeshModel, surf_doublet: np.ndarray, wake_doublet: np.ndarray) -> None:
        sol.doublet[:] = surf_doublet[:]
        sol.wake_doublet[:mesh.wake.nte] = wake_doublet[:]
        return

    # Init wake mesh
    mesh.init_wake(u_func, length, time_step)

    # Create solution
    sol = SolutionModel()
    sol.init(mesh.surface.nf, mesh.surface.nv, mesh.wake.nte, mesh.wake.nw)

    # Freestream
    freestream = np.asarray([u_func(x) for x in mesh.surface.p_avg], dtype=np.double)

    # Calculate source strength
    calculate_source_strength.main(freestream, mesh.surface.e3, sol.add_source)

    # Potential
    a_ij = np.empty((mesh.surface.nf, mesh.surface.nf), dtype=np.double)
    b_ij = np.empty((mesh.surface.nf, mesh.surface.nf), dtype=np.double)
    c_ik = np.zeros((mesh.surface.nf, mesh.wake.nte), dtype=np.double)
    d_i = np.zeros((mesh.surface.nf), dtype=np.double)

    # Initialize a_j, b_ij
    calculate_surface_coefs.main(mesh.surface.nf, mesh.surface.faces[:, 0], mesh.surface.p_avg, mesh.surface.p_ctrl_minus, mesh.surface.e1, mesh.surface.e2, mesh.surface.e3, mesh.surface.p1, mesh.surface.p2, mesh.surface.p3, mesh.surface.p4, a_ij, b_ij)

    # Solve using linear Kutta condition
    warnings.title('Initial condition', verbose)
    linear_solution.initial_condition(mesh.surface.nf, a_ij, b_ij, sol.source, sol.add_doublet)

    # Vertices connection
    vertices_connection = find_vertices_connection.main(mesh.surface.nv, mesh.surface.vertices, mesh.surface.faces)

    # Calculate surface parameters
    calculate_panels_parameters.main(mesh.surface.nf, mesh.surface.nv, a_ij, b_ij, mesh.surface.faces, sol.source, sol.doublet, mesh.surface.p1, mesh.surface.p2, mesh.surface.p3, mesh.surface.p4, mesh.surface.e1, mesh.surface.e2, mesh.surface.e3, freestream, vertices_connection, sol.add_multiple)
    
    # Add wake
    add_solution = ft.partial(add_surface_and_wake_doublet, sol, mesh)
    warnings.title('Creating free wake [{} sections]'.format(mesh.wake.nw), verbose)
    for section in range(1, mesh.wake.nw + 1):
        warnings.title('section {}'.format(section), verbose)
        
        # Add wake section
        add_wake_section.main(mesh.wake.nte, mesh.wake.nv_te, mesh.surface.nf, section, mesh.wake.vertices, mesh.wake.faces, mesh.wake.trailing_edge_ids, mesh.wake.areas, sol.wake_doublet, sol.vel_v, mesh.surface.faces[:, 0], mesh.surface.p_avg, mesh.surface.e1, mesh.surface.e2, mesh.surface.e3, mesh.surface.p1, mesh.surface.p2, mesh.surface.p3, mesh.surface.p4, sol.source, sol.doublet, time_step, u_func)

        calculate_wake_potential.main(mesh.surface.nf, mesh.wake.nte, mesh.wake.nw, mesh.wake.faces, mesh.wake.vertices, mesh.wake.areas, sol.wake_doublet, mesh.surface.p_ctrl_minus, c_ik, d_i)

        # Solve using non linear Kutta condition
        if kutta_condition == 'linear':
            linear_solution.main(mesh.surface.nf, a_ij, b_ij, c_ik, d_i, sol.source, mesh.surface.trailing_edge_faces, add_solution)
        elif kutta_condition == 'non-linear':
            non_linear_solution.main()
    
    # # Calculate panels parameters
    # calculate_panels_parameters.main()

    # # Calculate vertices parameters
    # calculate_vertices_parameters.main()

    return sol
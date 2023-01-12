import typing as tp
import numpy as np

from pypm3D.models.mesh_model import MeshModel
from pypm3D.models.solution_model import SolutionModel
from pypm3D.modules.aero.utils import (calculate_aij,
                                       calculate_bij,
                                       calculate_ckj,
                                       calculate_ckj_and_dj,
                                       linear_solution,
                                       non_linear_solution,
                                       add_wake_section,
                                       calculate_panels_parameters,
                                       calculate_vertices_parameters)


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
    
    # Init wake mesh
    mesh.init_wake(u_func, length, time_step)
    
    # Create solution
    sol = SolutionModel()
    sol.init(mesh.surface.nf, mesh.surface.nv, mesh.wake.nte, mesh.wake.nw)

    # Local parameters
    freestream = np.asarray([u_func(x) for x in mesh.surface.p_avg], dtype=np.double)

    a_ij = np.empty((mesh.surface.nf, mesh.surface.nf, 3), dtype=np.double)
    b_ij = np.empty((mesh.surface.nf, mesh.surface.nf, 3), dtype=np.double)
    c_kj = np.empty((mesh.surface.nf, mesh.surface.nte, 3), dtype=np.double)
    d_j = np.zeros((mesh.surface.nf, 3), dtype=np.double)

    # Initialize a_ij, b_ij, and c_kj
    calculate_aij.main()
    calculate_bij.main()
    calculate_ckj.main()

    # Solve using linear Kutta condition
    linear_solution.main()

    for section in range(mesh.wake.nw):

        # Add wake section
        add_wake_section.main()

        # Calcula c_kj and d_j
        calculate_ckj_and_dj.main()

        # Solve using non linear Kutta condition
        if kutta_condition == 'linear':
            linear_solution.main()
        elif kutta_condition == 'non-linear':
            non_linear_solution.main()
    
    # Calculate panels parameters
    calculate_panels_parameters.main()

    # Calculate vertices parameters
    calculate_vertices_parameters.main()

    return
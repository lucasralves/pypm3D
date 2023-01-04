from typing import List
from abc import ABC, abstractmethod
from numpy import ndarray, empty, double
from numpy.linalg import solve

from pypm3D.models.mesh_model import MeshModel

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
class PotentialFlowAbs(ABC):
    """
    Solve a source/doublet distribuition using a free wake model
    """

    def __init__(self, mesh: MeshModel) -> None:
        pass

    @abstractmethod
    def solve_source_doublet(self, freestream: ndarray, traspiration: ndarray) -> None:
        """
        It calculates the surface source and wake doublet based on the
        velocity transpiration and Kutta condition.
        """
        pass

    @abstractmethod
    def update_wake(self) -> None:
        """
        It updates the wake position by shedding new panels into the
        freestream by trailing edge.
        """
        pass
    
    @abstractmethod
    def calculate_surface_parameters(self) -> List[ndarray]:
        """
        It calculates the velocity, pressure coefficient and transpiration
        at each panel
        """
        pass

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
class PotentialFlow(PotentialFlowAbs):

    def __init__(self, mesh: MeshModel) -> None:

        # Input
        self.__mesh = mesh

        # Current wake id
        self.__wake_id = 0
        
        # Influence coefficients tensors
        self.__a_ij = empty((mesh.nf, mesh.nf, 3), dtype=double)     # surface source
        self.__b_kj = empty((mesh.nf, mesh.nte, 3), dtype=double)    # t.e. inner doublet
        self.__c_kj = empty((mesh.nf, mesh.nte, 3), dtype=double)    # t.e. external doublet
        self.__d_j = empty((mesh.nf, 3), dtype=double)               # remaining wake

        # Linear system
        self.__matrix = empty((mesh.nf + 2 * mesh.nte, mesh.nf + 2 * mesh.nte), dtype=double)
        self.__array = empty((mesh.nf + 2 * mesh.nte,), dtype=double)

        # Initialize constant coefficients
        self.__calculate_constant_coefs()

        # Singularities
        self.__source = empty((mesh.nf,), dtype=double)
        self.__wake_circulation = empty((mesh.nte, mesh.nw), dtype=double)
        self.__wake_initial_area = empty((mesh.nte, mesh.nw), dtype=double)

        return
    
    """Implementation"""
    def solve_source_doublet(self, freestream: ndarray, traspiration: ndarray) -> None:

        # Initialize non constant coefficients
        self.__calculate_non_constant_coefs()

        # Create rhs
        self.__array[:] = traspiration[:] - (freestream[0] * self.__mesh.e3[:, 0] + freestream[1] * self.__mesh.e3[:, 1] + freestream[2] * self.__mesh.e3[:, 2])

        # Create lhs
        self.__matrix[:, :] = 0.0

        # Solve
        sol = solve(self.__matrix, self.__array)

        # Save
        self.__source[:] = sol[:self.__mesh.nf]
        self.__wake_circulation[:, 0] = sol[self.__mesh.nf:self.__mesh.nf + self.__mesh.nte]
        self.__wake_circulation[:, 1] = sol[self.__mesh.nf + self.__mesh.nte:self.__mesh.nf + 2 * self.__mesh.nte]

        return
    
    """Implementation"""
    def update_wake(self) -> None:
        return
    
    """Implementation"""
    def calculate_surface_parameters(self) -> List[ndarray]:
        return
    
    def __calculate_constant_coefs(self) -> None:
        """Calculate a_ij and b_kj"""
        return
    
    def __calculate_non_constant_coefs(self) -> None:
        """Calculate c_kj and d_j"""
        return
    
    def __calculate_point_velocity(self) -> None:
        """Calculate the velocity induced by the surface and wake at a given point"""
        return
    
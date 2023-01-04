from typing import List
from numpy import ndarray, empty, double, int32

from pypm3D.models.mesh_model import MeshModel
from pypm3D.models.solution_model import SolutionModel

def process_mesh_data(nv: int,
                      nf: int,
                      nte: int,
                      nw: int,
                      vertices: ndarray,
                      faces: ndarray,
                      trailing_edge: ndarray) -> MeshModel:
    
    nv_wake = 0
    te_ids = []

    for i in range(nte):
        if trailing_edge[i, 0] not in te_ids:
            te_ids.append(trailing_edge[i, 0])
            nv_wake += 1
        
        if trailing_edge[i, 1] not in te_ids:
            te_ids.append(trailing_edge[i, 1])
            nv_wake += 1

    mesh = MeshModel(
        vertices=vertices,
        faces=faces,
        p_avg=empty((nf, 3), dtype=double),
        p_ctrl=empty((nf, 3), dtype=double),
        e1=empty((nf, 3), dtype=double),
        e2=empty((nf, 3), dtype=double),
        e3=empty((nf, 3), dtype=double),
        p1=empty((nf, 2), dtype=double),
        p2=empty((nf, 2), dtype=double),
        p3=empty((nf, 2), dtype=double),
        p4=empty((nf, 2), dtype=double),
        te_faces=empty((nte, 2), dtype=int32),
        wake_faces=empty((nte, nw, 5), dtype=int32),
        wake_vertices=empty((nv_wake * (nw + 1), 3), dtype=int32),
        wake_panel_scale=empty((nte, 2), dtype=double),
    )

    return

def calculate_a_b_coefs(mesh: MeshModel, a: ndarray, b: ndarray) -> None:
    return

def calculate_c_d_coefs(mesh: MeshModel, c: ndarray, d: ndarray) -> None:
    return

def create_and_solve_linear_system() -> ndarray:
    return

def update_wake() -> None:
    return

def solve(vertices: ndarray,
          faces: ndarray,
          trailing_edge: ndarray,
          freestream: float,
          alpha: float,
          beta: float,
          wake_length: float,
          wake_time_step: float) -> List[MeshModel, SolutionModel]:
    """
    Solves the potential flow around an arbitrary lifting surface

    Parameters:
    -----------
    - vertices (ndarray): it contains all panel's vertices with the shape
    (nv, 3), where nv is the number of vertices and each line correspond
    to a point in space.
    
    - faces (ndarray): it contains all faces with the shape (nf, 5), where
    nf is the number of faces. Each line corresponds to a panel, where the
    first element is the number of points that can be 3 or 4 (the mesh may
    contain only tri or quad panels), and the rest are the vertices ids that
    make the panel (positively oriented).

    - trailing_edge (ndarray): it contains the pair of vertices ids of an edge
    in the trailing edge with the shape (nte, 2), where nte is the number of
    edges in the trailing edge.

    - freestream: freestream velocity [m/s]

    - alpha: attack angle [degrees]

    - beta: side slip angle [degrees]

    - wake_length: wake length [m]

    - wake_time_step: time step used to generate the free wake [s]
    """

    # Parameters
    nv = vertices.shape[0]
    nf = faces.shape[0]
    nte = trailing_edge.shape[0]
    nw = int(wake_length / (freestream * wake_time_step)) + 1

    # a = empty((nf + 2 * nte, nf + 2 * nte, 3), dtype=double)        # surface surface coef
    # b = empty((nf, nte, 3), dtype=double)                           # te inner doublet coef
    # c = empty((nf, nte, 3), dtype=double)                           # te external doublet coef
    # d = empty((nf, 3), dtype=double)                                # wake coef (remaining)

    # source = empty((nf,), dtype=double)
    # doublet = empty((nte, nw), dtype=double)

    # Process mesh data
    mesh = process_mesh_data(vertices, faces, trailing_edge)

    # # Calculate constant coefs
    # calculate_a_b_coefs(mesh, a, b)

    # # Create wake
    # for i in range(nw):

    #     # Calculate coefs
    #     calculate_c_d_coefs(mesh, c, d)

    #     # Create system
    #     sol = create_and_solve_linear_system()

    #     # Save variables
    #     source[:] = sol[:nf]
    #     doublet[:nte, 0] = sol[nf:nf + nte]
    #     doublet[:nte, 1] = sol[nf + nte:nf + 2 * nte]

    #     # Update wake
    #     update_wake()

    return []
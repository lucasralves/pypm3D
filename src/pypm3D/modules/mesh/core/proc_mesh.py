import numpy as np
import functools as ft

from pypm3D.models.mesh_model import MeshModel
from pypm3D.modules.mesh.utils import (find_faces_at_trailing_edge,
                                       calculate_panels_parameters,
                                       correct_trailing_edge_faces)

def proc_mesh(vertices: np.ndarray,
              faces: np.ndarray,
              trailing_edge: np.ndarray) -> MeshModel:
    """
    Process mesh. The longitudinal axis must be align with the x axis,
    and the span must be align with the y axis.

    Parameters:
    -----------
    - vertices: points in space [x, y, z]
    - faces: number of sides and vertices ids [n, id1, id2, id3, id4]
    - trailing_edge: vertices that makes an edge [id1, id2]
    """

    # Wrappers
    def _add_trailing_edge_ids(mesh: MeshModel, index: int, id1: int, id2: int) -> None:
        mesh.surface.trailing_edge_faces[index, 0] = id1
        mesh.surface.trailing_edge_faces[index, 1] = id2
        return

    def _add_face_parameters(mesh: MeshModel, p_avg: np.ndarray, p_ctrl_plus: np.ndarray, p_ctrl_minus: np.ndarray, e1: np.ndarray, e2: np.ndarray, e3: np.ndarray, p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, p4: np.ndarray) -> None:
        mesh.surface.p_avg[:, :] = p_avg[:, :]
        mesh.surface.p_ctrl_plus[:, :] = p_ctrl_plus[:, :]
        mesh.surface.p_ctrl_minus[:, :] = p_ctrl_minus[:, :]
        mesh.surface.e1[:, :] = e1[:, :]
        mesh.surface.e2[:, :] = e2[:, :]
        mesh.surface.e3[:, :] = e3[:, :]
        mesh.surface.p1[:, :] = p1[:, :]
        mesh.surface.p2[:, :] = p2[:, :]
        mesh.surface.p3[:, :] = p3[:, :]
        mesh.surface.p4[:, :] = p4[:, :]
        return
    
    def _correct_trailing_edge_faces_ids(mesh: MeshModel, trailing_edge_faces: np.ndarray) -> None:
        mesh.surface.trailing_edge_faces[:, :] = trailing_edge_faces
        return

    # Create mesh
    mesh = MeshModel()

    # Init surface
    mesh.init_surface(vertices, faces, trailing_edge)

    # Process mesh
    func = ft.partial(_add_trailing_edge_ids, mesh)
    find_faces_at_trailing_edge.main(mesh.surface.nte, mesh.surface.faces, mesh.surface.trailing_edge, func)

    func = ft.partial(_add_face_parameters, mesh)
    calculate_panels_parameters.main(mesh.surface.nf, mesh.surface.vertices, mesh.surface.faces, mesh.surface.trailing_edge, mesh.surface.trailing_edge_faces, func)

    func = ft.partial(_correct_trailing_edge_faces_ids, mesh)
    correct_trailing_edge_faces.main(mesh.surface.trailing_edge_faces, mesh.surface.e3, func)

    return mesh
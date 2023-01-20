import typing as tp
import numpy as np


def main(nf: int,
         vertices: np.ndarray,
         faces: np.ndarray,
         trailing_edge: np.ndarray,
         trailing_edge_faces: np.ndarray,
         func: tp.Callable[[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray], None]) -> None:

    # Create parameters
    p_avg = np.empty((nf, 3), dtype=np.double)
    p_ctrl_plus = np.empty((nf, 3), dtype=np.double)
    p_ctrl_minus = np.empty((nf, 3), dtype=np.double)
    e1 = np.empty((nf, 3), dtype=np.double)
    e2 = np.empty((nf, 3), dtype=np.double)
    e3 = np.empty((nf, 3), dtype=np.double)
    p1 = np.empty((nf, 2), dtype=np.double)
    p2 = np.empty((nf, 2), dtype=np.double)
    p3 = np.empty((nf, 2), dtype=np.double)
    p4 = np.empty((nf, 2), dtype=np.double)

    # Panel geometry
    tri_panel_ids = faces[:, 0] == 3
    quad_panel_ids = faces[:, 0] == 4

    # Panel center
    p_avg[tri_panel_ids, :] = (1. / 3.) * (vertices[faces[tri_panel_ids, 1], :] + vertices[faces[tri_panel_ids, 2], :] + vertices[faces[tri_panel_ids, 3], :])
    p_avg[quad_panel_ids, :] = (1. / 4.) * (vertices[faces[quad_panel_ids, 1], :] + vertices[faces[quad_panel_ids, 2], :] + vertices[faces[quad_panel_ids, 3], :] + vertices[faces[quad_panel_ids, 4], :])
    
    # Orthogonal base (e3)
    e3[tri_panel_ids, :] = np.cross(vertices[faces[tri_panel_ids, 2], :] - vertices[faces[tri_panel_ids, 1], :], vertices[faces[tri_panel_ids, 3], :] - vertices[faces[tri_panel_ids, 1], :])
    e3[quad_panel_ids, :] = np.cross(vertices[faces[quad_panel_ids, 2], :] - vertices[faces[quad_panel_ids, 4], :], vertices[faces[quad_panel_ids, 3], :] - vertices[faces[quad_panel_ids, 1], :])
    
    # Orthogonal base (e1)
    e1[:, :] = 0.5 * (vertices[faces[:, 1], :] + vertices[faces[:, 2], :]) - p_avg

    # Orthogonal base (e2)
    e2[:, :] = np.cross(e3, e1)

    # Orthogonal base (e3 - correção)
    e3[:, :] = np.cross(e1, e2)

    # Correct trailing edge faces
    e1[trailing_edge_faces[:, 0], :] = vertices[trailing_edge[:, 0], :] - vertices[trailing_edge[:, 1], :]
    e1[trailing_edge_faces[:, 1], :] = e1[trailing_edge_faces[:, 0], :]

    # Upper panel
    e2[trailing_edge_faces[:, 0], :] = np.cross(e3[trailing_edge_faces[:, 0], :], e1[trailing_edge_faces[:, 0], :])
    e3[trailing_edge_faces[:, 0], :] = np.cross(e1[trailing_edge_faces[:, 0], :], e2[trailing_edge_faces[:, 0], :])

    # Lower panel
    e2[trailing_edge_faces[:, 1], :] = np.cross(e3[trailing_edge_faces[:, 1], :], e1[trailing_edge_faces[:, 1], :])
    e3[trailing_edge_faces[:, 1], :] = np.cross(e1[trailing_edge_faces[:, 1], :], e2[trailing_edge_faces[:, 1], :])

    # Normalize orthogonal base
    aux_norm = np.sqrt(e1[:, 0] * e1[:, 0] + e1[:, 1] * e1[:, 1] + e1[:, 2] * e1[:, 2])

    e1[:, 0] = e1[:, 0] / aux_norm
    e1[:, 1] = e1[:, 1] / aux_norm
    e1[:, 2] = e1[:, 2] / aux_norm

    aux_norm = np.sqrt(e2[:, 0] * e2[:, 0] + e2[:, 1] * e2[:, 1] + e2[:, 2] * e2[:, 2])

    e2[:, 0] = e2[:, 0] / aux_norm
    e2[:, 1] = e2[:, 1] / aux_norm
    e2[:, 2] = e2[:, 2] / aux_norm

    aux_norm = np.sqrt(e3[:, 0] * e3[:, 0] + e3[:, 1] * e3[:, 1] + e3[:, 2] * e3[:, 2])

    e3[:, 0] = e3[:, 0] / aux_norm
    e3[:, 1] = e3[:, 1] / aux_norm
    e3[:, 2] = e3[:, 2] / aux_norm

    # Control point
    p_ctrl_plus[:, :] = p_avg[:, :] + e3[:, :] * 1e-12
    p_ctrl_minus[:, :] = p_avg[:, :] - e3[:, :] * 1e-12

    # Local vertices
    aux_vec = vertices[faces[:, 1], :] - p_avg
    p1[:, 0] = e1[:, 0] * aux_vec[:, 0] + e1[:, 1] * aux_vec[:, 1] + e1[:, 2] * aux_vec[:, 2]
    p1[:, 1] = e2[:, 0] * aux_vec[:, 0] + e2[:, 1] * aux_vec[:, 1] + e2[:, 2] * aux_vec[:, 2]

    aux_vec = vertices[faces[:, 2], :] - p_avg
    p2[:, 0] = e1[:, 0] * aux_vec[:, 0] + e1[:, 1] * aux_vec[:, 1] + e1[:, 2] * aux_vec[:, 2]
    p2[:, 1] = e2[:, 0] * aux_vec[:, 0] + e2[:, 1] * aux_vec[:, 1] + e2[:, 2] * aux_vec[:, 2]

    aux_vec = vertices[faces[:, 3], :] - p_avg
    p3[:, 0] = e1[:, 0] * aux_vec[:, 0] + e1[:, 1] * aux_vec[:, 1] + e1[:, 2] * aux_vec[:, 2]
    p3[:, 1] = e2[:, 0] * aux_vec[:, 0] + e2[:, 1] * aux_vec[:, 1] + e2[:, 2] * aux_vec[:, 2]

    aux_vec = vertices[faces[quad_panel_ids, 4], :] - p_avg[quad_panel_ids, :]
    p4[quad_panel_ids, 0] = e1[quad_panel_ids, 0] * aux_vec[:, 0] + e1[quad_panel_ids, 1] * aux_vec[:, 1] + e1[quad_panel_ids, 2] * aux_vec[:, 2]
    p4[quad_panel_ids, 1] = e2[quad_panel_ids, 0] * aux_vec[:, 0] + e2[quad_panel_ids, 1] * aux_vec[:, 1] + e2[quad_panel_ids, 2] * aux_vec[:, 2]

    # Save
    func(p_avg, p_ctrl_plus, p_ctrl_minus, e1, e2, e3, p1, p2, p3, p4)

    return

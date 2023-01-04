from typing import Dict
from numpy import ndarray, argwhere, empty, double, int32, asarray, dot, cross, sqrt, zeros


def get_trailing_edge_faces(faces: ndarray,
                            trailing_edge: ndarray) -> Dict[str, ndarray]:
    """
    It calculates the upper and lower panel of the trailing edge.
    """

    trailing_edge_faces = empty((trailing_edge.shape[0], 2), dtype=int32)

    for i in range(trailing_edge.shape[0]):
        check1 = (faces[:, 1] == trailing_edge[i, 0]) | (faces[:, 2] == trailing_edge[i, 0]) | (faces[:, 3] == trailing_edge[i, 0]) | (faces[:, 4] == trailing_edge[i, 0])
        check2 = (faces[:, 1] == trailing_edge[i, 1]) | (faces[:, 2] == trailing_edge[i, 1]) | (faces[:, 3] == trailing_edge[i, 1]) | (faces[:, 4] == trailing_edge[i, 1])
        index = argwhere(check1 & check2)
        trailing_edge_faces[i, 0] = index[0][0]
        trailing_edge_faces[i, 1] = index[1][0]
    
    out = {'trailing_edge_faces': trailing_edge_faces}
    return out

def calculate_faces_parameters(vertices: ndarray,
                               faces: ndarray,
                               trailing_edge: ndarray,
                               trailing_edge_faces: ndarray) -> Dict[str, ndarray]:

    # Output
    p_avg = empty((faces.shape[0], 3), dtype=double)
    p_ctrl = empty((faces.shape[0], 3), dtype=double)
    e1 = empty((faces.shape[0], 3), dtype=double)
    e2 = empty((faces.shape[0], 3), dtype=double)
    e3 = empty((faces.shape[0], 3), dtype=double)
    p1 = empty((faces.shape[0], 2), dtype=double)
    p2 = empty((faces.shape[0], 2), dtype=double)
    p3 = empty((faces.shape[0], 2), dtype=double)
    p4 = empty((faces.shape[0], 2), dtype=double)

    # Panel type
    tri_panel_ids = faces[:, 0] == 3
    quad_panel_ids = faces[:, 0] == 4

    # Panel center
    p_avg[tri_panel_ids] = (1. / 3.) * (vertices[faces[tri_panel_ids, 1], :] + vertices[faces[tri_panel_ids, 2], :] + vertices[faces[tri_panel_ids, 3], :])
    p_avg[quad_panel_ids] = (1. / 4.) * (vertices[faces[quad_panel_ids, 1], :] + vertices[faces[quad_panel_ids, 2], :] + vertices[faces[quad_panel_ids, 3], :] + vertices[faces[quad_panel_ids, 4], :])

    # Orthogonal base (e1)
    te_faces_ids = []
    for i in range(faces.shape[0]):
        check1 = (faces[i, 1] in trailing_edge_faces[:, 0]) or (faces[i, 2] in trailing_edge_faces[:, 0]) or (faces[i, 3] in trailing_edge_faces[:, 0]) or (faces[i, 4] in trailing_edge_faces[:, 0])
        check2 = (faces[i, 1] in trailing_edge_faces[:, 1]) or (faces[i, 2] in trailing_edge_faces[:, 1]) or (faces[i, 3] in trailing_edge_faces[:, 1]) or (faces[i, 4] in trailing_edge_faces[:, 1])
        
        if check1 or check2:
            te_faces_ids.append([True, ])

    te_faces_ids = asarray([((face[1] in trailing_edge_faces[:, 0]) or (face[2] in trailing_edge_faces[:, 0]) or (face[3] in trailing_edge_faces[:, 0]) or (face[4] in trailing_edge_faces[:, 0])) and ((face[1] in trailing_edge_faces[:, 1]) or (face[2] in trailing_edge_faces[:, 1]) or (face[3] in trailing_edge_faces[:, 1]) or (face[4] in trailing_edge_faces[:, 1])) for face in faces])
    
    e1[te_faces_ids, :] = vertices[trailing_edge[te_faces_ids, 1], :] - vertices[trailing_edge[te_faces_ids, 0], :]
    e1[not te_faces_ids, :] = 0.5 * (vertices[faces[not te_faces_ids, 1], :] + vertices[faces[not te_faces_ids, 2], :]) - p_avg[not te_faces_ids, :]

    aux_norm = sqrt(e1[:, 0] * e1[:, 0] + e1[:, 1] * e1[:, 1] + e1[:, 2] * e1[:, 2])

    e1[:, 0] = e1[:, 0] / aux_norm
    e1[:, 1] = e1[:, 1] / aux_norm
    e1[:, 2] = e1[:, 2] / aux_norm

    # Orthogonal base (e3)
    e3[tri_panel_ids, :] = cross(vertices[faces[tri_panel_ids, 2], :] - vertices[faces[tri_panel_ids, 1], :], vertices[faces[tri_panel_ids, 3], :] - vertices[faces[tri_panel_ids, 1], :])
    e3[quad_panel_ids, :] = cross(vertices[faces[quad_panel_ids, 2], :] - vertices[faces[quad_panel_ids, 4], :], vertices[faces[quad_panel_ids, 3], :] - vertices[faces[quad_panel_ids, 1], :])

    aux_norm = sqrt(e3[:, 0] * e3[:, 0] + e3[:, 1] * e3[:, 1] + e3[:, 2] * e3[:, 2])

    e3[:, 0] = e3[:, 0] / aux_norm
    e3[:, 1] = e3[:, 1] / aux_norm
    e3[:, 2] = e3[:, 2] / aux_norm

    # Orthogonal base (e2)
    e2[:, :] = cross(e3, e1)

    aux_norm = sqrt(e2[:, 0] * e2[:, 0] + e2[:, 1] * e2[:, 1] + e2[:, 2] * e2[:, 2])

    e2[:, 0] = e2[:, 0] / aux_norm
    e2[:, 1] = e2[:, 1] / aux_norm
    e2[:, 2] = e2[:, 2] / aux_norm

    # Orthogonal base (e3 correction)
    e3[:, :] = cross(e1, e2)

    aux_norm = sqrt(e3[:, 0] * e3[:, 0] + e3[:, 1] * e3[:, 1] + e3[:, 2] * e3[:, 2])

    e3[:, 0] = e3[:, 0] / aux_norm
    e3[:, 1] = e3[:, 1] / aux_norm
    e3[:, 2] = e3[:, 2] / aux_norm

    # Control point
    p_ctrl[:, :] = p_avg + e3 * 1e-12

    # Local points
    p1[:, 0] = dot(vertices[faces[:, 1], :] - p_avg, e1.T)
    p1[:, 1] = dot(vertices[faces[:, 1], :] - p_avg, e2.T)

    p2[:, 0] = dot(vertices[faces[:, 2], :] - p_avg, e1.T)
    p2[:, 1] = dot(vertices[faces[:, 2], :] - p_avg, e2.T)

    p3[:, 0] = dot(vertices[faces[:, 3], :] - p_avg, e1.T)
    p3[:, 1] = dot(vertices[faces[:, 3], :] - p_avg, e2.T)

    p4[quad_panel_ids, 0] = dot(vertices[faces[quad_panel_ids, 3], :] - p_avg[quad_panel_ids, :], e1[quad_panel_ids, :].T)
    p4[quad_panel_ids, 1] = dot(vertices[faces[quad_panel_ids, 3], :] - p_avg[quad_panel_ids, :], e2[quad_panel_ids, :].T)

    # Output
    out = {'p_avg': p_avg, 'p_ctrl': p_ctrl, 'e1': e1, 'e2': e2, 'e3': e3, 'p1': p1, 'p2': p2, 'p3': p3, 'p4': p4}
    
    return out

def calculate_inner_doublet_pannels(vertices: ndarray,
                                    faces: ndarray,
                                    trailing_edge: ndarray,
                                    trailing_edge_faces: ndarray,
                                    e3: ndarray) -> Dict[str, ndarray]:

    def __find_ids(a: ndarray, b: ndarray, face: int) -> int:

        n = a.shape[0] - 1
        ids = None

        # Find the start id
        for i in range(n):
            for j in range(n):

                if a[i + 1] == b[j + 1]:
                    ids = [i + 1, j + 1]
                    break
                
            if ids is not None:
                break
            
        # Create new arrays
        flip_face = True if e3[trailing_edge_faces[face, 0], 2] < 0 else False

        if a[0] == 4:

            if ids[0] == 1:
                anew = [a[1], a[4], a[3], a[2]] if flip_face else [a[1], a[2], a[3], a[4]]
            elif ids[0] == 2:
                anew = [a[2], a[1], a[4], a[3]] if flip_face else [a[2], a[3], a[4], a[1]]
            elif ids[0] == 3:
                anew = [a[3], a[2], a[1], a[4]] if flip_face else [a[3], a[4], a[1], a[2]]
            elif ids[0] == 4:
                anew = [a[4], a[3], a[2], a[1]] if flip_face else [a[4], a[1], a[2], a[3]]
                
            flip_face = not flip_face
                
            if ids[1] == 1:
                bnew = [b[1], b[4], b[3], b[2]] if flip_face else [b[1], b[2], b[3], b[4]]
            elif ids[1] == 2:
                bnew = [b[2], b[1], b[4], b[3]] if flip_face else [b[2], b[3], b[4], b[1]]
            elif ids[1] == 3:
                bnew = [b[3], b[2], b[1], b[4]] if flip_face else [b[3], b[4], b[1], b[2]]
            elif ids[1] == 4:
                bnew = [b[4], b[3], b[2], b[1]] if flip_face else [b[4], b[1], b[2], b[3]]
            
        if a[0] == 3:

            if ids[0] == 1:
                anew = [a[1], a[3], a[2]] if flip_face else [a[1], a[2], a[3]]
            elif ids[0] == 2:
                anew = [a[2], a[1], a[3]] if flip_face else [a[2], a[3], a[1]]
            elif ids[0] == 3:
                anew = [a[3], a[2], a[1]] if flip_face else [a[3], a[1], a[2]]
            elif ids[0] == 4:
                anew = [a[3], a[2], a[1]] if flip_face else [a[1], a[2], a[3]]
                
            flip_face = not flip_face

            if ids[1] == 1:
                bnew = [b[1], b[3], b[2]] if flip_face else [b[1], b[2], b[3]]
            elif ids[1] == 2:
                bnew = [b[2], b[1], b[3]] if flip_face else [b[2], b[3], b[1]]
            elif ids[1] == 3:
                bnew = [b[3], b[2], b[1]] if flip_face else [b[3], b[1], b[2]]
            elif ids[1] == 4:
                bnew = [b[3], b[2], b[1]] if flip_face else [b[1], b[2], b[3]]
            
        return [anew, bnew]

    inner_faces = zeros((trailing_edge.shape[0], 5), dtype=int32)
    inner_vertices = []

    ids = []

    for i in range(trailing_edge.shape[0]):
            
        # Upper and lower faces
        face1 = faces[trailing_edge_faces[i, 0], :]
        face2 = faces[trailing_edge_faces[i, 1], :]

        # Correct order
        ids1, ids2 = __find_ids(face1, face2, i)

        # Add ids and vertices
        for j in range(face1[0]):
            if ids1[j] not in ids:
                ids.append(ids1[j])
                inner_vertices.append(0.5 * (vertices[ids1[j], :] + vertices[ids2[j], :]))
            
        # Add face
        inner_faces[i, 0] = face1[0]
        inner_faces[i, 1] = ids.index(ids1[0])
        inner_faces[i, 2] = ids.index(ids1[1])
        inner_faces[i, 3] = ids.index(ids1[2])
        if face1[0] == 4: inner_faces[i, 4] = ids.index(ids1[3])

    inner_vertices = asarray(inner_vertices)

    # Output
    out = { 'inner_faces': inner_faces, 'inner_vertices': inner_vertices }

    return out

from typing import List
from abc import ABC, abstractmethod
from numpy import ndarray, empty, double
from numpy.linalg import solve

from pypm3D.models.mesh_model import MeshModel

#-----------------------------------------------------------------------------#
#-----------------------------------------------------------------------------#
class MeshAbs(ABC):
    """
    Generate the mesh parameters
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
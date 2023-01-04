from abc import ABC, abstractmethod
from typing import List
from numpy import int32, ndarray, empty, double, cross, sqrt, fabs, dot, argwhere, asarray, zeros
from numpy.linalg import norm

from pypm3D.models.mesh_model import MeshModel

def __get_trailing_edge_faces(mesh: MeshModel) -> None:

    for i in range(mesh.nte):
        check1 = (mesh.faces[:, 1] == mesh.trailing_edge[i, 0]) | (mesh.faces[:, 2] == mesh.trailing_edge[i, 0]) | (mesh.faces[:, 3] == mesh.trailing_edge[i, 0]) | (mesh.faces[:, 4] == mesh.trailing_edge[i, 0])
        check2 = (mesh.faces[:, 1] == mesh.trailing_edge[i, 1]) | (mesh.faces[:, 2] == mesh.trailing_edge[i, 1]) | (mesh.faces[:, 3] == mesh.trailing_edge[i, 1]) | (mesh.faces[:, 4] == mesh.trailing_edge[i, 1])
        index = argwhere(check1 & check2)
        mesh.te_faces[i, 0] = index[0][0]
        mesh.te_faces[i, 1] = index[1][0]
    
    return

def __calculate_faces_parameters(mesh: MeshModel) -> None:

    # Panel type
    tri_panel_ids = mesh.faces[:, 0] == 3
    quad_panel_ids = mesh.faces[:, 0] == 4

    # Panel center
    mesh.p_avg[tri_panel_ids] = (1. / 3.) * (mesh.vertices[mesh.faces[tri_panel_ids, 1], :] + mesh.vertices[mesh.faces[tri_panel_ids, 2], :] + mesh.vertices[mesh.faces[tri_panel_ids, 3], :])
    mesh.p_avg[quad_panel_ids] = (1. / 4.) * (mesh.vertices[mesh.faces[quad_panel_ids, 1], :] + mesh.vertices[mesh.faces[quad_panel_ids, 2], :] + mesh.vertices[mesh.faces[quad_panel_ids, 3], :] + mesh.vertices[mesh.faces[quad_panel_ids, 4], :])

    # Orthogonal base (e1)
    te_faces_ids = asarray([((face[1] in mesh.te_faces[:, 0]) & (face[2] in mesh.te_faces[:, 0]) & (face[3] in mesh.te_faces[:, 0]) & (face[4] in mesh.te_faces[:, 0])) & ((face[1] in mesh.te_faces[:, 1]) & (face[2] in mesh.te_faces[:, 1]) & (face[3] in mesh.te_faces[:, 1]) & (face[4] in mesh.te_faces[:, 1])) for face in mesh.faces])
    
    mesh.e1[te_faces_ids, :] = mesh.vertices[mesh.trailing_edge[te_faces_ids, 1], :] - mesh.vertices[mesh.trailing_edge[te_faces_ids, 0], :]
    mesh.e1[not te_faces_ids, :] = 0.5 * (mesh.vertices[mesh.faces[not te_faces_ids, 1], :] + mesh.vertices[mesh.faces[not te_faces_ids, 2], :]) - mesh.p_avg[not te_faces_ids, :]

    aux_norm = sqrt(mesh.e1[:, 0] * mesh.e1[:, 0] + mesh.e1[:, 1] * mesh.e1[:, 1] + mesh.e1[:, 2] * mesh.e1[:, 2])

    mesh.e1[:, 0] = mesh.e1[:, 0] / aux_norm
    mesh.e1[:, 1] = mesh.e1[:, 1] / aux_norm
    mesh.e1[:, 2] = mesh.e1[:, 2] / aux_norm

    # Orthogonal base (e3)
    mesh.e3[tri_panel_ids, :] = cross(mesh.vertices[mesh.faces[tri_panel_ids, 2], :] - mesh.vertices[mesh.faces[tri_panel_ids, 1], :], mesh.vertices[mesh.faces[tri_panel_ids, 3], :] - mesh.vertices[mesh.faces[tri_panel_ids, 1], :])
    mesh.e3[quad_panel_ids, :] = cross(mesh.vertices[mesh.faces[quad_panel_ids, 2], :] - mesh.vertices[mesh.faces[quad_panel_ids, 4], :], mesh.vertices[mesh.faces[quad_panel_ids, 3], :] - mesh.vertices[mesh.faces[quad_panel_ids, 1], :])

    aux_norm = sqrt(mesh.e3[:, 0] * mesh.e3[:, 0] + mesh.e3[:, 1] * mesh.e3[:, 1] + mesh.e3[:, 2] * mesh.e3[:, 2])

    mesh.e3[:, 0] = mesh.e3[:, 0] / aux_norm
    mesh.e3[:, 1] = mesh.e3[:, 1] / aux_norm
    mesh.e3[:, 2] = mesh.e3[:, 2] / aux_norm

    # Orthogonal base (e2)
    mesh.e2[:, :] = cross(mesh.e3, mesh.e1)

    aux_norm = sqrt(mesh.e2[:, 0] * mesh.e2[:, 0] + mesh.e2[:, 1] * mesh.e2[:, 1] + mesh.e2[:, 2] * mesh.e2[:, 2])

    mesh.e2[:, 0] = mesh.e2[:, 0] / aux_norm
    mesh.e2[:, 1] = mesh.e2[:, 1] / aux_norm
    mesh.e2[:, 2] = mesh.e2[:, 2] / aux_norm

    # Orthogonal base (e3 correction)
    mesh.e3[:, :] = cross(mesh.e1, mesh.e2)

    aux_norm = sqrt(mesh.e3[:, 0] * mesh.e3[:, 0] + mesh.e3[:, 1] * mesh.e3[:, 1] + mesh.e3[:, 2] * mesh.e3[:, 2])

    mesh.e3[:, 0] = mesh.e3[:, 0] / aux_norm
    mesh.e3[:, 1] = mesh.e3[:, 1] / aux_norm
    mesh.e3[:, 2] = mesh.e3[:, 2] / aux_norm

    # Control point
    mesh.p_ctrl[:, :] = mesh.p_avg + mesh.e3 * 1e-12

    # Local points
    mesh.p1[:, 0] = dot(mesh.vertices[mesh.faces[:, 1], :] - mesh.p_avg, mesh.e1.T)
    mesh.p1[:, 1] = dot(mesh.vertices[mesh.faces[:, 1], :] - mesh.p_avg, mesh.e2.T)

    mesh.p2[:, 0] = dot(mesh.vertices[mesh.faces[:, 2], :] - mesh.p_avg, mesh.e1.T)
    mesh.p2[:, 1] = dot(mesh.vertices[mesh.faces[:, 2], :] - mesh.p_avg, mesh.e2.T)

    mesh.p3[:, 0] = dot(mesh.vertices[mesh.faces[:, 3], :] - mesh.p_avg, mesh.e1.T)
    mesh.p3[:, 1] = dot(mesh.vertices[mesh.faces[:, 3], :] - mesh.p_avg, mesh.e2.T)

    mesh.p4[quad_panel_ids, 0] = dot(mesh.vertices[mesh.faces[quad_panel_ids, 3], :] - mesh.p_avg[quad_panel_ids, :], mesh.e1[quad_panel_ids, :].T)
    mesh.p4[quad_panel_ids, 1] = dot(mesh.vertices[mesh.faces[quad_panel_ids, 3], :] - mesh.p_avg[quad_panel_ids, :], mesh.e2[quad_panel_ids, :].T)

    return

def __get_inner_doublet_pannels(mesh: MeshModel) -> List[ndarray]:

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
        flip_face = True if mesh.e3[mesh.te_faces[face, 0], 2] < 0 else False

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

    inner_faces = zeros((mesh.nte, 5), dtype=int32)
    inner_vertices = []

    ids = []

    for i in range(mesh.nte):
            
        # Upper and lower faces
        face1 = mesh.faces[mesh.te_faces[i, 0], :]
        face2 = mesh.faces[mesh.te_faces[i, 1], :]

        # Correct order
        ids1, ids2 = __find_ids(face1, face2, i)

        # Add ids and vertices
        for j in range(face1[0]):
            if ids1[j] not in ids:
                ids.append(ids1[j])
                inner_vertices.append(0.5 * (mesh.vertices[ids1[j], :] + mesh.vertices[ids2[j], :]))
            
        # Add face
        inner_faces[i, 0] = face1[0]
        inner_faces[i, 1] = ids.index(ids1[0])
        inner_faces[i, 2] = ids.index(ids1[1])
        inner_faces[i, 3] = ids.index(ids1[2])
        if face1[0] == 4: inner_faces[i, 4] = ids.index(ids1[3])

    inner_vertices = asarray(inner_vertices)

    return [inner_faces, inner_vertices]

def create_surface_mesh(vertices: ndarray,
                        faces: ndarray,
                        trailing_edge: ndarray,):
    
    """
    It calculates the surface mesh parameters and creates the inner
    wake panels.
    """

    # Size
    nv = vertices.shape[0]
    nf = faces.shape[0]
    nte = trailing_edge.shape[0]
    
    # Output
    mesh = MeshModel(nv=nv, nf=nf, nte=nte, vertices=vertices, faces=faces, trailing_edge=trailing_edge, p_avg=empty((nf, 3), dtype=double), p_ctrl=empty((nf, 3), dtype=double), e1=empty((nf, 3), dtype=double), e2=empty((nf, 3), dtype=double), e3=empty((nf, 3), dtype=double), p1=empty((nf, 2), dtype=double), p2=empty((nf, 2), dtype=double), p3=empty((nf, 2), dtype=double), p4=empty((nf, 2), dtype=double), te_faces=empty((nte, 2), dtype=double))

    # Correct last tri panel id
    mesh.faces[mesh.faces[:, 0] == 3, 4] = -1
    
    # Trailing edge faces
    __get_trailing_edge_faces(mesh)

    # Surface parameters
    __calculate_faces_parameters(mesh)

    return

def process_mesh(vertices: ndarray,
                 faces: ndarray,
                 trailing_edge: ndarray,
                 freestream: float,
                 wake_length: float,
                 wake_time_step: float) -> MeshModel:
    """
    It calculates the surface mesh parameters and creates the inner
    wake panels.
    """

    # Size
    nv = vertices.shape[0]
    nf = faces.shape[0]
    nte = trailing_edge.shape[0]
    nw = int(wake_length / (freestream / wake_time_step))
    nwv = 0

    te_faces_ids = []

    for i in range(nte):
        if trailing_edge[i, 0] not in te_faces_ids: te_faces_ids.append(trailing_edge[i, 0])
        if trailing_edge[i, 1] not in te_faces_ids: te_faces_ids.append(trailing_edge[i, 1])
        
    nwv = len(te_faces_ids)

    # Output
    mesh = MeshModel(nv=nv, nf=nf, nte=nte, nw=nw, nwv=nwv, vertices=vertices, faces=faces, trailing_edge=trailing_edge, p_avg=empty((nf, 3), dtype=double), p_ctrl=empty((nf, 3), dtype=double), e1=empty((nf, 3), dtype=double), e2=empty((nf, 3), dtype=double), e3=empty((nf, 3), dtype=double), p1=empty((nf, 2), dtype=double), p2=empty((nf, 2), dtype=double), p3=empty((nf, 2), dtype=double), p4=empty((nf, 2), dtype=double), te_faces=empty((nte, 2), dtype=double))

    # Correct last tri panel id
    for i in range(mesh.nf):
        if mesh.faces[i, 0] == 3:
            mesh.faces[i, 4] = -1

    # Trailing edge faces
    __get_trailing_edge_faces(mesh)

    # Surface parameters
    __calculate_faces_parameters(mesh)
    
    return mesh
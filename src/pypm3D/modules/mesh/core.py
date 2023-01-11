import numpy as np

from pypm3D import models


def find_faces_at_trailing_edge(mesh: models.MeshModel) -> None:

    for i in range(mesh.surface.nte):
        check1 = (mesh.surface.faces[:, 1] == mesh.surface.trailing_edge[i, 0]) | (mesh.surface.faces[:, 2] == mesh.surface.trailing_edge[i, 0]) | (mesh.surface.faces[:, 3] == mesh.surface.trailing_edge[i, 0]) | (mesh.surface.faces[:, 4] == mesh.surface.trailing_edge[i, 0])
        check2 = (mesh.surface.faces[:, 1] == mesh.surface.trailing_edge[i, 1]) | (mesh.surface.faces[:, 2] == mesh.surface.trailing_edge[i, 1]) | (mesh.surface.faces[:, 3] == mesh.surface.trailing_edge[i, 1]) | (mesh.surface.faces[:, 4] == mesh.surface.trailing_edge[i, 1])
        index = np.argwhere(check1 & check2)
        mesh.surface.trailing_edge_faces[i, 0] = index[0][0]
        mesh.surface.trailing_edge_faces[i, 1] = index[1][0]

    return

def calculate_panels_parameters(mesh: models.MeshModel) -> None:

    # Panel geometry
    tri_panel_ids = mesh.surface.faces[:, 0] == 3
    quad_panel_ids = mesh.surface.faces[:, 0] == 4

    # Panel center
    mesh.surface.p_avg[tri_panel_ids, :] = (1. / 3.) * (mesh.surface.vertices[mesh.surface.faces[tri_panel_ids, 1], :] + mesh.surface.vertices[mesh.surface.faces[tri_panel_ids, 2], :] + mesh.surface.vertices[mesh.surface.faces[tri_panel_ids, 3], :])
    mesh.surface.p_avg[quad_panel_ids, :] = (1. / 4.) * (mesh.surface.vertices[mesh.surface.faces[quad_panel_ids, 1], :] + mesh.surface.vertices[mesh.surface.faces[quad_panel_ids, 2], :] + mesh.surface.vertices[mesh.surface.faces[quad_panel_ids, 3], :] + mesh.surface.vertices[mesh.surface.faces[quad_panel_ids, 4], :])
    
    # Orthogonal base (e3)
    mesh.surface.e3[tri_panel_ids, :] = np.cross(mesh.surface.vertices[mesh.surface.faces[tri_panel_ids, 2], :] - mesh.surface.vertices[mesh.surface.faces[tri_panel_ids, 1], :], mesh.surface.vertices[mesh.surface.faces[tri_panel_ids, 3], :] - mesh.surface.vertices[mesh.surface.faces[tri_panel_ids, 1], :])
    mesh.surface.e3[quad_panel_ids, :] = np.cross(mesh.surface.vertices[mesh.surface.faces[quad_panel_ids, 2], :] - mesh.surface.vertices[mesh.surface.faces[quad_panel_ids, 4], :], mesh.surface.vertices[mesh.surface.faces[quad_panel_ids, 3], :] - mesh.surface.vertices[mesh.surface.faces[quad_panel_ids, 1], :])
    
    # Orthogonal base (e1)
    mesh.surface.e1[:, :] = 0.5 * (mesh.surface.vertices[mesh.surface.faces[:, 1], :] + mesh.surface.vertices[mesh.surface.faces[:, 2], :]) - mesh.surface.p_avg

    # Orthogonal base (e2)
    mesh.surface.e2[:, :] = np.cross(mesh.surface.e3, mesh.surface.e1)

    # Orthogonal base (e3 - correção)
    mesh.surface.e3[:, :] = np.cross(mesh.surface.e1, mesh.surface.e2)

    # Correct trailing edge faces
    mesh.surface.e1[mesh.surface.trailing_edge_faces[:, 0], :] = mesh.surface.vertices[mesh.surface.trailing_edge[:, 0], :] - mesh.surface.vertices[mesh.surface.trailing_edge[:, 1], :]
    mesh.surface.e1[mesh.surface.trailing_edge_faces[:, 1], :] = mesh.surface.e1[mesh.surface.trailing_edge_faces[:, 0], :]

    # Upper panel
    mesh.surface.e2[mesh.surface.trailing_edge_faces[:, 0], :] = np.cross(mesh.surface.e3[mesh.surface.trailing_edge_faces[:, 0], :], mesh.surface.e1[mesh.surface.trailing_edge_faces[:, 0], :])
    mesh.surface.e3[mesh.surface.trailing_edge_faces[:, 0], :] = np.cross(mesh.surface.e1[mesh.surface.trailing_edge_faces[:, 0], :], mesh.surface.e2[mesh.surface.trailing_edge_faces[:, 0], :])

    # Lower panel
    mesh.surface.e2[mesh.surface.trailing_edge_faces[:, 1], :] = np.cross(mesh.surface.e3[mesh.surface.trailing_edge_faces[:, 1], :], mesh.surface.e1[mesh.surface.trailing_edge_faces[:, 1], :])
    mesh.surface.e3[mesh.surface.trailing_edge_faces[:, 1], :] = np.cross(mesh.surface.e1[mesh.surface.trailing_edge_faces[:, 1], :], mesh.surface.e2[mesh.surface.trailing_edge_faces[:, 1], :])

    # Normalize orthogonal base
    aux_norm = np.sqrt(mesh.surface.e1[:, 0] * mesh.surface.e1[:, 0] + mesh.surface.e1[:, 1] * mesh.surface.e1[:, 1] + mesh.surface.e1[:, 2] * mesh.surface.e1[:, 2])

    mesh.surface.e1[:, 0] = mesh.surface.e1[:, 0] / aux_norm
    mesh.surface.e1[:, 1] = mesh.surface.e1[:, 1] / aux_norm
    mesh.surface.e1[:, 2] = mesh.surface.e1[:, 2] / aux_norm

    aux_norm = np.sqrt(mesh.surface.e2[:, 0] * mesh.surface.e2[:, 0] + mesh.surface.e2[:, 1] * mesh.surface.e2[:, 1] + mesh.surface.e2[:, 2] * mesh.surface.e2[:, 2])

    mesh.surface.e2[:, 0] = mesh.surface.e2[:, 0] / aux_norm
    mesh.surface.e2[:, 1] = mesh.surface.e2[:, 1] / aux_norm
    mesh.surface.e2[:, 2] = mesh.surface.e2[:, 2] / aux_norm

    aux_norm = np.sqrt(mesh.surface.e3[:, 0] * mesh.surface.e3[:, 0] + mesh.surface.e3[:, 1] * mesh.surface.e3[:, 1] + mesh.surface.e3[:, 2] * mesh.surface.e3[:, 2])

    mesh.surface.e3[:, 0] = mesh.surface.e3[:, 0] / aux_norm
    mesh.surface.e3[:, 1] = mesh.surface.e3[:, 1] / aux_norm
    mesh.surface.e3[:, 2] = mesh.surface.e3[:, 2] / aux_norm

    # Control point
    mesh.surface.p_ctrl[:, :] = mesh.surface.p_avg[:, :] + mesh.surface.e3[:, :] * 1e-8

    # Local vertices
    aux_vec = mesh.surface.vertices[mesh.surface.faces[:, 1], :] - mesh.surface.p_avg
    mesh.surface.p1[:, 0] = mesh.surface.e1[:, 0] * aux_vec[:, 0] + mesh.surface.e1[:, 1] * aux_vec[:, 1] + mesh.surface.e1[:, 2] * aux_vec[:, 2]
    mesh.surface.p1[:, 1] = mesh.surface.e2[:, 0] * aux_vec[:, 0] + mesh.surface.e2[:, 1] * aux_vec[:, 1] + mesh.surface.e2[:, 2] * aux_vec[:, 2]

    aux_vec = mesh.surface.vertices[mesh.surface.faces[:, 2], :] - mesh.surface.p_avg
    mesh.surface.p2[:, 0] = mesh.surface.e1[:, 0] * aux_vec[:, 0] + mesh.surface.e1[:, 1] * aux_vec[:, 1] + mesh.surface.e1[:, 2] * aux_vec[:, 2]
    mesh.surface.p2[:, 1] = mesh.surface.e2[:, 0] * aux_vec[:, 0] + mesh.surface.e2[:, 1] * aux_vec[:, 1] + mesh.surface.e2[:, 2] * aux_vec[:, 2]

    aux_vec = mesh.surface.vertices[mesh.surface.faces[:, 3], :] - mesh.surface.p_avg
    mesh.surface.p3[:, 0] = mesh.surface.e1[:, 0] * aux_vec[:, 0] + mesh.surface.e1[:, 1] * aux_vec[:, 1] + mesh.surface.e1[:, 2] * aux_vec[:, 2]
    mesh.surface.p3[:, 1] = mesh.surface.e2[:, 0] * aux_vec[:, 0] + mesh.surface.e2[:, 1] * aux_vec[:, 1] + mesh.surface.e2[:, 2] * aux_vec[:, 2]

    aux_vec = mesh.surface.vertices[mesh.surface.faces[quad_panel_ids, 4], :] - mesh.surface.p_avg[quad_panel_ids, :]
    mesh.surface.p4[quad_panel_ids, 0] = mesh.surface.e1[quad_panel_ids, 0] * aux_vec[:, 0] + mesh.surface.e1[quad_panel_ids, 1] * aux_vec[:, 1] + mesh.surface.e1[quad_panel_ids, 2] * aux_vec[:, 2]
    mesh.surface.p4[quad_panel_ids, 1] = mesh.surface.e2[quad_panel_ids, 0] * aux_vec[:, 0] + mesh.surface.e2[quad_panel_ids, 1] * aux_vec[:, 1] + mesh.surface.e2[quad_panel_ids, 2] * aux_vec[:, 2]

    return

def create_internal_panels(mesh: models.MeshModel) -> None:

    def __find_internal_ids(a: np.ndarray, b: np.ndarray, face: int) -> int:

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
        flip_face = True if mesh.surface.e3[mesh.surface.trailing_edge_faces[face, 0], 2] < 0 else False

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

    inner_faces = np.zeros((mesh.surface.nte, 5), dtype=np.int32)
    inner_vertices = []

    ids = []

    for i in range(mesh.surface.nte):
            
        # Upper and lower faces
        face1 = mesh.surface.faces[mesh.surface.trailing_edge_faces[i, 0], :]
        face2 = mesh.surface.faces[mesh.surface.trailing_edge_faces[i, 1], :]

        # Correct order
        ids1, ids2 = __find_internal_ids(face1, face2, i)

        # Add ids and vertices
        for j in range(face1[0]):
            if ids1[j] not in ids:
                ids.append(ids1[j])
                inner_vertices.append(0.5 * (mesh.surface.vertices[ids1[j], :] + mesh.surface.vertices[ids2[j], :]))
            
        # Add face
        inner_faces[i, 0] = face1[0]
        inner_faces[i, 1] = ids.index(ids1[0])
        inner_faces[i, 2] = ids.index(ids1[1])
        inner_faces[i, 3] = ids.index(ids1[2])
        if face1[0] == 4: inner_faces[i, 4] = ids.index(ids1[3])

    inner_vertices = np.asarray(inner_vertices)

    mesh.internal.nf = inner_faces.shape[0]
    mesh.internal.nv = inner_vertices.shape[0]

    mesh.internal.faces = np.empty((mesh.surface.nte, 5), dtype=np.int32)
    mesh.internal.vertices = np.empty((mesh.internal.nv, 3), dtype=np.double)
    
    mesh.internal.faces[:, :] = inner_faces.astype(np.int32)[:, :]
    mesh.internal.vertices[:, :] = inner_vertices.astype(np.double)[:, :]

    return

def calculate_scale_factor(mesh: models.MeshModel) -> None:

    def __point2line(a, b, p):
        return np.linalg.norm(np.cross(p - a, b - a)) / np.linalg.norm(b - a)
    
    mesh.internal.scale = np.empty((mesh.surface.nte,), dtype=np.double)

    for i in range(mesh.surface.nte):
        if mesh.internal.faces[i, 0] == 3:
            v1 = mesh.internal.vertices[mesh.internal.faces[i, 1], :]
            v2 = mesh.internal.vertices[mesh.internal.faces[i, 2], :]
            v3 = mesh.internal.vertices[mesh.internal.faces[i, 3], :]
            p_avg = (1. / 3.) * (v1 + v2 + v3)
            mesh.internal.scale[i] = min([__point2line(v1, v2, p_avg), __point2line(v2, v3, p_avg), __point2line(v3, v1, p_avg)])
        else:
            v1 = mesh.internal.vertices[mesh.internal.faces[i, 1], :]
            v2 = mesh.internal.vertices[mesh.internal.faces[i, 2], :]
            v3 = mesh.internal.vertices[mesh.internal.faces[i, 3], :]
            v4 = mesh.internal.vertices[mesh.internal.faces[i, 4], :]
            p_avg = (1. / 4.) * (v1 + v2 + v3 + v4)
            mesh.internal.scale[i] = min([__point2line(v1, v2, p_avg), __point2line(v2, v3, p_avg), __point2line(v3, v4, p_avg), __point2line(v4, v1, p_avg)])

    return
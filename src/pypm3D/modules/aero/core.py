import typing as tp
import numpy as np
import ctypes as ct
import math

from pypm3D import models


def create_wake_mesh(mesh: models.MeshModel,
                     u_func: tp.Callable[[np.ndarray], np.ndarray],
                     length: float,
                     time_step: float,
                     aero: models.AeroModel) -> None:
    
    # Approximate the mean velocity
    u_mean = np.linalg.norm(u_func(np.mean(mesh.surface.vertices, axis=0)))

    # Find the unique trailing edge ids
    trailing_edge_ids = []

    for edge in mesh.surface.trailing_edge:
        if edge[0] not in trailing_edge_ids: trailing_edge_ids.append(edge[0])
        if edge[1] not in trailing_edge_ids: trailing_edge_ids.append(edge[1])

    # Save and create parameters
    aero.wake.nw = math.ceil(length / (time_step * u_mean))
    aero.wake.nv = len(trailing_edge_ids)

    aero.wake.trailing_edge_ids = np.asarray(trailing_edge_ids, dtype=np.int32)
    aero.wake.vertices = np.empty((aero.wake.nv * (aero.wake.nw + 1), 3), dtype=np.double)
    aero.wake.faces = np.empty((mesh.surface.nte * aero.wake.nw, 5), dtype=np.int32)

    # Save trailing edge vertices
    for i in range(aero.wake.nv):
        aero.wake.vertices[i, :] = mesh.surface.vertices[aero.wake.trailing_edge_ids[i], :]
    
    # Erase here
    u = np.zeros((aero.wake.nv, 3), dtype=np.double)
    u[:, 0] = - 0.1
    for i in range(aero.wake.nw):
        u = u_func(aero.wake.vertices[i * aero.wake.nv:(i + 1) * aero.wake.nv, :])
        aero.wake.vertices[(i + 1) * aero.wake.nv:(i + 2) * aero.wake.nv, :] = aero.wake.vertices[i * aero.wake.nv:(i + 1) * aero.wake.nv, :] + u * time_step

    # Create faces ids
    for i in range(mesh.surface.nte):

        id1 = trailing_edge_ids.index(mesh.surface.trailing_edge[i, 0])
        id2 = trailing_edge_ids.index(mesh.surface.trailing_edge[i, 1])

        for j in range(aero.wake.nw):
            aero.wake.faces[i + j * mesh.surface.nte, 0] = 4
            aero.wake.faces[i + j * mesh.surface.nte, 1] = id1 + j * aero.wake.nv
            aero.wake.faces[i + j * mesh.surface.nte, 2] = id2 + j * aero.wake.nv
            aero.wake.faces[i + j * mesh.surface.nte, 3] = id2 + (j + 1) * aero.wake.nv
            aero.wake.faces[i + j * mesh.surface.nte, 4] = id1 + (j + 1) * aero.wake.nv
    
    # Create wake parameters
    aero.wake.circulations = np.empty((mesh.surface.nte, aero.wake.nw), dtype=np.double)
    aero.wake.areas = np.empty((mesh.surface.nte, aero.wake.nw), dtype=np.double)

    return

def create_vertices_connection(mesh: models.MeshModel) -> tp.List[models.VerticeConnectionModel]:

    vertices_connection = []

    for v_id in range(mesh.surface.nv):

        # Encontra as faces que se conectam com o v_id
        face_ids = [i[0] for i in np.argwhere((mesh.surface.faces[:, 1] == v_id) | (mesh.surface.faces[:, 2] == v_id) | (mesh.surface.faces[:, 3] == v_id) | (mesh.surface.faces[:, 4] == v_id))]

        # Encontra os ângulos
        angles = []

        for face_id in face_ids:
            
            if v_id == mesh.surface.faces[face_id, 1]:
                
                v1 = mesh.surface.vertices[mesh.surface.faces[face_id, 2], :] - mesh.surface.vertices[mesh.surface.faces[face_id, 1], :]
                
                if mesh.surface.faces[face_id, 0] == 3:
                    v2 = mesh.surface.vertices[mesh.surface.faces[face_id, 3], :] - mesh.surface.vertices[mesh.surface.faces[face_id, 1], :]
                else:
                    v2 = mesh.surface.vertices[mesh.surface.faces[face_id, 4], :] - mesh.surface.vertices[mesh.surface.faces[face_id, 1], :]

            elif v_id == mesh.surface.faces[face_id, 2]:
                v1 = mesh.surface.vertices[mesh.surface.faces[face_id, 3], :] - mesh.surface.vertices[mesh.surface.faces[face_id, 2], :]
                v2 = mesh.surface.vertices[mesh.surface.faces[face_id, 1], :] - mesh.surface.vertices[mesh.surface.faces[face_id, 2], :]

            elif v_id == mesh.surface.faces[face_id, 3]:

                v2 = mesh.surface.vertices[mesh.surface.faces[face_id, 2], :] - mesh.surface.vertices[mesh.surface.faces[face_id, 3], :]

                if mesh.surface.faces[face_id, 0] == 3:
                    v1 = mesh.surface.vertices[mesh.surface.faces[face_id, 1], :] - mesh.surface.vertices[mesh.surface.faces[face_id, 3], :]
                else:
                    v1 = mesh.surface.vertices[mesh.surface.faces[face_id, 4], :] - mesh.surface.vertices[mesh.surface.faces[face_id, 3], :]

            else:
                v1 = mesh.surface.vertices[mesh.surface.faces[face_id, 1], :] - mesh.surface.vertices[mesh.surface.faces[face_id, 4], :]
                v2 = mesh.surface.vertices[mesh.surface.faces[face_id, 3], :] - mesh.surface.vertices[mesh.surface.faces[face_id, 4], :]
            
            angles.append(np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))))

        sum_angles = sum(angles)

        # Adiciona um modelo
        vertices_connection.append(
            models.VerticeConnectionModel(
                n=len(face_ids),
                faces=face_ids,
                coefs=[a / sum_angles for a in angles],
            )
        )
    
    return vertices_connection

def set_a_ij(mesh: models.MeshModel, a_ij: np.ndarray) -> None:

    # a_ij
    a_ij_x = np.empty((mesh.surface.nf * mesh.surface.nf,), dtype=np.double)
    a_ij_y = np.empty((mesh.surface.nf * mesh.surface.nf,), dtype=np.double)
    a_ij_z = np.empty((mesh.surface.nf * mesh.surface.nf,), dtype=np.double)

    n_sides = np.empty((mesh.surface.nf,), dtype=np.int32)

    n_sides[:] = mesh.surface.faces[:, 0]

    lib = ct.CDLL('./src/pypm3D/modules/aero/bin/lib.so')

    lib.get_a_ij_coefs.argtypes = [
        ct.c_int,                         # nf
        ct.POINTER(ct.c_int),             # n_sides
        ct.POINTER(ct.c_double),          # p_avg
        ct.POINTER(ct.c_double),          # p_ctrl
        ct.POINTER(ct.c_double),          # e1
        ct.POINTER(ct.c_double),          # e2
        ct.POINTER(ct.c_double),          # e3
        ct.POINTER(ct.c_double),          # p1
        ct.POINTER(ct.c_double),          # p2
        ct.POINTER(ct.c_double),          # p3
        ct.POINTER(ct.c_double),          # p4
        ct.POINTER(ct.c_double),          # a_ij_x
        ct.POINTER(ct.c_double),          # a_ij_y
        ct.POINTER(ct.c_double),          # a_ij_z
    ]

    lib.get_a_ij_coefs.restype = None

    lib.get_a_ij_coefs(
        mesh.surface.nf,
        np.ctypeslib.as_ctypes(n_sides),
        np.ctypeslib.as_ctypes(mesh.surface.p_avg.reshape(mesh.surface.p_avg.size)),
        np.ctypeslib.as_ctypes(mesh.surface.p_ctrl.reshape(mesh.surface.p_ctrl.size)),
        np.ctypeslib.as_ctypes(mesh.surface.e1.reshape(mesh.surface.e1.size)),
        np.ctypeslib.as_ctypes(mesh.surface.e2.reshape(mesh.surface.e2.size)),
        np.ctypeslib.as_ctypes(mesh.surface.e3.reshape(mesh.surface.e3.size)),
        np.ctypeslib.as_ctypes(mesh.surface.p1.reshape(mesh.surface.p1.size)),
        np.ctypeslib.as_ctypes(mesh.surface.p2.reshape(mesh.surface.p2.size)),
        np.ctypeslib.as_ctypes(mesh.surface.p3.reshape(mesh.surface.p3.size)),
        np.ctypeslib.as_ctypes(mesh.surface.p4.reshape(mesh.surface.p4.size)),
        np.ctypeslib.as_ctypes(a_ij_x),
        np.ctypeslib.as_ctypes(a_ij_y),
        np.ctypeslib.as_ctypes(a_ij_z)
    )

    a_ij_x = a_ij_x.reshape((mesh.surface.nf, mesh.surface.nf))
    a_ij_y = a_ij_y.reshape((mesh.surface.nf, mesh.surface.nf))
    a_ij_z = a_ij_z.reshape((mesh.surface.nf, mesh.surface.nf))

    a_ij[:, :, 0] = a_ij_x[:, :]
    a_ij[:, :, 1] = a_ij_y[:, :]
    a_ij[:, :, 2] = a_ij_z[:, :]

    return

def set_b_kj(mesh: models.MeshModel, b_kj: np.ndarray) -> None:
    
    # b_kj
    b_kj_x = np.empty((mesh.surface.nf * mesh.surface.nte,), dtype=np.double)
    b_kj_y = np.empty((mesh.surface.nf * mesh.surface.nte,), dtype=np.double)
    b_kj_z = np.empty((mesh.surface.nf * mesh.surface.nte,), dtype=np.double)

    lib = ct.CDLL('./src/pypm3D/modules/aero/bin/lib.so')

    lib.get_b_kj_coefs.argtypes = [
        ct.c_int,                         # nf
        ct.c_int,                         # nte
        ct.POINTER(ct.c_int),             # faces
        ct.POINTER(ct.c_double),          # vertices
        ct.POINTER(ct.c_double),          # scale
        ct.POINTER(ct.c_double),          # p_ctrl
        ct.POINTER(ct.c_double),          # b_kj_x
        ct.POINTER(ct.c_double),          # b_kj_y
        ct.POINTER(ct.c_double),          # b_kj_z
    ]

    lib.get_b_kj_coefs.restype = None
    
    lib.get_b_kj_coefs(
        mesh.surface.nf,
        mesh.surface.nte,
        np.ctypeslib.as_ctypes(mesh.internal.faces.reshape(mesh.internal.faces.size)),
        np.ctypeslib.as_ctypes(mesh.internal.vertices.reshape(mesh.internal.vertices.size)),
        np.ctypeslib.as_ctypes(mesh.internal.scale.reshape(mesh.internal.scale.size)),
        np.ctypeslib.as_ctypes(mesh.surface.p_ctrl.reshape(mesh.surface.p_ctrl.size)),
        np.ctypeslib.as_ctypes(b_kj_x),
        np.ctypeslib.as_ctypes(b_kj_y),
        np.ctypeslib.as_ctypes(b_kj_z)
    )

    b_kj_x = b_kj_x.reshape((mesh.surface.nf, mesh.surface.nte))
    b_kj_y = b_kj_y.reshape((mesh.surface.nf, mesh.surface.nte))
    b_kj_z = b_kj_z.reshape((mesh.surface.nf, mesh.surface.nte))

    b_kj[:, :, 0] = b_kj_x[:, :]
    b_kj[:, :, 1] = b_kj_y[:, :]
    b_kj[:, :, 2] = b_kj_z[:, :]

    return

def set_c_kj(mesh: models.MeshModel, c_kj: np.ndarray) -> None:
    
    # a_ij
    c_kj_x = np.empty((mesh.surface.nf * mesh.surface.nte,), dtype=np.double)
    c_kj_y = np.empty((mesh.surface.nf * mesh.surface.nte,), dtype=np.double)
    c_kj_z = np.empty((mesh.surface.nf * mesh.surface.nte,), dtype=np.double)

    lib = ct.CDLL('./src/pypm3D/modules/aero/bin/lib.so')

    lib.get_c_kj_vortex_line_coefs.argtypes = [
        ct.c_int,                         # nf
        ct.c_int,                         # nte
        ct.POINTER(ct.c_int),             # te
        ct.POINTER(ct.c_double),          # vertices
        ct.POINTER(ct.c_double),          # scale
        ct.POINTER(ct.c_double),          # p_ctrl
        ct.POINTER(ct.c_double),          # c_kj_x
        ct.POINTER(ct.c_double),          # c_kj_y
        ct.POINTER(ct.c_double),          # c_kj_z
    ]

    lib.get_c_kj_vortex_line_coefs.restype = None

    lib.get_c_kj_vortex_line_coefs(
        mesh.surface.nf,
        mesh.surface.nte,
        np.ctypeslib.as_ctypes(mesh.surface.trailing_edge.reshape(mesh.surface.trailing_edge.size)),
        np.ctypeslib.as_ctypes(mesh.surface.vertices.reshape(mesh.surface.vertices.size)),
        np.ctypeslib.as_ctypes(mesh.internal.scale.reshape(mesh.internal.scale.size)),
        np.ctypeslib.as_ctypes(mesh.surface.p_ctrl.reshape(mesh.surface.p_ctrl.size)),
        np.ctypeslib.as_ctypes(c_kj_x),
        np.ctypeslib.as_ctypes(c_kj_y),
        np.ctypeslib.as_ctypes(c_kj_z)
    )

    c_kj_x = c_kj_x.reshape((mesh.surface.nf, mesh.surface.nte))
    c_kj_y = c_kj_y.reshape((mesh.surface.nf, mesh.surface.nte))
    c_kj_z = c_kj_z.reshape((mesh.surface.nf, mesh.surface.nte))

    c_kj[:, :, 0] = c_kj_x[:, :]
    c_kj[:, :, 1] = c_kj_y[:, :]
    c_kj[:, :, 2] = c_kj_z[:, :]

    return

def solve_initial_condition(a_ij: np.ndarray,
                            b_kj: np.ndarray,
                            c_kj: np.ndarray,
                            mesh: models.MeshModel,
                            freestream: np.ndarray,
                            aero: models.AeroModel) -> None:
    
    # Linear system
    matrix = np.empty((mesh.surface.nf + 2 * mesh.surface.nte, mesh.surface.nf + 2 * mesh.surface.nte), dtype=np.double)
    array = np.empty((mesh.surface.nf + 2 * mesh.surface.nte), dtype=np.double)

    # Create rhs
    array[:mesh.surface.nf] = - (freestream[:, 0] * mesh.surface.e3[:, 0] + freestream[:, 1] * mesh.surface.e3[:, 1] + freestream[:, 2] * mesh.surface.e3[:, 2])
    array[mesh.surface.nf:] = 0.0

    # Create lhs
    # Transpiration
    for i in range(mesh.surface.nf):
        matrix[i, :mesh.surface.nf] = a_ij[i, :, 0] * mesh.surface.e3[i, 0] + a_ij[i, :, 1] * mesh.surface.e3[i, 1] + a_ij[i, :, 2] * mesh.surface.e3[i, 2]
        matrix[i, mesh.surface.nf:mesh.surface.nf + mesh.surface.nte] = b_kj[i, :, 0] * mesh.surface.e3[i, 0] + b_kj[i, :, 1] * mesh.surface.e3[i, 1] + b_kj[i, :, 2] * mesh.surface.e3[i, 2]
        matrix[i, mesh.surface.nf + mesh.surface.nte:mesh.surface.nf + 2 * mesh.surface.nte] = c_kj[i, :, 0] * mesh.surface.e3[i, 0] + c_kj[i, :, 1] * mesh.surface.e3[i, 1] + c_kj[i, :, 2] * mesh.surface.e3[i, 2]
        
    # Kutta condition
    for i in range(mesh.surface.nte):

        face1 = mesh.surface.trailing_edge_faces[i, 0]
        face2 = mesh.surface.trailing_edge_faces[i, 1]

        # e1
        matrix[i + mesh.surface.nf, :mesh.surface.nf] = (a_ij[face1, :, 0] * mesh.surface.e1[face1, 0] + a_ij[face1, :, 1] * mesh.surface.e1[face1, 1] + a_ij[face1, :, 2] * mesh.surface.e1[face1, 2]) - (a_ij[face2, :, 0] * mesh.surface.e1[face2, 0] + a_ij[face2, :, 1] * mesh.surface.e1[face2, 1] + a_ij[face2, :, 2] * mesh.surface.e1[face2, 2])
        matrix[i + mesh.surface.nf, mesh.surface.nf:mesh.surface.nf + mesh.surface.nte] = (b_kj[face1, :, 0] * mesh.surface.e1[face1, 0] + b_kj[face1, :, 1] * mesh.surface.e1[face1, 1] + b_kj[face1, :, 2] * mesh.surface.e1[face1, 2]) - (b_kj[face2, :, 0] * mesh.surface.e1[face2, 0] + b_kj[face2, :, 1] * mesh.surface.e1[face2, 1] + b_kj[face2, :, 2] * mesh.surface.e1[face2, 2])
        matrix[i + mesh.surface.nf, mesh.surface.nf + mesh.surface.nte:mesh.surface.nf + 2 * mesh.surface.nte] = (c_kj[face1, :, 0] * mesh.surface.e1[face1, 0] + c_kj[face1, :, 1] * mesh.surface.e1[face1, 1] + c_kj[face1, :, 2] * mesh.surface.e1[face1, 2]) - (c_kj[face2, :, 0] * mesh.surface.e1[face2, 0] + c_kj[face2, :, 1] * mesh.surface.e1[face2, 1] + c_kj[face2, :, 2] * mesh.surface.e1[face2, 2])

        # e2
        matrix[i + mesh.surface.nf + mesh.surface.nte, :mesh.surface.nf] = (a_ij[face1, :, 0] * mesh.surface.e2[face1, 0] + a_ij[face1, :, 1] * mesh.surface.e2[face1, 1] + a_ij[face1, :, 2] * mesh.surface.e2[face1, 2]) + (a_ij[face2, :, 0] * mesh.surface.e2[face2, 0] + a_ij[face2, :, 1] * mesh.surface.e2[face2, 1] + a_ij[face2, :, 2] * mesh.surface.e2[face2, 2])
        matrix[i + mesh.surface.nf + mesh.surface.nte, mesh.surface.nf:mesh.surface.nf + mesh.surface.nte] = (b_kj[face1, :, 0] * mesh.surface.e2[face1, 0] + b_kj[face1, :, 1] * mesh.surface.e2[face1, 1] + b_kj[face1, :, 2] * mesh.surface.e2[face1, 2]) + (b_kj[face2, :, 0] * mesh.surface.e2[face2, 0] + b_kj[face2, :, 1] * mesh.surface.e2[face2, 1] + b_kj[face2, :, 2] * mesh.surface.e2[face2, 2])
        matrix[i + mesh.surface.nf + mesh.surface.nte, mesh.surface.nf + mesh.surface.nte:mesh.surface.nf + 2 * mesh.surface.nte] = (c_kj[face1, :, 0] * mesh.surface.e2[face1, 0] + c_kj[face1, :, 1] * mesh.surface.e2[face1, 1] + c_kj[face1, :, 2] * mesh.surface.e2[face1, 2]) + (c_kj[face2, :, 0] * mesh.surface.e2[face2, 0] + c_kj[face2, :, 1] * mesh.surface.e2[face2, 1] + c_kj[face2, :, 2] * mesh.surface.e2[face2, 2])

    # Solve
    sol = np.linalg.solve(matrix, array)

    # Save
    aero.surface.source[:] = sol[:mesh.surface.nf]
    aero.internal.doublet[:] = sol[mesh.surface.nf:mesh.surface.nf + mesh.surface.nte]
    aero.wake.circulations[:, 0] = sol[mesh.surface.nf + mesh.surface.nte:mesh.surface.nf + 2 * mesh.surface.nte]
    
    return

def calculate_surface_parameters(a_ij: np.ndarray,
                                 b_kj: np.ndarray,
                                 c_kj: np.ndarray,
                                 freestream: np.ndarray,
                                 mesh: models.MeshModel,
                                 aero: models.AeroModel) -> None:
    
    aero.surface.vel[:, 0] = freestream[:, 0] + np.dot(a_ij[:, :, 0], aero.surface.source) + np.dot(b_kj[:, :, 0], aero.internal.doublet) + np.dot(c_kj[:, :, 0], aero.wake.circulations[:, 0])
    aero.surface.vel[:, 1] = freestream[:, 1] + np.dot(a_ij[:, :, 1], aero.surface.source) + np.dot(b_kj[:, :, 1], aero.internal.doublet) + np.dot(c_kj[:, :, 1], aero.wake.circulations[:, 0])
    aero.surface.vel[:, 2] = freestream[:, 2] + np.dot(a_ij[:, :, 2], aero.surface.source) + np.dot(b_kj[:, :, 2], aero.internal.doublet) + np.dot(c_kj[:, :, 2], aero.wake.circulations[:, 0])

    aero.surface.cp[:] = 1 - (aero.surface.vel[:, 0] * aero.surface.vel[:, 0] + aero.surface.vel[:, 1] * aero.surface.vel[:, 1] + aero.surface.vel[:, 2] * aero.surface.vel[:, 2]) / (freestream[:, 0] * freestream[:, 0] + freestream[:, 1] * freestream[:, 1] + freestream[:, 2] * freestream[:, 2])

    aero.surface.transpiration[:] = aero.surface.vel[:, 0] * mesh.surface.e3[:, 0] + aero.surface.vel[:, 1] * mesh.surface.e3[:, 1] + aero.surface.vel[:, 2] * mesh.surface.e3[:, 2]

    return

def calculate_vertices_parameters(vertices_connection: tp.List[models.VerticeConnectionModel],
                                  mesh: models.MeshModel,
                                  aero: models.AeroModel) -> None:
    
    for i in range(mesh.surface.nv):

        aero.surface.source_v[i] = 0.0
        aero.surface.vel_v[i, :] = 0.0
        aero.surface.cp_v[i] = 0.0
        aero.surface.transpiration_v[i] = 0.0

        for j in range(vertices_connection[i].n):

            aero.surface.source_v[i] = aero.surface.source_v[i] + vertices_connection[i].coefs[j] * aero.surface.source[vertices_connection[i].faces[j]]
            aero.surface.vel_v[i, :] = aero.surface.vel_v[i, :] + vertices_connection[i].coefs[j] * aero.surface.vel[vertices_connection[i].faces[j], :]
            aero.surface.cp_v[i] = aero.surface.cp_v[i] + vertices_connection[i].coefs[j] * aero.surface.cp[vertices_connection[i].faces[j]]
            aero.surface.transpiration_v[i] = aero.surface.transpiration_v[i] + vertices_connection[i].coefs[j] * aero.surface.transpiration[vertices_connection[i].faces[j]]

    return
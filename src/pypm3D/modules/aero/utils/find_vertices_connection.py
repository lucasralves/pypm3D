import typing as tp
import numpy as np

from pypm3D.models.vertice_model import VerticeModel


def main(nv: int,
         vertices: np.ndarray,
         faces: np.ndarray) -> tp.List[VerticeModel]:
    
    vertices_connection = []

    for v_id in range(nv):

        # Encontra as faces que se conectam com o v_id
        face_ids = [i[0] for i in np.argwhere((faces[:, 1] == v_id) | (faces[:, 2] == v_id) | (faces[:, 3] == v_id) | (faces[:, 4] == v_id))]

        # Encontra os ângulos
        angles = []

        for face_id in face_ids:
            
            if v_id == faces[face_id, 1]:
                
                v1 = vertices[faces[face_id, 2], :] - vertices[faces[face_id, 1], :]
                
                if faces[face_id, 0] == 3:
                    v2 = vertices[faces[face_id, 3], :] - vertices[faces[face_id, 1], :]
                else:
                    v2 = vertices[faces[face_id, 4], :] - vertices[faces[face_id, 1], :]

            elif v_id == faces[face_id, 2]:
                v1 = vertices[faces[face_id, 3], :] - vertices[faces[face_id, 2], :]
                v2 = vertices[faces[face_id, 1], :] - vertices[faces[face_id, 2], :]

            elif v_id == faces[face_id, 3]:

                v2 = vertices[faces[face_id, 2], :] - vertices[faces[face_id, 3], :]

                if faces[face_id, 0] == 3:
                    v1 = vertices[faces[face_id, 1], :] - vertices[faces[face_id, 3], :]
                else:
                    v1 = vertices[faces[face_id, 4], :] - vertices[faces[face_id, 3], :]

            else:
                v1 = vertices[faces[face_id, 1], :] - vertices[faces[face_id, 4], :]
                v2 = vertices[faces[face_id, 3], :] - vertices[faces[face_id, 4], :]
            
            angles.append(np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))))

        sum_angles = sum(angles)

        # Adiciona um modelo
        vertices_connection.append(
            VerticeModel(
                n=len(face_ids),
                faces=face_ids,
                coefs=[a / sum_angles for a in angles],
            )
        )
        
    return vertices_connection
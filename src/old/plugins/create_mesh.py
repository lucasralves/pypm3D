from typing import Dict
from numpy import ndarray

from pypm3D.core.mesh import get_trailing_edge_faces, calculate_faces_parameters, calculate_inner_doublet_pannels

def create_mesh(vertices: ndarray, faces: ndarray, trailing_edge: ndarray) -> Dict[str, ndarray]:

    dict_1 = get_trailing_edge_faces(faces, trailing_edge)
    dict_2 = calculate_faces_parameters(vertices, faces, trailing_edge, dict_1['trailing_edge_faces'])
    dict_3 = calculate_inner_doublet_pannels(vertices, faces, trailing_edge, dict_1['trailing_edge_faces'], dict_2['e3'])
    
    out = {'vertices': vertices, 'faces': faces, 'trailing_edge': trailing_edge} | dict_1 | dict_2 | dict_3
    
    return out
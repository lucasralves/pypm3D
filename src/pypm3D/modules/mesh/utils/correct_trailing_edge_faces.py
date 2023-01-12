import typing as tp
import numpy as np


def main(trailing_edge_faces: np.ndarray,
         e3: np.ndarray,
         func: tp.Callable[[int, int, int], None]) -> None:
    

    flip = e3[trailing_edge_faces[:, 1], 2] > 0

    id1 = trailing_edge_faces[flip, 1]
    id2 = trailing_edge_faces[flip, 0]

    trailing_edge_faces[flip, 0] = id1[:]
    trailing_edge_faces[flip, 1] = id2[:]

    func(trailing_edge_faces)

    return
import typing as tp
import numpy as np


def main(nte: int,
         faces: np.ndarray,
         trailing_edge: np.ndarray,
         func: tp.Callable[[int, int, int], None]) -> None:

    for i in range(nte):
        check1 = (faces[:, 1] == trailing_edge[i, 0]) | (faces[:, 2] == trailing_edge[i, 0]) | (faces[:, 3] == trailing_edge[i, 0]) | (faces[:, 4] == trailing_edge[i, 0])
        check2 = (faces[:, 1] == trailing_edge[i, 1]) | (faces[:, 2] == trailing_edge[i, 1]) | (faces[:, 3] == trailing_edge[i, 1]) | (faces[:, 4] == trailing_edge[i, 1])
        index = np.argwhere(check1 & check2)
        func(i, index[0][0], index[1][0])

    return
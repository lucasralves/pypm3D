import typing as tp
import numpy as np


def main(freestream: np.ndarray,
         e3: np.ndarray,
         func: tp.Callable[[np.ndarray], None]) -> None:
    
    source = - (freestream[:, 0] * e3[:, 0] + freestream[:, 1] * e3[:, 1] + freestream[:, 2] * e3[:, 2])
    func(source)

    return
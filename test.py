import numpy as np
from typing import NamedTuple

class Test(NamedTuple):
    a: np.ndarray

a = Test(np.empty((5, 2), dtype=np.double))

print(a.a)

a.a[:, 0] = -1

print(a.a)
from typing import NamedTuple
from numpy import ndarray

class SolutionModel(NamedTuple):
    source_f: ndarray             # (nf,)
    transpiration_f: ndarray      # (nf,)
    vel_f: ndarray                # (nf, 3)
    cp_f: ndarray                 # (nf,)

    source_v: ndarray             # (nv,)
    transpiration_v: ndarray      # (nv,)
    vel_v: ndarray                # (nv, 3)
    cp_v: ndarray                 # (nv,)

    doublet: ndarray            # (nte, nw)
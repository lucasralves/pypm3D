from dataclasses import dataclass
from numpy import ndarray

@dataclass
class MeshModel:
    nv: int = None
    nf: int = None
    nte: int = None
    nw: int = None
    nwv: int = None
    vertices: ndarray = None           # (nv, 3)
    faces: ndarray = None              # (nf, 5)
    trailing_edge: ndarray = None      # (nte, 2)
    p_avg: ndarray = None              # (nf, 3)
    p_ctrl: ndarray = None             # (nf, 3)
    e1: ndarray = None                 # (nf, 3)
    e2: ndarray = None                 # (nf, 3)
    e3: ndarray = None                 # (nf, 3)
    p1: ndarray = None                 # (nf, 2)
    p2: ndarray = None                 # (nf, 2)
    p3: ndarray = None                 # (nf, 2)
    p4: ndarray = None                 # (nf, 2)
    te_faces: ndarray = None           # (nte, 2)
    inner_faces: ndarray = None        # (nte, 5)
    inner_vertices: ndarray = None     # (?, 3)
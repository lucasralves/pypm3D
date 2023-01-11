import typing as tp
import numpy as np
import math

class AeroWakeModel:

    def __init__(self) -> None:
        
        self.nv: int = None
        self.nw: int = None

        self.trailing_edge_ids: np.ndarray = None
        self.vertices: np.ndarray = None
        self.faces: np.ndarray = None

        self.circulations: np.ndarray = None
        self.areas: np.ndarray = None

        return

class AeroSurfaceModel:

    def __init__(self, nf: int,
                       nv: int) -> None:
        
        self.source: np.ndarray = np.empty((nf,), dtype=np.double)
        self.vel: np.ndarray = np.empty((nf, 3), dtype=np.double)
        self.cp: np.ndarray = np.empty((nf,), dtype=np.double)
        self.transpiration: np.ndarray = np.empty((nf,), dtype=np.double)

        self.source_v: np.ndarray = np.empty((nv,), dtype=np.double)
        self.vel_v: np.ndarray = np.empty((nv, 3), dtype=np.double)
        self.cp_v: np.ndarray = np.empty((nv,), dtype=np.double)
        self.transpiration_v: np.ndarray = np.empty((nv,), dtype=np.double)

class AeroInternalModel:

    def __init__(self, nte: int) -> None:

        self.doublet: np.ndarray = np.empty((nte,), dtype=np.double)

        return

class AeroModel:
    
    def __init__(self, nf: int,
                       nv: int,
                       nte: int) -> None:

        self.wake: AeroWakeModel = AeroWakeModel()
        self.surface: AeroSurfaceModel = AeroSurfaceModel(nf, nv)
        self.internal: AeroInternalModel = AeroInternalModel(nte)

        return
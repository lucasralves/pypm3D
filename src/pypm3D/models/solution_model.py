from dataclasses import dataclass
import numpy as np


@dataclass
class SolutionModel:
    source: np.ndarray = None
    doublet: np.ndarray = None
    vel: np.ndarray = None
    cp: np.ndarray = None
    transpiration: np.ndarray = None
    source_v: np.ndarray = None
    doublet_v: np.ndarray = None
    vel_v: np.ndarray = None
    cp_v: np.ndarray = None
    transpiration_v: np.ndarray = None
    wake_doublet: np.ndarray = None

    def init(self, nf: int, nv: int, nte: int, nw: int) -> None:

        self.source = np.empty((nf,), dtype=np.double)
        self.doublet = np.empty((nf,), dtype=np.double)
        self.vel = np.empty((nf, 3), dtype=np.double)
        self.cp = np.empty((nf,), dtype=np.double)
        self.transpiration = np.empty((nf,), dtype=np.double)

        self.source_v = np.empty((nv,), dtype=np.double)
        self.doublet_v = np.empty((nv,), dtype=np.double)
        self.vel_v = np.empty((nv, 3), dtype=np.double)
        self.cp_v = np.empty((nv,), dtype=np.double)
        self.transpiration_v = np.empty((nv,), dtype=np.double)

        self.wake_doublet = np.empty((nte, nw), dtype=np.double)

        return

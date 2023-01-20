from dataclasses import dataclass
import numpy as np


@dataclass
class SolutionModel:
    source: np.ndarray = None
    doublet: np.ndarray = None
    vel: np.ndarray = None
    cp: np.ndarray = None
    transpiration: np.ndarray = None
    potential: np.ndarray = None
    source_v: np.ndarray = None
    doublet_v: np.ndarray = None
    vel_v: np.ndarray = None
    cp_v: np.ndarray = None
    transpiration_v: np.ndarray = None
    potential_v: np.ndarray = None
    wake_doublet: np.ndarray = None

    def init(self, nf: int, nv: int, nte: int, nw: int) -> None:

        self.source = np.empty((nf,), dtype=np.double)
        self.doublet = np.empty((nf,), dtype=np.double)
        self.vel = np.empty((nf, 3), dtype=np.double)
        self.cp = np.empty((nf,), dtype=np.double)
        self.transpiration = np.empty((nf,), dtype=np.double)
        self.potential = np.empty((nf,), dtype=np.double)

        self.source_v = np.empty((nv,), dtype=np.double)
        self.doublet_v = np.empty((nv,), dtype=np.double)
        self.vel_v = np.empty((nv, 3), dtype=np.double)
        self.cp_v = np.empty((nv,), dtype=np.double)
        self.transpiration_v = np.empty((nv,), dtype=np.double)
        self.potential_v = np.empty((nv,), dtype=np.double)

        self.wake_doublet = np.empty((nte * nw), dtype=np.double)

        return

    def add_doublet(self, vals: np.ndarray) -> None:
        self.doublet[:] = vals[:]
        return
    
    def add_source(self, vals: np.ndarray) -> None:
        self.source[:] = vals[:]
        return
    
    def add_multiple(self, vel: np.ndarray, cp: np.ndarray, transpiration: np.ndarray, potential: np.ndarray, vel_v: np.ndarray, cp_v: np.ndarray, transpiration_v: np.ndarray, potential_v: np.ndarray, source_v: np.ndarray, doublet_v: np.ndarray) -> None:
        self.vel[:, :] = vel[:, :]
        self.cp[:] = cp[:]
        self.transpiration[:] = transpiration[:]
        self.potential[:] = potential[:]
        self.vel_v[:, :] = vel_v[:, :]
        self.cp_v[:] = cp_v[:]
        self.transpiration_v[:] = transpiration_v[:]
        self.potential_v[:] = potential_v[:]
        self.source_v[:] = source_v[:]
        self.doublet_v[:] = doublet_v[:]
        return
    
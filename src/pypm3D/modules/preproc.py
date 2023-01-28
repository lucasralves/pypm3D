import typing as tp
import numpy as np
import math


class SurfaceMesh:

    def __init__(self) -> None:
        self.nf: int = None                      # number of faces
        self.nv: int = None                      # number of vertices
        self.nte: int = None                     # number of edges at the trailing edge
        self.vt: np.ndarray = None               # vertices: [[x1, y1, z1], ..., [xnv, ynv, znv]]
        self.fc: np.ndarray = None               # faces: [[n1, id1, id2, id3, id4], ..., [nnf, idnf, idnf, idnf, idnf]]
        self.te: np.ndarray = None               # trailing edge vertices: [[id11, id12], ..., [idnte1, idnte2]]
        self.te_fc: np.ndarray = None            # faces that share the same trailing edge: [[id11, id12], ..., [idnte1, idnte2]]
        self.te_vt: np.ndarray = None            # ids of vertices at trailing edge
        self.p_avg: np.ndarray = None            # face's center: [[x1, y1, z1], ..., [xnf, ynf, znf]]
        self.p_ctrl: np.ndarray = None           # control point: [[x1, y1, z1], ..., [xnf, ynf, znf]]
        self.p1: np.ndarray = None               # face's vertice (id1) at local system: [[x1, y1], .., [xnf, ynf]]
        self.p2: np.ndarray = None               # face's vertice (id2) at local system: [[x1, y1], .., [xnf, ynf]]
        self.p3: np.ndarray = None               # face's vertice (id3) at local system: [[x1, y1], .., [xnf, ynf]]
        self.p4: np.ndarray = None               # face's vertice (id4) at local system: [[x1, y1], .., [xnf, ynf]]
        self.e1: np.ndarray = None               # face's orthogonal base: [[x1, y1, z1], .., [xnf, ynf, znf]]
        self.e2: np.ndarray = None               # face's orthogonal base: [[x1, y1, z1], .., [xnf, ynf, znf]]
        self.e3: np.ndarray = None               # face's orthogonal base: [[x1, y1, z1], .., [xnf, ynf, znf]]
        return
    
    def init(self, vt: np.ndarray, fc: np.ndarray, te: np.ndarray) -> None:
        """
        Creates the surface's parameters.
        
        The geometry must be oriented as:
            - The span must be in the y direction
            - The longitidunal axis must be pointing to the negative direction of x.
            - The upper side must be in the positive z direction

        Parameters:
        -----------
        - vt: vertices
        - fc: faces
        - te: trailing edge
        """
        
        # Parameters
        self.nf: int = fc.shape[0]
        self.nv: int = vt.shape[0]
        self.nte: int = te.shape[0]

        self.vt: np.ndarray = np.copy(vt).astype(np.double)
        self.fc: np.ndarray = np.copy(fc).astype(np.int32); self.fc[self.fc[:, 0] == 3, 4] = -1
        self.te: np.ndarray = np.copy(te).astype(np.int32)

        self.te_fc = np.empty((self.nte, 2), dtype=np.int32)
        self.p_avg = np.empty((self.nf, 3), dtype=np.double)
        self.p_ctrl = np.empty((self.nf, 3), dtype=np.double)
        self.p1 = np.empty((self.nf, 2), dtype=np.double)
        self.p2 = np.empty((self.nf, 2), dtype=np.double)
        self.p3 = np.empty((self.nf, 2), dtype=np.double)
        self.p4 = np.empty((self.nf, 2), dtype=np.double)
        self.e1 = np.empty((self.nf, 3), dtype=np.double)
        self.e2 = np.empty((self.nf, 3), dtype=np.double)
        self.e3 = np.empty((self.nf, 3), dtype=np.double)

        # Face geometry
        tri_panel_ids = self.fc[:, 0] == 3
        quad_panel_ids = self.fc[:, 0] == 4

        self.p_avg[quad_panel_ids, :] = 0.25 * (self.vt[self.fc[quad_panel_ids, 1], :] + self.vt[self.fc[quad_panel_ids, 2], :] + self.vt[self.fc[quad_panel_ids, 3], :] + self.vt[self.fc[quad_panel_ids, 4], :])
        self.p_avg[tri_panel_ids, :] = (1. / 3.) * (self.vt[self.fc[tri_panel_ids, 1], :] + self.vt[self.fc[tri_panel_ids, 2], :] + self.vt[self.fc[tri_panel_ids, 3], :])
        
        vec3 = np.empty((self.nf, 3), dtype=np.double)
        vec3[quad_panel_ids, :] = np.cross(self.vt[self.fc[quad_panel_ids, 2], :] - self.vt[self.fc[quad_panel_ids, 4], :], self.vt[self.fc[quad_panel_ids, 3], :] - self.vt[self.fc[quad_panel_ids, 1], :])
        vec3[tri_panel_ids, :] = np.cross(self.vt[self.fc[tri_panel_ids, 2], :] - self.vt[self.fc[tri_panel_ids, 1], :], self.vt[self.fc[tri_panel_ids, 3], :] - self.vt[self.fc[tri_panel_ids, 1], :])
        vec1 = 0.5 * (self.vt[self.fc[:, 2], :] + self.vt[self.fc[:, 3], :]) - self.p_avg[:, :]
        vec2 = np.cross(vec3, vec1)
        vec3 = np.cross(vec1, vec2)

        self.e1[:, 0] = vec1[:, 0] / np.sqrt(vec1[:, 0] * vec1[:, 0] + vec1[:, 1] * vec1[:, 1] + vec1[:, 2] * vec1[:, 2])
        self.e1[:, 1] = vec1[:, 1] / np.sqrt(vec1[:, 0] * vec1[:, 0] + vec1[:, 1] * vec1[:, 1] + vec1[:, 2] * vec1[:, 2])
        self.e1[:, 2] = vec1[:, 2] / np.sqrt(vec1[:, 0] * vec1[:, 0] + vec1[:, 1] * vec1[:, 1] + vec1[:, 2] * vec1[:, 2])

        self.e2[:, 0] = vec2[:, 0] / np.sqrt(vec2[:, 0] * vec2[:, 0] + vec2[:, 1] * vec2[:, 1] + vec2[:, 2] * vec2[:, 2])
        self.e2[:, 1] = vec2[:, 1] / np.sqrt(vec2[:, 0] * vec2[:, 0] + vec2[:, 1] * vec2[:, 1] + vec2[:, 2] * vec2[:, 2])
        self.e2[:, 2] = vec2[:, 2] / np.sqrt(vec2[:, 0] * vec2[:, 0] + vec2[:, 1] * vec2[:, 1] + vec2[:, 2] * vec2[:, 2])

        self.e3[:, 0] = vec3[:, 0] / np.sqrt(vec3[:, 0] * vec3[:, 0] + vec3[:, 1] * vec3[:, 1] + vec3[:, 2] * vec3[:, 2])
        self.e3[:, 1] = vec3[:, 1] / np.sqrt(vec3[:, 0] * vec3[:, 0] + vec3[:, 1] * vec3[:, 1] + vec3[:, 2] * vec3[:, 2])
        self.e3[:, 2] = vec3[:, 2] / np.sqrt(vec3[:, 0] * vec3[:, 0] + vec3[:, 1] * vec3[:, 1] + vec3[:, 2] * vec3[:, 2])

        self.p_ctrl[:, :] = self.p_avg[:, :] - self.e3[:, :] * 1e-8
        
        vec1 = self.vt[self.fc[:, 1], :] - self.p_avg[:, :]
        vec2 = self.vt[self.fc[:, 2], :] - self.p_avg[:, :]
        vec3 = self.vt[self.fc[:, 3], :] - self.p_avg[:, :]
        vec4 = self.vt[self.fc[quad_panel_ids, 4], :] - self.p_avg[quad_panel_ids, :]

        self.p1[:, 0] = self.e1[:, 0] * vec1[:, 0] + self.e1[:, 1] * vec1[:, 1] + self.e1[:, 2] * vec1[:, 2]
        self.p2[:, 0] = self.e1[:, 0] * vec2[:, 0] + self.e1[:, 1] * vec2[:, 1] + self.e1[:, 2] * vec2[:, 2]
        self.p3[:, 0] = self.e1[:, 0] * vec3[:, 0] + self.e1[:, 1] * vec3[:, 1] + self.e1[:, 2] * vec3[:, 2]
        self.p4[quad_panel_ids, 0] = self.e1[quad_panel_ids, 0] * vec4[:, 0] + self.e1[quad_panel_ids, 1] * vec4[:, 1] + self.e1[quad_panel_ids, 2] * vec4[:, 2]

        self.p1[:, 1] = self.e2[:, 0] * vec1[:, 0] + self.e2[:, 1] * vec1[:, 1] + self.e2[:, 2] * vec1[:, 2]
        self.p2[:, 1] = self.e2[:, 0] * vec2[:, 0] + self.e2[:, 1] * vec2[:, 1] + self.e2[:, 2] * vec2[:, 2]
        self.p3[:, 1] = self.e2[:, 0] * vec3[:, 0] + self.e2[:, 1] * vec3[:, 1] + self.e2[:, 2] * vec3[:, 2]
        self.p4[quad_panel_ids, 1] = self.e2[quad_panel_ids, 0] * vec4[:, 0] + self.e2[quad_panel_ids, 1] * vec4[:, 1] + self.e2[quad_panel_ids, 2] * vec4[:, 2]


        # Trailing edge

        # The vertices must be in the negative direction of y
        flip = self.vt[self.te[:, 0], 1] < self.vt[self.te[:, 1], 1]
        id1 = self.te[flip, 1]
        id2 = self.te[flip, 0]
        self.te[flip, 0] = id1[:]
        self.te[flip, 1] = id2[:]

        # Find the vertices at trailing edge
        te_vt = []
        for edge in self.te:
            if edge[0] not in te_vt: te_vt.append(edge[0])
            if edge[1] not in te_vt: te_vt.append(edge[1])
        self.te_vt = np.asarray(te_vt)

        # Find the trailing edge faces
        for i in range(self.nte):
            check1 = (self.fc[:, 1] == self.te[i, 0]) | (self.fc[:, 2] == self.te[i, 0]) | (self.fc[:, 3] == self.te[i, 0]) | (self.fc[:, 4] == self.te[i, 0])
            check2 = (self.fc[:, 1] == self.te[i, 1]) | (self.fc[:, 2] == self.te[i, 1]) | (self.fc[:, 3] == self.te[i, 1]) | (self.fc[:, 4] == self.te[i, 1])
            index = np.argwhere(check1 & check2)
            id1, id2 = index[0][0], index[1][0]
            if self.e3[id2, 2] > self.e3[id1, 2]:
                self.te_fc[i, 0] = id2
                self.te_fc[i, 1] = id1
            else:
                self.te_fc[i, 0] = id1
                self.te_fc[i, 1] = id2
        
        return

class WakeMesh:

    def __init__(self, sf: SurfaceMesh) -> None:
        # private
        self._sf: SurfaceMesh = sf

        # public
        self.nf: int = None                      # number of faces
        self.nv: int = None                      # number of vertices
        self.nw: int = None                      # number of sections in the freestream direction
        self.nv_te: int = None                   # number of vertices at the trailing edge
        self.fc: np.ndarray = None               # faces: [[n1, id1, id2, id3, id4], ..., [nnf, idnf, idnf, idnf, idnf]]
        self.ts: float = None                    # pseudo time step

        # The parameters below vary at every time step
        self.sec: int = None                     # number of wake sections (starts at 1)
        self.vt: np.ndarray = None               # vertices: [[x1, y1, z1], ..., [xnv, ynv, znv]]
        self.p_avg: np.ndarray = None            # face's center: [[x1, y1, z1], ..., [xnf, ynf, znf]]
        self.p1: np.ndarray = None               # face's vertice (id1) at local system: [[x1, y1], .., [xnf, ynf]]
        self.p2: np.ndarray = None               # face's vertice (id2) at local system: [[x1, y1], .., [xnf, ynf]]
        self.p3: np.ndarray = None               # face's vertice (id3) at local system: [[x1, y1], .., [xnf, ynf]]
        self.p4: np.ndarray = None               # face's vertice (id4) at local system: [[x1, y1], .., [xnf, ynf]]
        self.e1: np.ndarray = None               # face's orthogonal base: [[x1, y1, z1], .., [xnf, ynf, znf]]
        self.e2: np.ndarray = None               # face's orthogonal base: [[x1, y1, z1], .., [xnf, ynf, znf]]
        self.e3: np.ndarray = None               # face's orthogonal base: [[x1, y1, z1], .., [xnf, ynf, znf]]
        self.areas: np.ndarray = None            # face's areas: [a1, ..., anf]

        return

    @property
    def nte(self) -> int:
        return self._sf.nte

    def init(self, u_call: tp.Callable[[np.ndarray], np.ndarray], l: float, ts: float) -> None:
        """
        Creates the wake's parameters

        Parameters:
        -----------
        - u: freestream velocity
        - l: wake's length
        - ts: pseudo time step
        """
        
        u = np.linalg.norm(u_call(np.mean(self._sf.vt, axis=0)))

        self.nw = math.ceil(l / (u * ts))
        self.nv_te: int = self._sf.te_vt.shape[0]
        self.nv = self.nv_te * (self.nw + 1)
        self.nf = self.nte * self.nw
        
        self.ts = ts

        self.vt = np.empty((self.nv, 3), dtype=np.double)
        self.fc = np.empty((self.nf, 5), dtype=np.int32)

        for i in range(self.nv_te):
            self.vt[i, :] = self._sf.vt[self._sf.te_vt[i], :]
            
        for i in range(self.nw):
            for j in range(self.nte):
                id1 = np.argwhere(self._sf.te_vt == self._sf.te[j, 0])[0][0]
                id2 = np.argwhere(self._sf.te_vt == self._sf.te[j, 1])[0][0]

                self.fc[i * self.nte + j, 0] = 4
                self.fc[i * self.nte + j, 1] = id1 + i * self.nv_te
                self.fc[i * self.nte + j, 2] = id2 + i * self.nv_te
                self.fc[i * self.nte + j, 3] = self.fc[i * self.nte + j, 2] + self.nv_te
                self.fc[i * self.nte + j, 4] = self.fc[i * self.nte + j, 1] + self.nv_te

        self.sec = 0

        self.p_avg = np.empty((self.nf, 3), dtype=np.double)
        self.p1 = np.empty((self.nf, 2), dtype=np.double)
        self.p2 = np.empty((self.nf, 2), dtype=np.double)
        self.p3 = np.empty((self.nf, 2), dtype=np.double)
        self.p4 = np.empty((self.nf, 2), dtype=np.double)
        self.e1 = np.empty((self.nf, 3), dtype=np.double)
        self.e2 = np.empty((self.nf, 3), dtype=np.double)
        self.e3 = np.empty((self.nf, 3), dtype=np.double)
        self.areas = np.empty(self.nf, dtype=np.double)

        # Create one section
        v = np.empty((self.nv_te, 3), dtype=np.double)
        v_sum = np.empty(3, dtype=np.double)
        for i in range(self.nv_te):
            v_sum[:] = .0
            count = 0
            for edge_id in range(self.nte):
                if self._sf.te_vt[i] in self._sf.te[edge_id, :]:
                    count += 1
                    id1 = self._sf.te_fc[edge_id, 0]
                    id2 = self._sf.te_fc[edge_id, 1]
                    f1 = u_call(self._sf.p_avg[id1, :])
                    f2 = u_call(self._sf.p_avg[id2, :])
                    v1 = f1 - self._sf.e3[id1, :] * np.dot(f1, self._sf.e3[id1, :])
                    v2 = f2 - self._sf.e3[id2, :] * np.dot(f2, self._sf.e3[id2, :])
                    v_sum[:] = v_sum[:] + v1 / np.linalg.norm(v1) + v2 / np.linalg.norm(v2)
                if count == 2: break
            v_sum_norm = np.linalg.norm(v_sum)
            v[i, :] = v_sum[:] / v_sum_norm
        
        v[:, :] = u * v[:, :]

        self.update(v)

        return
    
    def update(self, v: np.ndarray) -> None:
        """
        Add one section to the wake

        Parameters:
        -----------
        - v: velocity field
        """
        
        if self.sec == self.nw:
            raise "Cannot add a section to the wake (self.sec == self.nw)"
        else:
            self.sec += 1
        
        # Update vertices
        nmax = v.shape[0]
        for i in range(nmax):
            old = nmax - i - 1
            new = old + self.nv_te
            self.vt[new, :] = self.vt[old, :] + v[old, :] * self.ts

        # Calculate parameters
        idmax = self.nte * self.sec

        # Face center
        self.p_avg[:idmax, :] = (1. / 4.) * (self.vt[self.fc[:idmax, 1], :] + self.vt[self.fc[:idmax, 2], :] + self.vt[self.fc[:idmax, 3], :] + self.vt[self.fc[:idmax, 4], :])
        
        # Orthogonal base (e3)
        self.e3[:idmax, :] = np.cross(self.vt[self.fc[:idmax, 2], :] - self.vt[self.fc[:idmax, 4], :], self.vt[self.fc[:idmax, 3], :] - self.vt[self.fc[:idmax, 1], :])
        
        # Orthogonal base (e1)
        self.e1[:idmax, :] = 0.5 * (self.vt[self.fc[:idmax, 1], :] + self.vt[self.fc[:idmax, 2], :]) - self.p_avg[:idmax, :]

        # Orthogonal base (e2)
        self.e2[:idmax, :] = np.cross(self.e3[:idmax, :], self.e1[:idmax, :])

        # Orthogonal base (e3 - correção)
        self.e3[:idmax, :] = np.cross(self.e1[:idmax, :], self.e2[:idmax, :])

        # Normalize orthogonal base
        aux_norm = np.sqrt(self.e1[:idmax, 0] * self.e1[:idmax, 0] + self.e1[:idmax, 1] * self.e1[:idmax, 1] + self.e1[:idmax, 2] * self.e1[:idmax, 2])

        self.e1[:idmax, 0] = self.e1[:idmax, 0] / aux_norm
        self.e1[:idmax, 1] = self.e1[:idmax, 1] / aux_norm
        self.e1[:idmax, 2] = self.e1[:idmax, 2] / aux_norm

        aux_norm = np.sqrt(self.e2[:idmax, 0] * self.e2[:idmax, 0] + self.e2[:idmax, 1] * self.e2[:idmax, 1] + self.e2[:idmax, 2] * self.e2[:idmax, 2])

        self.e2[:idmax, 0] = self.e2[:idmax, 0] / aux_norm
        self.e2[:idmax, 1] = self.e2[:idmax, 1] / aux_norm
        self.e2[:idmax, 2] = self.e2[:idmax, 2] / aux_norm

        aux_norm = np.sqrt(self.e3[:idmax, 0] * self.e3[:idmax, 0] + self.e3[:idmax, 1] * self.e3[:idmax, 1] + self.e3[:idmax, 2] * self.e3[:idmax, 2])

        self.e3[:idmax, 0] = self.e3[:idmax, 0] / aux_norm
        self.e3[:idmax, 1] = self.e3[:idmax, 1] / aux_norm
        self.e3[:idmax, 2] = self.e3[:idmax, 2] / aux_norm

        # Local vertices
        aux_vec = self.vt[self.fc[:idmax, 1], :] - self.p_avg[:idmax, :]
        self.p1[:idmax, 0] = self.e1[:idmax, 0] * aux_vec[:, 0] + self.e1[:idmax, 1] * aux_vec[:, 1] + self.e1[:idmax, 2] * aux_vec[:, 2]
        self.p1[:idmax, 1] = self.e2[:idmax, 0] * aux_vec[:, 0] + self.e2[:idmax, 1] * aux_vec[:, 1] + self.e2[:idmax, 2] * aux_vec[:, 2]

        aux_vec = self.vt[self.fc[:idmax, 2], :] - self.p_avg[:idmax, :]
        self.p2[:idmax, 0] = self.e1[:idmax, 0] * aux_vec[:, 0] + self.e1[:idmax, 1] * aux_vec[:, 1] + self.e1[:idmax, 2] * aux_vec[:, 2]
        self.p2[:idmax, 1] = self.e2[:idmax, 0] * aux_vec[:, 0] + self.e2[:idmax, 1] * aux_vec[:, 1] + self.e2[:idmax, 2] * aux_vec[:, 2]

        aux_vec = self.vt[self.fc[:idmax, 3], :] - self.p_avg[:idmax, :]
        self.p3[:idmax, 0] = self.e1[:idmax, 0] * aux_vec[:, 0] + self.e1[:idmax, 1] * aux_vec[:, 1] + self.e1[:idmax, 2] * aux_vec[:, 2]
        self.p3[:idmax, 1] = self.e2[:idmax, 0] * aux_vec[:, 0] + self.e2[:idmax, 1] * aux_vec[:, 1] + self.e2[:idmax, 2] * aux_vec[:, 2]

        aux_vec = self.vt[self.fc[:idmax, 4], :] - self.p_avg[:idmax, :]
        self.p4[:idmax, 0] = self.e1[:idmax, 0] * aux_vec[:, 0] + self.e1[:idmax, 1] * aux_vec[:, 1] + self.e1[:idmax, 2] * aux_vec[:, 2]
        self.p4[:idmax, 1] = self.e2[:idmax, 0] * aux_vec[:, 0] + self.e2[:idmax, 1] * aux_vec[:, 1] + self.e2[:idmax, 2] * aux_vec[:, 2]

        self.areas[:idmax] = 0.5 * np.abs(self.p1[:idmax, 0] * self.p2[:idmax, 1] + self.p2[:idmax, 0] * self.p3[:idmax, 1] + self.p3[:idmax, 0] * self.p4[:idmax, 1] + self.p4[:idmax, 0] * self.p1[:idmax, 1] - (self.p1[:idmax, 1] * self.p2[:idmax, 0] + self.p2[:idmax, 1] * self.p3[:idmax, 0] + self.p3[:idmax, 1] * self.p4[:idmax, 0] + self.p4[:idmax, 1] * self.p1[:idmax, 0]))

        return

class _Mesh:

    def __init__(self) -> None:
        self.sf: SurfaceMesh = SurfaceMesh()
        self.wk: WakeMesh = WakeMesh(self.sf)
        return

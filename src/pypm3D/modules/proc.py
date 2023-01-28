import typing as tp
import numpy as np

from pypm3D.modules.preproc import _Mesh
from pypm3D.bin import potential_wrapper, velocity_wrapper


class FacesConnection(tp.NamedTuple):
    n_sides: int
    c1: tp.List[int]
    c2: tp.List[int]
    c3: tp.List[int]
    c4: tp.List[int] | None

class _Solver():

    def __init__(self, mesh: _Mesh) -> None:
        # private
        self._mesh: _Mesh = mesh

        # public
        self.source_sf: np.ndarray = None            # face's source strength
        self.doublet_sf: np.ndarray = None           # face's doublet strength
        self.vel_sf: np.ndarray = None               # face's velocity
        self.cp_sf: np.ndarray = None                # face's pressure coefficient
        self.pt: np.ndarray = None                   # face's potential
        self.trp: np.ndarray = None                  # face's transpiration
        self.doublet_wk: np.ndarray = None           # wake's doublet strength

    @property
    def done(self) -> bool:
        return self.source_sf is not None

    def run(self, u_call: tp.Callable[[np.ndarray], np.ndarray]) -> None:
        """
        Solves a potential flow using dirichlet boundarycondition
        over an arbitrary lifting object.

        Parameters:
        -----------
        - u_call: callable that receives a position (x, y, z) and returns
          a velocity (ux, uy, uz).
        """

        # Initialize public parameters
        self.source_sf = np.empty(self._mesh.sf.nf, dtype=np.double)
        self.doublet_sf = np.empty(self._mesh.sf.nf, dtype=np.double)
        self.vel_sf = np.empty((self._mesh.sf.nf, 3), dtype=np.double)
        self.cp_sf = np.empty(self._mesh.sf.nf, dtype=np.double)
        self.pt = np.empty(self._mesh.sf.nf, dtype=np.double)
        self.trp = np.empty(self._mesh.sf.nf, dtype=np.double)
        self.doublet_wk = np.empty(self._mesh.wk.nf, dtype=np.double)

        # Initialize local parameters
        a_ij = np.empty((self._mesh.sf.nf, self._mesh.sf.nf), dtype=np.double)          # source coefs
        b_ij = np.empty((self._mesh.sf.nf, self._mesh.sf.nf), dtype=np.double)          # doublet coefs
        c_ik = np.empty((self._mesh.sf.nf, self._mesh.sf.nte), dtype=np.double)         # wake te coefs
        d_i = np.empty(self._mesh.sf.nf, dtype=np.double)                               # wake coefs
        lhs = np.empty((self._mesh.sf.nf, self._mesh.sf.nf), dtype=np.double)           # linear system lhs
        rhs = np.empty(self._mesh.sf.nf, dtype=np.double)                               # linear system rhs
        areas = np.empty(self._mesh.wk.nf, dtype=np.double)                             # initial area of the wake sections

        # Calculate source strngth
        self.source_sf[:] = np.asarray([np.dot(self._mesh.sf.e3[face, :], u_call(self._mesh.sf.p_avg[face, :])) for face in range(self._mesh.sf.nf)])

        # Calculate a_ij, b_ij
        for i in range(self._mesh.sf.nf):
            a_ij[:, i], b_ij[:, i] = potential_wrapper.main(self._mesh.sf.fc[i, 0], self._mesh.sf.p1[i, :], self._mesh.sf.p2[i, :], self._mesh.sf.p3[i, :], self._mesh.sf.p4[i, :], self._mesh.sf.e1[i, :], self._mesh.sf.e2[i, :], self._mesh.sf.e3[i, :], self._mesh.sf.p_avg[i, :], self._mesh.sf.p_ctrl[:, :])
        
        # Vertices connection
        fc_cn = self._faces_connection()

        # Create free wake
        while True:

            print('Generating free wak: {}/{}'.format(self._mesh.wk.sec, self._mesh.wk.nw), end='\r')
            
            # Calculate c_ik, and d_i
            d_i[:] = .0
            for i in range(self._mesh.wk.nte * self._mesh.wk.sec):
                if i < self._mesh.wk.nte:
                    _, c_ik[:, i] = potential_wrapper.main(self._mesh.wk.fc[i, 0], self._mesh.wk.p1[i, :], self._mesh.wk.p2[i, :], self._mesh.wk.p3[i, :], self._mesh.wk.p4[i, :], self._mesh.wk.e1[i, :], self._mesh.wk.e2[i, :], self._mesh.wk.e3[i, :], self._mesh.wk.p_avg[i, :], self._mesh.sf.p_ctrl[:, :])
                else:
                    _, d = potential_wrapper.main(self._mesh.wk.fc[i, 0], self._mesh.wk.p1[i, :], self._mesh.wk.p2[i, :], self._mesh.wk.p3[i, :], self._mesh.wk.p4[i, :], self._mesh.wk.e1[i, :], self._mesh.wk.e2[i, :], self._mesh.wk.e3[i, :], self._mesh.wk.p_avg[i, :], self._mesh.sf.p_ctrl[:, :])
                    d_i[:] = d_i[:] + d[:] * self.doublet_wk[i] * areas[i] / self._mesh.wk.areas[i]

            # Create linear system
            rhs[:] = - np.dot(a_ij[:, :], self.source_sf) - d_i[:]
            for i in range(self._mesh.sf.nf):
                lhs[i, :] = b_ij[i, :]
            
            for i in range(self._mesh.wk.nte):
                lhs[:, self._mesh.sf.te_fc[i, 0]] = lhs[:, self._mesh.sf.te_fc[i, 0]] + c_ik[:, i]
                lhs[:, self._mesh.sf.te_fc[i, 1]] = lhs[:, self._mesh.sf.te_fc[i, 1]] - c_ik[:, i]

            # Calculate doublet
            self.doublet_sf[:] = np.linalg.solve(lhs, rhs)
            self.doublet_wk[:self._mesh.wk.nte] = self.doublet_sf[self._mesh.sf.te_fc[:, 0]] - self.doublet_sf[self._mesh.sf.te_fc[:, 1]]

            self._surface_parameters(u_call, fc_cn)

            # Add layer
            if self._mesh.wk.sec < self._mesh.wk.nw:

                # Wake velocity
                v = np.zeros((self._mesh.wk.nv_te * (self._mesh.wk.sec + 1), 3), dtype=np.double)

                # Trailing edge and freestream
                v_sum = np.empty(3, dtype=np.double)
                for i in range(self._mesh.wk.nv_te * (self._mesh.wk.sec + 1)):
                    if i < self._mesh.wk.nv_te:
                        id_v = self._mesh.sf.te_vt[i]
                        v_sum[:] = .0
                        count = 0
                        for j in range(self._mesh.sf.nte):
                            if id_v in self._mesh.sf.te[j, :]:
                                count += 1
                                v_sum[:] = v_sum[:] + self.vel_sf[self._mesh.sf.te_fc[j, 0], :] + self.vel_sf[self._mesh.sf.te_fc[j, 1], :]
                            if count == 2:
                                break
                        v[i, :] = v_sum[:] / (2 * count)
                    else:
                        v[i, :] = u_call(self._mesh.wk.vt[i, :])
                        
                # Surface
                for i in range(self._mesh.sf.nf):
                    s, d = velocity_wrapper.main(self._mesh.sf.fc[i, 0], self._mesh.sf.p1[i, :], self._mesh.sf.p2[i, :], self._mesh.sf.p3[i, :], self._mesh.sf.p4[i, :], self._mesh.sf.e1[i, :], self._mesh.sf.e2[i, :], self._mesh.sf.e3[i, :], self._mesh.sf.p_avg[i, :], self._mesh.wk.vt[self._mesh.wk.nv_te:self._mesh.wk.nv_te * (self._mesh.wk.sec + 1), :])
                    v[self._mesh.wk.nv_te:, :] = v[self._mesh.wk.nv_te:, :] + s[:, :] * self.source_sf[i] + d[:, :] * self.doublet_sf[i]

                # Wake
                for i in range(self._mesh.wk.nte * self._mesh.wk.sec):
                    d = velocity_wrapper.quad_vortex_panel(self._mesh.wk.vt[self._mesh.wk.fc[i, 1], :], self._mesh.wk.vt[self._mesh.wk.fc[i, 2], :], self._mesh.wk.vt[self._mesh.wk.fc[i, 3], :], self._mesh.wk.vt[self._mesh.wk.fc[i, 4], :], self._mesh.wk.vt[self._mesh.wk.nv_te:self._mesh.wk.nv_te * (self._mesh.wk.sec + 1), :])
                    v[self._mesh.wk.nv_te:, :] = v[self._mesh.wk.nv_te:, :] + d[:] * self.doublet_wk[i] * areas[i] / self._mesh.wk.areas[i]

                # Update wake
                self._mesh.wk.update(v)
                areas[:self._mesh.wk.nte] = self._mesh.wk.areas[:self._mesh.wk.nte]
                self._update_wake(areas)
            
            else:
                break
        
        return

    def _update_wake(self, areas: np.ndarray) -> None:
        for i in range(self._mesh.wk.sec - 1):
            self.doublet_wk[self._mesh.wk.nte * (self._mesh.wk.sec - 1 - i):self._mesh.wk.nte * (self._mesh.wk.sec - i)] = self.doublet_wk[self._mesh.wk.nte * (self._mesh.wk.sec - 2 - i):self._mesh.wk.nte * (self._mesh.wk.sec - 1 - i)]
            areas[self._mesh.wk.nte * (self._mesh.wk.sec - 1 - i):self._mesh.wk.nte * (self._mesh.wk.sec - i)] = areas[self._mesh.wk.nte * (self._mesh.wk.sec - 2 - i):self._mesh.wk.nte * (self._mesh.wk.sec - 1 - i)]
        return

    def _faces_connection(self) -> tp.List[FacesConnection]:

        def _find_faces_with_vertices(face: int, id1: int, id2: int) -> tp.List[int]:
            if ((id1 in self._mesh.sf.te[:, 0]) or (id1 in self._mesh.sf.te[:, 1])) and ((id2 in self._mesh.sf.te[:, 0]) or (id2 in self._mesh.sf.te[:, 1])):
                return [face]
            else:
                check_1 = (id1 == self._mesh.sf.fc[:, 1]) | (id1 == self._mesh.sf.fc[:, 2]) | (id1 == self._mesh.sf.fc[:, 3]) | (id1 == self._mesh.sf.fc[:, 4])
                check_2 = (id2 == self._mesh.sf.fc[:, 1]) | (id2 == self._mesh.sf.fc[:, 2]) | (id2 == self._mesh.sf.fc[:, 3]) | (id2 == self._mesh.sf.fc[:, 4])
                ids = np.argwhere(check_1 & check_2)
                return[id[0] for id in ids]
        
        connections = []

        for face in range(self._mesh.sf.nf):

            is_quad = True if self._mesh.sf.fc[face, 0] == 4 else False

            c1 = _find_faces_with_vertices(face, self._mesh.sf.fc[face, 1], self._mesh.sf.fc[face, 2])
            c2 = _find_faces_with_vertices(face, self._mesh.sf.fc[face, 2], self._mesh.sf.fc[face, 3])
            c3 = _find_faces_with_vertices(face, self._mesh.sf.fc[face, 3], self._mesh.sf.fc[face, 4])
            c4 = _find_faces_with_vertices(face, self._mesh.sf.fc[face, 4], self._mesh.sf.fc[face, 1]) if is_quad else None

            connections.append(FacesConnection(self._mesh.sf.fc[face, 0], c1, c2,c3, c4))

        return connections

    def _surface_parameters(self, u_call: tp.Callable[[np.ndarray], np.ndarray], fc_cn: tp.List[FacesConnection]) -> None:

        us = np.asarray([u_call(x) for x in self._mesh.sf.p_avg])

        for face in range(self._mesh.sf.nf):
            
            self.vel_sf[face, :] = us[face] - self._mesh.sf.e3[face, :] * np.dot(us[face], self._mesh.sf.e3[face, :])

            d0 = self.doublet_sf[face]
            d1 = sum([self.doublet_sf[id] for id in fc_cn[face].c1]) / len(fc_cn[face].c1)
            d2 = sum([self.doublet_sf[id] for id in fc_cn[face].c2]) / len(fc_cn[face].c2)
            d3 = sum([self.doublet_sf[id] for id in fc_cn[face].c3]) / len(fc_cn[face].c3)
            d4 = sum([self.doublet_sf[id] for id in fc_cn[face].c4]) / len(fc_cn[face].c4)

            # Correct edges that are connect to only one face
            if len(fc_cn[face].c1) == 1: d1 = d0 + 1.2 * (d0 - d3)
            if len(fc_cn[face].c2) == 1: d2 = d0 + 1.2 * (d0 - d4)
            if len(fc_cn[face].c3) == 1: d3 = d0 + 1.2 * (d0 - d1)
            if len(fc_cn[face].c4) == 1: d4 = d0 + 1.2 * (d0 - d2)

            v1 = np.concatenate([0.5 * (self._mesh.sf.p1[face, :] + self._mesh.sf.p2[face, :]), [d1 - d0]])
            v2 = np.concatenate([0.5 * (self._mesh.sf.p2[face, :] + self._mesh.sf.p3[face, :]), [d2 - d0]])
            v3 = np.concatenate([0.5 * (self._mesh.sf.p3[face, :] + self._mesh.sf.p4[face, :]), [d3 - d0]])
            v4 = np.concatenate([0.5 * (self._mesh.sf.p4[face, :] + self._mesh.sf.p1[face, :]), [d4 - d0]])

            n1 = np.cross(v1, v2)
            n2 = np.cross(v2, v3)
            n3 = np.cross(v3, v4)
            n4 = np.cross(v4, v1)

            n = n1 + n2 + n3 + n4

            self.vel_sf[face, :] = self.vel_sf[face, :] + (n[0] / n[2]) * self._mesh.sf.e1[face, :] + (n[1] / n[2]) * self._mesh.sf.e2[face, :]

        self.cp_sf[:] = 1 - (self.vel_sf[:, 0] * self.vel_sf[:, 0] + self.vel_sf[:, 1] * self.vel_sf[:, 1] * self.vel_sf[:, 2] * self.vel_sf[:, 2]) / (us[:, 0] * us[:, 0] + us[:, 1] * us[:, 1] + us[:, 2] * us[:, 2])

        return

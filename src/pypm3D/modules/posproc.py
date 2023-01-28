import typing as tp
import numpy as np
import vtkmodules.all as vtk


class VerticeModel(tp.NamedTuple):
    n: int
    faces: tp.List[int]
    coefs: tp.List[float]

class _VtpFile:

    def __init__(self) -> None:
        self.nv = 0

        self.pd = vtk.vtkPolyData()
        
        self.points = vtk.vtkPoints()
        self.cells = vtk.vtkCellArray()
        self.lines = vtk.vtkCellArray()

        self.pd.SetPoints(self.points)
        self.pd.SetPolys(self.cells)
        self.pd.SetLines(self.lines)

        self.vt_cns = []

        return

    def add_surface(self,
                    vt: np.ndarray,
                    fc: np.ndarray) -> None:

        nv = self.nv

        for vertice in vt:
            self.points.InsertNextPoint(vertice[0], vertice[1], vertice[2])
            self.nv += 1

        for face in fc:
            if face[0] == 3:
                cell = vtk.vtkTriangle()
                cell.GetPointIds().SetId(0, nv + face[1])
                cell.GetPointIds().SetId(1, nv + face[2])
                cell.GetPointIds().SetId(2, nv + face[3])
            else:
                cell = vtk.vtkQuad()
                cell.GetPointIds().SetId(0, nv + face[1])
                cell.GetPointIds().SetId(1, nv + face[2])
                cell.GetPointIds().SetId(2, nv + face[3])
                cell.GetPointIds().SetId(3, nv + face[4])
            
            self.cells.InsertNextCell(cell) 

        self.vt_cns.append(self._calc_vertices_connection(vt, fc))
        
        return
    
    def add_grid(self,
                 vt: np.ndarray,
                 fc: np.ndarray) -> None:

        nv = self.nv

        for vertice in vt:
            self.points.InsertNextPoint(vertice[0], vertice[1], vertice[2])
            self.nv += 1

        for face in fc:

            polyLine = vtk.vtkPolyLine()

            if face[0] == 3:            
                polyLine.GetPointIds().SetNumberOfIds(4)
                polyLine.GetPointIds().SetId(0, nv + face[1])
                polyLine.GetPointIds().SetId(1, nv + face[2])
                polyLine.GetPointIds().SetId(2, nv + face[3])
                polyLine.GetPointIds().SetId(3, nv + face[1])
            else:
                polyLine.GetPointIds().SetNumberOfIds(5)
                polyLine.GetPointIds().SetId(0, nv + face[1])
                polyLine.GetPointIds().SetId(1, nv + face[2])
                polyLine.GetPointIds().SetId(2, nv + face[3])
                polyLine.GetPointIds().SetId(3, nv + face[4])
                polyLine.GetPointIds().SetId(4, nv + face[1])
            
            self.lines.InsertNextCell(polyLine)
        
        self.vt_cns.append(self._calc_vertices_connection(vt, fc))

        return
    
    def add_param(self, name: str, dim: int, vals: tp.List[np.ndarray]) -> None:

        assert len(vals) == len(self.vt_cns), 'Incorrect number of groups'

        def _scalar(vertices_connection: tp.List[VerticeModel],
                    vals: np.ndarray) -> np.ndarray:

            nv = len(vertices_connection)

            vals_v = np.empty((nv,), dtype=np.double)

            for i in range(nv):
                vals_v[i] = .0
                for j in range(vertices_connection[i].n):
                    vals_v[i] = vals_v[i] + vertices_connection[i].coefs[j] * vals[vertices_connection[i].faces[j]]
            
            return vals_v

        def _vec(vertices_connection: tp.List[VerticeModel],
                 vals: np.ndarray) -> np.ndarray:

            nv = len(vertices_connection)

            vals_v = np.empty((nv, 3), dtype=np.double)

            for i in range(nv):
                vals_v[i, :] = .0
                for j in range(vertices_connection[i].n):
                    vals_v[i, :] = vals_v[i, :] + vertices_connection[i].coefs[j] * vals[vertices_connection[i].faces[j], :]
            
            return vals_v
        
        vt_vals = []

        if dim == 1:
            for i in range(len(vals)):
                vt_vals.append(_scalar(self.vt_cns[i], vals[i]))
        
        if dim == 3:
            for i in range(len(vals)):
                vt_vals.append(_vec(self.vt_cns[i], vals[i]))
        
        
        values = vtk.vtkFloatArray()
        if dim == 3: values.SetNumberOfComponents(3)
        values.SetName(name)

        for vals in vt_vals:
            if dim == 1:
                for x in vals: values.InsertNextTuple1(x)
            
            if dim == 3:
                for x in vals: values.InsertNextTuple3(x[0], x[1], x[2])
        
        self.pd.GetPointData().AddArray(values)

        return

    def write(self, fn: str) -> None:
        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName('{}.vtp'.format(fn))
        writer.SetInputData(self.pd)
        writer.Write()
        return

    def _calc_vertices_connection(self, vt: np.ndarray, fc: np.ndarray) -> tp.List[VerticeModel]:

        vt_cn = []

        for v_id in range(vt.shape[0]):

            # Encontra as faces que se conectam com o v_id
            face_ids = [i[0] for i in np.argwhere((fc[:, 1] == v_id) | (fc[:, 2] == v_id) | (fc[:, 3] == v_id) | (fc[:, 4] == v_id))]

            # Encontra os ângulos
            angles = []

            for face_id in face_ids:
                
                if v_id == fc[face_id, 1]:
                    
                    v1 = vt[fc[face_id, 2], :] - vt[fc[face_id, 1], :]
                    
                    if fc[face_id, 0] == 3:
                        v2 = vt[fc[face_id, 3], :] - vt[fc[face_id, 1], :]
                    else:
                        v2 = vt[fc[face_id, 4], :] - vt[fc[face_id, 1], :]

                elif v_id == fc[face_id, 2]:
                    v1 = vt[fc[face_id, 3], :] - vt[fc[face_id, 2], :]
                    v2 = vt[fc[face_id, 1], :] - vt[fc[face_id, 2], :]

                elif v_id == fc[face_id, 3]:

                    v2 = vt[fc[face_id, 2], :] - vt[fc[face_id, 3], :]

                    if fc[face_id, 0] == 3:
                        v1 = vt[fc[face_id, 1], :] - vt[fc[face_id, 3], :]
                    else:
                        v1 = vt[fc[face_id, 4], :] - vt[fc[face_id, 3], :]

                else:
                    v1 = vt[fc[face_id, 1], :] - vt[fc[face_id, 4], :]
                    v2 = vt[fc[face_id, 3], :] - vt[fc[face_id, 4], :]
                
                angles.append(np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))))

            sum_angles = sum(angles)

            # Adiciona um modelo
            vt_cn.append(
                VerticeModel(
                    n=len(face_ids),
                    faces=face_ids,
                    coefs=[a / sum_angles for a in angles],
                )
            )
            
        return vt_cn

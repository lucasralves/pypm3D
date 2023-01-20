import numpy as np
import vtkmodules.all as vtk

from pypm3D.models.mesh_model import MeshModel
from pypm3D.models.solution_model import SolutionModel
from pypm3D.modules.posproc.utils import vtk_file

def create_vtp_file(filename: str,
                    mesh: MeshModel,
                    sol: SolutionModel = None,
                    show_wake: bool = None) -> None:
    
    pd = vtk.vtkPolyData()

    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()
    lines = vtk.vtkCellArray()

    pd.SetPoints(points)
    pd.SetPolys(cells)
    pd.SetLines(lines)

    vtk_file.add_surface(0, mesh.surface.vertices, mesh.surface.faces, points, cells)
    
    if show_wake is not None:
        check = show_wake # mesh.wake is not None
    else:
        check = mesh.wake is not None
    
    if check:
        vtk_file.add_grid(mesh.surface.nv, mesh.wake.vertices, mesh.wake.faces, points, lines)
    
    if sol is not None:
        source = vtk_file.add_scalar_value('source', np.concatenate([sol.source_v, np.mean(sol.source) * np.ones((mesh.wake.vertices.shape[0],), dtype=np.double)]) if check else sol.source_v)
        pd.GetPointData().AddArray(source)

        doublet = vtk_file.add_scalar_value('doublet', np.concatenate([sol.doublet_v, np.mean(sol.doublet) * np.ones((mesh.wake.vertices.shape[0],), dtype=np.double)]) if check else sol.doublet_v)
        pd.GetPointData().AddArray(doublet)

        cp = vtk_file.add_scalar_value('cp', np.concatenate([sol.cp_v, np.mean(sol.cp) * np.ones((mesh.wake.vertices.shape[0],), dtype=np.double)]) if check else sol.cp_v)
        pd.GetPointData().AddArray(cp)

        transpiration = vtk_file.add_scalar_value('transpiration', np.concatenate([sol.transpiration_v, np.mean(sol.transpiration) * np.ones((mesh.wake.vertices.shape[0],), dtype=np.double)]) if check else sol.transpiration_v)
        pd.GetPointData().AddArray(transpiration)

        potential = vtk_file.add_scalar_value('potential', np.concatenate([sol.potential_v, np.mean(sol.potential) * np.ones((mesh.wake.vertices.shape[0],), dtype=np.double)]) if check else sol.potential_v)
        pd.GetPointData().AddArray(potential)

        vel_mean = np.mean(sol.vel, axis=0)
        vel_aux = np.ones((mesh.wake.vertices.shape[0], 3), dtype=np.double)
        vel_aux[:, 0] = vel_aux[:, 0] * vel_mean[0]
        vel_aux[:, 1] = vel_aux[:, 1] * vel_mean[1]
        vel_aux[:, 2] = vel_aux[:, 2] * vel_mean[2]
        vel = vtk_file.add_vector3D_value('velocity', np.concatenate([sol.vel_v, vel_aux], axis=0) if check else sol.vel_v)
        pd.GetPointData().AddArray(vel)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName('{}.vtp'.format(filename))
    writer.SetInputData(pd)
    writer.Write()

    return
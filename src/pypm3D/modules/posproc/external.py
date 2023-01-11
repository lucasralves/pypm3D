from pypm3D import models
from pypm3D.modules.posproc import core

import vtkmodules.all as vtk

def create_vtp_file(filename: str, mesh: models.MeshModel, aero: models.AeroModel = None) -> None:

    pd = vtk.vtkPolyData()
    
    # Surface
    points, cells = core.vtk_create_opaque_surface(mesh.surface.vertices, mesh.surface.faces)
    lines = core.vtk_create_grid_surface(mesh.surface.nv, mesh.internal.vertices, mesh.internal.faces, points)

    # Wake
    if aero is not None:
        core.vtk_create_grid_surface(mesh.surface.nv + mesh.internal.nv, aero.wake.vertices, aero.wake.faces, points, lines=lines)
        arrays = core.vtk_add_vertices_values(mesh, aero)

        for array in arrays: pd.GetPointData().AddArray(array)

    pd.SetPoints(points)
    pd.SetPolys(cells)
    pd.SetLines(lines)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName('{}.vtp'.format(filename))
    writer.SetInputData(pd)
    writer.Write()

    return
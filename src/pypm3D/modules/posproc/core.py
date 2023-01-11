import typing as tp
import numpy as np
import vtkmodules.all as vtk

from pypm3D import models


def vtk_create_opaque_surface(vertices: np.ndarray,
                              faces: np.ndarray) -> tp.List[tp.Union[vtk.vtkPoints, vtk.vtkCellArray]]:
    
    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()
    
    for vertice in vertices:
        points.InsertNextPoint(vertice[0], vertice[1], vertice[2])

    for face in faces:
        if face[0] == 3:
            cell = vtk.vtkTriangle()
            cell.GetPointIds().SetId(0, face[1])
            cell.GetPointIds().SetId(1, face[2])
            cell.GetPointIds().SetId(2, face[3])
        else:
            cell = vtk.vtkQuad()
            cell.GetPointIds().SetId(0, face[1])
            cell.GetPointIds().SetId(1, face[2])
            cell.GetPointIds().SetId(2, face[3])
            cell.GetPointIds().SetId(3, face[4])
        
        cells.InsertNextCell(cell)
    
    return points, cells

def vtk_create_grid_surface(nv: int,
                            vertices: np.ndarray,
                            faces: np.ndarray,
                            points: vtk.vtkPoints,
                            lines: vtk.vtkCellArray = None) -> vtk.vtkCellArray | None:

    check = lines is None

    if check:
        lines = vtk.vtkCellArray()

    for vertice in vertices:
        points.InsertNextPoint(vertice[0], vertice[1], vertice[2])
    
    for face in faces:
        if face[0] == 3:

            polyLine = vtk.vtkPolyLine()
            polyLine.GetPointIds().SetNumberOfIds(4)
            polyLine.GetPointIds().SetId(0, nv + face[1])
            polyLine.GetPointIds().SetId(1, nv + face[2])
            polyLine.GetPointIds().SetId(2, nv + face[3])
            polyLine.GetPointIds().SetId(3, nv + face[1])
            lines.InsertNextCell(polyLine)
        else:
            polyLine = vtk.vtkPolyLine()
            polyLine.GetPointIds().SetNumberOfIds(5)
            polyLine.GetPointIds().SetId(0, nv + face[1])
            polyLine.GetPointIds().SetId(1, nv + face[2])
            polyLine.GetPointIds().SetId(2, nv + face[3])
            polyLine.GetPointIds().SetId(3, nv + face[4])
            polyLine.GetPointIds().SetId(4, nv + face[1])
            lines.InsertNextCell(polyLine)
    
    if check:
        return lines
    else:
        return None

def vtk_add_vertices_values(mesh: models.MeshModel,
                            aero: models.AeroModel) -> tp.List[vtk.vtkFloatArray]:

    source = vtk.vtkFloatArray(); source.SetName('source')
    vel = vtk.vtkFloatArray(); vel.SetNumberOfComponents(3); vel.SetName('vel')
    cp = vtk.vtkFloatArray(); cp.SetName('cp')
    transp = vtk.vtkFloatArray(); transp.SetName('transp')

    for i in range(mesh.surface.nv):
        source.InsertNextTuple1(aero.surface.source_v[i])
        vel.InsertNextTuple3(aero.surface.vel_v[i, 0], aero.surface.vel_v[i, 1], aero.surface.vel_v[i, 2])
        cp.InsertNextTuple1(aero.surface.cp_v[i])
        transp.InsertNextTuple1(aero.surface.transpiration_v[i])
    
    mean_source = np.mean(aero.surface.source)
    mean_vel = np.mean(aero.surface.vel, axis=0)
    mean_cp = np.mean(aero.surface.cp)
    mean_transpiration = np.mean(aero.surface.transpiration)

    for i in range(mesh.internal.nv + aero.wake.nv * aero.wake.nw + aero.wake.nv):
        source.InsertNextTuple1(mean_source)
        vel.InsertNextTuple3(mean_vel[0], mean_vel[1], mean_vel[2])
        cp.InsertNextTuple1(mean_cp)
        transp.InsertNextTuple1(mean_transpiration)

    return [source, vel, cp, transp]
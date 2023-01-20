import numpy as np
import vtkmodules.all as vtk


def add_surface(nv: int,
                vertices: np.ndarray,
                faces: np.ndarray,
                points: vtk.vtkPoints,
                cells: vtk.vtkCellArray) -> None:

    for vertice in vertices:
        points.InsertNextPoint(vertice[0], vertice[1], vertice[2])

    for face in faces:
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
        
        cells.InsertNextCell(cell)

    return

def add_grid(nv: int,
             vertices: np.ndarray,
             faces: np.ndarray,
             points: vtk.vtkPoints,
             lines: vtk.vtkCellArray) -> None:

    for vertice in vertices:
        points.InsertNextPoint(vertice[0], vertice[1], vertice[2])

    for face in faces:

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
        
        lines.InsertNextCell(polyLine)

    return

def add_scalar_value(name: str, vals: np.ndarray) -> None:

    values = vtk.vtkFloatArray()
    values.SetName(name)

    for x in vals:
        values.InsertNextTuple1(x)

    return values

def add_vector3D_value(name: str, vals: np.ndarray) -> None:

    values = vtk.vtkFloatArray()
    values.SetNumberOfComponents(3)
    values.SetName(name)

    for x in vals:
        values.InsertNextTuple3(x[0], x[1], x[2])

    return values
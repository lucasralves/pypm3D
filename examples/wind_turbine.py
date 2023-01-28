"""
This script solves the potential flow of a rectangular wing with a cross
section of a NACA 4 digits airfoil.
"""

import sys
sys.path.append('./src/')

import pypm3D

import numpy as np
import typing as tp
import gmsh as gm


#---------------------------------------#
#   Geometry
#---------------------------------------#
def wind_turbine(n_span: int, n_chord: int, show: bool = False):
    """
    It creates a wind turbine mesh with a symmetrical Naca 4 digits airfoil.

    Parameters:
    -----------
    - n_span: number of points in the span direction
    - n_chord: number of points in the chord

    Vertices:
    ---------

    t.e.    1e ------------------ 1 ------------------ 1d
            |                     |                     |
            4e ------------------ 4 ------------------ 1d
            |                     |                     |
    l.e.    3e ------------------ 3 ------------------ 1d
            |                     |                     |
            2e ------------------ 2 ------------------ 1d
            |                     |                     |
    t.e     1e ------------------ 1 ------------------ 1d
    """

    def correct_vertices_ids(vertices: np.ndarray, faces3: np.ndarray, faces4: np.ndarray, trailing_edge_list: tp.List[np.ndarray]):
    
        # Saída
        vertices_out = []
        faces_out = []
        trailing_edge_out = []

        # Encontre os ids dos vértices utilizados na malha
        vertices_ids = []

        for id in range(vertices.shape[0]):

            is_in_f3 = id in faces3[:, 0] or id in faces3[:, 1] or id in faces3[:, 2]
            is_in_f4 = id in faces4[:, 0] or id in faces4[:, 1] or id in faces4[:, 2] or id in faces4[:, 3]

            if is_in_f3 or is_in_f4:
                vertices_ids.append(id)
        
        vertices_ids = np.asarray(vertices_ids)

        # Corrige os valores dos vértices
        for id in vertices_ids:
            vertices_out.append(vertices[id, :])

        # Corrige os valores das faces
        for face in faces4:
            id1 = int(np.argwhere(face[0] == vertices_ids)[0])
            id2 = int(np.argwhere(face[1] == vertices_ids)[0])
            id3 = int(np.argwhere(face[2] == vertices_ids)[0])
            id4 = int(np.argwhere(face[3] == vertices_ids)[0])
            faces_out.append([4, id1, id2, id3, id4])
        
        for face in faces3:
            id1 = int(np.argwhere(face[0] == vertices_ids)[0])
            id2 = int(np.argwhere(face[1] == vertices_ids)[0])
            id3 = int(np.argwhere(face[2] == vertices_ids)[0])
            faces_out.append([3, id1, id2, id3, -1])
        
        # Corrige os valores do bordo de fuga
        for points in trailing_edge_list:

            d1 = np.linalg.norm(vertices[points[0], :] - vertices[points[2], :])
            d2 = np.linalg.norm(vertices[points[0], :] - vertices[points[-1], :])

            if d1 < d2:
                points_ordered = points[2:]
            else:
                points_ordered = np.flip(points[2:])
            
            new_points = [points[0]]
            for i in range(len(points_ordered)):
                new_points.append(points_ordered[i])
            new_points.append(points[1])
            
            for i in range(len(new_points) - 1):
                id1 = int(np.argwhere(new_points[i] == vertices_ids)[0])
                id2 = int(np.argwhere(new_points[i + 1] == vertices_ids)[0])
                trailing_edge_out.append([id1, id2])
        
        # Converte para numpy array
        vertices_out = np.asarray(vertices_out)
        faces_out = np.asarray(faces_out)
        trailing_edge_out = np.asarray(trailing_edge_out)

        # Flip trailing edge vertices if necessary
        for i in range(trailing_edge_out.shape[0]):
            if vertices_out[trailing_edge_out[i, 1], 1] - vertices_out[trailing_edge_out[i, 0], 1] < 0:
                id1 = trailing_edge_out[i, 0]
                id2 = trailing_edge_out[i, 1]
                trailing_edge_out[i, 0] = id2
                trailing_edge_out[i, 1] = id1
        
        return [vertices_out, faces_out, trailing_edge_out]

    def naca_point(thickness: float, x: float) -> np.ndarray:
        return 5 * thickness * (0.2969 * np.power(x, 0.5) - 0.1260 * x - 0.3516 * np.power(x, 2) + 0.2843 * np.power(x, 3) - 0.1015 * np.power(x, 4)) - 0.010499999999999815 * x * thickness
    
    th: float = 0.12
    span: float = 5.
    y_add = 1. + 0.5 * span

    root = 1.3
    tip = 0.7

    gm.initialize()

    gm.option.setNumber('General.Verbosity', 1)

    # Points
    p1 = gm.model.geo.add_point(1.0 - 0.5, y_add, naca_point(th, 1.0))
    p2 = gm.model.geo.add_point(0.5 - 0.5, y_add, naca_point(th, 0.5))
    p3 = gm.model.geo.add_point(0.0 - 0.5, y_add, naca_point(th, 0.0))
    p4 = gm.model.geo.add_point(0.5 - 0.5, y_add,-naca_point(th, 0.5))
    
    p1e = gm.model.geo.add_point(root * (1.0 - 0.5), y_add -0.5 * span, root * (naca_point(th, 1.0)))
    p2e = gm.model.geo.add_point(root * (0.5 - 0.5), y_add -0.5 * span, root * (naca_point(th, 0.5)))
    p3e = gm.model.geo.add_point(root * (0.0 - 0.5), y_add -0.5 * span, root * (naca_point(th, 0.0)))
    p4e = gm.model.geo.add_point(root * (0.5 - 0.5), y_add -0.5 * span,root * (-naca_point(th, 0.5)))

    p1d = gm.model.geo.add_point(tip * (1.0 - 0.5),y_add + 0.5 * span, tip * (naca_point(th, 1.0)))
    p2d = gm.model.geo.add_point(tip * (0.5 - 0.5),y_add + 0.5 * span, tip * (naca_point(th, 0.5)))
    p3d = gm.model.geo.add_point(tip * (0.0 - 0.5),y_add + 0.5 * span, tip * (naca_point(th, 0.0)))
    p4d = gm.model.geo.add_point(tip * (0.5 - 0.5),y_add + 0.5 * span,tip * (-naca_point(th, 0.5)))

    # Foils curves
    foil_1 = np.asarray([[x - 0.5, naca_point(th, x)] for x in np.linspace(1.0, 0.5, num=100)])[1:-1, :]
    foil_2 = np.asarray([[x - 0.5, naca_point(th, x)] for x in np.geomspace(0.6, 0.1, num=200) - 0.1])[1:-1, :]
    foil_3 = np.asarray([[x - 0.5,-naca_point(th, x)] for x in np.flip(np.geomspace(0.6, 0.1, num=200) - 0.1)])[1:-1, :]
    foil_4 = np.asarray([[x - 0.5,-naca_point(th, x)] for x in np.linspace(0.5, 1.0, num=100)])[1:-1, :]

    # Curves
    c1 = gm.model.geo.add_polyline([p1] + [gm.model.geo.add_point(x[0], y_add, x[1]) for x in foil_1] + [p2])
    c2 = gm.model.geo.add_polyline([p2] + [gm.model.geo.add_point(x[0], y_add, x[1]) for x in foil_2] + [p3])
    c3 = gm.model.geo.add_polyline([p3] + [gm.model.geo.add_point(x[0], y_add, x[1]) for x in foil_3] + [p4])
    c4 = gm.model.geo.add_polyline([p4] + [gm.model.geo.add_point(x[0], y_add, x[1]) for x in foil_4] + [p1])

    c1e = gm.model.geo.add_polyline([p1e] + [gm.model.geo.add_point(root * x[0],y_add-0.5 * span, root * x[1]) for x in foil_1] + [p2e])
    c2e = gm.model.geo.add_polyline([p2e] + [gm.model.geo.add_point(root * x[0],y_add-0.5 * span, root * x[1]) for x in foil_2] + [p3e])
    c3e = gm.model.geo.add_polyline([p3e] + [gm.model.geo.add_point(root * x[0],y_add-0.5 * span, root * x[1]) for x in foil_3] + [p4e])
    c4e = gm.model.geo.add_polyline([p4e] + [gm.model.geo.add_point(root * x[0],y_add-0.5 * span, root * x[1]) for x in foil_4] + [p1e])

    c1d = gm.model.geo.add_polyline([p1d] + [gm.model.geo.add_point(tip * x[0],y_add+0.5 * span, tip * x[1]) for x in foil_1] + [p2d])
    c2d = gm.model.geo.add_polyline([p2d] + [gm.model.geo.add_point(tip * x[0],y_add+0.5 * span, tip * x[1]) for x in foil_2] + [p3d])
    c3d = gm.model.geo.add_polyline([p3d] + [gm.model.geo.add_point(tip * x[0],y_add+0.5 * span, tip * x[1]) for x in foil_3] + [p4d])
    c4d = gm.model.geo.add_polyline([p4d] + [gm.model.geo.add_point(tip * x[0],y_add+0.5 * span, tip * x[1]) for x in foil_4] + [p1d])

    c11e = gm.model.geo.add_line(p1, p1e)
    c22e = gm.model.geo.add_line(p2, p2e)
    c33e = gm.model.geo.add_line(p3, p3e)
    c44e = gm.model.geo.add_line(p4, p4e)

    c11d = gm.model.geo.add_line(p1, p1d)
    c22d = gm.model.geo.add_line(p2, p2d)
    c33d = gm.model.geo.add_line(p3, p3d)
    c44d = gm.model.geo.add_line(p4, p4d)

    # Curve loops
    cl1e = gm.model.geo.add_curve_loop([c1, c22e, -c1e, -c11e])
    cl2e = gm.model.geo.add_curve_loop([c2, c33e, -c2e, -c22e])
    cl3e = gm.model.geo.add_curve_loop([c3, c44e, -c3e, -c33e])
    cl4e = gm.model.geo.add_curve_loop([c4, c11e, -c4e, -c44e])

    cl1d = gm.model.geo.add_curve_loop([-c1, c11d, c1d, -c22d])
    cl2d = gm.model.geo.add_curve_loop([-c2, c22d, c2d, -c33d])
    cl3d = gm.model.geo.add_curve_loop([-c3, c33d, c3d, -c44d])
    cl4d = gm.model.geo.add_curve_loop([-c4, c44d, c4d, -c11d])

    # Surfaces
    s1e = gm.model.geo.add_surface_filling([cl1e])
    s2e = gm.model.geo.add_surface_filling([cl2e])
    s3e = gm.model.geo.add_surface_filling([cl3e])
    s4e = gm.model.geo.add_surface_filling([cl4e])

    s1d = gm.model.geo.add_surface_filling([cl1d])
    s2d = gm.model.geo.add_surface_filling([cl2d])
    s3d = gm.model.geo.add_surface_filling([cl3d])
    s4d = gm.model.geo.add_surface_filling([cl4d])

    # Synchronize
    gm.model.geo.synchronize()

    # Set transfinite
    gm.model.mesh.set_transfinite_curve(c1, int(n_chord / 4), 'Progress', 1.2)
    gm.model.mesh.set_transfinite_curve(c2, int(n_chord / 4), 'Progress',-1.2)
    gm.model.mesh.set_transfinite_curve(c3, int(n_chord / 4), 'Progress', 1.2)
    gm.model.mesh.set_transfinite_curve(c4, int(n_chord / 4), 'Progress',-1.2)

    gm.model.mesh.set_transfinite_curve(c1e, int(n_chord / 4), 'Progress', 1.2)
    gm.model.mesh.set_transfinite_curve(c2e, int(n_chord / 4), 'Progress',-1.2)
    gm.model.mesh.set_transfinite_curve(c3e, int(n_chord / 4), 'Progress', 1.2)
    gm.model.mesh.set_transfinite_curve(c4e, int(n_chord / 4), 'Progress',-1.2)

    gm.model.mesh.set_transfinite_curve(c1d, int(n_chord / 4), 'Progress', 1.2)
    gm.model.mesh.set_transfinite_curve(c2d, int(n_chord / 4), 'Progress',-1.2)
    gm.model.mesh.set_transfinite_curve(c3d, int(n_chord / 4), 'Progress', 1.2)
    gm.model.mesh.set_transfinite_curve(c4d, int(n_chord / 4), 'Progress',-1.2)

    gm.model.mesh.set_transfinite_curve(c11e, int(n_span / 2), 'Progress',-1.2)
    gm.model.mesh.set_transfinite_curve(c22e, int(n_span / 2), 'Progress',-1.2)
    gm.model.mesh.set_transfinite_curve(c33e, int(n_span / 2), 'Progress',-1.2)
    gm.model.mesh.set_transfinite_curve(c44e, int(n_span / 2), 'Progress',-1.2)

    gm.model.mesh.set_transfinite_curve(c11d, int(n_span / 2), 'Progress',-1.2)
    gm.model.mesh.set_transfinite_curve(c22d, int(n_span / 2), 'Progress',-1.2)
    gm.model.mesh.set_transfinite_curve(c33d, int(n_span / 2), 'Progress',-1.2)
    gm.model.mesh.set_transfinite_curve(c44d, int(n_span / 2), 'Progress',-1.2)
    
    gm.model.mesh.set_transfinite_surface(s1e)
    gm.model.mesh.set_transfinite_surface(s2e)
    gm.model.mesh.set_transfinite_surface(s3e)
    gm.model.mesh.set_transfinite_surface(s4e)

    gm.model.mesh.set_transfinite_surface(s1d)
    gm.model.mesh.set_transfinite_surface(s2d)
    gm.model.mesh.set_transfinite_surface(s3d)
    gm.model.mesh.set_transfinite_surface(s4d)
    
    # Recombine
    gm.model.mesh.set_recombine(2, s1e)
    gm.model.mesh.set_recombine(2, s2e)
    gm.model.mesh.set_recombine(2, s3e)
    gm.model.mesh.set_recombine(2, s4e)

    gm.model.mesh.set_recombine(2, s1d)
    gm.model.mesh.set_recombine(2, s2d)
    gm.model.mesh.set_recombine(2, s3d)
    gm.model.mesh.set_recombine(2, s4d)

    # Create
    gm.model.mesh.generate(2)

    if '-nopopup' not in sys.argv and show:
        gm.fltk.run()
    
    # Vertices
    data = gm.model.mesh.get_nodes()
    vertices = data[1].reshape((data[0].size, 3)).astype(np.double)

    # Faces
    gm.model.mesh.create_faces()
    data = gm.model.mesh.get_all_faces(4)
    faces_4 = data[1].reshape((data[0].size, 4)).astype(np.int32) - 1
    data = gm.model.mesh.get_all_faces(3)
    faces_3 = data[1].reshape((data[0].size, 3)).astype(np.int32) - 1

    faces = np.empty((faces_3.shape[0] + faces_4.shape[0], 5), dtype=np.int32)

    for i in range(faces_4.shape[0]):
        faces[i, 0] = 4
        faces[i, 1] = faces_4[i, 0]
        faces[i, 2] = faces_4[i, 1]
        faces[i, 3] = faces_4[i, 2]
        faces[i, 4] = faces_4[i, 3]

    for i in range(faces_3.shape[0]):
        faces[faces_4.shape[0] + i, 0] = 3
        faces[faces_4.shape[0] + i, 1] = faces_3[i, 0]
        faces[faces_4.shape[0] + i, 2] = faces_3[i, 1]
        faces[faces_4.shape[0] + i, 3] = faces_3[i, 2]
        faces[faces_4.shape[0] + i, 4] = -1

    # Trailing edge
    trailing_edge_1_e = gm.model.add_physical_group(1, [c11e])
    trailing_edge_1_d = gm.model.add_physical_group(1, [c11d])

    trailing_edge_1_e = gm.model.mesh.get_nodes_for_physical_group(1, trailing_edge_1_e)[0] - 1
    trailing_edge_1_d = gm.model.mesh.get_nodes_for_physical_group(1, trailing_edge_1_d)[0] - 1
    
    vertices_out, faces_out, trailing_edge_out = correct_vertices_ids(vertices, faces_3, faces_4, [trailing_edge_1_e, trailing_edge_1_d])

    return [vertices_out, faces_out, trailing_edge_out]

#---------------------------------------#
#   Main
#---------------------------------------#
def u_call(x: np.ndarray) -> np.ndarray:
    """Freestream velocity at a point x"""
    omega = 8.
    return np.array([omega * x[1], -omega * x[0], 2.])

n_span = 20    # number of points
n_chord = 50   # number of points
l = 30.         # wake length
ts = 0.02       # pseudo time step

vt, fc, te = wind_turbine(n_span, n_chord, show=False)

pypm3D.init()
pypm3D.proc_mesh(vt, fc, te, u_call, l, ts)
pypm3D.solve(u_call)
pypm3D.gen_vtk('./examples/wind_turbine', True, u_call)

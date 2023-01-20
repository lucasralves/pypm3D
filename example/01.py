"""
Calculation of Cp distribution over the middle section of a rectangular
wing with a Naca 4 digits symmetrical airfoil.
"""

import sys
sys.path.append('./src/')

import matplotlib.pyplot as plt
import numpy as np
import pypm3D
import mesh

def u_func(x: np.ndarray) -> np.ndarray:
    """Freestream velocity at a point x"""
    alpha = 1.0
    return np.array([np.cos(np.deg2rad(alpha)), .0, np.sin(np.deg2rad(alpha))])

if __name__ == '__main__':

    # mesh
    th = .12       # Naca 4 digits thickness
    span = 5.      # wing's span
    n_span = 20     # number of points
    n_chord = 50    # number of points

    vertices, faces, trailing_edge = mesh.rectangular_wing(th, span, n_span, n_chord, show=False)

    # solution
    length = 2.
    time_step = .2

    msh = pypm3D.mesh.proc_mesh(vertices, faces, trailing_edge)
    sol = pypm3D.aero.solve(msh, u_func, length, time_step, verbose=True)
    pypm3D.posproc.create_vtp_file('./example/01', msh, sol)

    # cp distribution
    xz_vertices = mesh.find_vertices_at_xz_plane(vertices, 0.0)

    plt.figure()
    plt.scatter(msh.surface.vertices[xz_vertices, 0], sol.cp_v[xz_vertices])
    plt.gca().invert_yaxis()
    plt.grid()
    plt.show()

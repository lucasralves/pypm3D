"""
Calculate the forces over a rectangular wing with a Naca 4 digits symmetrical airfoil.
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
    n_chord = 120    # number of points

    vertices, faces, trailing_edge = mesh.rectangular_wing(th, span, n_span, n_chord, show=False)

    # solution
    length = 1.
    time_step = .1
    ref_area = span

    msh = pypm3D.mesh.proc_mesh(vertices, faces, trailing_edge)
    sol = pypm3D.aero.solve(msh, u_func, length, time_step)
    pypm3D.posproc.create_vtp_file('./example/02', msh, sol)
    force = pypm3D.posproc.forces_and_moments(msh, sol, ref_area)

    print(force)
"""
    cp distribution
"""

import sys
sys.path.append('./src/')

import numpy as np
import matplotlib.pyplot as plt
import pypm3D as pm

def get_geo_1():
    """
    vertices:
    0   6   12  18
    5   11  17  23
    4   10  16  22
    3   9   15  21
    2   8   14  20
    1   7   13  19
    0   6   12  18

    faces:
    5   11  17
    4   10  16
    3   9   15
    2   8   14
    1   7   13
    0   6   12
    """

    h0 = 0.
    h1 = 0.05
    h2 = 0.08
    h3 = 0.

    r0 = -0.5
    r1 = -0.166
    r2 = 0.166
    r3 = 0.5

    s0 = -0.5
    s1 = -0.16
    s2 = 0.16
    s3 = 0.5

    vertices = np.array([
        [r0, s0, h0], [r1, s0, h1], [r2, s0, h2], [r3, s0, h3], [r2, s0, -h2], [r1, s0, -h1],
        [r0, s1, h0], [r1, s1, h1], [r2, s1, h2], [r3, s1, h3], [r2, s1, -h2], [r1, s1, -h1],
        [r0, s2, h0], [r1, s2, h1], [r2, s2, h2], [r3, s2, h3], [r2, s2, -h2], [r1, s2, -h1],
        [r0, s3, h0], [r1, s3, h1], [r2, s3, h2], [r3, s3, h3], [r2, s3, -h2], [r1, s3, -h1],
    ], dtype=np.double)

    faces = np.array([
        [4, 0, 1, 7, 6], [4, 1, 2, 8, 7], [4, 2, 3, 9, 8], [4, 3, 4, 10, 9], [4, 4, 5, 11, 10], [4, 5, 0, 6, 11],
        [4, 6, 7, 13, 12], [4, 7, 8, 14, 13], [4, 8, 9, 15, 14], [4, 9, 10, 16, 15], [4, 10, 11, 17, 16], [4, 11, 6, 12, 17],
        [4, 12, 13, 19, 18], [4, 13, 14, 20, 19], [4, 14, 15, 21, 20], [4, 15, 16, 22, 21], [4, 16, 17, 23, 22], [4, 17, 12, 18, 23],
    ], dtype=np.int32)

    trailing_edge = np.array([
        [0, 6],
        [6, 12],
        [12, 18],
    ], dtype=np.int32)

    return vertices, faces, trailing_edge

def get_geo_2():
    """
    vertices:
        0   6   12  18
        5   11  17  23
        4   10  16  22
    24  3   9   15  21  25
        2   8   14  20
        1   7   13  19
        0   6   12  18

    faces:
    23    5   11  17    29
    22    4   10  16    28
    21    3   9   15    27
    20    2   8   14    26
    19    1   7   13    25
    18    0   6   12    24
    """

    h0 = 0.
    h1 = 0.05
    h2 = 0.08
    h3 = 0.

    r0 = -0.5
    r1 = -0.166
    r2 = 0.166
    r3 = 0.5

    s0 = -0.5
    s1 = -0.16
    s2 = 0.16
    s3 = 0.5

    vertices = np.array([
        [r0, s0, h0], [r1, s0, h1], [r2, s0, h2], [r3, s0, h3], [r2, s0, -h2], [r1, s0, -h1],
        [r0, s1, h0], [r1, s1, h1], [r2, s1, h2], [r3, s1, h3], [r2, s1, -h2], [r1, s1, -h1],
        [r0, s2, h0], [r1, s2, h1], [r2, s2, h2], [r3, s2, h3], [r2, s2, -h2], [r1, s2, -h1],
        [r0, s3, h0], [r1, s3, h1], [r2, s3, h2], [r3, s3, h3], [r2, s3, -h2], [r1, s3, -h1],
        [.0, -0.8, .0], [.0, 0.8, .0],
    ], dtype=np.double)

    faces = np.array([
        [4, 0, 1, 7, 6], [4, 1, 2, 8, 7], [4, 2, 3, 9, 8], [4, 3, 4, 10, 9], [4, 4, 5, 11, 10], [4, 5, 0, 6, 11],
        [4, 6, 7, 13, 12], [4, 7, 8, 14, 13], [4, 8, 9, 15, 14], [4, 9, 10, 16, 15], [4, 10, 11, 17, 16], [4, 11, 6, 12, 17],
        [4, 12, 13, 19, 18], [4, 13, 14, 20, 19], [4, 14, 15, 21, 20], [4, 15, 16, 22, 21], [4, 16, 17, 23, 22], [4, 17, 12, 18, 23],
        [3, 0, 24, 1, -1], [3, 1, 24, 2, -1], [3, 2, 24, 3, -1], [3, 3, 24, 4, -1], [3, 4, 24, 5, -1], [3, 5, 24, 0, -1],
        [3, 18, 19, 25, -1], [3, 19, 20, 25, -1], [3, 20, 21, 25, -1], [3, 21, 22, 25, -1], [3, 22, 23, 25, -1], [3, 23, 18, 25, -1],
    ], dtype=np.int32)

    trailing_edge = np.array([
        [24, 0],
        [0, 6],
        [6, 12],
        [12, 18],
        [18, 25],
    ], dtype=np.int32)
    
    return vertices, faces, trailing_edge

def get_geo_3(file: str, span: np.ndarray, section: int):

    foil = np.loadtxt(file)

    nv = foil.shape[0] * span.shape[0]
    nf = foil.shape[0] * (span.shape[0] - 1)
    nte = span.shape[0] - 1

    vertices = np.empty((nv, 3), dtype=np.double)
    faces = np.empty((nf, 5), dtype=np.int32)
    trailing_edge = np.empty((nte, 2), dtype=np.int32)
    faces_plot = np.empty((foil.shape[0]), dtype=np.int32)

    count = 0
    for i in range(span.shape[0]):
        for j in range(foil.shape[0]):
            vertices[count, :] = np.array([foil[j, 0], span[i], foil[j, 1]], dtype=np.double)
            count += 1
    
    count = 0
    for i in range(span.shape[0] - 1):
        for j in range(foil.shape[0]):
            if j == foil.shape[0] - 1:
                id1 = i * foil.shape[0] + j
                id2 = (i + 1) * foil.shape[0] + j
                id3 = (i + 1) * foil.shape[0]
                id4 = i * foil.shape[0]
            else:
                id1 = i * foil.shape[0] + j
                id2 = (i + 1) * foil.shape[0] + j
                id3 = (i + 1) * foil.shape[0] + j + 1
                id4 = i * foil.shape[0] + j + 1
            faces[count, :] = np.array([4, id1, id2, id3, id4], dtype=np.int32)
            count += 1
    
    count = 0
    for i in range(span.shape[0] - 1):
        id1 = i * foil.shape[0]
        id2 = (i + 1) * foil.shape[0]
        trailing_edge[i, :] = np.array([id1, id2], dtype=np.int32)
        count += 1
    
    count = 0
    for i in range(foil.shape[0]):
        faces_plot[i] = i + section * foil.shape[0]

    return vertices, faces, trailing_edge, faces_plot

def plot(faces: np.ndarray, mesh: pm.MeshModel, aero: pm.AeroModel) -> None:

    x = mesh.surface.p_avg[faces, 0]
    y = aero.surface.cp[faces]

    plt.figure()
    plt.title('cp')
    plt.xlabel('chord [m]')
    plt.ylabel('cp')
    plt.plot(x, y)
    plt.gca().invert_yaxis()
    plt.grid()

    plt.figure()
    plt.title('doublet')
    plt.plot(aero.internal.doublet, label='internal')
    plt.plot(aero.wake.circulations[:, 0], label='external')
    plt.grid()


    plt.show()

    return

if __name__ == '__main__':

    def u_func(x: np.ndarray) -> np.ndarray:
        U = 1.0
        alpha = 0.0
        beta = 0.0
        return np.array(
                [
                    U * np.cos(np.deg2rad(beta)) * np.cos(np.deg2rad(alpha)),
                    U * np.sin(np.deg2rad(beta)) * np.cos(np.deg2rad(alpha)),
                    U * np.sin(np.deg2rad(alpha))
                ]
            )

    vt, fc, te, fp = get_geo_3('./src/tests/data/naca0012-100.dat', np.linspace(-1, 1, num=6), 3)

    mesh = pm.mesh.proc_mesh(vt, fc, te)
    aero = pm.aero.solve(mesh, u_func, 1, 0.1)

    pm.posproc.create_vtp_file('./src/tests/test', mesh, aero)

    # faces = np.array([6, 7, 8, 9, 10, 11])

    plot(fp, mesh, aero)

"""
    find_faces_at_trailing_edge
"""

import sys
sys.path.append('./src/')

import numpy as np
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

def get_geo_3():
    """
    vertices:
    0   4   8
    3   7   11
    2   6   10
    1   5   9
    0   4   8

    faces:
    3 7
    2 6
    1 5
    0 4
    """

    h0 = 0.
    h1 = 0.05
    h2 = 0.

    r0 = -0.5
    r1 = 0.0
    r2 = 0.5

    s0 = -0.5
    s1 = 0.0
    s2 = 0.5

    vertices = np.array([
        [r0, s0, h0], [r1, s0, h1], [r2, s0, h2], [r1, s0, -h1],
        [r0, s1, h0], [r1, s1, h1], [r2, s1, h2], [r1, s1, -h1],
        [r0, s2, h0], [r1, s2, h1], [r2, s2, h2], [r1, s2, -h1],
    ], dtype=np.double)

    faces = np.array([
        [4, 0, 1, 5, 4], [4, 1, 2, 6, 5], [4, 2, 3, 7, 6], [4, 3, 0, 4, 7],
        [4, 4, 5, 9, 8], [4, 5, 6, 10, 9], [4, 6, 7, 11, 10], [4, 7, 4, 8, 11],
    ], dtype=np.int32)

    trailing_edge = np.array([
        [0, 4],
        [4, 8],
    ], dtype=np.int32)

    return vertices, faces, trailing_edge

if __name__ == '__main__':

    vt_1, fc_1, te_1 = get_geo_1()
    vt_2, fc_2, te_2 = get_geo_2()
    vt_3, fc_3, te_3 = get_geo_3()

    sol_1 = np.array([[0, 5], [6, 11], [12, 17]])
    sol_2 = np.array([[18, 23], [5, 0], [11, 6], [17, 12], [29, 24]])
    sol_3 = np.array([[0, 3], [7, 4]])

    mesh_1 = pm.mesh.proc_mesh(vt_1, fc_1, te_1)
    mesh_2 = pm.mesh.proc_mesh(vt_2, fc_2, te_2)
    mesh_3 = pm.mesh.proc_mesh(vt_3, fc_3, te_3)

    print(mesh_1.surface.trailing_edge_faces)
    print('')
    print(mesh_2.surface.trailing_edge_faces)
    print('')
    print(mesh_3.surface.trailing_edge_faces)
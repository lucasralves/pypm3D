"""
Process mesh
"""
import sys
sys.path.append('./src/')

import numpy as np
import pypm3D as pm

def _get_geo_1():
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

if __name__ == '__main__':
    vertices, faces, trailing_edge = _get_geo_1()
    mesh = pm.mesh.proc_mesh(vertices, faces, trailing_edge)
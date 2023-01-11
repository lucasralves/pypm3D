"""
    calculate_panels_parameters, create_internal_panels, and calculate_scale_factor
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

def view(mesh: pm.MeshModel) -> None:

    plt.figure()
    plt.plot(mesh.internal.scale)
    plt.grid()

    scale = 0.1
    max_size = max([max(np.abs(mesh.surface.vertices[:, 0])), max(np.abs(mesh.surface.vertices[:, 1])), max(np.abs(mesh.surface.vertices[:, 2]))])

    ax = plt.figure().add_subplot(projection='3d')

    # Surface
    for i in range(mesh.surface.nf):
        
        # faces
        if mesh.surface.faces[i, 0] == 4:
            v1 = mesh.surface.vertices[mesh.surface.faces[i, 1], :]
            v2 = mesh.surface.vertices[mesh.surface.faces[i, 2], :]
            v3 = mesh.surface.vertices[mesh.surface.faces[i, 3], :]
            v4 = mesh.surface.vertices[mesh.surface.faces[i, 4], :]
            ax.plot([v1[0], v2[0], v3[0], v4[0], v1[0]], [v1[1], v2[1], v3[1], v4[1], v1[1]], [v1[2], v2[2], v3[2], v4[2], v1[2]], 'k')
        else:
            v1 = mesh.surface.vertices[mesh.surface.faces[i, 1], :]
            v2 = mesh.surface.vertices[mesh.surface.faces[i, 2], :]
            v3 = mesh.surface.vertices[mesh.surface.faces[i, 3], :]
            ax.plot([v1[0], v2[0], v3[0], v1[0]], [v1[1], v2[1], v3[1], v1[1]], [v1[2], v2[2], v3[2], v1[2]], 'k')

        # Orthogonal base
        ax.plot([mesh.surface.p_avg[i, 0], mesh.surface.p_avg[i, 0] + scale * mesh.surface.e3[i, 0]], [mesh.surface.p_avg[i, 1], mesh.surface.p_avg[i, 1] + scale * mesh.surface.e3[i, 1]], [mesh.surface.p_avg[i, 2], mesh.surface.p_avg[i, 2] + scale * mesh.surface.e3[i, 2]], 'r')
        ax.plot([mesh.surface.p_avg[i, 0], mesh.surface.p_avg[i, 0] + scale * mesh.surface.e1[i, 0]], [mesh.surface.p_avg[i, 1], mesh.surface.p_avg[i, 1] + scale * mesh.surface.e1[i, 1]], [mesh.surface.p_avg[i, 2], mesh.surface.p_avg[i, 2] + scale * mesh.surface.e1[i, 2]], 'g')
        ax.plot([mesh.surface.p_avg[i, 0], mesh.surface.p_avg[i, 0] + scale * mesh.surface.e2[i, 0]], [mesh.surface.p_avg[i, 1], mesh.surface.p_avg[i, 1] + scale * mesh.surface.e2[i, 1]], [mesh.surface.p_avg[i, 2], mesh.surface.p_avg[i, 2] + scale * mesh.surface.e2[i, 2]], 'b')

        # p_avg
        ax.scatter(mesh.surface.p_avg[i, 0], mesh.surface.p_avg[i, 1], mesh.surface.p_avg[i, 2], color='r')

        # p_ctrl
        ax.scatter(mesh.surface.p_ctrl[i, 0], mesh.surface.p_ctrl[i, 1], mesh.surface.p_ctrl[i, 2], color='b')

    # Internal
    for i in range(mesh.internal.nf):
        
        # faces
        if mesh.internal.faces[i, 0] == 4:
            v1 = mesh.internal.vertices[mesh.internal.faces[i, 1], :]
            v2 = mesh.internal.vertices[mesh.internal.faces[i, 2], :]
            v3 = mesh.internal.vertices[mesh.internal.faces[i, 3], :]
            v4 = mesh.internal.vertices[mesh.internal.faces[i, 4], :]
            ax.plot([v1[0], v2[0], v3[0], v4[0], v1[0]], [v1[1], v2[1], v3[1], v4[1], v1[1]], [v1[2], v2[2], v3[2], v4[2], v1[2]], '--k')
        else:
            v1 = mesh.internal.vertices[mesh.internal.faces[i, 1], :]
            v2 = mesh.internal.vertices[mesh.internal.faces[i, 2], :]
            v3 = mesh.internal.vertices[mesh.internal.faces[i, 3], :]
            ax.plot([v1[0], v2[0], v3[0], v1[0]], [v1[1], v2[1], v3[1], v1[1]], [v1[2], v2[2], v3[2], v1[2]], '--k')

    ax.scatter(0, 0, 1.1 * max_size, color='w')
    ax.scatter(0, 0, -1.1 * max_size, color='w')

    ax.scatter(1.1 * max_size, 0, 0, color='w')
    ax.scatter(-1.1 * max_size, 0, 0, color='w')

    ax.scatter(0, 1.1 * max_size, 0, color='w')
    ax.scatter(0, -1.1 * max_size, 0, color='w')

    plt.show()

    return

if __name__ == '__main__':

    vt, fc, te = get_geo_3()

    mesh = pm.mesh.proc_mesh(vt, fc, te)
    
    print('Orthogonal base')
    for i in range(mesh.surface.nf):
        print('[{}, {}, {}]'.format(np.linalg.norm(mesh.surface.e1[i, :]), np.linalg.norm(mesh.surface.e2[i, :]), np.linalg.norm(mesh.surface.e3[i, :])))
    
    print('')

    print('Faces center')
    for i in range(mesh.surface.nf):
        if mesh.surface.faces[i, 0] == 4:
            v1 = mesh.surface.p1[i, :]
            v2 = mesh.surface.p2[i, :]
            v3 = mesh.surface.p3[i, :]
            v4 = mesh.surface.p4[i, :]
            p_avg = (1 / 4) * (v1 + v2 + v3 + v4)
        else:
            v1 = mesh.surface.p1[i, :]
            v2 = mesh.surface.p2[i, :]
            v3 = mesh.surface.p3[i, :]
            p_avg = (1 / 3) * (v1 + v2 + v3)
        print(p_avg)

    # view(mesh)
"""
Verifica se os parâmetros calculados na superfície estão corretos
"""

import sys
sys.path.append('./src/')

import numpy as np
import matplotlib.pyplot as plt
import pypm3D as pm


if __name__ == '__main__':

    """
    Geometria quadrangular:

    > ids:
        0   5   10
        4   9   14
        3   8   13
        2   7   12
        1   6   11
        0   5   10
    """

    vertices = np.array([ [0.0, 0.0, 0.0], [0.5, 0.0, 0.08], [1.0, 0.0, 0.02], [1.0, 0.0, -0.02], [0.5, 0.0, -0.08], [0.0, 0.5, 0.0], [0.5, 0.5, 0.08], [1.0, 0.5, 0.02], [1.0, 0.5, -0.02], [0.5, 0.5, -0.08], [0.0, 1.0, 0.0], [0.5, 1.0, 0.08], [1.0, 1.0, 0.02], [1.0, 1.0, -0.02], [0.5, 1.0, -0.08] ], dtype=np.double)
    faces = np.array([ [4, 0, 1, 6, 5], [4, 1, 2, 7, 6], [4, 2, 3, 8, 7], [4, 3, 4, 9, 8], [4, 4, 0, 5, 9], [4, 5, 6, 11, 10], [4, 6, 7, 12, 11], [4, 7, 8, 13, 12], [4, 8, 9, 14, 13], [4, 9, 5, 10, 14] ], dtype=np.int32)
    trailing_edge = np.array([[0, 5], [5, 10]], dtype=np.int32)

    data = pm.proc_mesh(vertices, faces, trailing_edge)

    ax = plt.figure().add_subplot(projection='3d')
    
    scale = 0.1

    for i in range(data['nf']):

        v1 = data['vt'][data['fc'][i, 1], :]
        v2 = data['vt'][data['fc'][i, 2], :]
        v3 = data['vt'][data['fc'][i, 3], :]
        v4 = data['vt'][data['fc'][i, 4], :]

        ax.plot([v1[0], v2[0], v3[0], v4[0], v1[0]], [v1[1], v2[1], v3[1], v4[1], v1[1]], [v1[2], v2[2], v3[2], v4[2], v1[2]], 'k')
        
        ax.plot([data['p_avg'][i, 0], data['p_avg'][i, 0] + scale * data['e3'][i, 0]], [data['p_avg'][i, 1], data['p_avg'][i, 1] + scale * data['e3'][i, 1]], [data['p_avg'][i, 2], data['p_avg'][i, 2] + scale * data['e3'][i, 2]], 'r')
        ax.plot([data['p_avg'][i, 0], data['p_avg'][i, 0] + scale * data['e1'][i, 0]], [data['p_avg'][i, 1], data['p_avg'][i, 1] + scale * data['e1'][i, 1]], [data['p_avg'][i, 2], data['p_avg'][i, 2] + scale * data['e1'][i, 2]], 'g')
        ax.plot([data['p_avg'][i, 0], data['p_avg'][i, 0] + scale * data['e2'][i, 0]], [data['p_avg'][i, 1], data['p_avg'][i, 1] + scale * data['e2'][i, 1]], [data['p_avg'][i, 2], data['p_avg'][i, 2] + scale * data['e2'][i, 2]], 'b')

    for i in range(data['nte']):

        v1 = data['vt_d'][data['fc_d'][i, 1], :]
        v2 = data['vt_d'][data['fc_d'][i, 2], :]
        v3 = data['vt_d'][data['fc_d'][i, 3], :]
        v4 = data['vt_d'][data['fc_d'][i, 4], :]

        ax.plot([v1[0], v2[0], v3[0], v4[0], v1[0]], [v1[1], v2[1], v3[1], v4[1], v1[1]], [v1[2], v2[2], v3[2], v4[2], v1[2]], '--k')
        

    ax.scatter(0, 0, 0.5, color='w')
    ax.scatter(0, 0, -0.5, color='w')

    plt.show()
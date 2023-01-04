import numpy as np

def proc_mesh(vt: np.ndarray,
              fc: np.ndarray,
              te: np.ndarray) -> None:
    """
    Processa os dados de malha.

    Parâmetros:
    -----------
    - vt: vértices da malha; (nv, 3)
    - fc: ids dos vértices, sendo o primeiro elemento o número de lados; (nf, 5)
    - te: ids dos vértices que compõe as arestas do bordo de fuga; (nte, 2)

    Saída:
    ------
    - nv: númerto de vértices
    - nf: número de faces
    - nte: número de arestas no bordo de fuga
    - nvd: número de vértices das faces internas
    - vt: vertices
    - fc: faces
    - te: arestas no borod de fuga
    - p_avg: ponto médio de cada face
    - p_ctrl: ponto de controle
    - e1: vetor da base ortogonal
    - e2: vetor da base ortogonal
    - e3: vetor da base ortogonal
    - p1: ponto 1 no sistema local
    - p2: ponto 2 no sistema local
    - p3: ponto 3 no sistema local
    - p4: ponto 4 no sistema local
    - te_fc: id das faces do bordo de fuga
    - vt_d: vertices das faces internas
    - fc_d: faces internas
    """

    # Saída
    out = {'nv': vt.shape[0], 'nf': fc.shape[0], 'nte': te.shape[0], 'vt': vt, 'fc': fc, 'te': te}

    # Calcula as faces no bordo de fuga
    out['te_fc'] = np.empty((out['nte'], 2), dtype=np.int32)

    for i in range(out['nte']):
        check1 = (fc[:, 1] == te[i, 0]) | (fc[:, 2] == te[i, 0]) | (fc[:, 3] == te[i, 0]) | (fc[:, 4] == te[i, 0])
        check2 = (fc[:, 1] == te[i, 1]) | (fc[:, 2] == te[i, 1]) | (fc[:, 3] == te[i, 1]) | (fc[:, 4] == te[i, 1])
        index = np.argwhere(check1 & check2)
        out['te_fc'][i, 0] = index[0][0]
        out['te_fc'][i, 1] = index[1][0]
    
    # Parâmetros das faces
    out['p_avg'] = np.empty((out['nf'], 3), dtype=np.double)
    out['p_ctrl'] = np.empty((out['nf'], 3), dtype=np.double)
    out['e1'] = np.empty((out['nf'], 3), dtype=np.double)
    out['e2'] = np.empty((out['nf'], 3), dtype=np.double)
    out['e3'] = np.empty((out['nf'], 3), dtype=np.double)
    out['p1'] = np.empty((out['nf'], 2), dtype=np.double)
    out['p2'] = np.empty((out['nf'], 2), dtype=np.double)
    out['p3'] = np.empty((out['nf'], 2), dtype=np.double)
    out['p4'] = np.empty((out['nf'], 2), dtype=np.double)

    # Geometria do painel
    tri_panel_ids = fc[:, 0] == 3
    quad_panel_ids = fc[:, 0] == 4

    # Centro do painel
    out['p_avg'][tri_panel_ids, :] = (1. / 3.) * (vt[fc[tri_panel_ids, 1], :] + vt[fc[tri_panel_ids, 2], :] + vt[fc[tri_panel_ids, 3], :])
    out['p_avg'][quad_panel_ids, :] = (1. / 4.) * (vt[fc[quad_panel_ids, 1], :] + vt[fc[quad_panel_ids, 2], :] + vt[fc[quad_panel_ids, 3], :] + vt[fc[quad_panel_ids, 4], :])
    
    # Base ortogonal (e3)
    out['e3'][tri_panel_ids, :] = np.cross(vt[fc[tri_panel_ids, 2], :] - vt[fc[tri_panel_ids, 1], :], vt[fc[tri_panel_ids, 3], :] - vt[fc[tri_panel_ids, 1], :])
    out['e3'][quad_panel_ids, :] = np.cross(vt[fc[quad_panel_ids, 2], :] - vt[fc[quad_panel_ids, 4], :], vt[fc[quad_panel_ids, 3], :] - vt[fc[quad_panel_ids, 1], :])

    aux_norm = np.sqrt(out['e3'][:, 0] * out['e3'][:, 0] + out['e3'][:, 1] * out['e3'][:, 1] + out['e3'][:, 2] * out['e3'][:, 2])

    out['e3'][:, 0] = out['e3'][:, 0] / aux_norm
    out['e3'][:, 1] = out['e3'][:, 1] / aux_norm
    out['e3'][:, 2] = out['e3'][:, 2] / aux_norm

    # Base ortogonal (e1)
    out['e1'][:, :] = 0.5 * (vt[fc[:, 1], :] + vt[fc[:, 2], :]) - out['p_avg']

    aux_norm = np.sqrt(out['e1'][:, 0] * out['e1'][:, 0] + out['e1'][:, 1] * out['e1'][:, 1] + out['e1'][:, 2] * out['e1'][:, 2])

    out['e1'][:, 0] = out['e1'][:, 0] / aux_norm
    out['e1'][:, 1] = out['e1'][:, 1] / aux_norm
    out['e1'][:, 2] = out['e1'][:, 2] / aux_norm

    # Base ortogonal (e2)
    out['e2'][:, :] = np.cross(out['e3'], out['e1'])

    aux_norm = np.sqrt(out['e2'][:, 0] * out['e2'][:, 0] + out['e2'][:, 1] * out['e2'][:, 1] + out['e2'][:, 2] * out['e2'][:, 2])

    out['e2'][:, 0] = out['e2'][:, 0] / aux_norm
    out['e2'][:, 1] = out['e2'][:, 1] / aux_norm
    out['e2'][:, 2] = out['e2'][:, 2] / aux_norm

    # Base ortogonal (e3 - correção)
    out['e3'][:, :] = np.cross(out['e1'], out['e2'])

    aux_norm = np.sqrt(out['e3'][:, 0] * out['e3'][:, 0] + out['e3'][:, 1] * out['e3'][:, 1] + out['e3'][:, 2] * out['e3'][:, 2])

    out['e3'][:, 0] = out['e3'][:, 0] / aux_norm
    out['e3'][:, 1] = out['e3'][:, 1] / aux_norm
    out['e3'][:, 2] = out['e3'][:, 2] / aux_norm

    # Corrige as faces do bordo de fuga
    for i in range(out['nte']):

        # Vetor em comum
        aux_vec = out['vt'][out['te'][i, 1], :] - out['vt'][out['te'][i, 0], :]
        aux_norm = np.linalg.norm(aux_vec)
        e1 = aux_vec / aux_norm

        # Face 1
        aux_vec = np.cross(out['e3'][out['te_fc'][i, 0], :], e1)
        aux_norm = np.linalg.norm(aux_vec)
        e2 = aux_vec / aux_norm

        aux_vec = np.cross(e1, e2)
        aux_norm = np.linalg.norm(aux_vec)
        e3 = aux_vec / aux_norm

        out['e1'][out['te_fc'][i, 0], :] = e1[:]
        out['e2'][out['te_fc'][i, 0], :] = e2[:]
        out['e3'][out['te_fc'][i, 0], :] = e3[:]

        # Face 2
        aux_vec = np.cross(out['e3'][out['te_fc'][i, 1], :], e1)
        aux_norm = np.linalg.norm(aux_vec)
        e2 = aux_vec / aux_norm

        aux_vec = np.cross(e1, e2)
        aux_norm = np.linalg.norm(aux_vec)
        e3 = aux_vec / aux_norm

        out['e1'][out['te_fc'][i, 1], :] = e1[:]
        out['e2'][out['te_fc'][i, 1], :] = e2[:]
        out['e3'][out['te_fc'][i, 1], :] = e3[:]

    # Ponto de controle
    out['p_ctrl'][:, :] = out['p_avg'][:, :] + out['e3'][:, :] * 1e-12

    # Vértices locais
    aux_vec = out['vt'][out['fc'][:, 1], :] - out['p_avg']
    out['p1'][:, 0] = out['e1'][:, 0] * aux_vec[:, 0] + out['e1'][:, 1] * aux_vec[:, 1] + out['e1'][:, 2] * aux_vec[:, 2]
    out['p1'][:, 1] = out['e2'][:, 0] * aux_vec[:, 0] + out['e2'][:, 1] * aux_vec[:, 1] + out['e2'][:, 2] * aux_vec[:, 2]

    aux_vec = out['vt'][out['fc'][:, 2], :] - out['p_avg']
    out['p2'][:, 0] = out['e1'][:, 0] * aux_vec[:, 0] + out['e1'][:, 1] * aux_vec[:, 1] + out['e1'][:, 2] * aux_vec[:, 2]
    out['p2'][:, 1] = out['e2'][:, 0] * aux_vec[:, 0] + out['e2'][:, 1] * aux_vec[:, 1] + out['e2'][:, 2] * aux_vec[:, 2]

    aux_vec = out['vt'][out['fc'][:, 3], :] - out['p_avg']
    out['p3'][:, 0] = out['e1'][:, 0] * aux_vec[:, 0] + out['e1'][:, 1] * aux_vec[:, 1] + out['e1'][:, 2] * aux_vec[:, 2]
    out['p3'][:, 1] = out['e2'][:, 0] * aux_vec[:, 0] + out['e2'][:, 1] * aux_vec[:, 1] + out['e2'][:, 2] * aux_vec[:, 2]

    aux_vec = out['vt'][out['fc'][quad_panel_ids, 3], :] - out['p_avg'][quad_panel_ids, :]
    out['p4'][:, 0] = out['e1'][quad_panel_ids, 0] * aux_vec[:, 0] + out['e1'][quad_panel_ids, 1] * aux_vec[:, 1] + out['e1'][quad_panel_ids, 2] * aux_vec[:, 2]
    out['p4'][:, 1] = out['e2'][quad_panel_ids, 0] * aux_vec[:, 0] + out['e2'][quad_panel_ids, 1] * aux_vec[:, 1] + out['e2'][quad_panel_ids, 2] * aux_vec[:, 2]

    # Painéis internos
    def __find_internal_ids(a: np.ndarray, b: np.ndarray, face: int) -> int:

        n = a.shape[0] - 1
        ids = None

        # Find the start id
        for i in range(n):
            for j in range(n):

                if a[i + 1] == b[j + 1]:
                    ids = [i + 1, j + 1]
                    break
                
            if ids is not None:
                break
            
        # Create new arrays
        flip_face = True if out['e3'][out['te_fc'][face, 0], 2] < 0 else False

        if a[0] == 4:

            if ids[0] == 1:
                anew = [a[1], a[4], a[3], a[2]] if flip_face else [a[1], a[2], a[3], a[4]]
            elif ids[0] == 2:
                anew = [a[2], a[1], a[4], a[3]] if flip_face else [a[2], a[3], a[4], a[1]]
            elif ids[0] == 3:
                anew = [a[3], a[2], a[1], a[4]] if flip_face else [a[3], a[4], a[1], a[2]]
            elif ids[0] == 4:
                anew = [a[4], a[3], a[2], a[1]] if flip_face else [a[4], a[1], a[2], a[3]]
                
            flip_face = not flip_face
                
            if ids[1] == 1:
                bnew = [b[1], b[4], b[3], b[2]] if flip_face else [b[1], b[2], b[3], b[4]]
            elif ids[1] == 2:
                bnew = [b[2], b[1], b[4], b[3]] if flip_face else [b[2], b[3], b[4], b[1]]
            elif ids[1] == 3:
                bnew = [b[3], b[2], b[1], b[4]] if flip_face else [b[3], b[4], b[1], b[2]]
            elif ids[1] == 4:
                bnew = [b[4], b[3], b[2], b[1]] if flip_face else [b[4], b[1], b[2], b[3]]
            
        if a[0] == 3:

            if ids[0] == 1:
                anew = [a[1], a[3], a[2]] if flip_face else [a[1], a[2], a[3]]
            elif ids[0] == 2:
                anew = [a[2], a[1], a[3]] if flip_face else [a[2], a[3], a[1]]
            elif ids[0] == 3:
                anew = [a[3], a[2], a[1]] if flip_face else [a[3], a[1], a[2]]
            elif ids[0] == 4:
                anew = [a[3], a[2], a[1]] if flip_face else [a[1], a[2], a[3]]
                
            flip_face = not flip_face

            if ids[1] == 1:
                bnew = [b[1], b[3], b[2]] if flip_face else [b[1], b[2], b[3]]
            elif ids[1] == 2:
                bnew = [b[2], b[1], b[3]] if flip_face else [b[2], b[3], b[1]]
            elif ids[1] == 3:
                bnew = [b[3], b[2], b[1]] if flip_face else [b[3], b[1], b[2]]
            elif ids[1] == 4:
                bnew = [b[3], b[2], b[1]] if flip_face else [b[1], b[2], b[3]]
            
        return [anew, bnew]

    inner_faces = np.zeros((out['nte'], 5), dtype=np.int32)
    inner_vertices = []

    ids = []

    for i in range(out['nte']):
            
        # Upper and lower faces
        face1 = out['fc'][out['te_fc'][i, 0], :]
        face2 = out['fc'][out['te_fc'][i, 1], :]

        # Correct order
        ids1, ids2 = __find_internal_ids(face1, face2, i)

        # Add ids and vertices
        for j in range(face1[0]):
            if ids1[j] not in ids:
                ids.append(ids1[j])
                inner_vertices.append(0.5 * (out['vt'][ids1[j], :] + out['vt'][ids2[j], :]))
            
        # Add face
        inner_faces[i, 0] = face1[0]
        inner_faces[i, 1] = ids.index(ids1[0])
        inner_faces[i, 2] = ids.index(ids1[1])
        inner_faces[i, 3] = ids.index(ids1[2])
        if face1[0] == 4: inner_faces[i, 4] = ids.index(ids1[3])

    inner_vertices = np.asarray(inner_vertices)

    out['vt_d'] = inner_vertices
    out['fc_d'] = inner_faces

    return out
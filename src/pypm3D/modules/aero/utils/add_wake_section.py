import typing as tp
import numpy as np


from pypm3D.modules.aero.utils import calculate_wake_section_velocity, calculate_surface_section_velocity


def main(nte: int,
         nv_te: int,
         nf: int,
         section: int,
         vertices_w: np.ndarray,
         faces_w: np.ndarray,
         trailing_edge_ids: np.ndarray,
         areas: np.ndarray,
         circulations: np.ndarray,
         vel_surface: np.ndarray,
         n_sides_s: np.ndarray,
         p_avg_s: np.ndarray,
         e1_s: np.ndarray,
         e2_s: np.ndarray,
         e3_s: np.ndarray,
         p1_s: np.ndarray,
         p2_s: np.ndarray,
         p3_s: np.ndarray,
         p4_s: np.ndarray,
         source_s: np.ndarray,
         doublet_s: np.ndarray,
         time_step: float,
         u_func: tp.Callable[[np.ndarray], np.ndarray]):

    # Parameters
    vel = np.empty((nv_te, section, 3), dtype=np.double)
    freestream = np.empty((nv_te, 3), dtype=np.double)
    ind_vel_wake = np.empty((nv_te, 3), dtype=np.double)
    ind_vel_surf = np.empty((nv_te, 3), dtype=np.double)

    # Calculate velocity field
    for i in range(section):
        if i == 0:
            vel[:, i, :] = vel_surface[trailing_edge_ids, :]
        else:
            freestream[:] = np.asarray([u_func(x) for x in vertices_w[nv_te * i:nv_te * (i + 1)]])
            ind_vel_wake[:] = calculate_wake_section_velocity.main(nv_te, nte, section, vertices_w, faces_w, areas, circulations, vertices_w[nv_te * i:nv_te * (i + 1), :])
            ind_vel_surf[:] = calculate_surface_section_velocity.main(nv_te, nf, n_sides_s, p_avg_s, e1_s, e2_s, e3_s, p1_s, p2_s, p3_s, p4_s, source_s, doublet_s, vertices_w[nv_te * i:nv_te * (i + 1), :])
            vel[:, i] = freestream[:] + ind_vel_wake[:] + ind_vel_surf[:]

    # Update position, area and circulation
    for i in range(section):
        vertices_w[nv_te * (section - i):nv_te * ((section + 1) - i), :] = vertices_w[nv_te * ((section - 1) - i):nv_te * (section - i), :] + time_step * vel[:, section - 1 - i, :]
    
    # Update area and circulation
    if section > 1:
        for i in range(section - 1, 0, -1):
            areas[i * nte:(i + 1) * nte] = areas[(i - 1) * nte:i * nte]
            circulations[i * nte:(i + 1) * nte] = circulations[(i - 1) * nte:i * nte]

    # Calculate area
    for i in range(nte):
        v1 = vertices_w[faces_w[i, 2], :] - vertices_w[faces_w[i, 1], :]
        v2 = vertices_w[faces_w[i, 3], :] - vertices_w[faces_w[i, 1], :]
        v3 = vertices_w[faces_w[i, 4], :] - vertices_w[faces_w[i, 1], :]
        areas[i] = 0.5 * (np.linalg.norm(np.cross(v1, v2)) + np.linalg.norm(np.cross(v2, v3)))
    
    return

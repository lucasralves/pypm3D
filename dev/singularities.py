import numpy as np
import plotly.graph_objects as go

from potential import potential_wrapper
from velocity import velocity_wrapper


def _test_quad_panel(singularity: str = 'doublet'):

    # Panel vertices
    p1 = np.array([ 1.0, 1.0, 0.0])
    p2 = np.array([-1.0, 1.0, 0.0])
    p3 = np.array([-1.0,-1.0, 0.0])
    p4 = np.array([ 1.0,-1.0, 0.0])

    p_avg = 0.25 * (p1 + p2 + p3 + p4)

    p1 = p1 - p_avg
    p2 = p2 - p_avg
    p3 = p3 - p_avg
    p4 = p4 - p_avg

    p_avg = 0.25 * (p1 + p2 + p3 + p4)

    e1 = np.array([ 1.0, 0.0, 0.0])
    e2 = np.array([ 0.0, 1.0, 0.0])
    e3 = np.array([ 0.0, 0.0, 1.0])

    # Grid parameters
    side = 2.0
    n = 30
    z0 = -0.5
    x0 = 0.0

    # Grid points
    xs = np.linspace(-side, side, num=n)
    ys = np.linspace(-side, side, num=n)
    zs = np.linspace(-side, side, num=n)

    xy, yx = np.meshgrid(xs, ys)
    yz, zy = np.meshgrid(ys, zs)

    # Control points
    p_ctrl_xy = np.empty((n * n, 3), dtype=np.double)

    p_ctrl_xy[:, 0] = xy.flatten()
    p_ctrl_xy[:, 1] = yx.flatten()
    p_ctrl_xy[:, 2] = z0

    p_ctrl_yz = np.empty((n * n, 3), dtype=np.double)

    p_ctrl_yz[:, 0] = x0
    p_ctrl_yz[:, 1] = yz.flatten()
    p_ctrl_yz[:, 2] = zy.flatten()

    # Potentials
    source_pot_xy = np.empty((n * n), dtype=np.double)
    source_pot_yz = np.empty((n * n), dtype=np.double)

    doublet_pot_xy = np.empty((n * n), dtype=np.double)
    doublet_pot_yz = np.empty((n * n), dtype=np.double)

    # Velocities
    source_vel_xy = np.empty((n * n, 3), dtype=np.double)
    source_vel_yz = np.empty((n * n, 3), dtype=np.double)

    doublet_vel_xy = np.empty((n * n, 3), dtype=np.double)
    doublet_vel_yz = np.empty((n * n, 3), dtype=np.double)

    # Call wrapper
    source_pot_xy, doublet_pot_xy = potential_wrapper.main(4, p1, p2, p3, p4, e1, e2, e3, p_avg, p_ctrl_xy)
    source_pot_yz, doublet_pot_yz = potential_wrapper.main(4, p1, p2, p3, p4, e1, e2, e3, p_avg, p_ctrl_yz)

    source_vel_xy, doublet_vel_xy = velocity_wrapper.main(4, p1, p2, p3, p4, e1, e2, e3, p_avg, p_ctrl_xy)
    source_vel_yz, doublet_vel_yz = velocity_wrapper.main(4, p1, p2, p3, p4, e1, e2, e3, p_avg, p_ctrl_yz)
    
    # Reshape
    source_pot_xy = source_pot_xy.reshape((n, n))
    source_pot_yz = source_pot_yz.reshape((n, n))
    doublet_pot_xy = doublet_pot_xy.reshape((n, n))
    doublet_pot_yz = doublet_pot_yz.reshape((n, n))

    if singularity == 'doublet':
        
        vmax1 = max([np.max(doublet_pot_xy), np.max(doublet_pot_yz)])
        vmin1 = min([np.min(doublet_pot_xy), np.min(doublet_pot_yz)])

        vmax2 = max([np.max(np.sqrt(doublet_vel_xy[:, 0] * doublet_vel_xy[:, 0] + doublet_vel_xy[:, 1] * doublet_vel_xy[:, 1] + doublet_vel_xy[:, 2] * doublet_vel_xy[:, 2])), np.max(np.sqrt(doublet_vel_yz[:, 0] * doublet_vel_yz[:, 0] + doublet_vel_yz[:, 1] * doublet_vel_yz[:, 1] + doublet_vel_yz[:, 2] * doublet_vel_yz[:, 2]))])
        vmin2 = min([np.min(np.sqrt(doublet_vel_xy[:, 0] * doublet_vel_xy[:, 0] + doublet_vel_xy[:, 1] * doublet_vel_xy[:, 1] + doublet_vel_xy[:, 2] * doublet_vel_xy[:, 2])), np.min(np.sqrt(doublet_vel_yz[:, 0] * doublet_vel_yz[:, 0] + doublet_vel_yz[:, 1] * doublet_vel_yz[:, 1] + doublet_vel_yz[:, 2] * doublet_vel_yz[:, 2]))])

        fig = go.Figure(
            data=[
                go.Surface(x=xy, y=yx, z=z0 * np.ones_like(doublet_pot_xy), surfacecolor=doublet_pot_xy, cmin=vmin1, cmax=vmax1, showscale=False),
                go.Surface(x=x0 * np.ones_like(doublet_pot_yz), y=yz, z=zy, surfacecolor=doublet_pot_yz, cmin=vmin1, cmax=vmax1, showscale=False),
                
                go.Scatter3d(x=[p1[0], p2[0], p3[0], p4[0], p1[0]], y=[p1[1], p2[1], p3[1], p4[1], p1[1]], z=[p1[2], p2[2], p3[2], p4[2], p1[2]], mode='lines', line=dict(color='black', width=2)),

                go.Cone(x=p_ctrl_xy[:, 0], y=p_ctrl_xy[:, 1], z=p_ctrl_xy[:, 2], u=doublet_vel_xy[:, 0], v=doublet_vel_xy[:, 1], w=doublet_vel_xy[:, 2], cmin=vmin2, cmax=vmax2, showscale=False),
                go.Cone(x=p_ctrl_yz[:, 0], y=p_ctrl_yz[:, 1], z=p_ctrl_yz[:, 2], u=doublet_vel_yz[:, 0], v=doublet_vel_yz[:, 1], w=doublet_vel_yz[:, 2], cmin=vmin2, cmax=vmax2, showscale=True),
            ]
        )
        fig.show()
    
    if singularity == 'source':

        vmax1 = max([np.max(source_pot_xy), np.max(source_pot_yz)])
        vmin1 = min([np.min(source_pot_xy), np.min(source_pot_yz)])

        vmax2 = max([np.max(np.sqrt(source_vel_xy[:, 0] * source_vel_xy[:, 0] + source_vel_xy[:, 1] * source_vel_xy[:, 1] + source_vel_xy[:, 2] * source_vel_xy[:, 2])), np.max(np.sqrt(source_vel_yz[:, 0] * source_vel_yz[:, 0] + source_vel_yz[:, 1] * source_vel_yz[:, 1] + source_vel_yz[:, 2] * source_vel_yz[:, 2]))])
        vmin2 = min([np.min(np.sqrt(source_vel_xy[:, 0] * source_vel_xy[:, 0] + source_vel_xy[:, 1] * source_vel_xy[:, 1] + source_vel_xy[:, 2] * source_vel_xy[:, 2])), np.min(np.sqrt(source_vel_yz[:, 0] * source_vel_yz[:, 0] + source_vel_yz[:, 1] * source_vel_yz[:, 1] + source_vel_yz[:, 2] * source_vel_yz[:, 2]))])

        fig = go.Figure(
            data=[
                go.Surface(x=xy, y=yx, z=z0 * np.ones_like(source_pot_xy), surfacecolor=source_pot_xy, cmin=vmin1, cmax=vmax1, showscale=False),
                go.Surface(x=x0 * np.ones_like(source_pot_yz), y=yz, z=zy, surfacecolor=source_pot_yz, cmin=vmin1, cmax=vmax1, showscale=False),
                
                go.Scatter3d(x=[p1[0], p2[0], p3[0], p4[0], p1[0]], y=[p1[1], p2[1], p3[1], p4[1], p1[1]], z=[p1[2], p2[2], p3[2], p4[2], p1[2]], mode='lines', line=dict(color='black', width=2)),

                go.Cone(x=p_ctrl_xy[:, 0], y=p_ctrl_xy[:, 1], z=p_ctrl_xy[:, 2], u=source_vel_xy[:, 0], v=source_vel_xy[:, 1], w=source_vel_xy[:, 2], cmin=vmin2, cmax=vmax2, showscale=False),
                go.Cone(x=p_ctrl_yz[:, 0], y=p_ctrl_yz[:, 1], z=p_ctrl_yz[:, 2], u=source_vel_yz[:, 0], v=source_vel_yz[:, 1], w=source_vel_yz[:, 2], cmin=vmin2, cmax=vmax2, showscale=True),
            ]
        )
        fig.show()

    return

def _test_tri_panel(singularity: str = 'doublet'):

    # Panel vertices
    p1 = np.array([ 1.0, 1.0, 0.0])
    p2 = np.array([-1.0, 1.0, 0.0])
    p3 = np.array([-1.0,-1.0, 0.0])
    p4 = np.empty(3)

    p_avg = (1.0 / 3.0) * (p1 + p2 + p3)

    p1 = p1 - p_avg
    p2 = p2 - p_avg
    p3 = p3 - p_avg

    p_avg = (1.0 / 3.0) * (p1 + p2 + p3)

    e1 = np.array([ 1.0, 0.0, 0.0])
    e2 = np.array([ 0.0, 1.0, 0.0])
    e3 = np.array([ 0.0, 0.0, 1.0])

    # Grid parameters
    side = 2.0
    n = 30
    z0 = -0.5
    x0 = 0.0

    # Grid points
    xs = np.linspace(-side, side, num=n)
    ys = np.linspace(-side, side, num=n)
    zs = np.linspace(-side, side, num=n)

    xy, yx = np.meshgrid(xs, ys)
    yz, zy = np.meshgrid(ys, zs)

    # Control points
    p_ctrl_xy = np.empty((n * n, 3), dtype=np.double)

    p_ctrl_xy[:, 0] = xy.flatten()
    p_ctrl_xy[:, 1] = yx.flatten()
    p_ctrl_xy[:, 2] = z0

    p_ctrl_yz = np.empty((n * n, 3), dtype=np.double)

    p_ctrl_yz[:, 0] = x0
    p_ctrl_yz[:, 1] = yz.flatten()
    p_ctrl_yz[:, 2] = zy.flatten()

    # Potentials
    source_pot_xy = np.empty((n * n), dtype=np.double)
    source_pot_yz = np.empty((n * n), dtype=np.double)

    doublet_pot_xy = np.empty((n * n), dtype=np.double)
    doublet_pot_yz = np.empty((n * n), dtype=np.double)

    # Velocities
    source_vel_xy = np.empty((n * n, 3), dtype=np.double)
    source_vel_yz = np.empty((n * n, 3), dtype=np.double)

    doublet_vel_xy = np.empty((n * n, 3), dtype=np.double)
    doublet_vel_yz = np.empty((n * n, 3), dtype=np.double)

    # Call wrapper
    source_pot_xy, doublet_pot_xy = potential_wrapper.main(3, p1, p2, p3, p4, e1, e2, e3, p_avg, p_ctrl_xy)
    source_pot_yz, doublet_pot_yz = potential_wrapper.main(3, p1, p2, p3, p4, e1, e2, e3, p_avg, p_ctrl_yz)

    source_vel_xy, doublet_vel_xy = velocity_wrapper.main(3, p1, p2, p3, p4, e1, e2, e3, p_avg, p_ctrl_xy)
    source_vel_yz, doublet_vel_yz = velocity_wrapper.main(3, p1, p2, p3, p4, e1, e2, e3, p_avg, p_ctrl_yz)
    
    # Reshape
    source_pot_xy = source_pot_xy.reshape((n, n))
    source_pot_yz = source_pot_yz.reshape((n, n))
    doublet_pot_xy = doublet_pot_xy.reshape((n, n))
    doublet_pot_yz = doublet_pot_yz.reshape((n, n))

    if singularity == 'doublet':
        
        vmax1 = max([np.max(doublet_pot_xy), np.max(doublet_pot_yz)])
        vmin1 = min([np.min(doublet_pot_xy), np.min(doublet_pot_yz)])

        vmax2 = max([np.max(np.sqrt(doublet_vel_xy[:, 0] * doublet_vel_xy[:, 0] + doublet_vel_xy[:, 1] * doublet_vel_xy[:, 1] + doublet_vel_xy[:, 2] * doublet_vel_xy[:, 2])), np.max(np.sqrt(doublet_vel_yz[:, 0] * doublet_vel_yz[:, 0] + doublet_vel_yz[:, 1] * doublet_vel_yz[:, 1] + doublet_vel_yz[:, 2] * doublet_vel_yz[:, 2]))])
        vmin2 = min([np.min(np.sqrt(doublet_vel_xy[:, 0] * doublet_vel_xy[:, 0] + doublet_vel_xy[:, 1] * doublet_vel_xy[:, 1] + doublet_vel_xy[:, 2] * doublet_vel_xy[:, 2])), np.min(np.sqrt(doublet_vel_yz[:, 0] * doublet_vel_yz[:, 0] + doublet_vel_yz[:, 1] * doublet_vel_yz[:, 1] + doublet_vel_yz[:, 2] * doublet_vel_yz[:, 2]))])

        fig = go.Figure(
            data=[
                go.Surface(x=xy, y=yx, z=z0 * np.ones_like(doublet_pot_xy), surfacecolor=doublet_pot_xy, cmin=vmin1, cmax=vmax1, showscale=False),
                go.Surface(x=x0 * np.ones_like(doublet_pot_yz), y=yz, z=zy, surfacecolor=doublet_pot_yz, cmin=vmin1, cmax=vmax1, showscale=False),
                
                go.Scatter3d(x=[p1[0], p2[0], p3[0], p1[0]], y=[p1[1], p2[1], p3[1], p1[1]], z=[p1[2], p2[2], p3[2], p1[2]], mode='lines', line=dict(color='black', width=2)),

                go.Cone(x=p_ctrl_xy[:, 0], y=p_ctrl_xy[:, 1], z=p_ctrl_xy[:, 2], u=doublet_vel_xy[:, 0], v=doublet_vel_xy[:, 1], w=doublet_vel_xy[:, 2], cmin=vmin2, cmax=vmax2, showscale=False),
                go.Cone(x=p_ctrl_yz[:, 0], y=p_ctrl_yz[:, 1], z=p_ctrl_yz[:, 2], u=doublet_vel_yz[:, 0], v=doublet_vel_yz[:, 1], w=doublet_vel_yz[:, 2], cmin=vmin2, cmax=vmax2, showscale=True),
            ]
        )
        fig.show()
    
    if singularity == 'source':

        vmax1 = max([np.max(source_pot_xy), np.max(source_pot_yz)])
        vmin1 = min([np.min(source_pot_xy), np.min(source_pot_yz)])

        vmax2 = max([np.max(np.sqrt(source_vel_xy[:, 0] * source_vel_xy[:, 0] + source_vel_xy[:, 1] * source_vel_xy[:, 1] + source_vel_xy[:, 2] * source_vel_xy[:, 2])), np.max(np.sqrt(source_vel_yz[:, 0] * source_vel_yz[:, 0] + source_vel_yz[:, 1] * source_vel_yz[:, 1] + source_vel_yz[:, 2] * source_vel_yz[:, 2]))])
        vmin2 = min([np.min(np.sqrt(source_vel_xy[:, 0] * source_vel_xy[:, 0] + source_vel_xy[:, 1] * source_vel_xy[:, 1] + source_vel_xy[:, 2] * source_vel_xy[:, 2])), np.min(np.sqrt(source_vel_yz[:, 0] * source_vel_yz[:, 0] + source_vel_yz[:, 1] * source_vel_yz[:, 1] + source_vel_yz[:, 2] * source_vel_yz[:, 2]))])

        fig = go.Figure(
            data=[
                go.Surface(x=xy, y=yx, z=z0 * np.ones_like(source_pot_xy), surfacecolor=source_pot_xy, cmin=vmin1, cmax=vmax1, showscale=False),
                go.Surface(x=x0 * np.ones_like(source_pot_yz), y=yz, z=zy, surfacecolor=source_pot_yz, cmin=vmin1, cmax=vmax1, showscale=False),
                
                go.Scatter3d(x=[p1[0], p2[0], p3[0], p1[0]], y=[p1[1], p2[1], p3[1], p1[1]], z=[p1[2], p2[2], p3[2], p1[2]], mode='lines', line=dict(color='black', width=2)),

                go.Cone(x=p_ctrl_xy[:, 0], y=p_ctrl_xy[:, 1], z=p_ctrl_xy[:, 2], u=source_vel_xy[:, 0], v=source_vel_xy[:, 1], w=source_vel_xy[:, 2], cmin=vmin2, cmax=vmax2, showscale=False),
                go.Cone(x=p_ctrl_yz[:, 0], y=p_ctrl_yz[:, 1], z=p_ctrl_yz[:, 2], u=source_vel_yz[:, 0], v=source_vel_yz[:, 1], w=source_vel_yz[:, 2], cmin=vmin2, cmax=vmax2, showscale=True),
            ]
        )
        fig.show()

    return


if __name__ == '__main__':
    # _test_tri_panel('doublet') # 'source'
    _test_quad_panel('source') # 'doublet'
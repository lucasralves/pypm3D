import typing as tp
import numpy as np

from pypm3D.modules.preproc import _Mesh
from pypm3D.modules.proc import _Solver
from pypm3D.modules.posproc import _VtpFile


def init() -> None:
    """
    pypm3D.init()

    Initialize _Case class. This must be called before any call to the other
    functions.
    """
    
    global _mesh
    global _solver

    _mesh = _Mesh()
    _solver = _Solver(_mesh)

    return

def proc_mesh(vt: np.ndarray, fc: np.ndarray, te: np.ndarray,
              u_call: tp.Callable[[np.ndarray], np.ndarray], l: float, ts: float) -> None:
    """
    pypm3D.proc_mesh(vt, fc, te, u_call, l, ts)

    The surface's parameters are created.
        
    The geometry must be oriented as:
     - The span must be in the y direction
     - The longitidunal axis must be pointing to the negative direction of x.
     - The upper side must be in the positive z direction

    Parameters:
    -----------
    - vertices: points in space [x, y, z]
    - faces: number of sides and vertices ids positive oriented [n, id1, id2, id3, id4]
    - trailing_edge: vertices that makes an edge [id1, id2]
    - u_call([x, y, z]): calculate the freestream velocity [u, v, w] based at a point [x, y, z]
    - l: wake's length
    - ts: pseudo time step
    """

    _mesh.sf.init(vt, fc, te)
    _mesh.wk.init(u_call, l, ts)

    return

def solve(u_call: tp.Callable[[np.ndarray], np.ndarray]) -> None:
    """
    pypm3D.solve(u_call)

    Calculate a potential flow using dirichlet boundarycondition over
    an arbitrary lifting object.

    Parameters:
    -----------
    - u_call: callable that receives a position (x, y, z) and returns
      a velocity (ux, uy, uz).
    """

    _solver.run(u_call)

    return

def gen_vtk(fn: str, wake: bool = False, u_call: tp.Callable[[np.ndarray], np.ndarray] = None, type: str = 'grid') -> None:
    """
    Create a vtp file.

    Parameters:
    -----------
    - fn: file name
    - wake: show the wake if True.
    - u_call: callable that receives a position (x, y, z) and returns
      a velocity (ux, uy, uz). Cannot be None if wake is True.
    - type: wake type. Must be 'grid' or 'surface'
    """

    assert type in ['grid', 'surface'], 'Incorrect wake_type'
    if wake:
        assert u_call is not None, 'u_call cannot be None if wake is True'
    
    vtp = _VtpFile()

    vtp.add_surface(_mesh.sf.vt, _mesh.sf.fc)

    if wake:
        if type == 'grid':
            vtp.add_grid(_mesh.wk.vt, _mesh.wk.fc)
        else:
            vtp.add_surface(_mesh.wk.vt, _mesh.wk.fc)
    
    if _solver.done:
        if wake:
            vtp.add_param('source', 1, [_solver.source_sf, np.zeros(_mesh.wk.nf)])
            vtp.add_param('doublet', 1, [_solver.doublet_sf, _solver.doublet_wk])
            vtp.add_param('cp', 1, [_solver.cp_sf, np.zeros(_mesh.wk.nf)])
            vtp.add_param('vel', 3, [_solver.vel_sf, np.asarray([u_call(x) for x in _mesh.wk.p_avg])])
        else:
            vtp.add_param('source', 1, [_solver.source_sf])
            vtp.add_param('doublet', 1, [_solver.source_sf])
            vtp.add_param('cp', 1, [_solver.cp_sf])
            vtp.add_param('vel', 3, [_solver.vel_sf])

    vtp.write(fn)

    return

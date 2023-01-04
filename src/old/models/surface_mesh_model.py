from __future__ import annotations
from dataclasses import dataclass
from numpy import ndarray, empty, double, int32

@dataclass
class SurfaceMeshModel:

    # Size
    nv: int                     # number of vertices
    nf: int                     # number of faces
    nte: int                    # number of edges at the trailing edge

    # Surface panels
    vertices: ndarray           # (nv, 3) [double - x, y and z positions]
    faces: ndarray              # (nf, 5) [int - vertices ids]

    # Trailing edge
    trailing_edge: ndarray      # (nte, 2) - [int - vertices ids]

    # Surface parameters
    p_avg: ndarray              # (nf, 3) [double - x, y and z positions]
    p_ctrl: ndarray             # (nf, 3) [double - x, y and z positions]
    e1: ndarray                 # (nf, 3) [double - x, y and z positions]
    e2: ndarray                 # (nf, 3) [double - x, y and z positions]
    e3: ndarray                 # (nf, 3) [double - x, y and z positions]
    p1: ndarray                 # (nf, 2) [double - x and y positions]
    p2: ndarray                 # (nf, 2) [double - x and y positions]
    p3: ndarray                 # (nf, 2) [double - x and y positions]
    p4: ndarray                 # (nf, 2) [double - x and y positions]

    # Trailing edge parameters
    te_faces: ndarray           # (nte, 2) [int - faces ids]
    inner_faces: ndarray        # (nte, 5) [int - inner vertices ids]
    inner_vertices: ndarray     # (?, 3)   [double - x, y and z positions]
    inner_scale_factor: ndarray # (nte,)   [double - scale factor]

    def surface_mesh_from_sizes(vertices: ndarray, faces: ndarray, trailing_edge: ndarray) -> SurfaceMeshModel:

        nv = vertices.shape[0]
        nf = faces.shape[0]
        nte = trailing_edge.shape[0]

        mesh = SurfaceMeshModel(
            nv=nv,
            nf=nf,
            nte=nte,
            vertices=vertices,
            faces=faces,
            trailing_edge=trailing_edge,
            p_avg=empty((nf, 3), dtype=double),
            p_ctrl=empty((nf, 3), dtype=double),
            e1=empty((nf, 3), dtype=double),
            e2=empty((nf, 3), dtype=double),
            e3=empty((nf, 3), dtype=double),
            p1=empty((nf, 2), dtype=double),
            p2=empty((nf, 2), dtype=double),
            p3=empty((nf, 2), dtype=double),
            p4=empty((nf, 2), dtype=double),
            te_faces=empty((nte, 2), dtype=int32),
            inner_faces=empty((nte, 2), dtype=int32),
            inner_vertices=empty((nte, 2), dtype=int32),
            inner_scale_factor=empty((nte, 2), dtype=int32),
        )

        return mesh
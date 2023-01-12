from dataclasses import dataclass
import math
import typing as tp
import numpy as np


@dataclass
class SurfaceMeshModel:
    nv: int = None                          # number of vertices
    nf: int = None                          # number of faces
    nte: int = None                         # number of edges at the trailing edge
    vertices: np.ndarray = None             # points in space [x, y, z]
    faces: np.ndarray = None                # number of sides and vertices ids [n, id1, id2, id3, id4]
    trailing_edge: np.ndarray = None        # vertices that makes an edge [id1, id2]
    p_avg: np.ndarray = None                # center of panel
    p_ctrl: np.ndarray = None               # control point
    e1: np.ndarray = None                   # orthogonal base
    e2: np.ndarray = None                   # orthogonal base
    e3: np.ndarray = None                   # orthogonal base
    p1: np.ndarray = None                   # panel vertice at local system
    p2: np.ndarray = None                   # panel vertice at local system
    p3: np.ndarray = None                   # panel vertice at local system
    p4: np.ndarray = None                   # panel vertice at local system
    trailing_edge_faces: np.ndarray = None  # faces ids that contain the same edge

    def init(self, vertices: np.ndarray,
                   faces: np.ndarray,
                   trailing_edge: np.ndarray) -> None:
        
        self.nv: int = vertices.shape[0]
        self.nf: int = faces.shape[0]
        self.nte: int = trailing_edge.shape[0]

        self.vertices: np.ndarray = np.copy(vertices.astype(np.double))
        self.faces: np.ndarray = np.copy(faces.astype(np.int32))
        self.trailing_edge: np.ndarray = np.copy(trailing_edge.astype(np.int32))

        self.p_avg = np.empty((self.nf, 3), dtype=np.double)
        self.p_ctrl = np.empty((self.nf, 3), dtype=np.double)
        self.e1 = np.empty((self.nf, 3), dtype=np.double)
        self.e2 = np.empty((self.nf, 3), dtype=np.double)
        self.e3 = np.empty((self.nf, 3), dtype=np.double)
        self.p1 = np.empty((self.nf, 2), dtype=np.double)
        self.p2 = np.empty((self.nf, 2), dtype=np.double)
        self.p3 = np.empty((self.nf, 2), dtype=np.double)
        self.p4 = np.empty((self.nf, 2), dtype=np.double)

        self.trailing_edge_faces = np.empty((self.nte, 2), dtype=np.int32)

        return

@dataclass
class WakeMeshModel:
    nv_te: int = None                       # number of vertices at the trailing edge
    nte: int = None                         # number of edges at the trailing edge
    nw: int = None                          # number of sections on the wake
    vertices: np.ndarray = None             # points in space [x, y, z]
    faces: np.ndarray = None                # number of sides and vertices ids [n, id1, id2, id3, id4]
    trailing_edge_ids: np.ndarray = None    # ids of the surface vertices at the trailing edge
    areas: np.ndarray = None                # areas of each panel when they are shedded into the wake
    time_step: float = None                 # time step used to create the wake

    def init(self, trailing_edge: np.ndarray,
                   surface_vertices: np.ndarray,
                   u_ref: np.ndarray,
                   length: float,
                   time_step: float) -> None:
        
        # Find the unique trailing edge ids
        trailing_edge_ids = []

        for edge in trailing_edge:
            if edge[0] not in trailing_edge_ids: trailing_edge_ids.append(edge[0])
            if edge[1] not in trailing_edge_ids: trailing_edge_ids.append(edge[1])

        # Save and create parameters
        self.nv_te = len(trailing_edge_ids)
        self.nte = trailing_edge.shape[0]
        self.nw = math.ceil(length / (time_step * np.linalg.norm(u_ref)))

        self.vertices = np.empty((self.nv_te * (self.nw + 1), 3), dtype=np.double)
        self.faces = np.empty((self.nte * self.nw, 5), dtype=np.int32)
        self.trailing_edge_ids = np.asarray(trailing_edge_ids, dtype=np.int32)

        # Save trailing edge vertices
        for i in range(self.nv_te):
            self.vertices[i, :] = surface_vertices[self.trailing_edge_ids[i], :]
        
        # Erase here
        for i in range(self.nw):
            self.vertices[(i + 1) * self.nv_te:(i + 2) * self.nv_te, :] = self.vertices[i * self.nv_te:(i + 1) * self.nv_te, :] + u_ref * time_step

        # Create faces ids
        for i in range(self.nte):

            id1 = trailing_edge_ids.index(trailing_edge[i, 0])
            id2 = trailing_edge_ids.index(trailing_edge[i, 1])

            for j in range(self.nw):
                self.faces[i + j * self.nte, 0] = 4
                self.faces[i + j * self.nte, 1] = id1 + j * self.nv_te
                self.faces[i + j * self.nte, 2] = id2 + j * self.nv_te
                self.faces[i + j * self.nte, 3] = id2 + (j + 1) * self.nv_te
                self.faces[i + j * self.nte, 4] = id1 + (j + 1) * self.nv_te
        
        # Create wake parameters
        self.areas = np.empty((self.nte, self.nw), dtype=np.double)

        # Time step
        self.time_step = time_step

        return

@dataclass
class MeshModel:
    surface: SurfaceMeshModel = None
    wake: WakeMeshModel = None

    def init_surface(self, vertices: np.ndarray,
                           faces: np.ndarray,
                           trailing_edge: np.ndarray) -> None:

        self.surface = SurfaceMeshModel()
        self.surface.init(vertices, faces, trailing_edge)

        return
    
    def init_wake(self, u_func: tp.Callable[[np.ndarray], np.ndarray],
                        length: float,
                        time_step: float) -> None:

        self.wake = WakeMeshModel()
        self.wake.init(self.surface.trailing_edge, self.surface.vertices, u_func(np.mean(self.surface.vertices, axis=0)), length, time_step)

        return
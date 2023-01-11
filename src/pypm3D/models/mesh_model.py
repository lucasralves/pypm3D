import numpy as np

class SurfaceMeshModel:

    def __init__(self, vertices: np.ndarray,
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

class InternalMeshModel:
    
    def __init__(self) -> None:

        self.nv: int = None
        self.nf: int = None

        self.vertices: np.ndarray = None
        self.faces: np.ndarray = None

        self.scale: np.ndarray = None

        return


class MeshModel:
    
    def __init__(self, vertices: np.ndarray,
                       faces: np.ndarray,
                       trailing_edge: np.ndarray) -> None:
        
        self.surface: SurfaceMeshModel = SurfaceMeshModel(vertices, faces, trailing_edge)
        self.internal: InternalMeshModel = InternalMeshModel()
        
        return
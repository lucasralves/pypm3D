import typing as tp

class VerticeConnectionModel:

    def __init__(self, n: int, faces: tp.List[int], coefs: tp.List[float]) -> None:

        self.n: int = n
        self.faces: tp.List[int] = faces
        self.coefs: tp.List[float] = coefs

        return
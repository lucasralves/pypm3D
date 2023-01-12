from dataclasses import dataclass
import typing as tp


@dataclass
class VerticeModel:
    n: int
    faces: tp.List[int]
    coefs: tp.List[float]
from dataclasses import dataclass
import numpy as np

@dataclass
class SystemState:
    N: int                     # number of particles
    xs: np.ndarray             # positions (N, 3)
    ps: np.ndarray             # momenta (N, 3)
    ms: np.ndarray             # masses (N, 1)
    labels: list[str]          
    box: np.ndarray | None     # (3,), None = no box (infinite space)
    
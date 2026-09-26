from abc import ABC, abstractmethod

import numpy as np
import numba as nb

from .state import SystemState

class Interaction(ABC):
    """Base class for interactions returning (forces, potential_energy)."""

    @abstractmethod
    def compute(self, state: SystemState) -> tuple[np.ndarray, float]:
        """Compute forces and potential energy for the given system state."""
        pass


class HarmonicPotential(Interaction):
    """Independent harmonic coordinates about the origin: U = k/2 * sum(x**2)."""

    def __init__(self, k: float = 1.0):
        if not np.isfinite(k) or k < 0:
            raise ValueError("k must be finite and non-negative.")
        self.k = k

    def compute(self, state: SystemState) -> tuple[np.ndarray, float]:
        
        if state.box is not None:
            raise ValueError("HarmonicPotential requires box=None.")
        
        U = self.k * 0.5 * (state.xs**2).sum()
        dU_dx = self.k * state.xs
        
        return -dU_dx, float(U)


@nb.jit(nopython=True)
def interactions_LJ_numba(x: np.ndarray, box: np.ndarray | None, sigma: float = 1.0, eps: float = 1.0, cutoff: float | None = None, shift_energy: bool = True):
    """Lennard-Jones energy and derivatives, from ../simplemd/interactions.py.

    The pair calculation follows the original interactions_LJ_numba.
    Added support for box=None, an optional cutoff and an energy shift.
    """
    N, d = x.shape

    if d != 3:
        raise ValueError("Positions must have 3 spatial dimensions.")

    sigma2 = sigma**2
    eps4 = 4.0 * eps
    eps24 = 24.0 * eps

    boxinv = None
    if box is not None:
        boxinv = 1.0 / box

    xij = np.zeros(3)
    
    U = 0.0
    dU_dx = np.zeros_like(x)

    r_c2 = np.inf if cutoff is None else cutoff**2
    u_shift = 0.0
    if cutoff is not None and shift_energy:
        sdinv6_c = (sigma / cutoff)**6
        u_shift = eps4 * (sdinv6_c**2 - sdinv6_c)

    for i in range(N-1):
        for j in range(i+1, N):

            d2 = 0.0
            for k in range(3):
                dxk = x[j, k] - x[i, k]
                if box is not None:
                    dxk -= box[k] * round(dxk * boxinv[k])
                xij[k] = dxk
                d2 += dxk * dxk

            if d2 == 0.0:
                raise ValueError("Distinct particles overlap; the LJ potential is singular.")
            if d2 >= r_c2:
                continue

            dinv2 = 1.0 / d2
            sdinv6 = (sigma2 * dinv2)**3
            sdinv12 = sdinv6**2
            u = eps4 * (sdinv12 - sdinv6)
            du_ddd = -eps24 * (2*sdinv12 - sdinv6) * dinv2
            
            U += u - u_shift

            for k in range(3):
                dU_dxk = du_ddd * xij[k]
                dU_dx[i, k] -= dU_dxk
                dU_dx[j, k] += dU_dxk

    return U, dU_dx


class LennardJonesPotential(Interaction):
    """Lennard-Jones pair potential, optionally with periodic boundaries.

    Adapted from ../simplemd/interactions.py: interactions_LJ_numba.
    Uses a Numba-compiled pair loop, O(N**2) work and O(N) memory.
    Without cutoff, reproduces the original all-pairs calculation.

    cutoff is an absolute distance in the same units as positions and sigma.
    If shift_energy=True, U(cutoff) is subtracted from each included pair.
    This makes the energy continuous at cutoff, but does not shift forces.
    """

    def __init__(self, epsilon: float = 1.0, sigma: float = 1.0, cutoff: float | None = None, shift_energy: bool = True):
        
        if not np.isfinite(epsilon) or epsilon <= 0:
            raise ValueError("epsilon must be finite and positive.")
        
        if not np.isfinite(sigma) or sigma <= 0:
            raise ValueError("sigma must be finite and positive.")
        
        if cutoff is not None and (not np.isfinite(cutoff) or cutoff <= 0):
            raise ValueError("cutoff must be finite and positive, or None.")
        
        self.epsilon = epsilon
        self.sigma = sigma
        self.cutoff = cutoff
        self.shift_energy = shift_energy

    def compute(self, state: SystemState) -> tuple[np.ndarray, float]:
        
        x = np.asarray(state.xs, dtype=float)
        if x.ndim != 2 or x.shape[1] != 3 or not np.isfinite(x).all():
            raise ValueError("Positions must be a finite array of shape (N, 3).")

        box = None
        if state.box is not None:
            box = np.asarray(state.box, dtype=float)
            
            if box.shape != (3,) or not np.isfinite(box).all() or np.any(box <= 0):
                raise ValueError("Box must contain three finite positive lengths.")
            
            if self.cutoff is not None and self.cutoff >= 0.5 * box.min():
                raise ValueError("cutoff must be smaller than half the shortest box length.")

        U, dU_dx = interactions_LJ_numba(x, box, sigma=self.sigma, eps=self.epsilon, cutoff=self.cutoff, shift_energy=self.shift_energy)
        return -dU_dx, float(U)
    
class SumInteraction(Interaction):
    
    def __init__(self, interactions):
        self.interactions = interactions
        
    def compute(self, state: SystemState) -> tuple[np.ndarray, float]:
        
        forces = np.zeros_like(state.xs)
        potential_energy = 0.0
        
        for interaction in self.interactions:
            f, U = interaction.compute(state)
            forces += f
            potential_energy += U
            
        return forces, potential_energy

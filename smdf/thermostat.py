from abc import ABC, abstractmethod
import numpy as np

from .state import SystemState

class Thermostat(ABC):
    """Base class for thermostats."""
    
    @abstractmethod
    def apply(self, state: SystemState) -> None:
        pass
    
class NoThermostat(Thermostat):
    """A thermostat that does nothing."""
    
    def apply(self, state: SystemState) -> None:
        pass
    
class AndersenThermostat(Thermostat):
    """Andersen thermostat implementation."""
    
    def __init__(self, T: float, tau: float, dt: float, ms: np.ndarray, rng=None):
        self.rng = rng if rng is not None else np.random.default_rng()
        self.T = T
        self.tau = tau
        self.dt = dt
        self.ms = ms
        
        self.p_resample = dt / tau
        self.scale_p = np.sqrt(T * np.repeat(ms, 3, axis=1))
    
    def apply(self, state: SystemState) -> None:
        """Apply the Andersen thermostat to the system state."""
        
        idx = np.where(self.rng.random(state.N) < self.p_resample)
        state.ps[idx] = self.rng.normal(0, self.scale_p[idx])
        
class LangevinThermostat(Thermostat):
    """Langevin thermostat implementation."""
    
    def __init__(self, T: float, tau: float, dt: float, ms: np.ndarray, rng=None):
        
        self.rng = rng if rng is not None else np.random.default_rng()
        self.T = T
        self.tau = tau
        self.dt = dt
        
        self.gamma = 2 / tau
        self.A = np.exp(-self.gamma * dt)
        self.B = np.sqrt((1-self.A**2) * T * ms)
    
    def apply(self, state: SystemState) -> None:
        """Apply the Langevin thermostat to the system state."""
        state.ps *= self.A
        state.ps += self.B * self.rng.normal(0, 1, size=state.ps.shape)
        
class BerendsenThermostat(Thermostat):
    """Berendsen thermostat implementation."""
    
    def __init__(self, T: float, tau: float, dt: float):
        raise ValueError("Nope, don't do it. It does not sample the canonical ensemble, avoid it.")
    
    def apply(self, state: SystemState) -> None:
        """Apply the Berendsen thermostat to the system state."""
        raise ValueError("Nope, don't do it. It does not sample the canonical ensemble, avoid it.")
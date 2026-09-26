from abc import ABC, abstractmethod

import numpy as np

from .interactions import Interaction
from .state import SystemState
from .thermostat import NoThermostat, Thermostat

class Integrator(ABC):
    """Base class for integrators that update a SystemState in place."""
    
    def __init__(self, dt, interaction: Interaction, thermostat: Thermostat | None = None):
        
        if not np.isfinite(dt) or dt <= 0:
            raise ValueError("dt must be finite and positive.")
        if thermostat is not None and getattr(thermostat, "dt", dt) != dt:
            raise ValueError("Integrator and thermostat must use the same dt.")
        
        self.dt = dt
        self.interaction = interaction
        self.thermostat = thermostat if thermostat is not None else NoThermostat()
        self.forces, self.potential_energy = None, None
        
    def initialize(self, state: SystemState) -> None:
        self.forces, self.potential_energy = self.interaction.compute(state)

    @abstractmethod
    def step(self, state: SystemState) -> None:
        """Advance state using its current forces"""
        pass


class BAOAB(Integrator):
    """B/2 -> A/2 -> O -> A/2 -> B/2 splitting."""

    def step(self, state: SystemState) -> None:
        
        if self.forces is None:
            self.initialize(state)
        
        state.ps += 0.5 * self.dt * self.forces
        state.xs += 0.5 * self.dt * state.ps / state.ms

        self.thermostat.apply(state)

        state.xs += 0.5 * self.dt * state.ps / state.ms

        self.forces, self.potential_energy = self.interaction.compute(state)
        state.ps += 0.5 * self.dt * self.forces

class VelocityVerlet(Integrator):
    """Velocity Verlet followed by an optional full-dt thermostat step."""

    def step(self, state: SystemState) -> None:
        if self.forces is None:
            self.initialize(state)

        state.ps += 0.5 * self.dt * self.forces
        state.xs += self.dt * state.ps / state.ms

        self.forces, self.potential_energy = self.interaction.compute(state)
        state.ps += 0.5 * self.dt * self.forces

        self.thermostat.apply(state)


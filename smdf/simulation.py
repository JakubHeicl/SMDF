from .io import Output
from .state import SystemState
from .integrator import Integrator


class Simulation:
    """A state and its integrator, the integrator owns interaction and thermostat."""

    def __init__(self, state: SystemState, integrator: Integrator):
        self.state = state
        self.integrator = integrator
        self.step_count = 0
        self.integrator.initialize(state)

    @property
    def time(self):
        return self.step_count * self.integrator.dt

    def advance(self):
        self.integrator.step(self.state)
        self.step_count += 1

    def run(self, n_steps: int, output: Output | None = None):
        """Run additional steps. Always output the start and final state once."""
        
        if type(n_steps) is not int or n_steps < 0:
            raise ValueError("n_steps must be a non-negative integer.")
        
        if output is not None:
            output.write(self.state, self.step_count, self.time, self.integrator.potential_energy, force=True)
            
        for _ in range(n_steps):
            self.advance()
            if output is not None:
                output.write(self.state, self.step_count, self.time, self.integrator.potential_energy)
                
        if output is not None:
            output.write(self.state, self.step_count, self.time, self.integrator.potential_energy, force=True)
            output.write_final(self.state)

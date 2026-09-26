import csv
import json
import shutil
import tomllib
from pathlib import Path
from typing import TextIO

import numpy as np

from .state import SystemState


def read_toml(path: str | Path) -> dict:
    
    with Path(path).open("rb") as f_in:
        return tomllib.load(f_in)

def read_xyz_frame(f_in: TextIO) -> tuple[list[str], np.ndarray] | None:
    """Read one ordinary XYZ frame, return None at the end of the file."""
    
    line = f_in.readline()
    
    if line == "":
        return None
    
    try:
        N = int(line)
    except ValueError:
        raise ValueError("XYZ line 1: expected the number of particles.") from None
    
    if N <= 0:
        raise ValueError("XYZ line 1: particle count must be positive.")
    
    if f_in.readline() == "":
        raise ValueError("XYZ line 2: missing comment line.")

    labels = []
    xs = np.empty((N, 3))
    
    for i in range(N):
        
        items = f_in.readline().split()
        
        if len(items) != 4:
            raise ValueError(f"XYZ frame line {i + 3}: expected label and three coordinates.")
        
        labels.append(items[0])
        
        try:
            xs[i] = [float(value) for value in items[1:]]
            
        except ValueError:
            raise ValueError(f"XYZ frame line {i + 3}: invalid coordinates.") from None
        
        if not np.isfinite(xs[i]).all():
            raise ValueError(f"XYZ frame line {i + 3}: coordinates must be finite.")
        
    return labels, xs

def read_xyz(path: str | Path) -> tuple[list[str], np.ndarray]:
    """Read the first frame, in reduced units. Further frames are ignored."""
    
    with Path(path).open(encoding="utf-8") as f_in:
        frame = read_xyz_frame(f_in)
        
    if frame is None:
        raise ValueError(f"Empty XYZ file: {path}")
    
    return frame

def write_xyz_frame(f_out: TextIO, labels: list[str], xs: np.ndarray, comment: str = ""):
    """Write one ordinary XYZ frame."""

    if np.shape(xs) != (len(labels), 3):
        raise ValueError("XYZ coordinates must have shape (number of labels, 3).")
    
    f_out.write(f"{len(labels)}\n{comment}\n")
    
    for label, x in zip(labels, xs):
        f_out.write(f"{label} {x[0]:.16g} {x[1]:.16g} {x[2]:.16g}\n")


class Output:
    """Write energies and XYZ frames without retaining the trajectory in memory."""

    def __init__(self, directory: str | Path, thermo_stride: int = 10, trajectory_stride: int = 100):
        
        for stride in (thermo_stride, trajectory_stride):
            if type(stride) is not int or stride <= 0:
                raise ValueError("Output strides must be positive integers.")
            
        self.directory = Path(directory)
        self.thermo_stride = thermo_stride
        self.trajectory_stride = trajectory_stride
        
        self.last_thermo = None
        self.last_trajectory = None
        self.fs_out = {}

    def __enter__(self):
        """Create the output directory and open files for writing."""
        
        self.directory.mkdir(parents=True, exist_ok=False)
        
        try:
            self.fs_out['thermo'] = (self.directory / 'thermo.csv').open('w', newline='', encoding='utf-8')
            self.fs_out['trajectory'] = (self.directory / 'trajectory.xyz').open('w', encoding='utf-8')
            self.writer = csv.writer(self.fs_out['thermo'])
            self.writer.writerow(['step', 'time', 'kinetic', 'potential', 'total', 'temperature'])
        except Exception:
            
            self.__exit__(None, None, None)
            raise
        
        self.last_thermo = self.last_trajectory = None
        
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        for f_out in self.fs_out.values():
            f_out.close()
        self.fs_out = {}

    def save_setup(self, input_path: str | Path, settings: dict, state: SystemState):
        
        """Store original input and resolved settings, including generated box."""
        
        shutil.copyfile(input_path, self.directory / 'input.toml')
        
        metadata = dict(settings = settings, N = state.N,
                        box = None if state.box is None else state.box.tolist(),
                        degrees_of_freedom = 3 * state.N, units='reduced')
        
        with (self.directory / 'metadata.json').open('w', encoding='utf-8') as f_out:
            json.dump(metadata, f_out, indent=2, allow_nan=False)
            
        with (self.directory / 'initial.xyz').open('w', encoding='utf-8') as f_out:
            write_xyz_frame(f_out, state.labels, state.xs, 'step=0 time=0 units=reduced')

    def write(self, state: SystemState, step: int, time: float, potential_energy: float, force=False):
        
        """Write the current state to the output files."""
        
        if not self.fs_out:
            raise RuntimeError("Use Output inside a with block.")
        
        if self.last_thermo != step and (force or step % self.thermo_stride == 0):
            
            kinetic = float(0.5 * np.sum(state.ps**2 / state.ms))
            temperature = 2 * kinetic / (3 * state.N)
            self.writer.writerow([step, time, kinetic, potential_energy, kinetic + potential_energy, temperature])
            self.last_thermo = step
            
        if self.last_trajectory != step and (force or step % self.trajectory_stride == 0):
            
            comment = f'step={step} time={time:.16g} units=reduced'
            write_xyz_frame(self.fs_out['trajectory'], state.labels, state.xs, comment)
            self.last_trajectory = step

    def write_final(self, state: SystemState):
        
        with (self.directory / 'final.xyz').open('w', encoding='utf-8') as f_out:
            write_xyz_frame(f_out, state.labels, state.xs, 'units=reduced')

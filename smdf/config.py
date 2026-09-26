from copy import deepcopy
from pathlib import Path

import numpy as np

from .io import read_toml, read_xyz
from .state import SystemState
from .interactions import HarmonicPotential, LennardJonesPotential, SumInteraction
from .integrator import BAOAB, VelocityVerlet
from .thermostat import NoThermostat, AndersenThermostat, LangevinThermostat, BerendsenThermostat
from .simulation import Simulation


DEFAULTS = {
    'simulation': dict(units='reduced', dt=0.005, steps=1000, seed=42),
    'initial': dict(source='xyz', temperature=1.0),
    'system': dict(mass=1.0),
    'interaction': dict(type='lennard_jones'),
    'integrator': dict(type='velocity_verlet'),
    'thermostat': dict(type='none'),
    'output': dict(directory='results/run', thermo_stride=10, trajectory_stride=100),
}

def check_keys(settings, allowed, section):
    
    unknown = set(settings) - set(allowed)
    if unknown:
        raise ValueError(f"Unknown settings in {section}: {', '.join(sorted(unknown))}")


def number(value, name, minimum=0, allow_equal=False):
    
    if type(value) not in (int, float) or not np.isfinite(value):
        raise ValueError(f"{name} must be a finite number.")
    
    if value < minimum or (value == minimum and not allow_equal):
        comparison = '>=' if allow_equal else '>'
        raise ValueError(f"{name} must be {comparison} {minimum}.")


def integer(value, name, minimum=0):
    
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}.")

def prepare_interaction(potential, name='interaction'):
    """Fill defaults and validate one interaction configuration."""
    
    if not isinstance(potential, dict):
        raise ValueError(f"{name} must be a TOML table.")
    
    if 'type' not in potential:
        raise ValueError(f"{name}.type is required.")
    
    if potential['type'] == 'lennard_jones':
        
        defaults = dict(epsilon=1.0, sigma=1.0, cutoff=None, shift_energy=True)
        
        for key, value in defaults.items():
            potential.setdefault(key, value)
            
        check_keys(potential, ['type', *defaults], name)
        number(potential['epsilon'], f'{name}.epsilon')
        number(potential['sigma'], f'{name}.sigma')
        
        if potential['cutoff'] is not None:
            number(potential['cutoff'], f'{name}.cutoff')
            
        if type(potential['shift_energy']) is not bool:
            raise ValueError(f"{name}.shift_energy must be true or false.")
        
    elif potential['type'] == 'harmonic':
        potential.setdefault('k', 1.0)
        check_keys(potential, ['type', 'k'], name)
        number(potential['k'], f'{name}.k', allow_equal=True)
        
    elif potential['type'] == 'sum':
        check_keys(potential, ['type', 'terms'], name)
        terms = potential.get('terms')
        
        if not isinstance(terms, list) or len(terms) == 0:
            raise ValueError(f"{name}.terms must be a non-empty array of tables.")
        
        for i, term in enumerate(terms):
            prepare_interaction(term, f'{name}.terms[{i}]')
            
    else:
        raise ValueError(f"{name}.type must be lennard_jones, harmonic, or sum.")


def contains_interaction(potential, interaction_type):
    if potential['type'] == interaction_type:
        return True
    if potential['type'] == 'sum':
        return any(contains_interaction(term, interaction_type) for term in potential['terms'])
    return False

def build_interaction(potential):
    if potential['type'] == 'lennard_jones':
        return LennardJonesPotential(epsilon=potential['epsilon'], sigma=potential['sigma'], cutoff=potential['cutoff'], shift_energy=potential['shift_energy'])
    
    if potential['type'] == 'harmonic':
        return HarmonicPotential(k=potential['k'])
    
    return SumInteraction([build_interaction(term) for term in potential['terms']])


def load_config(path: str | Path) -> dict[str, object]:
    """Read settings, fill defaults, and resolve paths relative to the TOML."""
    
    path = Path(path).resolve()
    raw = read_toml(path)
    
    check_keys(raw, DEFAULTS, 'input')
    
    settings = deepcopy(DEFAULTS)
    
    for section, values in raw.items():
        
        if not isinstance(values, dict):
            raise ValueError(f"{section} must be a TOML table.")
        
        settings[section].update(values)

    initial = settings['initial']
    source = initial['source']
    
    if source == 'xyz':
        initial_defaults = dict(path=None)
        
    elif source == 'lattice':
        initial_defaults = dict(n_side=5, density=0.5, label='Ar')
        
    else:
        raise ValueError("initial.source must be xyz or lattice.")
    
    for key, value in initial_defaults.items():
        initial.setdefault(key, value)

    potential = settings['interaction']
    prepare_interaction(potential)

    thermostat = settings['thermostat']
    
    if thermostat['type'] not in ('none', 'andersen', 'langevin', 'berendsen'):
        raise ValueError("thermostat.type must be none, andersen, langevin, or berendsen.")
    
    if thermostat['type'] != 'none':
        thermostat.setdefault('temperature', 1.0)
        thermostat.setdefault('tau', 1.0)

    allowed = {
        'simulation': DEFAULTS['simulation'],
        'initial': ['source', 'temperature', *initial_defaults],
        'system': ['mass', 'box'],
        'integrator': ['type'],
        'thermostat': ['type'] if thermostat['type'] == 'none' else ['type', 'temperature', 'tau'],
        'output': DEFAULTS['output'],
    }
    
    for section in allowed:
        check_keys(settings[section], allowed[section], section)

    sim = settings['simulation']
    
    if sim['units'] != 'reduced':
        raise ValueError("Only reduced units are supported.")
    
    number(sim['dt'], 'simulation.dt')
    integer(sim['steps'], 'simulation.steps')
    integer(sim['seed'], 'simulation.seed')
    number(initial['temperature'], 'initial.temperature', allow_equal=True)
    number(settings['system']['mass'], 'system.mass')
    
    if settings['integrator']['type'] not in ('baoab', 'velocity_verlet'):
        raise ValueError("integrator.type must be baoab or velocity_verlet.")
    
    if thermostat['type'] != 'none':
        number(thermostat['temperature'], 'thermostat.temperature', allow_equal=True)
        number(thermostat['tau'], 'thermostat.tau')
        
        if thermostat['type'] == 'andersen' and sim['dt'] > thermostat['tau']:
            raise ValueError("Andersen requires dt <= tau (collision probability is dt/tau).")
        
    output = settings['output']
    
    integer(output['thermo_stride'], 'output.thermo_stride', 1)
    integer(output['trajectory_stride'], 'output.trajectory_stride', 1)
    
    if source == 'lattice':
        integer(initial['n_side'], 'initial.n_side', 1)
        number(initial['density'], 'initial.density')
        
        if not isinstance(initial['label'], str) or len(initial['label'].split()) != 1:
            raise ValueError("initial.label must be a non-empty word.")
        
        if 'box' in settings['system']:
            raise ValueError("A lattice derives its box from density; omit system.box.")
        
    if 'box' in settings['system']:
        box = np.asarray(settings['system']['box'], dtype=float)
        
        if box.shape != (3,) or not np.isfinite(box).all() or np.any(box <= 0):
            raise ValueError("system.box must contain three positive finite lengths.")
        
    if contains_interaction(potential, 'harmonic') and (source == 'lattice' or 'box' in settings['system']):
        raise ValueError("The harmonic potential requires XYZ input without a periodic box.")

    paths = [(output, 'directory')]
    
    if source == 'xyz':
        paths.append((initial, 'path'))
        
    for section, key in paths:
        if not isinstance(section[key], str) or not section[key].strip():
            raise ValueError(f"{key} must be a non-empty path.")
        section[key] = str((path.parent / section[key]).resolve())
    return settings

def build_state(settings, rng):
    """Build the initial system state."""
    
    initial = settings['initial']
    system = settings['system']
    
    if initial['source'] == 'xyz':
        labels, xs = read_xyz(initial['path'])
        box = None if 'box' not in system else np.array(system['box'], dtype=float)
        
    else:
        n = initial['n_side']
        L = (n**3 / initial['density'])**(1 / 3)
        axis = (np.arange(n, dtype=float) + 0.5) * L / n
        xs = np.stack(np.meshgrid(axis, axis, axis, indexing='ij'), axis=-1).reshape(-1, 3)
        labels = [initial['label']] * n**3
        box = np.full(3, L)
        
    N = len(labels)
    ms = np.full((N, 1), system['mass'], dtype=float)
    ps = rng.normal(size=(N, 3)) * np.sqrt(initial['temperature'] * ms)
    return SystemState(N=N, xs=xs, ps=ps, ms=ms, labels=labels, box=box)

def build_simulation(settings):
    
    rng = np.random.default_rng(settings['simulation']['seed'])
    state = build_state(settings, rng)
    interaction = build_interaction(settings['interaction'])
    dt = settings['simulation']['dt']
    t = settings['thermostat']
    if t['type'] == 'none':
        thermostat = NoThermostat()
    elif t['type'] == 'andersen':
        thermostat = AndersenThermostat(t['temperature'], t['tau'], dt, state.ms, rng=rng)
    elif t['type'] == 'langevin':
        thermostat = LangevinThermostat(t['temperature'], t['tau'], dt, state.ms, rng=rng)
    else:
        thermostat = BerendsenThermostat(t['temperature'], t['tau'], dt)
    if settings['integrator']['type'] == 'baoab':
        integrator = BAOAB(dt, interaction, thermostat)
    else:
        integrator = VelocityVerlet(dt, interaction, thermostat)
    return Simulation(state, integrator)

import argparse
from pathlib import Path

from .config import load_config, build_simulation
from .io import Output

def main():
    
    parser = argparse.ArgumentParser(description='Simple molecular dynamics framework')
    parser.add_argument('command', choices=['run', 'check'])
    parser.add_argument('input', type=Path, help='TOML configuration file')
    parser.add_argument('--output', type=Path, help='Override output directory (relative to current directory)')
    args = parser.parse_args()
    
    try:
        settings = load_config(args.input)
        
        if args.output is not None:
            settings['output']['directory'] = str(args.output.resolve())
            
        simulation = build_simulation(settings)
        
        print(f'Particles: {simulation.state.N}')
        print(f'Initial potential energy: {simulation.integrator.potential_energy:.10g}')
        
        if args.command == 'check':
            print('Input and initial interactions are valid. No files written.')
            return
        
        out = settings['output']
        
        with Output(out['directory'], out['thermo_stride'], out['trajectory_stride']) as output:
            output.save_setup(args.input, settings, simulation.state)
            simulation.run(settings['simulation']['steps'], output)
            
        print(f'Finished {simulation.step_count} steps, time = {simulation.time:g}')
        print(f'Results: {out["directory"]}')
        
    except (ValueError, OSError, TypeError, OverflowError) as error:
        parser.exit(1, f'Error: {error}\n')

if __name__ == '__main__':
    main()

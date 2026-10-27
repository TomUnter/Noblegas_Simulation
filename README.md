# Noble-Gas Particle Simulation

An educational Python project developed during a practical session
to explore molecular dynamics and the numerical simulation of
interacting noble-gas particles.

## Overview

The simulation models particles in three dimensions using
Lennard–Jones-type pair interactions and repulsive wall forces.
It includes:

- Random initial positions selected from a discrete grid.
- Normally distributed initial velocity components.
- Pair-interaction calculations with a distance cutoff.
- A Verlet-style integration routine.
- Berendsen-style velocity rescaling and gradual cooling.
- Export of particle coordinates and velocities.

## Files

- `AtomSimulationTomUnterleiter.py` — main simulation script,
  including an experimental, currently unused pressure function.
- `vimscript.py` — an alternative simulation variant.

## Requirements and execution

Dependencies: Python, NumPy, and SciPy.

```bash
pip install numpy scipy
python AtomSimulationTomUnterleiter.py
```

The main script starts a simulation with 150 particles.
Parameters can be changed in the final `AtomSim(...)` call.
Additional settings, including wall positions, iteration count,
and cooling schedule, are defined inside the code.

## Output

- `coords.pdb` — multi-frame, PDB-style coordinate output.
- `vels.txt` — particle velocity components for saved frames.

Output files are written to the current working directory and
overwritten on subsequent runs.

## Project status

This repository preserves a practical-course implementation.
It is not a validated molecular-dynamics package.

The force convention, integration sequence, cached-distance updates,
unit consistency, and boundary handling require review before
quantitative physical interpretation. The pressure function is
unfinished and is not used in the simulation.

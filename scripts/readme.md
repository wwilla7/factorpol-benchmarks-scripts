This directory contains scripts to setup benchmark calculations with factorpol electrostatic model and OpenFF Sage valence and Lennard-Jones (LJ) parameters.

1. `system.py`: a script to parameterize systems and set up liquid boxes with factorpol electrostatic model, OpenFF Sage valence parameters, and OpenFF Sage LJ parameters.

2. `run.py`: a script to simulate liquid boxes with OpenMM and MPIDOpenMMPlugin.

3. `setup.py`: a simple script to set up all systems included in this benchmark.

4. `analysis.py`: a script to post-process simulation trajectories for densities, dielectric constants, and simulation speeds.

5. `system.csv`: a csv file to specify systems to be simulated.
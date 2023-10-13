#!/usr/bin/env python
# coding: utf-8

import os
import pickle
from sys import stdout

import mdtraj
import mpidplugin
import numpy as np
import pandas as pd
import pint
from mdtraj.reporters import NetCDFReporter
from openmm import *
from openmm.app import *
from openmm.unit import *
from simtk.unit import *

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity


cwd = os.getcwd()

system_path = os.path.join(cwd, "system")
simulation_path = os.path.join(cwd, "simulation")

os.makedirs(simulation_path, exist_ok=True)

pdb = PDBFile(os.path.join(system_path, "liquid.pdb"))

forcefield = ForceField(os.path.join(system_path, "forcefield.xml"))
system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=8 * angstrom,
    constraints=HBonds,
    defaultTholeWidth=8,
    polarization="direct",
)

temp = TEMP

temperature = float(temp) * kelvin
integrator = LangevinIntegrator(temperature, 1 / picosecond, 2 * femtoseconds)
system.addForce(MonteCarloBarostat(1 * atmosphere, temperature, 25))

os.chdir(simulation_path)

with open("system.xml", "w") as f:
    f.write(XmlSerializer.serializeSystem(system))

num_particles = system.getNumParticles()

forces = {
    system.getForce(i).__class__.__name__: system.getForce(i)
    for i in range(system.getNumForces())
}

if forces["Force"]:
    pass
else:
    print("I don't see MPID Force here!")
    exit()

myplatform = Platform.getPlatformByName("CUDA")
deviceid = "CUDAID" 
myproperties = {"DeviceIndex": deviceid, "Precision": "mixed"}

simulation = Simulation(pdb.topology, system, integrator, myplatform, myproperties)

context = simulation.context
if pdb.topology.getPeriodicBoxVectors():
    context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

# Initialize
context.setPositions(pdb.positions)
simulation.minimizeEnergy()
# run 5ns of equilibration
nsteps = 2500000
# Dump simulation info every 10ps
simulation.reporters.append(
    StateDataReporter(
        "equilibration.log",
        5000,
        progress=True,
        totalSteps=nsteps,
        volume=True,
        step=True,
        potentialEnergy=True,
        totalEnergy=True,
        temperature=True,
        density=True,
        speed=True,
    )
)
simulation.reporters.append(PDBReporter("equilibrated.pdb", nsteps))
simulation.step(nsteps)

# clear simulation reporter
simulation.reporters.clear()

# # run 15ns of equilibration
# nsteps = 7500000
# run 5ns of equilibration
nsteps = 2500000 
# Dump simulation info every 10ps
simulation.reporters.append(
    StateDataReporter(
        "production.log",
        5000,
        progress=True,
        totalSteps=nsteps,
        volume=True,
        step=True,
        potentialEnergy=True,
        totalEnergy=True,
        temperature=True,
        density=True,
        speed=True,
    )
)

# Dump trajectory info every 10ps
simulation.reporters.append(NetCDFReporter("trajectory.nc", 5000))
simulation.reporters.append(PDBReporter("production.pdb", nsteps))
simulation.step(nsteps)
simulation.saveState("restart.xml")

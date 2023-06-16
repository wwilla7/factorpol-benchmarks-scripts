#!/usr/bin/env python
# coding: utf-8

import os
import json
import pickle
import shutil

import pint
from sys import argv

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

import numpy as np
import pandas as pd
import parmed as pmd
import mdtraj
from openmm import *
from openmm.app import *
from openmm.unit import *
from simtk.unit import *
import mpidplugin

cwd = os.getcwd()

simulation_path = os.path.join(cwd, "simulation")

smiles = SMILE
temp = TEMP
temperature = float(temp) * kelvin

os.chdir(simulation_path)

with open("system.xml", "r") as f:
    dt = f.read()
system = XmlSerializer.deserializeSystem(dt)

pdb = PDBFile("production.pdb")

num_residues = pdb.topology.getNumResidues()
num_particles = system.getNumParticles()

integrator = LangevinIntegrator(temperature, 1 / picosecond, 2 * femtoseconds)
myplatform = Platform.getPlatformByName("CUDA")
deviceid = "0"
myproperties = {"DeviceIndex": deviceid, "Precision": "mixed"}
simulation = Simulation(pdb.topology, system, integrator, myplatform, myproperties)

context = simulation.context
if pdb.topology.getPeriodicBoxVectors():
    context.setPeriodicBoxVectors(*pdb.topology.getPeriodicBoxVectors())

simulation.loadState("restart.xml")

forces = {
    system.getForce(i).__class__.__name__: system.getForce(i)
    for i in range(system.getNumForces())
}

if mpidplugin.MPIDForce.isinstance(forces["Force"]):
    mpidforce = mpidplugin.MPIDForce.cast(forces["Force"])
alphas = [
    mpidforce.getMultipoleParameters(i)[-1][0] for i in range(system.getNumParticles())
]
charges = [
    mpidforce.getMultipoleParameters(i)[0] for i in range(system.getNumParticles())
]

# calculate dielectric constant
df = pd.read_csv("production.log")
volume = Q_(df["Box Volume (nm^3)"].values, "nm**3").to("a0**3")
average_speed = np.mean(df["Speed (ns/day)"].values)

# $\epsilon_{infty}$
eps0 = Q_(8.854187812e-12, "F/m").to("e**2/a0/hartree")
sum_alphas = Q_(np.sum(alphas), "nm**3").to("a0**3")
prefactor00 = 1 / eps0
mdtraj_topology = mdtraj.load_frame("production.pdb", 0).topology

average_volumes = np.mean(volume)
average_density = Q_(np.mean(df["Density (g/mL)"].values), "g/mL")

dipole_moments = []

for idx2, chunk in enumerate(
    mdtraj.iterload("trajectory.nc", top=mdtraj_topology, chunk=1)
):
    context.setPeriodicBoxVectors(*chunk.openmm_boxes(0))
    context.setPositions(chunk.openmm_positions(0))

    dipole_moment = (
        Q_(mpidforce.getSystemMultipoleMoments(context)[1:4], "debye")
        .to("e*nm")
        .magnitude
    )
    dipole_moments.append(dipole_moment)

### fluctuation of dipoles
avg_sqr_dipole_moments = Q_(
    np.sum(np.mean(np.square(dipole_moments), axis=0)), "e**2 * nm**2"
).to("e**2 * a0**2")
avg_dipole_moment_sqr = Q_(
    np.sum(np.square(np.mean(dipole_moments, axis=0))), "e**2 * nm**2"
).to("e**2 * a0**2")
dipole_variance = avg_sqr_dipole_moments - avg_dipole_moment_sqr

prefactor = (
    4
    * np.pi
    / (
        3.0
        * Q_(1, "boltzmann_constant").to("J/kelvin")
        * Q_(temperature / kelvin, "kelvin")
    ).to("hartree")
)

eps = prefactor * (dipole_variance / average_volumes)
eps_infty = prefactor00 * (sum_alphas / average_volumes)

dipole_per_frame = (
    Q_(np.linalg.norm(dipole_moments, axis=1), "e*nm").to("debye").magnitude
)
standard_deviation_of_dipoles = np.std(dipole_per_frame)

data = {
    "molecule": smiles,
    "temperature (K)": temp,
    "density(g/mL)": average_density.magnitude,
    "high_frequency_epsilon": eps_infty.magnitude + 1,
    "epsilon": eps.magnitude + 1 + eps_infty.magnitude,
    "speed": average_speed,
    "num_particles": num_particles,
    "std_of_u": standard_deviation_of_dipoles,
}

json.dump(data, open("results.json", "w"), indent=2)

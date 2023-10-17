#!/usr/bin/env python

import os
import shutil

import parmed as pmd
from lxml import etree
from openff.evaluator.protocols.coordinates import BuildCoordinatesPackmol
from openff.evaluator.substances import Substance
from openff.toolkit.topology import Molecule as off_Molecule
from openff.toolkit.topology import Topology as off_Topology
from openff.toolkit.typing.engines.smirnoff import ForceField as off_ForceField
from openmm.app import ForceField
from openmmforcefields.generators import (GAFFTemplateGenerator,
                                          SMIRNOFFTemplateGenerator)
from parmed.openmm.parameters import OpenMMParameterSet

smile = SMILE

cwd = os.getcwd()

# sage
off_forcefield = off_ForceField("openff-2.0.0.offxml")

output_path = os.path.join(cwd, "system")

ff_name = "sage"
off_mol = off_Molecule.from_smiles(smile)
off_mol.generate_conformers(n_conformers=1)
off_mol.assign_partial_charges("am1bcc")
off_topology = off_Topology.from_molecules(off_mol)
omm_vac_topology = off_topology.to_openmm()

if os.path.exists(output_path):
    pass
else:
    os.makedirs(output_path)

# Create Openmm system
forcefield = ForceField()

if ff_name.lower() in ["gaff"]:
    gaff = GAFFTemplateGenerator(molecules=[off_mol])
    forcefield.registerTemplateGenerator(gaff.generator)
    omm_vac_system = forcefield.createSystem(omm_vac_topology)

elif ff_name.lower() in ["sage"]:
    offff = SMIRNOFFTemplateGenerator(molecules=[off_mol])
    forcefield.registerTemplateGenerator(offff.generator)
    omm_vac_system = forcefield.createSystem(omm_vac_topology)

pmd_structure = pmd.openmm.load_topology(omm_vac_topology, system=omm_vac_system)

openmm_params = OpenMMParameterSet()
openmm_xml = openmm_params.from_structure(pmd_structure)
get_residues = pmd.modeller.ResidueTemplateContainer.from_structure(
    pmd_structure
).to_library()
openmm_xml.residues.update(get_residues)
openmm_xml.write(os.path.join(output_path, "forcefield.xml"))

with open(os.path.join(output_path, "forcefield.xml"), "r") as f:
    ff_text = f.read()

root = etree.fromstring(ff_text)

sh = root.find("NonbondedForce")

sh.set("coulomb14scale", "0.83333")
sh.set("lj14scale", "0.5")

# organize XML file
tree = etree.ElementTree(root)
xml = etree.tostring(tree, encoding="utf-8", pretty_print=True).decode("utf-8")
xml = xml.replace("><", ">\n\t<")
xml = xml.replace("\t</ForceField>", "</ForceField>")
with open(os.path.join(output_path, "forcefield.xml"), "w") as f:
    f.write(xml)

off_mol.to_file(os.path.join(output_path, "molecule.mol2"), file_format="mol2")
off_mol.to_file(os.path.join(output_path, "molecule.pdb"), file_format="pdb")

bcrd = BuildCoordinatesPackmol("")
bcrd.substance = Substance.from_components(smile)
bcrd.max_molecules = 256
bcrd.execute(output_path)
shutil.copy(f"{output_path}/output.pdb", f"{output_path}/liquid.pdb")

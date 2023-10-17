#!/usr/bin/env python

import os
import shutil

import numpy as np
from openff.evaluator.protocols.coordinates import BuildCoordinatesPackmol
from openff.evaluator.substances import Substance
from openff.toolkit import ForceField
from pkg_resources import resource_filename

from factorpol.parameterization import parameterize_molecule
from factorpol.utilities import BondChargeCorrections, Polarizability

smile = SMILE

cwd = os.getcwd()


# # element polarizabilities
# off_forcefield = ForceField(
#     resource_filename("factorpol", os.path.join("data", "off_examples.offxml"))
# )
# polarizability = Polarizability(data_source=resource_filename("factorpol", os.path.join("data", "alphas", "alphas_element.csv")))
# bcc_dpol_library = BondChargeCorrections(data_source=resource_filename("factorpol", os.path.join("data", "bccs", "bcc_element.csv")))

# # sage polarizabilities
off_forcefield = ForceField("openff-2.0.0.offxml")
polarizability = Polarizability(
    data_source=resource_filename(
        "factorpol", os.path.join("data", "alphas", "alphas_sageLJ.csv")
    )
)
bcc_dpol_library = BondChargeCorrections(
    data_source=resource_filename(
        "factorpol", os.path.join("data", "bccs", "bcc_sageLJ.csv")
    )
)

output_path = os.path.join(cwd, "system")

parm = parameterize_molecule(
    smile=smile,
    ff_name="sage",
    polarizability=polarizability,
    BCCLibrary=bcc_dpol_library,
    output_path=output_path,
    off_forcefield=off_forcefield,
)


bcrd = BuildCoordinatesPackmol("")
bcrd.substance = Substance.from_components(smile)
bcrd.max_molecules = 256
bcrd.execute(output_path)
shutil.copy(f"{output_path}/output.pdb", f"{output_path}/liquid.pdb")

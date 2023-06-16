import pandas as pd
import copy
import os, shutil

dt = pd.read_csv("system.csv", sep=",", index_col="SMILES")

smiles = dt.index
cwd = os.getcwd()

with open("system.py", "r") as f:
    csys = f.read()

with open("run.py", "r") as f:
    runt = f.read()

with open("analysis.py", "r") as f:
    anly = f.read()

for i, smile in enumerate(smiles):
    temp = dt.loc[smile, "Temperature (K)"]
    system_path = os.path.join(cwd, "systems", f"molecule{i:02d}")
    os.makedirs(system_path, exist_ok=True)
    os.chdir(system_path)
    tcsys = csys.replace("SMILE", f"'{smile}'")
    trunt = runt.replace("TEMP", f"{temp}")
    tanly = anly.replace("TEMP", f"{temp}")
    tanly = tanly.replace("SMILE", f"'{smile}'")

    with open("system.py", "w") as f:
        f.write(tcsys)

    with open("run.py", "w") as f:
        f.write(trunt)

    with open("analysis.py", "w") as f:
        f.write(tanly)

    os.chdir(cwd)

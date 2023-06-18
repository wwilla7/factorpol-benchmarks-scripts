import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from openff.toolkit import ForceField, Molecule


def molecular_polarizability(
    smiles: str, forcefield: ForceField, pol_handler: str = "MPIDPolarizability"
):
    """
    Calculate the molecular polarizability of a molecule using typed atomic polarizabilities.

    Parameters
    ----------
    smiles : str
    forcefield : ForceField
    pol_handler : str, optional

    Returns
    -------
    float

    """

    offmol = Molecule.from_smiles(smiles, allow_undefined_stereo=True)
    mol_params = forcefield.label_molecules(offmol.to_topology())[0][pol_handler]
    pol_values = [p.polarizabilityXX for p in mol_params.values()]

    mol_pol = sum(pol_values)
    return mol_pol.to("angstrom ** 3").magnitude


if __name__ == "__main__":
    import os

    import pandas as pd
    from sklearn.linear_model import LinearRegression

    cwd = os.getcwd()
    data = pd.read_csv("expt_molpol.csv", sep=",", skiprows=1)

    types = "element"
    # types = "sagelj"

    ff = ForceField(f"../../forcefields/{types}.pols.offxml", load_plugins=True)

    pol_handler = "MPIDPolarizability"

    data["calc."] = data["smiles"].apply(
        lambda x: molecular_polarizability(
            smiles=x, forcefield=ff, pol_handler=pol_handler
        )
    )

    plot_data = data.dropna()

    x = plot_data["expt."].values.reshape(-1, 1)
    y = plot_data["calc."].values.reshape(-1, 1)

    reg = LinearRegression()
    reg.fit(x, y)
    r_sq = reg.score(x, y)
    rrms = np.sqrt(1 / len(x) * np.sum((y - x) ** 2) / np.sum(x**2))

    linear_regression_fit = f"y = {reg.coef_[0][0]:.2f}x + {reg.intercept_[0]:.2f},\n$r^2$ = {r_sq:.2f}, RRMS = {rrms:.2%}"

    print(f"Stats for {types} polarizabilities:")
    print(f"Linear regression fit: {linear_regression_fit}")
    print(f"Saving calculated data to expt_vs_calc_pol_{types}.csv")

    plot_data.to_csv(f"expt_vs_calc_pol_{types}.csv", index=True)

    fig = plt.figure(figsize=(5, 5))
    sns.scatterplot(data=plot_data, x="expt.", y="calc.", c="royalblue")
    plt.plot([0, 25], [0, 25], label="y=x", c="k")
    x0 = np.linspace(0, 25, 100)
    y0 = x0 * reg.coef_[0][0] + reg.intercept_[0]
    plt.plot(x0, y0, c="b", ls="--", label=linear_regression_fit)
    plt.xlabel("Experimental ($\AA^3$)")
    plt.ylabel("Calculated ($\AA^3$)")
    plt.title(f"Calculated vs. experimental \nmolecular polarizability ({types})")
    plt.tight_layout()
    plt.legend(frameon=False)
    plt.savefig(f"expt_vs_calc_pol_{types}.png", dpi=300)

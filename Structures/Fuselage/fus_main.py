#%%
""" Main file for analizing the fuselage structure"""
# Importing
from fus_idealization import *
from fus_loads import *
from fus_stress_calc import *
import pandas as pd
from matplotlib import pyplot as plt

#%%
# Plotting loads
def plot_shear(data: pd.DataFrame):
    """Plot the stresses in each of the booms

    Args:
        data (pd.DataFrame): _description_
    """
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    max_stress = data["shearstress"].max()
    min_stress = data["shearstress"].min()
    q = ax.scatter(
        data["xcoor"],
        data["ycoor"],
        data["zcoor"],
        c=data["shearstress"],
        cmap="jet",
        s=7,
    )
    ax.set_ylim(-810, 810)
    ax.set_xlim(3500, 5000)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z [mm]")
    fig.colorbar(q, label="Shearstress [MPa]")
    ax.set_facecolor("grey")
    plt.show()


def plot_moment(data):
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    q = ax.scatter(
        data["xcoor"],
        data["ycoor"],
        data["zcoor"],
        c=data["bendingstress"],
        cmap="hsv",
        s=7,
    )
    ax.set_ylim(-810, 810)
    ax.set_xlim(3000, 5000)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z [mm]")
    fig.colorbar(q, label="Bendingstress [MPa]")
    ax.set_facecolor("grey")
    plt.show()


def critical_values(data, shear_threshold, cripling_comp, cripling_tens):
    shear_crit = data.query("abs(shearstress) >@shear_threshold")

    comp_crit = data.query("bendingstress< -@cripling_comp")
    tens_crit = data.query("bendingstress > @cripling_tens")

    return shear_crit, tens_crit, comp_crit


#%%
if __name__ == "__main__":
    fuselage = Fuselage_idealized(
        spacing_3d=100,
        stringer_spacing=500,
        stringer_thickness=1.2,
        stringer_width=15,
        stringer_height=25,
        skin_thickness=0.96,
        density=1480,
        start_x=3696.0,
        end_x=4940,
    )
    fuselage.create_fuselage()
    results = Fuselage_apply_loads.apply_forces(fuselage, 12)
    shear, comp, tens = critical_values(
        results,
        fuselage.shear,
        fuselage.stringerCripplingComp,
        fuselage.stringerCripplingTens,
    )

    print(shear)
    print(f"Crititcal crippling compression = {fuselage.stringerCripplingComp} [MPa]")
    print(f"Crititcal crippling tension = {fuselage.stringerCripplingTens} [MPa]")
    print(f"Critical shear stress = {fuselage.shear} [MPa]")

    plot_moment(results)

# %%

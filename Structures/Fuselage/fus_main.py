#%%
""" Main file for analizing the fuselage structure
Niels Prins"""
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
        cmap="jet",
        s=7,
    )
    ax.set_ylim(-810, 810)
    ax.set_xlim(3000, 6000)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z [mm]")
    fig.colorbar(q, label="Bendingstress [MPa]")
    ax.set_facecolor("grey")
    plt.show()


def critical_values(data, shear_threshold, cripling_comp, cripling_tens):
    shear_crit = data.query("abs(shearstress) >@shear_threshold").drop(
        columns=["area", "totalinertia", "inertia", "neutralY"]
    )
    comp_crit = data.query("bendingstress > @cripling_comp").drop(
        columns=["area", "totalinertia", "inertia", "neutralY"]
    )
    tens_crit = data.query("bendingstress < -@cripling_tens").drop(
        columns=["area", "totalinertia", "inertia", "neutralY"]
    )

    return shear_crit, tens_crit, comp_crit


def plot_margin(data, max_compres, max_crippling, max_shear):
    max_bending = data[data["zcoor"] == 1000]
    min_bending = data[data["zcoor"] == -620]
    shear = data[data["ycoor"] == 400]
    print(shear)

    plt.plot(
        max_bending["xcoor"],
        [max_compres / x for x in max_bending["bendingstress"]],
        label="Bending Compression",
    )
    plt.plot(
        max_bending["xcoor"],
        [-max_crippling / x for x in min_bending["bendingstress"]],
        label="Bending Tensile",
    )
    plt.plot(
        shear["xcoor"], [max_shear / x for x in shear["shearstress"]], label="Shear"
    )
    plt.grid(visible=True, which="major", color="#666666", linestyle="-")
    plt.minorticks_on()
    plt.grid(visible=True, which="minor", color="#999999", linestyle="-", alpha=0.2)

    plt.xlabel("X [mm]")
    plt.ylabel("Stress ratio [-]")
    plt.legend()
    plt.title("Safety margins of structure")

    plt.show()


#%%
if __name__ == "__main__":
    # Manual y coordinate placement
    bot = [0.0, 250, 400]  # From LEFT to Right in y
    side = [250]  # From RIGHT to LEFT in y
    side2 = [750]  ## This one are the z coordinates, Bottom to TOP!
    top = [250, 0.0]  # From right to left in y

    fuselage1 = Fuselage_idealized(
        spacing_3d=20,
        stringer_spacing=500,
        flange_thickness=0.96,
        web_thickness=0.96,
        stringer_width=10,
        stringer_height=15,
        skin_thickness=0.96,
        density=1480,
        start_x=3696,
        end_x=400,
        top=top,
        bot=bot,
        side=side,
        side2=side2,
        manual_place=True,
    )

    fuselage2 = Fuselage_idealized(
        spacing_3d=20,
        stringer_spacing=500,
        flange_thickness=0.96,
        web_thickness=0.96,
        stringer_width=10,
        stringer_height=15,
        skin_thickness=0.96,
        density=1480,
        start_x=3696,
        end_x=400,
        top=top,
        bot=bot,
        side=side,
        side2=side2,
        manual_place=True,
    )
    fuselage1.create_fuselage()
    results1 = Fuselage_apply_loads.apply_forces(fuselage1, 8)

    shear1, comp1, tens1 = critical_values(
        results1,
        fuselage1.shear,
        fuselage1.stringerCripplingComp,
        fuselage1.stringerCripplingTens,
    )

    plot_margin(
        results,
        fuselage1.stringerCripplingComp,
        fuselage1.stringerCripplingTens,
        fuselage1.shear,
    )

    print(f"- Critical Compression points: {comp1.head(100)}\n")
    print(f"- Critical Tensile points: {tens1.head(100)}\n")
    print(f"- Critical shear points: {shear1.head(100)}\n")
    print(
        f"Crititcal crippling compression = {round(fuselage1.stringerCripplingComp, 2)} [MPa]"
    )
    print(
        f"Crititcal crippling tension = {round(fuselage1.stringerCripplingTens,2)} [MPa]"
    )
    print(f"Critical shear stress = {fuselage1.shear} [MPa]\n")

    print(f"The fuselage weight = {round(fuselage1.weight,2)} [kg]")

    plot_shear(results)
    plot_moment(results)


# %%

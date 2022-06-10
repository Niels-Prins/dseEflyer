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


def plot_margin(data1, objects):
    compression = []
    tensile = []
    shear = []
    x = []

    for idx, data in enumerate(data1):
        max_bending1 = data[data["zcoor"] == 1000]
        min_bending1 = data[data["zcoor"] == -620]
        shear1 = data[data["ycoor"] == 400.0]
        print(shear1)

        compression.extend(
            objects[idx].stringerCripplingComp / x
            for x in max_bending1["bendingstress"]
        )
        tensile.extend(
            -objects[idx].stringerCripplingTens / x
            for x in min_bending1["bendingstress"]
        )
        
        shear.extend(objects[idx].shear / x for x in shear1["shearstress"])
        x.extend(max_bending1["xcoor"].tolist())
    print(shear)
    plt.plot(
        x,
        compression,
        label="Bending Compression",
        color="deepskyblue",
    )
    plt.plot(
        x,
        tensile,
        label="Bending Tensile",
        color="firebrick",
    )
    plt.plot(
        x,
        shear,
        label="Shear",
        color="forestgreen",
    )

    plt.grid(visible=True, which="major", color="#666666", linestyle="-")
    plt.minorticks_on()
    plt.grid(visible=True, which="minor", color="#999999", linestyle="-", alpha=0.2)

    plt.xlabel("X [mm]")
    plt.ylabel("Allowable stress / Actual stress [-]")
    plt.ylim(1, 4.2)
    plt.hlines(y =1.5, xmin=x[0], xmax=x[-1], ls="--", label="Limit", color="darkorange")
    plt.legend(loc="upper right")
    plt.title("Safety margins of structure")

    plt.show()
    
def design_fuse():
    # Manual y coordinate placement
    bot = [0.0, 250, 400]  # From LEFT to Right in y
    side = [250]  # From RIGHT to LEFT in y
    side2 = [750]  ## This one are the z coordinates, Bottom to TOP!
    top = [250, 0.0]  # From right to left in y
    
    # for thickness in range

    fuselage = Fuselage_idealized(
        flange_thickness=0.96,
        web_thickness=0.64,
        stringer_width=10,
        stringer_height=15,
        skin_thickness=0.80,
        start_x= 4200,
        end_x=4450,
        top=top,
        bot=bot,
        side=side,
        side2=side2,
        manual_place=True,
    )
    
    fuselage.create_fuselage()
    results = Fuselage_apply_loads.apply_forces(fuselage, 8)
    shear, comp, tens = critical_values(
    results,
    fuselage.shear,
    fuselage.stringerCripplingComp,
    fuselage.stringerCripplingTens,
)

    print(f"- Critical Compression points: {comp.head(100)}\n")
    print(f"- Critical Tensile points: {tens.head(100)}\n")
    print(f"- Critical shear points: {shear.head(100)}\n")
    print(
        f"Crititcal crippling compression = {round(fuselage.stringerCripplingComp, 2)} [MPa]"
    )
    print(
        f"Crititcal crippling tension = {round(fuselage.stringerCripplingTens,2)} [MPa]"
    )
    print(f"Critical shear stress = {fuselage.shear} [MPa]\n")

    print(f"The fuselage weight = {round(fuselage.weight,2)} [kg]")

    # plot_shear(results)
    # plot_moment(results)
    
    plot_margin([results], [fuselage])
    
    return results
    
if __name__ == "__main__":
    results = design_fuse()

#%%
def final_design():
    # Manual y coordinate placement
    bot = [0.0, 250, 400]  # From LEFT to Right in y
    side = [250]  # From RIGHT to LEFT in y
    side2 = [750]  ## This one are the z coordinates, Bottom to TOP!
    top = [250, 0.0]  # From right to left in y

    fuselage1 = Fuselage_idealized(
        flange_thickness=0.96,
        web_thickness=0.96,
        stringer_width=10,
        stringer_height=15,
        skin_thickness=0.96,
        start_x=3696,
        end_x=4200,
        top=top,
        bot=bot,
        side=side,
        side2=side2,
        manual_place=True,
    )
    
    fuselage2 = Fuselage_idealized(
        flange_thickness=0.96,
        web_thickness=0.64,
        stringer_width=10,
        stringer_height=15,
        skin_thickness=0.80,
        start_x= 4200,
        end_x=4450,
        top=top,
        bot=bot,
        side=side,
        side2=side2,
        manual_place=True,
    )
  
    fuselage3 = Fuselage_idealized(
        flange_thickness=0.96,
        web_thickness=0.64,
        stringer_width=10,
        stringer_height=15,
        skin_thickness=0.64,
        start_x=4450,
        end_x=4800,
        top=top,
        bot=bot,
        side=side,
        side2=side2,
        manual_place=True,
    )
    
    fuselage4 = Fuselage_idealized(
        flange_thickness=0.96,
        web_thickness=0.64,
        stringer_width=10,
        stringer_height=15,
        skin_thickness=0.48,
        start_x=4800,
        end_x=4996,
        top=top,
        bot=bot,
        side=side,
        side2=side2,
        manual_place=True,
    )
    fuselage1.create_fuselage()
    fuselage2.create_fuselage()
    fuselage3.create_fuselage()
    fuselage4.create_fuselage()
    
    results1 = Fuselage_apply_loads.apply_forces(fuselage1, 8)
    results2 = Fuselage_apply_loads.apply_forces(fuselage2, 8)
    results3 = Fuselage_apply_loads.apply_forces(fuselage3,8)
    results4 = Fuselage_apply_loads.apply_forces(fuselage4, 8)
    
    plot_margin([results1, results2, results3, results4], [fuselage1, fuselage2,fuselage3,fuselage4])

    print(fuselage4.weight +fuselage1.weight +fuselage2.weight+ fuselage3.weight)
#%%
final_design()
    


# %%

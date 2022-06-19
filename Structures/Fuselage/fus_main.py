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
def plot_shear(data: pd.DataFrame, multiplication=1):
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
        c=data["shearstress"] * multiplication,
        cmap="jet",
        s=4,
    )

    q = ax.scatter(
        data["xcoor"],
        -data["ycoor"],
        data["zcoor"],
        c=data["shearstress"] * multiplication,
        cmap="jet",
        s=4,
    )
    ax.set_ylim(-810, 810)
    ax.set_xlim(3000, 6000)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z [mm]")
    fig.colorbar(q, label="Shearstress [MPa]")
    plt.show()


def plot_moment(data, multiplication=1):
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    q = ax.scatter(
        data["xcoor"],
        data["ycoor"],
        data["zcoor"],
        c=data["bendingstress"] * multiplication,
        cmap="jet",
        s=4,
    )

    q = ax.scatter(
        data["xcoor"],
        -data["ycoor"],
        data["zcoor"],
        c=data["bendingstress"] * multiplication,
        cmap="jet",
        s=4,
    )

    ax.set_ylim(-810, 810)
    ax.set_xlim(3000, 6000)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z [mm]")
    fig.colorbar(q, label="Bendingstress [MPa]")
    plt.legend()
    plt.show()


def critical_values(data, shear_threshold, cripling_comp, cripling_tens):
    shear_crit = data.query("abs(shearstress) >@shear_threshold/1.5").drop(
        columns=["area", "totalinertia", "inertia", "neutralY"]
    )
    comp_crit = data.query("bendingstress > @cripling_comp/1.5").drop(
        columns=["area", "totalinertia", "inertia", "neutralY"]
    )
    tens_crit = data.query("bendingstress < -@cripling_tens/1.5").drop(
        columns=["area", "totalinertia", "inertia", "neutralY"]
    )

    return shear_crit, tens_crit, comp_crit


def plot_margin(data1, objects):
    compression = []
    tensile = []
    shear = []
    x = []
    sections = [3696]

    for idx, data in enumerate(data1):
        max_bending1 = data[data["zcoor"] == 1000]
        min_bending1 = data[data["zcoor"] == -620]
        shear1 = data[data["ycoor"] == 400.0]

        compression.extend(
            abs(objects[idx].stringerCripplingComp / x)
            for x in max_bending1["bendingstress"]
        )
        tensile.extend(
            abs(-objects[idx].stringerCripplingTens / x)
            for x in min_bending1["bendingstress"]
        )

        shear.extend(abs(objects[idx].shear / x) for x in shear1["shearstress"])
        Xcoor = max_bending1["xcoor"].tolist()
        sections.append(max(Xcoor))
        x.extend(Xcoor)

    plt.plot(
        x,
        compression,
        label="Bending\ncompression",
        color="deepskyblue",
    )
    plt.plot(
        x,
        tensile,
        label="Bending\ntensile",
        color="firebrick",
    )
    plt.plot(
        x,
        shear,
        label="Shear",
        color="forestgreen",
    )

    plt.vlines(sections, ymin=1, ymax=8, color="mediumblue", ls="--", label="Section")

    plt.grid(visible=True, which="major", color="#666666", linestyle="-")
    plt.minorticks_on()
    plt.grid(visible=True, which="minor", color="#999999", linestyle="-", alpha=0.2)

    plt.xlabel("x [mm]")
    plt.ylabel("Allowable stress / Actual stress [-]")
    plt.ylim(1, 6)
    plt.hlines(y=1.5, xmin=x[0], xmax=x[-1], ls="--", label="Limit", color="darkorange")
    plt.legend(loc="upper left")
    # plt.title("Safety margins of structure")
    plt.savefig("marginsFuselage6")

    plt.show()


def design_fuse():
    # Manual y coordinate placement
    bot = [0.0, 250, 400]  # From LEFT to Right in y
    side = [300, 250]  # From RIGHT to LEFT in y
    side2 = [750]  ## This one are the z coordinates, Bottom to TOP!
    top = [250, 0.0]  # From right to left in y

    # for thickness in range

    fuselage = Fuselage_idealized(
        stringer_spacing=40,
        flange_thickness=0.80,
        web_thickness=0.80,
        stringer_width=10,
        stringer_height=12,
        skin_thickness=0.64,
        start_x=3696,
        end_x=4080,
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

    plot_shear(results, multiplication=1.2)
    plot_moment(results, multiplication=3.2)

    # plot_margin([results], [fuselage])

    return results


# design_fuse()
#%%
def final_design():
    # Manual y coordinate placement
    bot = [0.0, 250, 400]  # From LEFT to Right in y
    side = [300, 250]  # From RIGHT to LEFT in y
    side2 = [750]  ## This one are the z coordinates, Bottom to TOP!
    top = [0.0]  # From right to left in y

    fuselage1 = Fuselage_idealized(
        flange_thickness=0.96,
        web_thickness=0.96,
        stringer_width=12,
        stringer_height=15,
        skin_thickness=0.48,
        start_x=3696,
        end_x=3815,
        top=top,
        bot=bot,
        side=side,
        side2=side2,
        manual_place=True,
    )

    fuselage2 = Fuselage_idealized(
        flange_thickness=0.96,
        web_thickness=0.96,
        stringer_width=12,
        stringer_height=15,
        skin_thickness=0.64,
        start_x=3815,
        end_x=4100,
        top=top,
        bot=bot,
        side=side,
        side2=side2,
        manual_place=True,
    )
    
    fuselage3 = Fuselage_idealized(
        flange_thickness=0.96,
        web_thickness=0.96,
        stringer_width=12,
        stringer_height=15,
        skin_thickness=0.48,
        start_x=4100,
        end_x=4400,
        top=top,
        bot=bot,
        side=side,
        side2=side2,
        manual_place=True,
    )
    
    fuselage4 = Fuselage_idealized(
        flange_thickness=0.96,
        web_thickness=0.96,
        stringer_width=12,
        stringer_height=15,
        skin_thickness=0.32,
        start_x=4400,
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
    results3 = Fuselage_apply_loads.apply_forces(fuselage3, 8)
    results4 = Fuselage_apply_loads.apply_forces(fuselage4, 8)

    test = pd.concat([ results1, results2, results3, results4]) 
    plot_moment(test)
    plot_shear(test)

    plot_margin(
        [results1,results2, results3, results4], # 
        [fuselage1, fuselage2, fuselage3, fuselage4]) # 

    # fuselage1.plot_cross_sec()
    print(fuselage1.stringerArea)

    print(fuselage1.weight, fuselage2.weight)
    print(
        f"The weight of the fuselage = {round(fuselage1.weight +fuselage2.weight +fuselage3.weight + fuselage4.weight,2)} [kg]"
    )
    

# %%
def plot_final_design():
    bot = [0.0, 250, 400]  # From LEFT to Right in y
    side = [300, 250]  # From RIGHT to LEFT in y
    side2 = [750]  ## This one are the z coordinates, Bottom to TOP!
    top = [0.0]  # From right to left in y

    fuselage1 = Fuselage_idealized(
        stringer_spacing =30,
        spacing_3d=20,
        flange_thickness=0.96,
        web_thickness=0.96,
        stringer_width=12,
        stringer_height=15,
        skin_thickness=0.48,
        start_x=3696,
        end_x=3815,
        top=top,
        bot=bot,
        side=side,
        side2=side2,
        # manual_place=True,
    )

    fuselage2 = Fuselage_idealized(
        stringer_spacing =30,
        spacing_3d=20,
        flange_thickness=0.96,
        web_thickness=0.96,
        stringer_width=12,
        stringer_height=15,
        skin_thickness=0.64,
        start_x=3815,
        end_x=4100,
        top=top,
        bot=bot,
        side=side,
        side2=side2,
        # manual_place=True,
    )
    
    fuselage3 = Fuselage_idealized(
        stringer_spacing =30,
        spacing_3d=20,
        flange_thickness=0.96,
        web_thickness=0.96,
        stringer_width=12,
        stringer_height=15,
        skin_thickness=0.48,
        start_x=4100,
        end_x=4400,
        top=top,
        bot=bot,
        side=side,
        side2=side2,
        # manual_place=True,
    )
    
    fuselage4 = Fuselage_idealized(
        stringer_spacing =30,
        spacing_3d=20,
        flange_thickness=0.96,
        web_thickness=0.96,
        stringer_width=12,
        stringer_height=15,
        skin_thickness=0.32,
        start_x=4400,
        end_x=4996,
        top=top,
        bot=bot,
        side=side,
        side2=side2,
        # manual_place=True,
    )

    fuselage1.create_fuselage()
    fuselage2.create_fuselage()
    fuselage3.create_fuselage()
    fuselage4.create_fuselage()

    results1 = Fuselage_apply_loads.apply_forces(fuselage1, 8)
    results2 = Fuselage_apply_loads.apply_forces(fuselage2, 8)
    results3 = Fuselage_apply_loads.apply_forces(fuselage3, 8)
    results4 = Fuselage_apply_loads.apply_forces(fuselage4, 8)

    test = pd.concat([ results1, results2, results3, results4]) 
    plot_moment(test, multiplication=5)
    plot_shear(test, multiplication=0.9)

    plot_margin(
        [results1, results2], [fuselage1, fuselage2]
    )  # , fuselage3, fuselage4],)

    # print(
    #     f"The weight of the fuselage = {round(fuselage4.weight +fuselage1.weight +fuselage2.weight+ fuselage3.weight,2)} [kg]"
    # )


# %%
if __name__ == "__main__":
    final_design()
    plot_final_design()

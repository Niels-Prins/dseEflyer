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
def plot_stress(data: pd.DataFrame):
    """ Plot the stresses in each of the booms

    Args:
        data (pd.DataFrame): _description_
    """    
    fig = plt.figure()
    ax = plt.axes(projection="3d")
    max_stress = data["shearstress"].max()
    min_stress = data["shearstress"].min()
    q = ax.scatter(
        data["xcoor"], data["ycoor"], data["zcoor"], c=data["shearstress"], cmap="YlOrRd", s=2.5
    )
    ax.set_ylim(-810, 810)
    ax.set_xlim(3000, 5000)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z [mm]")
    fig.colorbar(q, label="Shearstress [N/mm2]")
    plt.show()


def critical_values(data, threshold):
    highest = data.query(
        "abs(shearstress) >@threshold | abs(bendingstress) > @threshold "
    )
    return highest


#%%
if __name__ == "__main__":
    fuselage = Fuselage_idealized(
        spacing_3d=100,
        stringer_spacing=140,
        stringer_thickness=2.5,
        stringer_width=15,
        stringer_height=30,
        skin_thickness=1.92,
        density=1480,
        start_x =3696.0,
        end_x =4940 
    )
    fuselage.create_fuselage()
    load_results = Fuselage_apply_loads.apply_forces(fuselage.x, fuselage.data, 8)
    plot_stress(load_results)

    criticalshit = critical_values(load_results, 50)


# %%

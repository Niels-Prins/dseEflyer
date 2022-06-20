#%%
from audioop import avg
from re import M
from tkinter import X
from unittest.loader import VALID_MODULE_NAME

import numpy as np
from matplotlib import pyplot as plt

g0 = 9.80665


aircraftOEW = 567  # aircraft mass, [kg]
bat_mass = 276
motor_mass = 33
full_fuse = 7.6  # fuselage length

cg_oew = 4.787
wing_pos = 4.112  # wing position
bat_pos = 1.955
motor_pos = 5.14

n = 8  # SET G-LOADING HERE
OEW = aircraftOEW * g0  # aircraft weight, [N]
BAT = bat_mass * g0
MOT = motor_mass * g0
Wing = (OEW + BAT + MOT) * n
Shear_FW = []
Bending_FW = []
x_pos_FW = np.linspace(0, cg_oew, 10000)

Shear_AFT = []
Bending_AFT = []
x_pos_AFT = np.linspace(0, full_fuse - cg_oew, 10000)

Vx = []
Mx = []

distributedL_weightFW = cg_oew  # fuselage span in front of the CG
distributedL_weightAFT = full_fuse - cg_oew  # fuselage span aft of the CG


qFW = (OEW * 2) / (distributedL_weightFW ** 2)
qAFT = (OEW * 2) / (distributedL_weightAFT ** 2)


def VM(x, n: int = 12):

    distributedL_weightFW = cg_oew #fuselage span in front of the CG
    distributedL_weightAFT = full_fuse - cg_oew #fuselage span aft of the CG


    qFW = (OEW*2)/(distributedL_weightFW**2)*n
    qAFT = (OEW*2)/(distributedL_weightAFT**2)*n

    if 0 <= x <= bat_pos:
        v = ((qFW/(2))*(x)**2)
        m = ((qFW/(6))*(x)**3)
        
    if bat_pos <= x <= cg_oew:
        v =  BAT + ((qFW/(2))*(x - bat_pos)**2)
        m =  BAT * (x - bat_pos) + ((qFW/(6))*(x - bat_pos)**3)
        
    if x > 3.81:
        x = x -3.81
        if 0 <= x <= (distributedL_weightAFT-(full_fuse-motor_pos)):
            v = ((qAFT/(2))*(distributedL_weightAFT - x)**2) + MOT - Wing
            m = ((qAFT/(6))*(distributedL_weightAFT - x)**3) + MOT*(distributedL_weightAFT-(full_fuse-motor_pos) - x) - Wing*(distributedL_weightAFT-(full_fuse-wing_pos) - x)

        if distributedL_weightAFT-(full_fuse-motor_pos) <= x <= distributedL_weightAFT:
            v = ((qAFT/(2))*(distributedL_weightAFT - x)**2)
            m = ((qAFT/(6))*(distributedL_weightAFT - x)**3)

    return v , m 

if __name__ == "__main__":
    Vx = []
    Mx = []

    distributedL_weightFW = cg_oew #fuselage span in front of the CG
    distributedL_weightAFT = full_fuse - cg_oew #fuselage span aft of the CG


    qFW = (OEW*2)/(distributedL_weightFW**2)*n
    qAFT = (OEW*2)/(distributedL_weightAFT**2)*n

    for i in range(len(x_pos_FW)):
        if 0 <= x_pos_FW[i] <= bat_pos:
            v1 = ((qFW/(2))*(x_pos_FW[i])**2)
            m1 = ((qFW/(6))*(x_pos_FW[i])**3)
            Shear_FW.append(v1/1000)
            Bending_FW.append(m1/1000)
        if bat_pos <= x_pos_FW[i] <= cg_oew:
            v2 =  BAT + ((qFW/(2))*(x_pos_FW[i] - bat_pos)**2)
            m2 =  BAT * (x_pos_FW[i] - bat_pos) + ((qFW/(6))*(x_pos_FW[i] - bat_pos)**3)
            Shear_FW.append(v2/1000)
            Bending_FW.append(m2/1000)
            
    for k in range(len(x_pos_AFT)):
        if 0 <= x_pos_AFT[k] <= (distributedL_weightAFT-(full_fuse-motor_pos)):
            v3 = ((qAFT/(2))*(distributedL_weightAFT - x_pos_AFT[k])**2) + MOT - Wing
            m3 = ((qAFT/(6))*(distributedL_weightAFT - x_pos_AFT[k])**3) + MOT*(distributedL_weightAFT-(full_fuse-motor_pos) - x_pos_AFT[k]) - Wing*(distributedL_weightAFT-(full_fuse-wing_pos) - x_pos_AFT[k])
            Shear_AFT.append(v3/1000)
            Bending_AFT.append(m3/1000)
        if distributedL_weightAFT-(full_fuse-motor_pos) <= x_pos_AFT[k] <= distributedL_weightAFT:
            v4 = ((qAFT/(2))*(distributedL_weightAFT - x_pos_AFT[k])**2)
            m4 = ((qAFT/(6))*(distributedL_weightAFT - x_pos_AFT[k])**3)
            Shear_AFT.append(v4/1000)
            Bending_AFT.append(m4/1000)

    x_pos_whole = np.linspace(0, full_fuse, (len(x_pos_FW) + len(x_pos_AFT)))
    Vx = [*Shear_FW, *Shear_AFT]
    Mx = [*Bending_FW, *Bending_AFT]



def plot_load_distribution():
    print(
        "Triangular distributed load on the front half of the fuselage:", qFW, "[N/m^2]"
    )
    print(
        "Triangular distributed load on the aft half of the fuselage:", qAFT, "[N/m^2]"
    )

    fig, axis = plt.subplots(nrows=2, ncols=1)
    axis[0].plot((x_pos_whole), (Vx), color="red")
    axis[0].set_ylabel("Shear Force [kN]")
    axis[0].set_ylim(min(Vx), max(Vx))
    axis[0].minorticks_on()
    axis[0].grid(True, which="minor", color="#999999", linestyle="-", alpha=0.2)
    axis[0].grid(True, which="major", color="#666666", linestyle="-")
    axis[0].vlines(3.696, ymin=min(Mx), ymax=max(Mx), color="g", ls="--")
    axis[0].vlines(4.996, ymin=min(Mx), ymax=max(Mx), color="g", ls="--")

    axis[1].plot((x_pos_whole), (Mx), color="blue")
    axis[1].set_xlabel("fuselage [mm]")
    axis[1].set_ylabel("Bending moment [kNm]")
    axis[1].set_ylim(min(Mx), max(Mx))
    axis[1].minorticks_on()
    axis[1].grid(True, which="minor", color="#999999", linestyle="-", alpha=0.2)
    axis[1].grid(True, which="major", color="#666666", linestyle="-")
    axis[1].vlines(3.696, ymin=min(Mx), ymax=max(Mx), color="g", ls="--")
    axis[1].vlines(4.996, ymin=min(Mx), ymax=max(Mx), color="g", ls="--")

    plt.show()

    SF = 1.5
    print(
        f"Maximum shear force including a safety factor of {SF} is {max(Vx)/1000 * SF}[kN]"
    )
    print(
        f"Maximum bending moment including a safety factor of {SF} is {max(Mx)/1000 * SF}[kNm]"
    )

# plot_load_distribution()
# %%

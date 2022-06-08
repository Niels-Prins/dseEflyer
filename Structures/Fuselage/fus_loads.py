#%%
from audioop import avg
from re import M
from tkinter import X
import numpy as np
from matplotlib import pyplot as plt


g0 = 9.81

aircraftMass = 975  # aircraft mass, [kg]
full_fuse = 7.6  # fuselage length
fw_CG = 3.8  # forward CG position
aft_CG = 4.0  # aft CG position
wing_pos = 4.112  # wing position
tailH_pos = 6.63  # horizontal tail position
avg_CG = (fw_CG + aft_CG) / 2  # average CG position

n = 1  # G-Force
W = aircraftMass * g0  # aircraft weight, [N]

Shear_FW = []
Bending_FW = []
x_pos_FW = np.linspace(0, avg_CG, 100)

Shear_AFT = []
Bending_AFT = []
x_pos_AFT = np.linspace(0, full_fuse - avg_CG, 100)

Vx = []
Mx = []

distributedL_weightFW = avg_CG  # fuselage span in front of the CG
distributedL_weightAFT = full_fuse - avg_CG  # fuselage span aft of the CG

qFW = ((W * n) / (distributedL_weightFW)) / 2
qAFT = ((W * n) / (distributedL_weightAFT)) / 2

distributedL_weightFW = avg_CG  # fuselage span in front of the CG
distributedL_weightAFT = full_fuse - avg_CG  # fuselage span aft of the CG

qFW = ((W * n) / (distributedL_weightFW)) / 2
qAFT = ((W * n) / (distributedL_weightAFT)) / 2

def VM(x, n: int=8):
    distributedL_weightFW = avg_CG  # fuselage span in front of the CG
    distributedL_weightAFT = full_fuse - avg_CG  # fuselage span aft of the CG
    
    qFW = ((W * n) / (distributedL_weightFW)) / 2
    qAFT = ((W * n) / (distributedL_weightAFT)) / 2
    

    if x <= avg_CG:
        v = (qFW / 2) * (x) ** 2
        m = (qFW / 6) * (x) ** 3
    if x >= avg_CG:
        v = (qAFT / 2) * (full_fuse -x) ** 2
        m = (qAFT / 6) * (full_fuse -x) ** 3

    return v, m


for i in range(len(x_pos_FW)):
    v1 = (qFW / 2) * (x_pos_FW[i]) ** 2
    m1 = (qFW / 6) * (x_pos_FW[i] ** 3)
    Shear_FW.append(v1)
    Bending_FW.append(m1)

for j in range(len(x_pos_AFT)):
    v2 = (qAFT / 2) * (distributedL_weightAFT - x_pos_AFT[j]) ** 2
    m2 = (qAFT / 6) * (distributedL_weightAFT - x_pos_AFT[j]) ** 3
    Shear_AFT.append(v2)
    Bending_AFT.append(m2)
    
x_pos_whole = np.linspace(0, full_fuse, (len(x_pos_FW) + len(x_pos_AFT)))
Vx = [*Shear_FW, *Shear_AFT]
Mx = [*Bending_FW, *Bending_AFT]

if __name__ == "__main__":
    print(
        "Triangular distributed load on the front half of the fuselage:", qFW, "[N/m^2]"
    )
    print(
        "Triangular distributed load on the aft half of the fuselage:", qAFT, "[N/m^2]"
    )

    fig, axis = plt.subplots(nrows=2, ncols=1)
    axis[0].plot(x_pos_whole, Vx, color="red")

    axis[1].plot(x_pos_whole, Mx, color="blue")
    plt.show()

    loc = 3.5
    print(f"Shear force at {loc}[m] is {VM(loc)[0]}[N]")
    print(f"Bending moment at {loc}[m] is {VM(loc)[1]}[Nm]")
#%%
# fig,axis=plt.subplots(nrows = 2,ncols = 2)
# axis[0,0].plot(x_pos_FW, Shear_FW)
# axis[1,0].plot(x_pos_FW, Bending_FW)
# axis[0,0].plot(x_pos_FW, np.zeros(len(x_pos_FW)), color='black')
# axis[0,0].set_ylim(-4900, 5000)
# axis[0,0].grid(True)
# axis[1,0].plot(x_pos_FW, np.zeros(len(x_pos_FW)), color='black')
# axis[1,0].grid(True)
# axis[1,0].set_ylim(-100, 7000)


# axis[0,1].plot(x_pos_AFT, Shear_AFT)
# axis[1,1].plot(x_pos_AFT, Bending_AFT)
# axis[0,1].plot(x_pos_AFT, np.zeros(len(x_pos_AFT)), color='black')
# axis[0,1].grid(True)
# axis[0,1].set_ylim(-4900, 5000)
# axis[1,1].plot(x_pos_AFT, np.zeros(len(x_pos_AFT)), color='black')
# axis[1,1].grid(True)
# axis[1,1].set_ylim(-100, 7000)

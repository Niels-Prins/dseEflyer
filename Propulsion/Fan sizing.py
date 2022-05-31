import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import quad

""""###############################################################################################################"""
""""Variables"""
""""###############################################################################################################"""
# constants
kt = 1.032914
kq = 0.1110977401

# velocities
V_cl = 77.17  # [m/s]
V_to = 28.29
V_la = 33.438

# Power values
P_required = 229000  # [W]
P_a = 260000

# aircraft characteristics
AR = 5.8
S = 12.3
e = 0.8
g = 9.81
MTOM = 979  # kg
MTOM_3g = 3 * MTOM
Cd0 = 0.025

# take-off
Cd_to = 0.202
Cl_to = 1.59
mu = 0.015
T = 3200
W = g * MTOM

# altitudes
alt_0 = 0
alt_6000 = 1828

# desired D for calculations
D = 0.63

# efficiencies
eff_engine = 0.95
eff_fan = 0.885
""""###############################################################################################################"""
""""Definitions"""
""""###############################################################################################################"""

#Calculator for the ISA conditions up to 11000 meters
def ISA(h):
    T0 = 288.15
    p0 = 101325.0
    g0 = 9.80665
    R = 287.0
    a1 = -0.0065
    T = T0 + a1 * (h - 0)
    p = p0 * (T / T0) ** (-g0 / (R * a1))
    rho1 = (p / (R * T))
    return [T, p, rho1]

#Mssflow formula
def Massflow(thr, V4, V):
    return thr / (V4 - V)

#Thrust formula (kt is based on the data available of the UL-39 ducted fan aircraft)
def Thrust(kt, n, D, alt):
    return kt * n ** 2 * D ** 4 * ISA(alt)[2]

#Torque formula (kq is based on the data available of the UL-39 ducted fan aircraft)
def Torque(kq, n, D, alt):
    return kq * n ** 2 * D ** 5 * ISA(alt)[2]

#formula for change in speed due to the spinning propeller
def V_induced(e_d, V, D, Thr, alt):
    return (0.5 * e_d - 1) * V + np.sqrt((e_d * V / 2) ** 2 + (e_d * Thr) / (ISA(alt)[2] * np.pi * D ** 2 / 4))

#formula for the torque power based on the induced airspeed
def P_torque(w, Q):
    return w * Q

#formula for calculating the power required for certain rpm, thrust and diameter, velocity
def Preq(alt, Thr, D, V, e_d):
    return 0.75 * Thr * V + np.sqrt(((Thr ** 2 * V ** 2) / 16) + ((Thr ** 3) / (ISA(alt)[2] * np.pi * e_d * D ** 2)))

#Velocity just after the propeller including the induced speed
def V3(V, w):
    return V + w

#Exit speed at the end of the duct
def V4(V, e_d):
    return V / e_d

#Definition for drag as a function of time
def Drag(V, MTOM, alt):
    Cl = (MTOM * g) / (0.5 * ISA(alt)[2] * V ** 2 * S)
    Cd = Cd0 + Cl ** 2 / (np.pi * AR * e)
    return 0.5 * ISA(alt)[2] * V ** 2 * S * Cd

#formula for Preq as a function of drag
def Preq_drag(V, MTOM, alt, D, e_d):
    return 0.75 * Drag(V, MTOM, alt) * V + np.sqrt(
        ((Drag(V, MTOM, alt) ** 2 * V ** 2) / 16) + ((Drag(V, MTOM, alt) ** 3) / (ISA(alt)[2] * np.pi * e_d * D ** 2)))

# function we want to integrate
def f(X):
    return 1 / ((g / W * T) - (g / W * mu * W) - (g / W * 0.5 * 1.225 * X * S * (Cd_to - mu * Cl_to)))


""""###############################################################################################################"""
""""Diameter and RPM calculation"""
""""###############################################################################################################"""
# Constants for loop specific
lst = []
thr_thres_cl_0 = 2001
thr_thres_cl_6000 = 1967
thr_thres_to_0 = 1344
thr_thres_to_6000 = 1321
thr_thres_la_0 = 1593
thr_thres_la_6000 = 1566
rpmlower = 4000
rpmupper = 8000
rpmspacing = 100
marginlower = 1
marginupper = 1.05


# for loop
for D_cl in np.arange(0.4, 0.75, 0.01):
    for e_d in np.arange(0.8, 1.21, 0.01):
        for n_cl in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
            Thr_cl_0 = Thrust(kt, n_cl, D_cl, alt_0)
            Torque_cl_0 = Torque(kq, n_cl, D_cl, alt_0)
            Pr_cl_0 = Preq(alt_0, Thr_cl_0, D_cl, V_cl, e_d)
            w_cl_0 = V_induced(e_d, V_cl, D_cl, Thr_cl_0, alt_0)
            Torque_power_cl_0 = P_torque(w_cl_0, Torque_cl_0)
            V_3_cl_0 = V3(V_cl, w_cl_0)
            V_4_cl_0 = V4(V_3_cl_0, e_d)
            mflow_cl_0 = Massflow(Thr_cl_0, V_4_cl_0, V_cl)
            if marginlower * thr_thres_cl_0 <= Thr_cl_0 <= marginupper * thr_thres_cl_0:
                for n_cl_2 in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                    Thr_cl_6000 = Thrust(kt, n_cl_2, D_cl, alt_6000)
                    Torque_cl_6000 = Torque(kq, n_cl_2, D_cl, alt_6000)
                    Pr_cl_6000 = Preq(alt_6000, Thr_cl_6000, D_cl, V_cl, e_d)
                    w_cl_6000 = V_induced(e_d, V_cl, D_cl, Thr_cl_6000, alt_6000)
                    Torque_power_cl_6000 = P_torque(w_cl_6000, Torque_cl_6000)
                    V_3_cl_6000 = V3(V_cl, w_cl_6000)
                    V_4_cl_6000 = V4(V_3_cl_6000, e_d)
                    mflow_cl_6000 = Massflow(Thr_cl_6000, V_4_cl_6000, V_cl)
                    if marginlower * thr_thres_cl_6000 <= Thr_cl_6000 <= marginupper * thr_thres_cl_6000:
                        for n_to in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                            Thr_to_0 = Thrust(kt, n_to, D_cl, alt_0)
                            Torque_to_0 = Torque(kq, n_to, D_cl, alt_0)
                            Pr_to_0 = Preq(alt_0, Thr_to_0, D_cl, V_to, e_d)
                            w_to_0 = V_induced(e_d, V_to, D_cl, Thr_to_0, alt_0)
                            Torque_power_to_0 = P_torque(w_to_0, Torque_to_0)
                            V_3_to_0 = V3(V_to, w_to_0)
                            V_4_to_0 = V4(V_3_to_0, e_d)
                            mflow_to_0 = Massflow(Thr_to_0, V_4_to_0, V_to)
                            if marginlower * thr_thres_to_0 <= Thr_to_0 <= marginupper * thr_thres_to_0:
                                for n_to_2 in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                                    Thr_to_6000 = Thrust(kt, n_to_2, D_cl, alt_6000)
                                    Torque_to_6000 = Torque(kq, n_to_2, D_cl, alt_6000)
                                    Pr_to_6000 = Preq(alt_6000, Thr_to_6000, D_cl, V_cl, e_d)
                                    w_to_6000 = V_induced(e_d, V_to, D_cl, Thr_to_6000, alt_6000)
                                    Torque_power_to_6000 = P_torque(w_to_6000, Torque_to_6000)
                                    V_3_to_6000 = V3(V_to, w_to_6000)
                                    V_4_to_6000 = V4(V_3_to_6000, e_d)
                                    mflow_to_6000 = Massflow(Thr_to_6000, V_4_to_6000, V_to)
                                    if marginlower * thr_thres_to_6000 <= Thr_to_6000 <= marginupper * thr_thres_to_6000:
                                        for n_la in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                                            Thr_la_0 = Thrust(kt, n_la, D_cl, alt_0)
                                            Torque_la_0 = Torque(kq, n_la, D_cl, alt_0)
                                            Pr_la_0 = Preq(alt_0, Thr_la_0, D_cl, V_cl, e_d)
                                            w_la_0 = V_induced(e_d, V_la, D_cl, Thr_to_0, alt_0)
                                            Torque_power_la_0 = P_torque(w_la_0, Torque_la_0)
                                            V_3_la_0 = V3(V_la, w_la_0)
                                            V_4_la_0 = V4(V_3_la_0, e_d)
                                            mflow_la_0 = Massflow(Thr_la_0, V_4_la_0, V_la)
                                            if marginlower * thr_thres_la_0 <= Thr_la_0 <= marginupper * thr_thres_la_0:
                                                for n_la_2 in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                                                    Thr_la_6000 = Thrust(kt, n_la_2, D_cl, alt_6000)
                                                    Torque_la_6000 = Torque(kq, n_la_2, D_cl, alt_6000)
                                                    Pr_la_6000 = Preq(alt_6000, Thr_la_6000, D_cl, V_la, e_d)
                                                    w_la_6000 = V_induced(e_d, V_la, D_cl, Thr_to_6000, alt_6000)
                                                    Torque_power_la_6000 = P_torque(w_la_6000, Torque_la_6000)
                                                    V_3_la_6000 = V3(V_la, w_la_6000)
                                                    V_4_la_6000 = V4(V_3_la_6000, e_d)
                                                    mflow_la_6000 = Massflow(Thr_la_6000, V_4_la_6000, V_la)
                                                    if marginlower * thr_thres_la_6000 <= Thr_la_6000 <= marginupper * thr_thres_la_6000:
                                                        lst.append([D_cl, D_cl ** 2 / 4 * np.pi, e_d,
                                                                    max(n_cl * 60, n_cl_2 * 60, n_to * 60,
                                                                        n_to_2 * 60, n_la * 60, n_la_2 * 60),
                                                                    min(n_cl * 60, n_cl_2 * 60, n_to * 60,
                                                                        n_to_2 * 60, n_la * 60, n_la_2 * 60),
                                                                    max(Pr_cl_0, Pr_cl_6000, Pr_to_0, Pr_to_6000,
                                                                        Pr_la_0,
                                                                        Pr_la_6000),
                                                                    max(Pr_cl_0, Pr_cl_6000, Pr_to_0, Pr_to_6000,
                                                                        Pr_la_0, Pr_la_6000) / eff_fan / eff_engine,
                                                                    n_cl * 60,
                                                                    n_cl_2 * 60, n_to * 60,
                                                                    n_to_2 * 60, n_la * 60, n_la_2 * 60,
                                                                    Pr_cl_0, Pr_cl_6000, Pr_to_0, Pr_to_6000, Pr_la_0,
                                                                    Pr_la_6000,
                                                                    Thr_cl_0, Thr_cl_6000, Thr_to_0, Thr_to_6000,
                                                                    Thr_la_0,
                                                                    Thr_la_6000,
                                                                    Torque_cl_0, Torque_cl_6000, Torque_to_0,
                                                                    Torque_to_6000,
                                                                    Torque_la_0, Torque_la_6000
                                                                       , mflow_cl_0, mflow_cl_6000, mflow_to_0,
                                                                    mflow_to_6000,
                                                                    mflow_la_0, mflow_la_6000, w_cl_0, w_cl_6000
                                                                       , w_to_0, w_to_6000, w_la_0, w_la_6000,
                                                                    V_3_cl_0, V_3_cl_6000, V_3_to_0, V_3_to_6000,
                                                                    V_3_la_0, V_3_la_6000,
                                                                    V_4_cl_0, V_4_cl_6000, V_4_to_0, V_4_to_6000,
                                                                    V_4_la_0, V_4_la_6000,
                                                                    Torque_power_cl_0, Torque_power_cl_6000,
                                                                    Torque_power_to_0,
                                                                    Torque_power_to_6000, Torque_power_la_0,
                                                                    Torque_power_la_6000])

lst_header = ["Diameter", "Area", "e_d", "Max RPM", 'Min RPM', "Max P prop", "Max P engine",
              "RPM_cl_0", "RPM_cl_6000", "RPM_to_0", "RPM_to_6000", "RPM_la_0", "RPM_la_6000",
              "P_cl_0", "P_cl_6000", "P_to_0", "P_to_6000", "P_la_0", "P_la_6000",
              "Thrust_cl_0", "Thrust_cl_6000", "Thrust_to_0", "Thrust_to_6000",
              "Thrust_la_0", "Thrust_la_6000", "Torque_cl_0", "Torque_cl_6000",
              "Torque_to_0", "Torque_to_6000", "Torque_la_0", "Torque_la_6000",
              "Mass_flow_cl_0", "Mass_flow_cl_6000", "Mass_flow_to_0",
              "Mass_flow_to_6000", "Mass_flow_la_0", "Mass_flow_la_6000",
              "Induced_airspeed_cl_0", "Induced_airspeed_cl_6000", "Induced_airspeed_to_0",
              "Induced_airspeed_to_6000", "Induced_airspeed_la_0", "Induced_airspeed_la_0",
              "Rotor_speed_V3_cl_0", "Rotor_speed_V3_cl_6000", "Rotor_speed_V3_to_0",
              "Rotor_speed_V3_to_6000", "Rotor_speed_V3_la_0", "Rotor_speed_V3_la_6000",
              "Exit_speed_V4_cl_0", "Exit_speed_V4_cl_6000", "Exit_speed_V4_to_0",
              "Exit_speed_V4_to_6000", "Exit_speed_V4_la_0", "Exit_speed_V4_la_6000",
              "Torque_power_cl_0", "Torque_power_cl_6000", "Torque_power_to_0", "Torque_power_to_6000",
              "Torque_power_la_0", "Torque_power_la_0"]
df = pd.DataFrame(lst, columns=lst_header)
df.to_excel("Prop_sizing_final.xlsx", index=False)

new_df = df[(round(df["Diameter"], 2) == round(D, 2)) & (round(df["e_d"], 2) == round(1, 2))]
print(new_df.iloc[0])

id = 0
df_choice = pd.DataFrame(np.array([new_df.iloc[id]]).transpose(), columns=["Propeller"], index=lst_header)
df_choice.to_excel("Propeller_choice.xlsx", index=True)

""""###############################################################################################################"""
""""Power curves"""
""""###############################################################################################################"""
alt_power = 0
lst_Preq = []
for V in np.arange(20, 120, 0.01):
    lst_Preq.append([Preq_drag(V, MTOM, alt_power, new_df.iloc[id]["Diameter"], new_df.iloc[id]["e_d"]), V])

# plot the curve for normal flight
df_curve = pd.DataFrame(lst_Preq, columns=["Power required", "Velocity"])
plt.plot(df_curve["Velocity"], df_curve["Power required"], label="Power required")
plt.plot(df_curve["Velocity"], P_a * np.ones(len(lst_Preq)), label="Power available")
plt.vlines(x=V_cl, ymin=0, ymax=df_curve[round(df_curve["Velocity"], 2) == round(V_cl, 2)]["Power required"].values[0],
           linestyle='dotted', color='green')
plt.xlabel("Velocity")
plt.ylabel("Power required")
plt.legend()
plt.title("Power curve")
#plt.show()

alt_power = 0
lst_Preq_g = []
g_load = 3
for V in np.arange(60, 90, 0.01):
    lst_Preq_g.append([Preq_drag(V, g_load * MTOM, alt_power, new_df.iloc[id]["Diameter"], new_df.iloc[id]["e_d"]), V])

df_curve_g = pd.DataFrame(lst_Preq_g, columns=["Power required", "Velocity"])
plt.plot(df_curve_g["Velocity"], df_curve_g["Power required"], label="Power required")
plt.plot(df_curve_g["Velocity"], P_a * np.ones(len(lst_Preq_g)), label="Power available")
plt.vlines(x=V_cl, ymin=0,
           ymax=df_curve_g[round(df_curve_g["Velocity"], 2) == round(V_cl, 2)]["Power required"].values[0],
           linestyle='dotted', color='green')
plt.xlabel("Velocity")
plt.ylabel("Power required")
plt.ylim(0.9 * min(df_curve_g["Power required"]), 1.1 * max(df_curve_g["Power required"]))
plt.legend()
plt.title("Power curve with g_load")
#plt.show()

""""###############################################################################################################"""
""""Take-off"""
""""###############################################################################################################"""

#integration for the determination of the landing distance forr a certain output of thrust
res, err = quad(f, 0.001, V_to ** 2)

print("The numerical result is {:f} (+-{:g})".format(res, err))

""""###############################################################################################################"""
""""Landing"""
""""###############################################################################################################"""


# print(Preq(0, 3200, 0.63, 28.29, 1) / eff_fan / eff_engine)

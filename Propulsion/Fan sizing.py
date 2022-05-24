import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# constants
eff_fan = 0.885  # [-] this value is taken from the snorri book, should be updated later on should reference
eff_engine = 0.95  # [-]
volt = 400  # [v]
e_d = 1  # [-] A4/AR, can be changed later
rho0 = 1.225  # [kg/m3]
h = 0  # [m]
AR = 5.84  # [-]
kt = 1.032914
kq = 0.1110977401
Q = 113

# definitions
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

def thrust(eff_fan, eff_engine, P_req, V):
    return (eff_fan * eff_engine * P_req) / V

def area_prop(Thr, V, P_req, rho, e_d):
    return (Thr ** (1.5) / P_req) ** 2 / (4 * rho * e_d)

def RPM(P_req, e_d, kt, h, D):
    T0 = 288.15
    p0 = 101325.0
    g0 = 9.80665
    R = 287.0
    a1 = -0.0065
    T = T0 + a1 * (h - 0)
    p = p0 * (T / T0) ** (-g0 / (R * a1))
    rho1 = (p / (R * T))
    return ((np.pi*P_req**2*e_d)/(kt**3*rho1**2*D**10))**(1/6)*60

def Preq(h,Thr, D,V,e_d):
    T0 = 288.15
    p0 = 101325.0
    g0 = 9.80665
    R = 287.0
    a1 = -0.0065
    T = T0 + a1 * (h - 0)
    p = p0 * (T / T0) ** (-g0 / (R * a1))
    rho1 = (p / (R * T))
    return 0.75 * Thr * V + np.sqrt(((Thr**2 * V**2)/16) + ((Thr**3) / (rho1 *np.pi*e_d*D**2)))

V_cl = 77.17  # [m/s]
V_to = 28.29
V_la = 33.438

P_required = 323000  # [W]
P_req_eff = P_required * eff_engine * eff_fan
lowerbound = P_req_eff - 100000
upperbound = P_req_eff + 100000
lst_cl_0= []
lst_cl_6000= []
lst_to_0= []
lst_to_6000= []
lst_la_0= []
lst_la_6000= []

alt_0 = 0
alt_6000= 1828

thr_thres_cl_0 = 2864
thr_thres_cl_6000 = 2394
thr_thres_to_0 = 1613
thr_thres_to_6000 = 1348
thr_thres_la_0 = 1812
thr_thres_la_6000 = 1598
rpmlower = 4000
rpmupper = 8000
rpmspacing = 100
marginlower = 1
marginupper = 1.1


for D_cl in np.arange(0.4, 1, 0.01):
    for e_d_cl in np.arange(0.8, 1.01, 0.01):
        for n_cl in np.arange(rpmlower/60, rpmupper/60, rpmspacing/60):
            mflow_cl_0 = ISA(alt_0)[2] * np.pi * D_cl ** 2 / 4 * V_cl
            Thr_cl_0 = kt * n_cl ** 2 * D_cl ** 4 * ISA(alt_0)[2]
            Torque_cl_0 = kq * n_cl ** 2 * D_cl ** 5 * ISA(alt_0)[2]
            Pr_cl_0 = Preq(alt_0, Thr_cl_0, D_cl, V_cl, e_d_cl)
            w_cl_0 = (0.5*e_d_cl-1)*V_cl+np.sqrt((e_d_cl*V_cl/2)**2+(e_d_cl*Thr_cl_0)/(ISA(alt_0)[2] * np.pi  *D_cl**2/4))
            Torque_power_cl_0 =  w_cl_0 * Torque_cl_0
            V_3_cl_0 = V_cl + w_cl_0
            V_4_cl_0 = V_3_cl_0 / e_d_cl
            if marginlower * thr_thres_cl_0 <= Thr_cl_0 <= marginupper * thr_thres_cl_0 and  P_required <= Pr_cl_0 <= 1.002 * P_required:
                for n_cl_2 in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                    mflow_cl_6000 = ISA(alt_6000)[2] * np.pi * D_cl ** 2 / 4 * V_cl
                    Thr_cl_6000 = kt * n_cl_2 ** 2 * D_cl ** 4 * ISA(alt_6000)[2]
                    Torque_cl_6000 = kq * n_cl_2 ** 2 * D_cl ** 5 * ISA(alt_6000)[2]
                    Pr_cl_6000 = Preq(alt_6000, Thr_cl_6000, D_cl, V_cl, e_d_cl)
                    w_cl_6000 = (0.5 * e_d_cl - 1) * V_cl + np.sqrt(
                        (e_d_cl * V_cl / 2) ** 2 + (e_d_cl * Thr_cl_6000) / (ISA(alt_6000)[2] * np.pi  *D_cl**2/4))
                    Torque_power_cl_6000 = w_cl_6000 * Torque_cl_6000
                    V_3_cl_6000 = V_cl + w_cl_6000
                    V_4_cl_6000 = V_3_cl_6000 / e_d_cl
                    if marginlower *thr_thres_cl_6000 <= Thr_cl_6000 <= marginupper * thr_thres_cl_6000:
                        for n_to in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                            mflow_to_0 = ISA(alt_0)[2] * np.pi * D_cl ** 2 / 4 * V_to
                            Thr_to_0 = kt * n_to ** 2 * D_cl ** 4 * ISA(alt_0)[2]
                            Torque_to_0 = kq * n_to ** 2 * D_cl ** 5 * ISA(alt_0)[2]
                            Pr_to_0 = Preq(alt_0,Thr_to_0, D_cl, V_to, e_d_cl)
                            w_to_0 = (0.5 * e_d_cl - 1) * V_to + np.sqrt(
                                (e_d_cl * V_to / 2) ** 2 + (e_d_cl * Thr_to_0) / (
                                            ISA(alt_0)[2] * np.pi  *D_cl**2/4))
                            Torque_power_to_0 = w_to_0 * Torque_to_0
                            V_3_to_0 = V_to + w_to_0
                            V_4_to_0 = V_3_to_0 / e_d_cl
                            if marginlower *thr_thres_to_0 <= Thr_to_0 <= marginupper * thr_thres_to_0:
                                for n_to_2 in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                                    mflow_to_6000 = ISA(alt_6000)[2] * np.pi * D_cl ** 2 / 4 * V_to
                                    Thr_to_6000 = kt * n_to_2 ** 2 * D_cl ** 4 * ISA(alt_6000)[2]
                                    Torque_to_6000 = kq * n_to_2 ** 2 * D_cl ** 5 * ISA(alt_6000)[2]
                                    Pr_to_6000 = Preq(alt_6000,Thr_to_6000, D_cl, V_cl, e_d_cl)
                                    w_to_6000 = (0.5 * e_d_cl - 1) * V_to + np.sqrt(
                                        (e_d_cl * V_to / 2) ** 2 + (e_d_cl * Thr_to_6000) / (
                                                ISA(alt_6000)[2] * np.pi  *D_cl**2/4))
                                    Torque_power_to_6000 = w_to_6000 * Torque_to_6000
                                    V_3_to_6000 = V_to + w_to_6000
                                    V_4_to_6000 = V_3_to_6000 / e_d_cl
                                    if marginlower *thr_thres_to_6000 <= Thr_to_6000 <= marginupper * thr_thres_to_6000:
                                        for n_la in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                                            mflow_la_0 = ISA(alt_0)[2] * np.pi * D_cl ** 2 / 4 * V_la
                                            Thr_la_0 = kt * n_la ** 2 * D_cl ** 4 * ISA(alt_0)[2]
                                            Torque_la_0 = kq * n_la ** 2 * D_cl ** 5 * ISA(alt_0)[2]
                                            Pr_la_0 = Preq(alt_0, Thr_la_0, D_cl, V_cl, e_d_cl)
                                            w_la_0 = (0.5 * e_d_cl - 1) * V_la + np.sqrt(
                                                (e_d_cl * V_la / 2) ** 2 + (e_d_cl * Thr_la_0) / (
                                                        ISA(alt_0)[2] * np.pi  *D_cl**2/4))
                                            Torque_power_la_0 = w_la_0 * Torque_la_0
                                            V_3_la_0 = V_la + w_la_0
                                            V_4_la_0 = V_3_la_0 / e_d_cl
                                            if marginlower *  thr_thres_la_0 <= Thr_la_0 <= marginupper *thr_thres_la_0:
                                                for n_la_2 in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                                                    mflow_la_6000 = ISA(alt_6000)[2] * np.pi * D_cl ** 2 / 4 * V_la
                                                    Thr_la_6000 = kt * n_la_2 ** 2 * D_cl ** 4 * ISA(alt_6000)[2]
                                                    Torque_la_6000 = kq * n_la_2 ** 2 * D_cl ** 5 * ISA(alt_6000)[2]
                                                    Pr_la_6000 = Preq(alt_6000, Thr_la_6000, D_cl, V_la, e_d_cl)
                                                    w_la_6000 = (0.5 * e_d_cl - 1) * V_la + np.sqrt((e_d_cl * V_la / 2) ** 2 + (e_d_cl * Thr_la_6000) / (
                                                                ISA(alt_6000)[2] * np.pi  *D_cl**2/4))
                                                    Torque_power_la_6000 = w_la_6000 * Torque_la_6000
                                                    V_3_la_6000 = V_la + w_la_6000
                                                    V_4_la_6000 = V_3_la_6000 / e_d_cl
                                                    if  marginlower * thr_thres_la_6000 <= Thr_la_6000 <= marginupper * thr_thres_la_6000:
                                                        lst_cl_0.append([D_cl,D_cl**2/4*np.pi, e_d_cl, max(n_cl * 60, n_cl_2 * 60, n_to * 60,
                                                                 n_to_2 * 60, n_la * 60, n_la_2 * 60), min(n_cl * 60, n_cl_2 * 60, n_to * 60,
                                                                 n_to_2 * 60, n_la * 60, n_la_2 * 60), n_cl * 60, n_cl_2 * 60, n_to * 60,
                                                                 n_to_2 * 60, n_la * 60, n_la_2 * 60,
                                                                 Pr_cl_0, Pr_cl_6000, Pr_to_0, Pr_to_6000, Pr_la_0,
                                                                 Pr_la_6000,
                                                                 Thr_cl_0, Thr_cl_6000, Thr_to_0, Thr_to_6000, Thr_la_0,
                                                                 Thr_la_6000,
                                                                 Torque_cl_0, Torque_cl_6000, Torque_to_0, Torque_to_6000,
                                                                 Torque_la_0, Torque_la_6000
                                                                    , mflow_cl_0, mflow_cl_6000, mflow_to_0, mflow_to_6000,
                                                                 mflow_la_0, mflow_la_6000, w_cl_0, w_cl_6000
                                                                            , w_to_0, w_to_6000, w_la_0, w_la_6000,
                                                                         V_3_cl_0,V_3_cl_6000,V_3_to_0,V_3_to_6000,V_3_la_0,V_3_la_6000,
                                                                         V_4_cl_0, V_4_cl_6000, V_4_to_0, V_4_to_6000,
                                                                         V_4_la_0, V_4_la_6000,
                                                                         Torque_power_cl_0,Torque_power_cl_6000, Torque_power_to_0,
                                                                         Torque_power_to_6000, Torque_power_la_0,
                                                                         Torque_power_la_6000])

df = pd.DataFrame(lst_cl_0, columns=[ "Diameter", "Area","e_d", "Max RPM", 'Min RPM',
                                        "RPM_cl_0", "RPM_cl_6000","RPM_to_0", "RPM_to_6000","RPM_la_0", "RPM_la_6000",
                                        "P_cl_0", "P_cl_6000","P_to_0", "P_to_6000","P_la_0", "P_la_6000",
                                        "Thrust_cl_0", "Thrust_cl_6000","Thrust_to_0", "Thrust_to_6000",
                                        "Thrust_la_0", "Thrust_la_6000", "Torque_cl_0","Torque_cl_6000",
                                        "Torque_to_0","Torque_to_6000","Torque_la_0","Torque_la_6000",
                                        "Mass_flow_cl_0", "Mass_flow_cl_6000","Mass_flow_to_0",
                                        "Mass_flow_to_6000","Mass_flow_la_0", "Mass_flow_la_6000",
                                     "Induced_airspeed_cl_0", "Induced_airspeed_cl_6000", "Induced_airspeed_to_0",
                                     "Induced_airspeed_to_6000", "Induced_airspeed_la_0","Induced_airspeed_la_0",
                                      "Rotor_speed_V3_cl_0", "Rotor_speed_V3_cl_6000","Rotor_speed_V3_to_0",
                                      "Rotor_speed_V3_to_6000","Rotor_speed_V3_la_0","Rotor_speed_V3_la_6000",
                                      "Exit_speed_V4_cl_0", "Exit_speed_V4_cl_6000", "Exit_speed_V4_to_0",
                                      "Exit_speed_V4_to_6000", "Exit_speed_V4_la_0", "Exit_speed_V4_la_6000",
                                      "Torque_power_cl_0", "Torque_power_cl_6000", "Torque_power_to_0", "Torque_power_to_6000",
                                     "Torque_power_la_0","Torque_power_la_0"])
df.to_excel("Prop_sizing_final.xlsx", index=False)

print(df)


plt.plot(df["Diameter"], df["Max RPM"], label = "Max RPM")
plt.plot(df["Diameter"], df["Min RPM"], label = "Min RPM")
plt.ylabel("RPM")
plt.xlabel("Diameter")
plt.legend()
plt.show()

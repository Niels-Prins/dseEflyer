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

P_required = 337000  # [W]
P_req_eff = P_required * eff_engine * eff_fan
lst_cl_0= []
lst_cl_6000= []
lst_to_0= []
lst_to_6000= []
lst_la_0= []
lst_la_6000= []
lowerbound = P_req_eff - 100000
upperbound = P_req_eff + 100000
alt_0 = 0
alt_6000= 1828
threshold = 240000

thr_thres_cl_0 = 2864
thr_thres_cl_6000 = 2394
thr_thres_to_0 = 1613
thr_thres_to_6000 = 1348
thr_thres_la_0 = 1812
thr_thres_la_6000 = 1598
rpmlower = 4000
rpmupper = 8000
rpmspacing = 50
margin = 1.03


for D_cl in np.arange(0.4, 1, 0.01):
    for e_d_cl in np.arange(0.8, 1.01, 0.01):
        for n_cl in np.arange(rpmlower/60, rpmupper/60, rpmspacing/60):
            mflow_cl_0 = ISA(alt_0)[2] * np.pi * D_cl ** 2 / 4 * V_cl
            Thr_cl_0 = kt * n_cl ** 2 * D_cl ** 4 * ISA(alt_0)[2]
            Torque_cl_0 = kq * n_cl ** 2 * D_cl ** 5 * ISA(alt_0)[2]
            Pr_cl_0 = Preq(alt_0, Thr_cl_0, D_cl, V_cl, e_d_cl)
            w_cl_0 = (0.5*e_d_cl-1)*V_cl+np.sqrt((e_d_cl*V_cl/2)**2+(e_d_cl*Thr_cl_0)/(ISA(alt_0)[2] *np.pi*e_d_cl*D_cl**2))
            if thr_thres_cl_0 <= Thr_cl_0 <= margin * thr_thres_cl_0:
                for n_cl_2 in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                    mflow_cl_6000 = ISA(alt_6000)[2] * np.pi * D_cl ** 2 / 4 * V_cl
                    Thr_cl_6000 = kt * n_cl_2 ** 2 * D_cl ** 4 * ISA(alt_6000)[2]
                    Torque_cl_6000 = kq * n_cl_2 ** 2 * D_cl ** 5 * ISA(alt_6000)[2]
                    Pr_cl_6000 = Preq(alt_6000, Thr_cl_6000, D_cl, V_cl, e_d_cl)
                    w_cl_6000 = (0.5 * e_d_cl - 1) * V_cl + np.sqrt(
                        (e_d_cl * V_cl / 2) ** 2 + (e_d_cl * Thr_cl_6000) / (ISA(alt_6000)[2] * np.pi * e_d_cl * D_cl ** 2))
                    if thr_thres_cl_6000 <= Thr_cl_6000 <= margin * thr_thres_cl_6000:
                        for n_to in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                            mflow_to_0 = ISA(alt_0)[2] * np.pi * D_cl ** 2 / 4 * V_to
                            Thr_to_0 = kt * n_to ** 2 * D_cl ** 4 * ISA(alt_0)[2]
                            Torque_to_0 = kq * n_to ** 2 * D_cl ** 5 * ISA(alt_0)[2]
                            Pr_to_0 = Preq(alt_0,Thr_to_0, D_cl, V_to, e_d_cl)
                            w_to_0 = (0.5 * e_d_cl - 1) * V_to + np.sqrt(
                                (e_d_cl * V_to / 2) ** 2 + (e_d_cl * Thr_to_0) / (
                                            ISA(alt_0)[2] * np.pi * e_d_cl * D_cl ** 2))
                            if thr_thres_to_0 <= Thr_to_0 <= margin * thr_thres_to_0:
                                for n_to_2 in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                                    mflow_to_6000 = ISA(alt_6000)[2] * np.pi * D_cl ** 2 / 4 * V_to
                                    Thr_to_6000 = kt * n_to_2 ** 2 * D_cl ** 4 * ISA(alt_6000)[2]
                                    Torque_to_6000 = kq * n_to_2 ** 2 * D_cl ** 5 * ISA(alt_6000)[2]
                                    Pr_to_6000 = Preq(alt_6000,Thr_to_6000, D_cl, V_cl, e_d_cl)
                                    w_to_6000 = (0.5 * e_d_cl - 1) * V_to + np.sqrt(
                                        (e_d_cl * V_to / 2) ** 2 + (e_d_cl * Thr_to_6000) / (
                                                ISA(alt_6000)[2] * np.pi * e_d_cl * D_cl ** 2))
                                    if thr_thres_to_6000 <= Thr_to_6000 <= margin * thr_thres_to_6000:
                                        for n_la in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                                            mflow_la_0 = ISA(alt_0)[2] * np.pi * D_cl ** 2 / 4 * V_la
                                            Thr_la_0 = kt * n_la ** 2 * D_cl ** 4 * ISA(alt_0)[2]
                                            Torque_la_0 = kq * n_la ** 2 * D_cl ** 5 * ISA(alt_0)[2]
                                            Pr_la_0 = Preq(alt_0, Thr_la_0, D_cl, V_cl, e_d_cl)
                                            w_la_0 = (0.5 * e_d_cl - 1) * V_la + np.sqrt(
                                                (e_d_cl * V_la / 2) ** 2 + (e_d_cl * Thr_la_0) / (
                                                        ISA(alt_0)[2] * np.pi * e_d_cl * D_cl ** 2))
                                            if  thr_thres_la_0 <= Thr_la_0 <= margin *thr_thres_la_0:
                                                for n_la_2 in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                                                    mflow_la_6000 = ISA(alt_6000)[2] * np.pi * D_cl ** 2 / 4 * V_la
                                                    Thr_la_6000 = kt * n_la_2 ** 2 * D_cl ** 4 * ISA(alt_6000)[2]
                                                    Torque_la_6000 = kq * n_la_2 ** 2 * D_cl ** 5 * ISA(alt_6000)[2]
                                                    Pr_la_6000 = Preq(alt_6000, Thr_la_6000, D_cl, V_la, e_d_cl)
                                                    w_la_6000 = (0.5 * e_d_cl - 1) * V_la + np.sqrt((e_d_cl * V_la / 2) ** 2 + (e_d_cl * Thr_la_6000) / (
                                                                ISA(alt_6000)[2] * np.pi * e_d_cl * D_cl ** 2))
                                                    if thr_thres_la_6000 <= Thr_la_6000 <= margin * thr_thres_la_6000:
                                                        lst_cl_0.append([D_cl, e_d_cl, max(n_cl * 60, n_cl_2 * 60, n_to * 60,
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
                                                                            , w_to_0, w_to_6000, w_la_0, w_la_6000])

df = pd.DataFrame(lst_cl_0, columns=["Diameter","e_d", "Max RPM", 'Min RPM',
                                        "RPM_cl_0", "RPM_cl_6000","RPM_to_0", "RPM_to_6000","RPM_la_0", "RPM_la_6000",
                                        "P_cl_0", "P_cl_6000","P_to_0", "P_to_6000","P_la_0", "P_la_6000",
                                        "Thrust_cl_0", "Thrust_cl_6000","Thrust_to_0", "Thrust_to_6000",
                                        "Thrust_la_0", "Thrust_la_6000", "Torque_cl_0","Torque_cl_6000",
                                        "Torque_to_0","Torque_to_6000","Torque_la_0","Torque_la_6000",
                                        "Mass_flow_cl_0", "Mass_flow_cl_6000","Mass_flow_to_0",
                                        "Mass_flow_to_6000","Mass_flow_la_0", "Mass_flow_la_6000",
                                     "Induced_airspeed_cl_0", "Induced_airspeed_cl_6000", "Induced_airspeed_to_0",
                                     "Induced_airspeed_to_6000", "Induced_airspeed_la_0","Induced_airspeed_la_0"])
df.to_excel("Prop_sizing_final.xlsx", index=False)
print(df)

plt.plot(df["Diameter"], df["Max RPM"], label = "Max RPM")
plt.plot(df["Diameter"], df["Min RPM"], label = "Min RPM")
plt.ylabel("RPM")
plt.xlabel("Diameter")
plt.legend()
plt.show()
#
# for D_cl in np.arange(0.1, 1.5, 0.01):
#     for n_cl in np.arange(5000 / 60, 8000 / 60, 50 / 60):
#         for e_d_cl in np.arange(0.8, 1.01, 0.01):
#             mflow_cl_6000 = ISA(alt_6000)[2] * np.pi * D_cl ** 2 / 4 * V_cl
#             Thr_cl_6000 = kt * n_cl ** 2 * D_cl ** 4 * ISA(alt_6000)[2]
#             Torque_cl_6000 = kq * n_cl ** 2 * D_cl ** 5 * ISA(alt_6000)[2]
#             Pr_cl_6000 = Preq(alt_6000, Thr_cl_6000, D_cl, V_cl, e_d_cl)
#             if Thr_cl_6000 >= thr_thres_cl_6000:
#                 lst_cl_6000.append(
#                     [ Pr_cl_6000,  Thr_cl_6000,  Torque_cl_6000, e_d_cl, D_cl, n_cl * 60,
#                      mflow_cl_6000])
#
# df_cl2 = pd.DataFrame(lst_cl_6000, columns=["P_6000", "Thrust_6000", "Torque_6000", "e_d", "Diameter", "RPM", "Mass_flow_0" ])
# df_cl2.to_excel("Prop_sizing_cl2.xlsx", index=False)
# print(df_cl2)




#
# for n_to in np.arange(5000/60, 8000/60, 50/60):
#     for D_to in np.arange(0.1, 1.5, 0.01):
#         for e_d_to in np.arange(0.8, 1.01, 0.01):
#             mflow_to_0 = ISA(alt_0)[2] * np.pi * D_to ** 2 / 4 * V_to
#             mflow_to_6000 = ISA(alt_6000)[2] * np.pi * D_to ** 2 / 4 * V_to
#             Thr_to_0 = kt * n_to ** 2 * D_to ** 4 * ISA(alt_0)[2]
#             Thr_to_6000 = kt * n_to ** 2 * D_to ** 4 * ISA(alt_6000)[2]
#             Torque_to_0 = kq * n_to ** 2 * D_to ** 5 * ISA(alt_0)[2]
#             Torque_to_6000 = kq * n_to ** 2 * D_to ** 5 * ISA(alt_6000)[2]
#             Pr_to_0 = Preq(alt_0,Thr_to_0, D_to, V_to, e_d_to)
#             Pr_to_6000 = Preq(alt_6000,Thr_to_6000, D_to, V_cl, e_d_to)
#             lst_to.append([Pr_to_0, Pr_to_6000, Thr_to_0, Thr_to_6000, Torque_to_0, Torque_to_6000, e_d_to, D_to, n_to * 60, mflow_to_0,mflow_to_6000])
#
# df_to = pd.DataFrame(lst_to, columns=["P_0","P_6000", "Thrust_0", "Thrust_6000", "Torque_0","Torque_6000", "e_d", "Diameter", "RPM", "Mass_flow_0", "Mass_flow_0" ])
# df_to.to_excel("Prop_sizing_to.xlsx", index=False)
# print(df_to)
#
# for n_la in np.arange(5000/60, 8000/60, 50/60):
#     for D_la in np.arange(0.1, 1.5, 0.01):
#         for e_d_la in np.arange(0.8, 1.01, 0.01):
#             mflow_la_0 = ISA(alt_0)[2] * np.pi * D_la ** 2 / 4 * V_la
#             mflow_la_6000 = ISA(alt_6000)[2] * np.pi * D_la ** 2 / 4 * V_la
#             Thr_la_0 = kt * n_la ** 2 * D_la ** 4 * ISA(alt_0)[2]
#             Thr_la_6000 = kt * n_la ** 2 * D_la ** 4 * ISA(alt_6000)[2]
#             Torque_la_0 = kq * n_la ** 2 * D_la ** 5 * ISA(alt_0)[2]
#             Torque_la_6000 = kq * n_la ** 2 * D_la ** 5 * ISA(alt_6000)[2]
#             Pr_la_0 = Preq(alt_0,Thr_la_0, D_la, V_cl, e_d_la)
#             Pr_la_6000 = Preq(alt_6000,Thr_la_6000, D_la, V_la, e_d_la)
#             lst_la.append([Pr_la_0, Pr_la_6000, Thr_la_0, Thr_la_6000, Torque_la_0, Torque_la_6000, e_d_la, D_la, n_la * 60,mflow_la_0, mflow_la_6000])
#
# df_la = pd.DataFrame(lst_la, columns=["P_0","P_6000", "Thrust_0", "Thrust_6000", "Torque_0","Torque_6000", "e_d", "Diameter", "RPM", "Mass_flow_0", "Mass_flow_0" ])
# df_la.to_excel("Prop_sizing_la.xlsx", index=False)
# print(df_la)
#
# Q = 113
# V = 77.17  # [m/s]
# P_req = 300000  # [W]
# P_req_eff = P_req * eff_engine * eff_fan
# lst = []
# lowerbound = P_req_eff - 1000
# upperbound = P_req_eff + 1000
# alt = 0
# for n in np.arange(50,9000, 50):
#     n = n/60
#     for D in np.arange(0.1, 1.5, 0.01):
#         for e_d in np.arange(0.8, 1, 0.01):
#             Thr = kt * n ** 2 * D ** 4 * ISA(alt)[2]
#             Torque = kq * n ** 2 * D ** 5 * ISA(alt)[2]
#             Pr = Preq(h, D, V, e_d)
#             if Pr >= lowerbound and Pr <= upperbound:
#                 lst.append([round(Pr), round(Thr), round(Torque), e_d, D, n * 60,ISA(alt)[2]*np.pi*D**2/4*V ])
#
#
# df = pd.DataFrame(lst, columns=["Pr", "Thrust","Torque", "e_d", "Diameter", "RPM", "Mass flow"])
# df.to_excel("Area vs Diameter graph 1.xlsx", index=False)
# print(df)

""""Important calculation of the propeller"""
# Q = 113
# V = 77.17  # [m/s]
# P_req = 300000  # [W]
# P_req_eff = P_req * eff_engine * eff_fan
# lst = []
# lowerbound = P_req_eff - 100000
# upperbound = P_req_eff + 100000
# alt = 1828
# for n in np.arange(5000,8000, 50):
#     n = n/60
#     for D in np.arange(0.1, 1.5, 0.01):
#         for e_d in np.arange(0.8, 1.01, 0.01):
#             Thr = kt * n ** 2 * D ** 4 * ISA(alt)[2]
#             Torque = kq * n ** 2 * D ** 5 * ISA(alt)[2]
#             Pr = Preq(h, D, V, e_d)
#             if Pr >= lowerbound and Pr <= upperbound:
#                 lst.append([(Pr), (Thr), (Torque), e_d, D, n * 60,ISA(alt)[2]*np.pi*D**2/4*V ])
#
# df = pd.DataFrame(lst, columns=["Pr", "Thrust","Torque", "e_d", "Diameter", "RPM", "Mass flow"])
# df.to_excel("Area vs Diameter graph 1.xlsx", index=False)
# print(df)

# #Take-off
# V = 28.29  # [m/s]
# P_req = 300000  # [W]
# P_req_eff = P_req * eff_engine * eff_fan
# lst = []
# lowerbound = P_req_eff - 100000
# upperbound = P_req_eff + 100000
# alt = 0
# for n in np.arange(5000, 8000, 50):
#     n = n / 60
#     for D in np.arange(0.1, 1.5, 0.01):
#         for e_d in np.arange(0.8, 1.01, 0.01):
#             Thr = kt * n ** 2 * D ** 4 * ISA(alt)[2]
#             Torque = kq * n ** 2 * D ** 5 * ISA(alt)[2]
#             Pr = Preq(h, D, V, e_d)
#             if Pr >= lowerbound and Pr <= upperbound:
#                 lst.append([Pr, Thr, Torque, e_d, D, n * 60, ISA(alt)[2] * np.pi * D ** 2 / 4 * V])



#
# for D_cl in np.arange(0.1, 1.5, 0.01):
#     for e_d_cl in np.arange(0.8, 1.01, 0.01):
#         for n_cl in np.arange(5000/60, 8000/60, 50/60):
#             mflow_cl_0 = ISA(alt_0)[2] * np.pi * D_cl ** 2 / 4 * V_cl
#             Thr_cl_0 = kt * n_cl ** 2 * D_cl ** 4 * ISA(alt_0)[2]
#             Torque_cl_0 = kq * n_cl ** 2 * D_cl ** 5 * ISA(alt_0)[2]
#             Pr_cl_0 = Preq(alt_0, Thr_cl_0, D_cl, V_cl, e_d_cl)
#             for n_cl_2 in np.arange(5000 / 60, 8000 / 60, 50 / 60):
#                 mflow_cl_6000 = ISA(alt_6000)[2] * np.pi * D_cl ** 2 / 4 * V_cl
#                 Thr_cl_6000 = kt * n_cl_2 ** 2 * D_cl ** 4 * ISA(alt_6000)[2]
#                 Torque_cl_6000 = kq * n_cl_2 ** 2 * D_cl ** 5 * ISA(alt_6000)[2]
#                 Pr_cl_6000 = Preq(alt_6000, Thr_cl_6000, D_cl, V_cl, e_d_cl)
#                 for n_to in np.arange(5000 / 60, 8000 / 60, 50 / 60):
#                     mflow_to_0 = ISA(alt_0)[2] * np.pi * D_cl ** 2 / 4 * V_to
#                     Thr_to_0 = kt * n_to ** 2 * D_cl ** 4 * ISA(alt_0)[2]
#                     Torque_to_0 = kq * n_to ** 2 * D_cl ** 5 * ISA(alt_0)[2]
#                     Pr_to_0 = Preq(alt_0,Thr_to_0, D_cl, V_to, e_d_cl)
#                     for n_to_2 in np.arange(5000 / 60, 8000 / 60, 50 / 60):
#                         mflow_to_6000 = ISA(alt_6000)[2] * np.pi * D_cl ** 2 / 4 * V_to
#                         Thr_to_6000 = kt * n_to_2 ** 2 * D_cl ** 4 * ISA(alt_6000)[2]
#                         Torque_to_6000 = kq * n_to_2 ** 2 * D_cl ** 5 * ISA(alt_6000)[2]
#                         Pr_to_6000 = Preq(alt_6000,Thr_to_6000, D_cl, V_cl, e_d_cl)
#                         for n_la in np.arange(5000 / 60, 8000 / 60, 50 / 60):
#                             mflow_la_0 = ISA(alt_0)[2] * np.pi * D_cl ** 2 / 4 * V_la
#                             Thr_la_0 = kt * n_la ** 2 * D_cl ** 4 * ISA(alt_0)[2]
#                             Torque_la_0 = kq * n_la ** 2 * D_cl ** 5 * ISA(alt_0)[2]
#                             Pr_la_0 = Preq(alt_0,Thr_la_0, D_cl, V_cl, e_d_cl)
#                             for n_la_2 in np.arange(5000 / 60, 8000 / 60, 50 / 60):
#                                 mflow_la_6000 = ISA(alt_6000)[2] * np.pi * D_cl ** 2 / 4 * V_la
#                                 Thr_la_6000 = kt * n_la_2 ** 2 * D_cl ** 4 * ISA(alt_6000)[2]
#                                 Torque_la_6000 = kq * n_la_2 ** 2 * D_cl ** 5 * ISA(alt_6000)[2]
#                                 Pr_la_6000 = Preq(alt_6000,Thr_la_6000, D_cl, V_la, e_d_cl)
#
#                 if thr_thres_cl_0 <= Thr_cl_0 <= 1.1 *thr_thres_cl_0 and thr_thres_cl_6000 <= Thr_cl_6000 <= 1.1* thr_thres_cl_6000 \
#                         and thr_thres_to_0 <= Thr_to_0 <= 1.1 *thr_thres_to_0 and thr_thres_to_6000 <= Thr_to_6000 <= 1.1* thr_thres_to_6000 \
#                         and thr_thres_la_0 <= Thr_la_0 <= 1.1 *thr_thres_la_0 and thr_thres_la_6000 <= Thr_la_6000 <= 1.1* thr_thres_la_6000:
#                     lst_cl_0.append([D_cl,e_d_cl,n_cl * 60, n_cl_2 * 60,n_to * 60, n_to_2 * 60, n_la * 60, n_la_2 * 60,
#                                      Pr_cl_0,Pr_cl_6000,Pr_to_0,Pr_to_6000, Pr_la_0,Pr_la_6000,
#                                      Thr_cl_0,Thr_cl_6000,Thr_to_0,Thr_to_6000,Thr_la_0,Thr_la_6000,
#                                      Torque_cl_0, Torque_cl_6000,Torque_to_0, Torque_to_6000,Torque_la_0, Torque_la_6000
#                                         , mflow_cl_0, mflow_cl_6000, mflow_to_0, mflow_to_6000, mflow_la_0, mflow_la_6000])
#
# df_cl = pd.DataFrame(lst_cl_0, columns=["Diameter","e_d",
#                                         "RPM_cl_0", "RPM_cl_6000","RPM_to_0", "RPM_to_6000","RPM_la_0", "RPM_la_6000",
#                                         "P_cl_0", "P_cl_6000","P_to_0", "P_to_6000","P_la_0", "P_la_6000",
#                                         "Thrust_cl_0", "Thrust_cl_6000","Thrust_to_0", "Thrust_to_6000",
#                                         "Thrust_la_0", "Thrust_la_6000", "Torque_cl_0","Torque_cl_6000",
#                                         "Torque_to_0","Torque_to_6000","Torque_la_0","Torque_la_6000",
#                                         "Mass_flow_cl_0", "Mass_flow_cl_6000","Mass_flow_to_0",
#                                         "Mass_flow_to_6000","Mass_flow_la_0", "Mass_flow_la_6000"])
# print(df_cl)

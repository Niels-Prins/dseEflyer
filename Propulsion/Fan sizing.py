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


#clean
V_cl = 77.17  # [m/s]
V_to = 28.29
V_la = 33.438

P_req = 300000  # [W]
P_req_eff = P_req * eff_engine * eff_fan
lst_cl= []
lst_to= []
lst_la= []
lowerbound = P_req_eff - 100000
upperbound = P_req_eff + 100000
alt_0 = 0
alt_6000= 1828
threshold = 240000/eff_fan/eff_engine


for n_cl in np.arange(5000/60, 8000/60, 50/60):
    for D_cl in np.arange(0.1, 1.5, 0.01):
        for e_d_cl in np.arange(0.8, 1.01, 0.01):
            mflow_cl_0 = ISA(alt_0)[2] * np.pi * D_cl ** 2 / 4 * V_cl
            mflow_cl_6000 = ISA(alt_6000)[2] * np.pi * D_cl ** 2 / 4 * V_cl
            Thr_cl_0 = kt * n_cl ** 2 * D_cl ** 4 * ISA(alt_0)[2]
            Thr_cl_6000 = kt * n_cl ** 2 * D_cl ** 4 * ISA(alt_6000)[2]
            Torque_cl_0 = kq * n_cl ** 2 * D_cl ** 5 * ISA(alt_0)[2]
            Torque_cl_6000 = kq * n_cl ** 2 * D_cl ** 5 * ISA(alt_6000)[2]
            Pr_cl_0 = Preq(alt_0, Thr_cl_0, D_cl, V_cl, e_d_cl)
            Pr_cl_6000 = Preq(alt_6000, Thr_cl_6000, D_cl, V_cl, e_d_cl)
            lst_cl.append([Pr_cl_0, Pr_cl_6000, Thr_cl_0, Thr_cl_6000, Torque_cl_0, Torque_cl_6000, e_d_cl, D_cl, n_cl * 60, mflow_cl_0, mflow_cl_6000])

df_cl = pd.DataFrame(lst_cl, columns=["P_0","P_6000", "Thrust_0", "Thrust_6000", "Torque_0","Torque_6000", "e_d", "Diameter", "RPM", "Mass_flow_0", "Mass_flow_0" ])
df_cl.to_excel("Prop_sizing_cl.xlsx", index=False)
print(df_cl)
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

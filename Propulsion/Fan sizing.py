import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
""""###############################################################################################################"""
""""Constants"""
""""###############################################################################################################"""
rho0 = 1.225  # [kg/m3]
kt = 1.032914
kq = 0.1110977401
Q = 113
V_cl = 77.17  # [m/s]
V_to = 28.29
V_la = 33.438
P_required = 323000  # [W]

lst= []

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

""""###############################################################################################################"""
""""Definitions"""
""""###############################################################################################################"""
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

def Massflow(alt,D,V):
    return ISA(alt)[2] * np.pi * D ** 2 / 4 * V

def Thrust(kt,n,D,alt):
    return kt * n ** 2 * D ** 4 * ISA(alt)[2]

def Torque(kq,n,D,alt):
    return kq * n ** 2 * D ** 5 * ISA(alt)[2]

def V_induced(e_d,V,D,Thr, alt):
    return (0.5 * e_d - 1) * V + np.sqrt( (e_d*V/2)**2 + (e_d*Thr) / (ISA(alt)[2] * np.pi *D**2/4))

def P_torque(w,Q):
    return w*Q

def Preq(alt,Thr, D,V,e_d):
    return 0.75 * Thr * V + np.sqrt(((Thr**2 * V**2)/16) + ((Thr**3) / (ISA(alt)[2] *np.pi*e_d*D**2)))

def V3(V,w):
    return V + w

def V4(V, e_d):
    return V / e_d

""""###############################################################################################################"""
""""Loop for calculation"""
""""###############################################################################################################"""
for D_cl in np.arange(0.4, 1, 0.01):
    for e_d in np.arange(0.8, 1.01, 0.01):
        for n_cl in np.arange(rpmlower/60, rpmupper/60, rpmspacing/60):
            mflow_cl_0 = Massflow(alt_0, D_cl, V_cl)
            Thr_cl_0 = Thrust(kt,n_cl,D_cl,alt_0)
            Torque_cl_0 = Torque(kq,n_cl,D_cl,alt_0)
            Pr_cl_0 = Preq(alt_0, Thr_cl_0, D_cl, V_cl, e_d)
            w_cl_0 = V_induced(e_d,V_cl,D_cl,Thr_cl_0, alt_0)
            Torque_power_cl_0 = P_torque(w_cl_0, Torque_cl_0)
            V_3_cl_0 = V3(V_cl,w_cl_0)
            V_4_cl_0 = V4(V_3_cl_0, e_d)
            if marginlower * thr_thres_cl_0 <= Thr_cl_0 <= marginupper * thr_thres_cl_0 and  P_required <= Pr_cl_0 <= 1.002 * P_required:
                for n_cl_2 in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                    mflow_cl_6000 = Massflow(alt_6000, D_cl, V_cl)
                    Thr_cl_6000 = Thrust(kt,n_cl_2,D_cl,alt_6000)
                    Torque_cl_6000 = Torque(kq,n_cl_2,D_cl,alt_6000)
                    Pr_cl_6000 = Preq(alt_6000, Thr_cl_6000, D_cl, V_cl, e_d)
                    w_cl_6000 = V_induced(e_d,V_cl,D_cl,Thr_cl_6000, alt_6000)
                    Torque_power_cl_6000 = P_torque(w_cl_6000, Torque_cl_6000)
                    V_3_cl_6000 = V3(V_cl,w_cl_6000)
                    V_4_cl_6000 = V4(V_3_cl_6000, e_d)
                    if marginlower *thr_thres_cl_6000 <= Thr_cl_6000 <= marginupper * thr_thres_cl_6000:
                        for n_to in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                            mflow_to_0 = Massflow(alt_0, D_cl, V_to)
                            Thr_to_0 = Thrust(kt,n_to,D_cl,alt_0)
                            Torque_to_0 = Torque(kq,n_to,D_cl,alt_0)
                            Pr_to_0 = Preq(alt_0,Thr_to_0, D_cl, V_to, e_d)
                            w_to_0 = V_induced(e_d,V_to,D_cl,Thr_to_0, alt_0)
                            Torque_power_to_0 = P_torque(w_to_0, Torque_to_0)
                            V_3_to_0 = V3(V_to,w_to_0)
                            V_4_to_0 = V4(V_3_to_0, e_d)
                            if marginlower *thr_thres_to_0 <= Thr_to_0 <= marginupper * thr_thres_to_0:
                                for n_to_2 in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                                    mflow_to_6000 = Massflow(alt_6000, D_cl, V_to)
                                    Thr_to_6000 = Thrust(kt,n_to_2,D_cl,alt_6000)
                                    Torque_to_6000 = Torque(kq,n_to_2,D_cl,alt_6000)
                                    Pr_to_6000 = Preq(alt_6000,Thr_to_6000, D_cl, V_cl, e_d)
                                    w_to_6000 = V_induced(e_d,V_to,D_cl,Thr_to_6000, alt_6000)
                                    Torque_power_to_6000 = P_torque(w_to_6000, Torque_to_6000)
                                    V_3_to_6000 = V3(V_to,w_to_6000)
                                    V_4_to_6000 = V4(V_3_to_6000, e_d)
                                    if marginlower *thr_thres_to_6000 <= Thr_to_6000 <= marginupper * thr_thres_to_6000:
                                        for n_la in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                                            mflow_la_0 = Massflow(alt_0, D_cl, V_la)
                                            Thr_la_0 = Thrust(kt,n_la,D_cl,alt_0)
                                            Torque_la_0 = Torque(kq,n_la,D_cl,alt_0)
                                            Pr_la_0 = Preq(alt_0, Thr_la_0, D_cl, V_cl, e_d)
                                            w_la_0 = V_induced(e_d,V_la,D_cl,Thr_to_0, alt_0)
                                            Torque_power_la_0 = P_torque(w_la_0, Torque_la_0)
                                            V_3_la_0 = V3(V_la,w_la_0)
                                            V_4_la_0 = V4(V_3_la_0, e_d)
                                            if marginlower *  thr_thres_la_0 <= Thr_la_0 <= marginupper *thr_thres_la_0:
                                                for n_la_2 in np.arange(rpmlower / 60, rpmupper / 60, rpmspacing / 60):
                                                    mflow_la_6000 = Massflow(alt_6000, D_cl, V_la)
                                                    Thr_la_6000 = Thrust(kt,n_la_2,D_cl,alt_6000)
                                                    Torque_la_6000 = Torque(kq,n_la_2,D_cl,alt_6000)
                                                    Pr_la_6000 = Preq(alt_6000, Thr_la_6000, D_cl, V_la, e_d)
                                                    w_la_6000 = V_induced(e_d,V_la,D_cl,Thr_to_6000, alt_6000)
                                                    Torque_power_la_6000 = P_torque(w_la_6000, Torque_la_6000)
                                                    V_3_la_6000 = V3(V_la,w_la_6000)
                                                    V_4_la_6000 = V4(V_3_la_6000, e_d)
                                                    if  marginlower * thr_thres_la_6000 <= Thr_la_6000 <= marginupper * thr_thres_la_6000:
                                                        lst.append([D_cl,D_cl**2/4*np.pi, e_d, max(n_cl * 60, n_cl_2 * 60, n_to * 60,
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

df = pd.DataFrame(lst, columns=[ "Diameter", "Area","e_d", "Max RPM", 'Min RPM',
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
df_new = df.loc[326]
print(df_new)
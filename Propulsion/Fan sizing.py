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

def Preq(h,D,V,e_d):
    T0 = 288.15
    p0 = 101325.0
    g0 = 9.80665
    R = 287.0
    a1 = -0.0065
    T = T0 + a1 * (h - 0)
    p = p0 * (T / T0) ** (-g0 / (R * a1))
    rho1 = (p / (R * T))
    return 0.75 * Thr * V + np.sqrt(((Thr**2 * V**2)/16) + ((Thr**3) / (rho1 *np.pi*e_d*D**2)))

Q = 113
V = 77.17  # [m/s]
P_req = 260000  # [W]
P_req_eff = P_req * eff_engine * eff_fan
lst = []
lowerbound = P_req_eff - 1000
upperbound = P_req_eff + 1000
alt = 0
for n in np.arange(50,9000, 50):
    n = n/60
    for D in np.arange(0.1, 1.5, 0.01):
        Thr = kt*n**2*D**4*ISA(alt)[2]
        Torque = kq*n**2*D**5*ISA(alt)[2]
        Pr = Preq(h,D,V,e_d)
        if Pr >= lowerbound and Pr <= upperbound:
            lst.append([round(Pr), round(Thr), round(Torque), D, n * 60])

df = pd.DataFrame(lst, columns=["Pr", "Thrust","Torque", "Diameter", "RPM"])
df.to_excel("Area vs Diameter graph 1.xlsx", index=False)
print(df)
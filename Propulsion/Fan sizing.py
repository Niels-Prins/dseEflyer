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

# def Preq(kt,h,n,D,V,e_d):
#     T0 = 288.15
#     p0 = 101325.0
#     g0 = 9.80665
#     R = 287.0
#     a1 = -0.0065
#     T = T0 + a1 * (h - 0)
#     p = p0 * (T / T0) ** (-g0 / (R * a1))
#     rho1 = (p / (R * T))
#     return 0.75 * kt * rho1 * n**2 * D**4 + np.sqrt(((kt**2 * rho1**2 * n**4 * D**8 * V**2)/16) + ((kt**3 * rho1**2 * n**6 * D**10) / (np.pi*e_d)))

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


V = 77.17  # [m/s]
P_req = 320000  # [W]
lst = []
lowerbound = P_req - 1000
upperbound = P_req + 1000
alt = 0
for n in np.arange(50,9000, 50):
    n = n/60
    for D in np.arange(0.1, 1.5, 0.01):
        Thr = kt*n**2*D**4*ISA(alt)[2]
        Pr = Preq(h,D,V,e_d)
        if Pr >= lowerbound and Pr <= upperbound:
            lst.append([round(Pr), round(Thr), D, n * 60])

df = pd.DataFrame(lst, columns=["Pr", "Thrust", "Diameter", "RPM"])
df.to_excel("Area vs Diameter graph 1.xlsx", index=False)

print(df)







# n = 6000
# D = 0.6
# Preq_df = Preq(kt, 0, n, D, V, e_d)
# print(Preq_df)

#
# lst = []
# for n in np.arange(1000,9000, 50):
#     for D in np.arange(0.30, 1.50, 0.01):
#         Preq_df = Preq(kt, 0, n, D, V, e_d)
#         lst.append([Preq_df, D, n*60])
#         # df = df.append([{'Pr': Preq}, {'D': D}, {'RPM': n*60}], ignore_index=True)
#
# df = pd.DataFrame(lst, columns=["Pr", "Diamter", "RPM"])
# df.to_excel("Area vs Diameter graph.xlsx", index=False)





# df = pd.DataFrame(columns =["Pr", "D", "RPM"])
# df = df.append(1, 1, 1, ignore_index=False)
# print(df)
#
# for n in np.arange(0,9000, 50):
#     for D in np.arange(0.01, 1.5, 0.01):
#         Preq_df = Preq(kt, 0, n, D, V, e_d)
#         df = df.append([{'Pr': Preq}, {'D': D}, {'RPM': n*60}], ignore_index=True)
# df.to_excel("Area vs Diameter graph.xlsx", index=False)


# D = np.arange(0.5,1,0.001)
#
# Thr = thrust(eff_fan, eff_engine, P_req, V)
# rho = ISA(h)[2]
# area_prop = area_prop(Thr, V, P_req, rho, e_d)
#
# df = pd.DataFrame()
# df["Diameter"] = D
# df["Area"] = np.pi*(D/2)**2
# df["RMP_0"] = RPM(P_req, e_d, kt, 0, D)
# df["RMP_500"] = RPM(P_req, e_d, kt, 500, D)
# df["RMP_1000"] = RPM(P_req, e_d, kt, 1000, D)
# df["RMP_1500"] = RPM(P_req, e_d, kt, 1500, D)
#
# plt.plot(df["RMP_0"],df["Diameter"], label = "0 meters altitude")
# plt.plot(df["RMP_500"],df["Diameter"], label = "500 meters altitude")
# plt.plot(df["RMP_1000"],df["Diameter"], label = "1000 meters altitude")
# plt.plot(df["RMP_1500"],df["Diameter"], label = "1500 meters altitude")
# plt.xlabel("RPM")
# plt.ylabel("Diameter")
# plt.title("Diameter vs RPM")
# plt.legend()
# plt.show()
#
# print("Thrust:", Thr, "N")
# print("Propeller area:", area_prop, "m2")
# print("Propeller diameter:", 2 * np.sqrt(area_prop / (np.pi)), "m")
# print("volumetric flow:", area_prop * V)



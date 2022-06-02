import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

Tau_c = 22.43*10**6
Tau_al = 207*10**6
#Tau-steel =
Tor = 800 #nm
sf = 2 #-
shaft_length = 0.953 #m
rho_c = 1800
Rho_al = 2710 #kg/m3
lst = []

#
# def Polarmom(c1,c2):
#     # global J
#     J = 0.5*np.pi*((c2**4))-((c1**4))
#     return J
#
# def Torque(c2,J, Tau):
#     # global T
#     T = Tau*J/(c2)
#     return T

for c1 in np.arange(0.01,0.05,0.0001):
    for c2 in np.arange(c1+0.0001,0.05,0.0001):
        J = 0.5*np.pi*(c2**4-c1**4)
        T = Tau_al * J / c2
        if Tor * sf < T < Tor*1.05*sf:
            Area = np.pi*(c2**2-c1**2)
            Volume = Area*shaft_length
            mass = Volume*Rho_al
            lst.append([c1*100, c2*100, Area, mass, Volume, T])
            # Torque(c2,J,Tau_al)
            # if Tor * sf <= Torque(c2,X,Tau_al) <= Tor*1.05 * sf:
            #     Area = np.pi*(c2**2-c1**2)
            #     Volume = Area*shaft_length
            #     mass = Volume*rho
            #     lst.append([c1,c2,Area, mass, Volume, Torque(c2,X,Tau_al)])

df = pd.DataFrame(lst, columns = ["c1", "c2", "Area","mass", "Volume", "Torque"])
df.to_excel("Driveshaft.xlsx", index=False)
print(df)


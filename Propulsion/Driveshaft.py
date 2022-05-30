import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

Tau_al = 206*10**6
#Tau-steel =
Tor = 800 #nm
sf = 2 #-
shaft_length = 0.928 #m
rho = 2710 #kg/m3
lst = []


def Polarmom(c1,c2):
    global J
    J = (np.pi/2)*(c2**4-c1**4)
    return J

def Torque(c2,J, Tau):
    global T
    T = Tau*J/(c2)
    return T

for c1 in np.arange(0,0.03,0.0001):
    for c2 in np.arange(0.01,0.05,0.0001):
        if c1 < c2:
            Polarmom(c1,c2)
            Torque(c2,J,Tau_al)
            if Tor*0.99 * sf <= T <= Tor*1.01 * sf:
                Area = np.pi*(c2**2-c1**2)
                Volume = Area*shaft_length
                mass = Volume*rho
                lst.append([c1,c2,Area, mass, Volume])

df = pd.DataFrame(lst, columns = ["c1", "c2", "Area","mass", "Volume"])
df.to_excel("Driveshaft.xlsx", index=False)
print(df)


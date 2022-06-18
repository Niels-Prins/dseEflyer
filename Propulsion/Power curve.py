import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

MTOM = 1153 #kg
MTOM_3g = 3* MTOM
g = 9.81
rho = 1.225 # kg/m^3
S = 14.5 #m^2
Ar = 5.82
e = 0.8 # clean config
Cd0 = 0.025
D = 0.63 #m
e_d = 0.84
n = 7400/60
kt = 1.032914
lstV = []
lstpreq = []
lstpower = []

def Cl(V, MTOM):
    return (MTOM*g)/(0.5*rho*V**2*S)

def Cd(V):
    return Cd0 + Cl(V)**2/(np.pi*Ar*e)


def Drag(V):
    return 0.5*rho*V**2*S*Cd(V)

def Preq(V):
    return 0.75 * Drag(V) * V + np.sqrt(((Drag(V)**2 * V**2)/16) + ((Drag(V)**3) / (rho *np.pi*e_d*D**2)))

def Pav(V):
    #global T
    T = kt * n ** 2 * D ** 4 * rho
    return 0.75 * T * V + np.sqrt(((T ** 2 * V ** 2) / 16) + ((T ** 3) / (rho * np.pi * e_d * D ** 2)))

for V in range(70,90):
    #lift_c = Cl(V)
    #Drag_c = Cd(lift_c)
    #D = Drag(V)
    Prequired = Preq(V)
    Poweravailable = Pav(V)
    lstV.append(V)
    lstpower.append(Poweravailable)
    lstpreq.append(Prequired)

plt.plot(lstV, lstpreq)
plt.xlabel("V")
plt.ylabel("P")
plt.show()

print(Drag(77.17))
print(Cd(77.17))
print(Cl(77.17))
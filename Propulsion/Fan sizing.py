import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

"""This file will size the ducted fan"""
"""-----------------------------------------------------------------------------------"""

#constants
s_runway = 700                  #[m/s]
V = 77.17                       #[m/s]
P_req = 119000                  #[W]
eff_fan = 0.885                 #[-] this value is taken from the snorri book, should be updated later on
eff_engine = 0.95               #[-]
volt = 400                      #[v]
e_d = 1                         #[-] A4/AR, can be changed later
rho0 = 1.225                    #[kg/m3]
h = 3000


#definitions
def ISA(h):
    T0 = 288.15
    p0 = 101325.0
    g0 = 9.80665
    R = 287.0
    a1= -0.0065
    T = T0 + a1 *(h-0)
    p = p0 * (T/T0)**(-g0/(R*a1))
    rho1 = (p/(R*T))
    return [T,p,rho1]

def thrust(eff_fan, eff_engine, P_req, V):
    return (eff_fan*eff_engine*P_req)/V

def area_prop(Thr, V, P_req, rho, e_d):
    return ((Thr**3)/((P_req-0.75*Thr*V)**2-((Thr**2*V**2)/16)))/(4*rho*e_d)


Thr = thrust(eff_fan, eff_engine, P_req, V)
rho = ISA(h)[2]
area_prop = area_prop(Thr, V, P_req, rho, e_d)

print("Propeller area", area_prop, "m2")
print("Propeller radius", np.sqrt(area_prop/(2*np.pi)), "m")
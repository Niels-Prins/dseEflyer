import numpy as np
import math as mt
from scipy import interpolate
import scipy as sc
import matplotlib.pyplot as plt

#PROPULSION UNIT PARAMETERS
D_fan = 0.63
pi_fan = 1.0648 #taken from UL-39
shaft_rad = 0.5 #assumed
a = 0.9         #assumed fuselage height behind the cockpit (consider clearance for FOD and wing placement)
kt = 1.032914   #Wessel computed this
e_d = 1         #Wessel decided this
thrust = 3200   #Wessel computed this
v3 = 121.55     #Wessel computed this

#FLIGHT CONDITIONS
v = 77.17 #operating velocity (freestream)

h = 0 * 0.3048
kin_vis = 1.460 * (10**-5)
kappa = 1.4
cp = 1000
T0 = 288.15
p0 = 101325.0
g0 = 9.80665
R = 287.0
a1 = -0.0065
T = T0 + a1 * (h - 0)
p = p0 * (T / T0) ** (-g0 / (R * a1))
rho = (p / (R * T))


def diam(area):
    diameter = np.sqrt((area * 4) / np.pi)
    return diameter

def velocities(freestream, aft_fan):
    v_0 = freestream
    v_3 = aft_fan
    v_4 = aft_fan/e_d
    return [v_0, v_3, v_4]

def massflow(thrust, V_4, freestram):
    massflow = thrust/(V_4 - freestram)
    return massflow

def W(Re, velocity, length, diam):
    inv_fd = 1.93/mt.log(10) * sc.special.lambertw(10**(-0.537/1.930) * mt.log(10)/1.930 * Re)
    f_d = np.real((inv_fd **-1) ** 2)
    head_loss = f_d * ((rho * velocity**2)/2) * (length / diam)
    return [f_d, head_loss]

def Reynolds(velocity, diameter):
    r = (velocity * diameter) / kin_vis
    return r

def total_conditions(v, temp, pres):
    tot_temp = temp + (v ** 2 / (2 * cp))
    tot_pres = pres * (tot_temp/temp)**(kappa/(kappa-1))
    return [tot_temp, tot_pres]

def inlet_conditions(v, temp, pres, in_eff):
    tot_temp_in = temp + (v ** 2 / (2 * cp))
    tot_pres_in = pres * (1 + in_eff*(v**2/(2*cp*tot_temp_in)))**(kappa/(kappa - 1))
    return [tot_temp_in, tot_pres_in]

def DynP(velocity):
    q = 0.5 * rho * (velocity**2)
    return q



fan_area = ((D_fan**2)*np.pi)/4
exhaust_area = fan_area*e_d
exhaust_diam = diam(exhaust_area)
v1 = velocities(v, v3)[0]
v4 = velocities(v, v3)[2]
mdot = massflow(thrust, v4, v)


Re3 = Reynolds(v3, D_fan)
Re4 = Reynolds(v4, D_fan*e_d)

if e_d < 1:
    mean_velocity34 = (v4-v3)/2
    mean_HyDiam34 = (D_fan - exhaust_diam)/2
    p_loss34 = W(Re3, mean_velocity34, 1.4, mean_HyDiam34)[1]
if e_d == 1:
    p_loss34 = W(Re3, v3, 1.4, D_fan)[1]

#stagnation point
T_0 = total_conditions(v, T, p)[0]
p_0 = total_conditions(v, T, p)[1]


p4 = p #ambient pressure (same method as UL39)
ptot_4 = p4 + DynP(v4)
ptot_3 = ptot_4 + p_loss34
p3 = ptot_3 - DynP(v3)

ptot_2 = ptot_3 / pi_fan
noz_eff = 0.85
fan_eff = 0.885
ram_eff = (ptot_2 - p)/(p_0 - p)

T_tot4 = (((ptot_4/p4) ** (kappa-1/kappa)) * T)
T_tot3 = T_tot4 / (1 + ((1/noz_eff) * ((ptot_4/ptot_3)**(kappa-1/kappa) - 1)))
T_tot2 = T_tot3 / (1 + ((1/fan_eff) * ((ptot_3/ptot_2)**(kappa-1/kappa) - 1)))


duct_acc_v = 80 - 62.5 #velocity increase due to the duct for the UL-39 (inlet velocity - freestream)
A2 = np.arange(fan_area, 2*fan_area, 0.01)
fanfaceV = []
fanfaceA = []
for i in range(len(A2)):
    vel2 = (fan_area*v3)/A2[i]
    if mt.isclose(vel2, v + duct_acc_v, rel_tol=0.05):
        fanfaceV.append(vel2)
        fanfaceA.append(A2[i])

v2 = np.max(fanfaceV)
fanface_area = np.min(fanfaceA)

Re2 = Reynolds(v2, diam(fanface_area))
ptot_2 = ptot_3 / pi_fan
p_loss12 = W(Re2, v2, 1.4, diam(fanface_area))[1]
p2 = ptot_2 - DynP(v2)
ptot_1 = ptot_2 - p_loss12
p1 = ptot_1 - DynP(v)

mdot_per_duct = mdot/2
ductend_area = mdot_per_duct/(v2 * rho)
inlet_area = (v2*ductend_area)/v

b = (2 * inlet_area) / (np.pi * a)
base_track = np.linspace(-1*(a/2), 1*(a/2), 3)
cowl_points = [0, b, 0]
cowl = interpolate.interp1d(base_track, cowl_points, kind='quadratic')
x = np.arange(-1*(a/2), 1*(a/2), 0.0001)
plt.xlim(-0.5, 0.5)
plt.ylim(-0.05, 0.5)
plt.plot(x, (cowl(x)), color='blue')
plt.plot([-1*(a/2), 1*(a/2)], [np.min(cowl(base_track)), np.min(cowl(base_track))], color='blue')
plt.grid(True)
plt.show()

print("Local Air Density:", rho)
print("Local Air Pressure:", p)
print("Local Air Temperature:", T)

print("Reynold's Number at Inlet:", Reynolds(v, diam(inlet_area)))
print("Reynold's Number Before Fan:", Re2)
print("Reynold's Number After Fan:", Re3)
print("Reynold's Number at Exhaust:", Re4)

print("Pressure loss between INLET & FAN:", p_loss12)
print("Pressure loss between FAN & EXHAUST:", p_loss34)

#print("Total Temperature at SP:", T_0)
#print("Total Pressure at SP:", p_0)

#print("Total Temperature at Inlet:", T_tot1)
print("Total pressure at Inlet:", ptot_1)
#print("Static Pressure Before Fan:", p2)

#print("Total Temperature Before Fan:", T_tot2)
print("Total pressure Before Fan:", ptot_2)
#print("Static Pressure Before Fan:", p2)

#print("Total Temperature After Fan:", T_tot3)
print("Total Pressure After Fan:", ptot_3)
#print("Static Pressure After Fan:", p3)

#print("Total Temperature at Exhaust:", T_tot4)
print("Total Pressure at Exhaust:", ptot_4)
#print("Static Pressure at Exhaust:", p4)

#print("Inlet duct pressure ratio:", ptot_2/ptot_1)
print("Fan pressure ratio:", pi_fan)
print("Intake duct pressure ratio:", ptot_2/ptot_1)
print("Exhaust duct pressure ratio:", ptot_4/ptot_3)

print("V_1 (assumed freestream):", v)
print("V_2:", v2)
print("Velocity increase due to duct convergence (taken from UL-39):", duct_acc_v)
print("V_3:", v3)
print("Fan induced velocity, w:", v3-v2)
print("V_4:", v4)
#print("delta V:", v4 - v)
print("Mass Flow:", mdot)

print("Inlet Area per side:", inlet_area)
print("Duct End Area per side:", ductend_area)
print("Fan Face Area:", fanface_area)
print("Fan Area:", fan_area)
print("Exhaust area:", exhaust_area)

print("Exhaust diameter:", diam(exhaust_area))
print("Fan face diameter:", diam(fanface_area))
print("Cowl max width:", b)
print("Cowl max height:", a)
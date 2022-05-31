import numpy as np
import math as mt
from scipy.optimize import fsolve
import scipy as sc
import matplotlib.pyplot as plt


D_fan = 0.63
pi_fan = 1.0648
shaft_rad = 0.5
a = 0.6 #fuselage heigh behind the cockpit
kt = 1.032914
e_d = 0.84
duct_taper = 0.95

h = 0 * 0.3048
thrust = 2810
v = 77.17
kin_vis = 1.460 * (10**-5)
w = 47.222
v3 = 124.39

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


def diam(area): #a3
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
    inv_fd = 1.93/mt.log(10) * sc.special.lambertw(10**(-0.537/1.930) * mt.log(10)/1.930 *Re)
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

def SoS(T):
    SoS = np.sqrt(287 * kappa * T)
    return SoS

fan_area = ((D_fan**2)*np.pi)/4
exhaust_area = fan_area*e_d
exhaust_diam = diam(exhaust_area)
v1 = velocities(v, v3)[0]
v4 = velocities(v, v3)[2]
mdot = massflow(thrust, v4, v)

#pressure loss stuff
Re3 = Reynolds(v3, D_fan)
Re4 = Reynolds(v4, D_fan*e_d)

mean_velocity34 = (v4-v3)/2
mean_HyDiam34 = (D_fan - exhaust_diam)/2

p_loss34 = W(Re3, mean_velocity34, 1.4, mean_HyDiam34)[1]

#stagnation point
T_0 = total_conditions(v, T, p)[0]
p_0 = total_conditions(v, T, p)[1]


p4 = p #ambient pressure (same method as UL39)
ptot_4 = p4 + DynP(v4)
ptot_3 = ptot_4 + p_loss34
p3 = ptot_3 - DynP(v3)

#UL39 numbers

ptot_2 = ptot_3 / pi_fan
ram_eff = (ptot_2 - p)/(p_0 - p)
p2 = ptot_2 / (1 + (ram_eff * (v**2/(2*cp*T_0))))
#v2 = np.sqrt(2*cp*(ptot_2-p2))

#inlet conditions
T_tot1 = inlet_conditions(v, T_0, p, ram_eff)[0]
ptot_1 = inlet_conditions(v, T_0, p, ram_eff)[1]

isen_eff = 0.85
fan_eff = 0.885

T_tot4 = (((ptot_4/p4) ** (kappa-1/kappa)) * T)
T_tot3 = T_tot4 / (1 + ((1/isen_eff) * ((ptot_4/ptot_3)**(kappa-1/kappa) - 1)))
T_tot2 = T_tot3 / (1 + ((1/fan_eff) * ((ptot_3/ptot_2)**(kappa-1/kappa) - 1)))

v2 = SoS(T_tot2) * np.sqrt(((ptot_2/p2)**((kappa-1)/kappa) - 1) / ((kappa-1)/2))



'''
duct_acc_v = 23.5 #velocity increase due to the duct for the UL-39
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

mdot_per_duct = mdot/2
ductend_area = mdot_per_duct/(v2 * rho)
inlet_area = (v2*ductend_area)/v'''

print("Local Air Density:", rho)
print("Local Air Pressure:", p)
print("Local Air Temperature:", T)

#print("Reynold's Number at Inlet:", Re1)
print("Reynold's Number at Fan:", Re3)
print("Reynold's Number at Exhaust:", Re4)

#print("Pressure loss between INLET & FAN:", p_loss12)
print("Pressure loss between FAN & EXHAUST:", p_loss34)

print("Total Temperature at SP:", T_0)
print("Total Pressure at SP:", p_0)

#print("Total Temperature at Inlet:", T_tot1)
#print("Total pressure at Inlet:", ptot_1)
#print("Static Pressure Before Fan:", p2)

print("Total Temperature Before Fan:", T_tot2)
print("Total pressure Before Fan:", ptot_2)
print("Static Pressure Before Fan:", p2)

#print("Total Temperature After Fan:", T_tot3)
print("Total Pressure After Fan:", ptot_3)
#print("Static Pressure After Fan:", p3)

#print("Total Temperature at Exhaust:", T_tot4)
print("Total Pressure at Exhaust:", ptot_4)
#print("Static Pressure at Exhaust:", p4)

#print("Inlet duct pressure ratio:", ptot_2/ptot_1)
print("Fan pressure ratio:", pi_fan)
print("Exhaust duct pressure ratio:", ptot_4/ptot_3)




#print("Total Temperature at Duct:", T_02)
#print("Total Pressure at Duct:", p_02)
print("V_0:", v)
print("V_2:", v2)
#print("Induced velocity, w:", v3-v2)
print("V_3:", v3)
print("V_4:", v4)
print("delta V:", v4 - v)
print("Mass Flow:", mdot)
#print("Inlet Area per side:", inlet_area)
#print("Duct End Area per side:", ductend_area)
#print("Fan Face Area:", fanface_area)
print("Fan Area:", fan_area)
#print("Fan diameter:", fan_diameter)
print("Exhaust area:", exhaust_area)
#print("Exhaust diameter:", exhaust_diam)


'''
inlet_diam = np.sqrt((areas(D_fan)[0] * 4)/np.pi)
fan_diam = np.sqrt((areas(D_fan)[1] * 4)/np.pi)
exit_diam = np.sqrt((areas(D_fan)[2] * 4)/np.pi)'''


'''
plt.plot([(-1*(b/2)), (1*(b/2))], [0, 0])
base_track = [-1*(b/2), 0, 1*(b/2)]
base_track_new = np.linspace(-1*(b/2), 1*(b/2), 1000)
cowl_points = [0, b, 0]
cowl = (base_track, cowl_points, base_track_new, order=2)
plt.plot(base_track_new, cowl)
plt.show()'''
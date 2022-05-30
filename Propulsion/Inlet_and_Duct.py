#%%
import numpy as np
import math as mt
import scipy.interpolate
import scipy as sc
import matplotlib.pyplot as plt
import pandas as pd

D_fan = 0.63
shaft_rad = 0.5
a = 0.6 #fuselage heigh behind the cockpit
kt = 1.032914
e_d = 0.84

h = 0
thrust = 2921.911
v = 77.17
v3 = 118.171

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


def fan_diam(thr, velocity): #a3
    area = thr / (rho * v3 * ((v3/e_d) - velocity))
    diameter = np.sqrt((area * 4) / np.pi)
    return [diameter, area]

def velocities(freestream, aft_fan):
    v_0 = freestream
    v_3 = aft_fan
    v_4 = aft_fan/e_d
    return [v_0, v_3, v_4]


def massflow(thrust, V_4, freestram):
    massflow = thrust/(V_4 - freestram)
    return massflow

def total_conditions(v, temp, pres):
    tot_temp = temp + (v ** 2 / (2 * cp))
    tot_pres = pres * (tot_temp/temp)**(kappa/(kappa-1))
    return [tot_temp, tot_pres]

def inlet_conditions(v, temp, pres, in_eff):
    tot_temp_in = temp + (v ** 2 / (2 * cp))
    tot_pres_in = pres * (1 + in_eff*(v**2/(2*cp*tot_temp_in)))**(kappa/(kappa - 1))
    return [tot_temp_in, tot_pres_in]

def smooth(Re):
    inv_fd = 1.93/mt.log(10) * sc.special.lambertw(10**(-0.537/1.930) * mt.log(10)/1.930 *Re)
    return inv_fd **-1



'''
w = ((0.5*e_d - 1)*v) + np.sqrt( ((e_d*v)/2)**2 + ((e_d*fan_thr_Ar_Vel(rho, rpm, D_fan, v)[0])/rho*fan_thr_Ar_Vel(rho, rpm, D_fan, v)[1]))
m_dot = fan_thr_Ar_Vel(rho, (rpm), D_fan, v)[2] * fan_thr_Ar_Vel(rho, rpm, D_fan, v)[1] * rho
exhaust_area = e_d * fan_thr_Ar_Vel(rho, rpm, D_fan, v)[1]
exhaust_vel = ((v + w)/e_d)'''

'''
inlet_area = (1.1 * fan_area(D_fan))
b = (2 * inlet_area) / (np.pi * a)
AR_inlet = a/b
inlet_efficiency = 0.98
duct_efficiency = 0.85
fan_pressure_ratio = 1.0648 #from UL-39, p2/p3
p_01 = total_conditions(v, T, p)[1]
p_02 = inlet_conditions(v, total_conditions(v, T, p)[0], total_conditions(v, T, p)[1], inlet_efficiency)[1]
T_02 = inlet_conditions(v, T, p, inlet_efficiency)[0]
p_03 = fan_pressure_ratio * p_02
T_03 = T_02 + (T_02/duct_efficiency)* (((fan_pressure_ratio) ** ((kappa - 1)/kappa)) - 1)
'''



#fan_face_area = (v3 * fan_area(D_fan))/v2

#fan_face_diam = np.sqrt(((fan_face_area/2) * 4)/np.pi)
# mdot = massflow(thrust, v4)
# fan_diameter = fan_diam(thrust, v)[0]
# fan_area = fan_diam(thrust, v)[1]
# exhaust_area = fan_diam(thrust, v)[1]*e_d
# exhaust_diam = np.sqrt((exhaust_area * 4)/np.pi)
# v4 = velocities(v, v3)[2]


# #print("V_2:", v2)'''
# print("V_3:", v3)
# print("V_4:", v4)
# print("delta V:", v4 - v)
# print("Mass Flow:", mdot)
# print("Fan diameter:", fan_diameter)
# print("Fan Area:", fan_area)
# print("Exhaust area:", exhaust_area)
# print("Exhaust diameter:", exhaust_diam)


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
# %%

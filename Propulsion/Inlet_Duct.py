Python 3.8.5 (v3.8.5:580fbb018f, Jul 20 2020, 12:11:27) 
[Clang 6.0 (clang-600.0.57)] on darwin
Type "help", "copyright", "credits" or "license()" for more information.
>>> import numpy as np
import math as mt
from scipy import interpolate
import scipy as sc
import matplotlib.pyplot as plt

h = 0 * 0.3048 # ft
kin_vis = 1.460 * (10**-5)
kappa = 1.4
cp = 1000
T0 = 288.15 # K
p0 = 101325.0 # Pa
g0 = 9.80665
R = 287.0
a1 = -0.0065
T = T0 + a1 * (h - 0) # K
p = p0 * (T / T0) ** (-g0 / (R * a1)) # Pa
rho = (p / (R * T)) # kg/m3

T0 = T
p0 = p

def inv_ptot(tot_pres, tot_temp, temp):
    pres = tot_pres / (tot_temp/temp)**(kappa/(kappa-1))
    return pres

def total_temp(v, temp):
    tot_temp = temp + (v ** 2 / (2 * cp))
    return tot_temp

def total_pres(pres, tot_temp, temp):
    tot_pres = pres * (tot_temp/temp)**(kappa/(kappa-1))
    return tot_pres

def W(Re, velocity, length, diam):
    inv_fd = 1.93/mt.log(10) * sc.special.lambertw(10**(-0.537/1.930) * mt.log(10)/1.930 * Re)
    f_d = np.real((inv_fd **-1) ** 2)
    head_loss = f_d * ((rho * velocity**2)/2) * (length / diam)
    return head_loss

def Reynolds(velocity, diameter):
    r = (velocity * diameter) / kin_vis
    return r

Pi_fan = 1.0614 # ..
kt = 1.032914 #
RPM = 6000 # RPM fixed
P_eng = 260000 # W fixed
thrust = 2060# N 
v0 = 77.17 # m/s

# Estimate of Stage and Duct Efficiency
fan_eff = 0.885
duct_eff = 0.9 
internal_eff = fan_eff * duct_eff

P_thrust = thrust * v0
P_av = P_eng * internal_eff # Power available W

# Engine Rating 
Omega = np.pi * RPM /30 # angular velocity rad/sec
Tau = 550 * (P_av * 1.341) / Omega * 1.3558 # required engine torque (Nm)

# Preliminary fan sizing
n = RPM/60

D_fan = 0.633

r_tip = D_fan/2 # tip radius fan
r_hub = 0.2 * r_tip # hub radius fan 
A2 = np.pi * (r_tip**2 - r_hub**2) # fan area

# Mass flow - first estimate
m_dot = (thrust**2)/(2*(P_av-P_thrust))

# Capture area
A0 = 0.466 # m2
  

def A2r(A):
    r = np.sqrt(A/np.pi)
    return r

# POWER AND PROPULSION METHOD
print('P&P Method')

Ttot0 = total_temp(v0, T0)
ptot0 = total_pres(p0, Ttot0, T0)
A0b = A0/2
# bifurcation intersection 1
taper01 = 0.9
A1b = A0b * taper01
A1 = A1b * 2
v1 = A0b * v0 /A1b
Ttot1 = Ttot0
T1 = Ttot1 - (v1**2/(2*cp))
p1 = p0 + 1/2 * rho * (v0**2 - v1**2)
ptot1 = total_pres(p1, Ttot1, T1)

v_avg01 = (v0+v1)/2
D_avg01 = A2r(A0b) + A2r(A1b)
ploss1b = W(Reynolds(v_avg01, D_avg01 ),v_avg01, 1, D_avg01)
ploss1 = 2 * ploss1b
ptot1bloss = ptot1 - ploss1b
ptot1loss = ptot1 - ploss1

deltap1 = ptot1loss - ptot0 


# fan face 2
v2 = A1 * v1 /A2
Ttot2 = Ttot1
T2 = Ttot2 - (v2 ** 2 / (2 * cp))
p2 = p0 + 1/2 * rho * (v0**2 - v2**2)
ptot2 = total_pres(p2, Ttot2, T2)

v_avg12 = (v2+v1)/2
D_avg12 = A2r(A0b) + A2r(A1)
ploss2 = W(Reynolds(v_avg12, D_avg12 ),v_avg12, 0.4, D_avg12)
ptot2loss = ptot2 - ploss2

deltap2 =  ptot2loss -ptot0


vnet_tip = np.sqrt((RPM* np.pi * np.sqrt(A2/np.pi)/30)**2+ v2**2)
c2 = np.sqrt(kappa * p0/rho)
M_tip = vnet_tip/c2

# after fan face 3
ptot3 = Pi_fan * ptot2loss
Ttot3 = Ttot2 + Ttot2/fan_eff * (Pi_fan ** ((kappa-1)/kappa)-1)
T3 = (Pi_fan * T2** (kappa/(kappa -1)))**((kappa-1)/kappa)
v3 = np.sqrt((Ttot3 - T3)* 2* cp)
p3 = ptot3 / (Ttot3/T3)**(kappa/(kappa -1))
A3 = A2

# exhaust 4
taper34 = 0.98
A4 = A3 * taper34
Ttot4 = Ttot3
v4 = A3 * v3 / A4
T4 = Ttot4 - (v4 ** 2 / (2 * cp))
p4 = p3 + 1/2 * rho * (v3**2 - v4**2)
ptot4 = total_pres(p4, Ttot4, T4)

v_avg34 = (v3+v4)/2
D_avg34 = A2r(A3) + A2r(A4)
ploss4 = W(Reynolds(v_avg34, D_avg34),v_avg34, 1.3, D_avg34)
ptot4loss = ptot4 - ploss4

deltap4 = ptot4loss - ptot0
print(f"Delta P 1 : {int(deltap1)}[Pa]")
print(f"Delta P 2: {int(deltap2)}[Pa]")
print(f"Delta P 4: {int(deltap4)}[Pa]\n")
print(f"Mach tip : {round(M_tip,2)}[M]\n")

T_fan = m_dot * (v4 - v0)
T_res = A4 * (p4 - p0)
T = T_fan + T_res

print(f"Mass Flow: {(m_dot)}[kg/s]\n")
print(f"Fan Thrust: {(thrust)}, Exhaust Thrust: {(T_fan)}[N], Residual Thrust: {(T_res)}[N]\n")
print(f"Total Thrust: {(T)}[N]\n")

total_ptotloss = ploss1 + ploss2 + ploss4
print(f"Complete total pressure loss: {int(total_ptotloss)}[Pa]\n")
print(f"Pressure loss in inlet duct: {(ploss2)}[Pa]\n")

print(f"A0: {round(A0,3)}[m2] \nA1: {round(A1,3)}[m2] \nA2: {round(A2,3)}[m2] \nA3: {round(A3,3)}[m2] \nA4: {round(A4,3)}[m2]\n")
print(f"p0: {(p0)}[Pa] \np1: {(p1)}[Pa] \np2: {(p2)}[Pa] \np3: {(p3)}[Pa] \np4: {(p4)}[Pa]\n")
print(f"p00: {(ptot0)}[Pa] \np01: {(ptot1)}[Pa] \np02: {(ptot2)}[Pa] \np03: {(ptot3)}[Pa] \np04: {(ptot4)}[Pa]\n")
print(f"T0: {(T0)}[K] \nT1: {(T1)}[K] \nT2: {(T2)}[K] \nT3: {(T3)}[K] \nT4: {(T4)}[K]\n")
print(f"T00: {(Ttot0)}[K] \nT01: {(Ttot1)}[K] \nT02: {(Ttot2)}[K] \nT03: {(Ttot3)}[K] \nT04: {(Ttot4)}[K]\n")
print(f"v0: {(v0)}[m/s] \nv1: {(v1)}[m/s] \nv2: {(v2)}[m/s] \nv3: {(v3)}[m/s] \nv4: {(v4)}[m/s]\n")

r0, r1, r2, r3, r4 = A2r(A0), A2r(A1), A2r(A2), A2r(A3), A2r(A4)

plt.plot([0,1,1.4,1.5,2.8], [r0, r1, r2, r3, r4], linestyle='solid', color='k', label='Single')

plt.axvline(x=0, linestyle='-.', label='0')
plt.axvline(x=1, linestyle='-.', label='1')
plt.axvline(x=1.4, linestyle='-.', label='2')
plt.axvline(x=1.5, linestyle='-.', label='3')
plt.axvline(x=2.8, linestyle='-.', label='4')
plt.xlabel("Length of the propulsion system, inlet to exhaust[m]")
plt.ylabel("Radius of ducting[m]")
plt.minorticks_on()
plt.grid(True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.grid(True, which='major', color='#666666', linestyle='-')
plt.ylim(0, 0.4)
plt.show()

plt.plot([0,1,1.4,1.5,2.8], [ptot0/1000, ptot1/1000, ptot2/1000, ptot3/1000, ptot4/1000], linestyle='solid', color='r', label='Nom.')
plt.plot([0,1,1.4,1.5,2.8], [ptot0, ptot1loss, ptot2loss, ptot3, ptot4loss], linestyle='solid', color='b', label='Loss')
plt.minorticks_on()
plt.grid(True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.grid(True, which='major', color='#666666', linestyle='-')
plt.xlabel("Length of the propulsion system, inlet to exhaust[m]")
plt.ylabel("Total Pressures[kPa]")
plt.show()

a = 0.8
b = (A0b*2) / (np.pi * (a / 2))
base_track = np.linspace(-1*(a/2), 1*(a/2), 3)
cowl_points = [0, b, 0]
cowl = interpolate.interp1d(base_track, cowl_points, kind='quadratic')
x = np.arange(-1*(a/2), 1*(a/2), 0.0001)
plt.ylim(-(a/2)-0.1, (a/2)+0.1)
plt.xlim(-0.05, 1.2)
plt.plot((cowl(x)), x, color='blue')
plt.plot([np.min(cowl(base_track)), np.min(cowl(base_track))], [-1*(a/2), 1*(a/2)], color='blue')
plt.xlabel("Aircraft Y-axis[m]")
plt.ylabel("Aircraft Z-axis[m]")
plt.minorticks_on()
plt.grid(True, which='minor', color='#999999', linestyle='-', alpha=0.2)
plt.grid(True, which='major', color='#666666', linestyle='-')
plt.show()


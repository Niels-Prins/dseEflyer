import numpy as np
import matplotlib.pyplot as plt

n_crew = 2  # Requirements
weight_crew = 80 * n_crew
weight_batteries = 424  # [kg] Class 1
n_batteries = int(weight_batteries / 27.7) + 1  # ELEO series 35
volume_batteries = n_batteries * (0.864 * 0.303 * 0.08)  # [m^3]
weight_engines = 60  # [kg] Class 1
n_engines = 2
volume_1_engine = np.pi * (0.15 ** 2) * 0.30  # [m^2] Saluqi
diameter_addition = 0.0508  # [m] Roskam
weight_landing_gear = 1  # TODO Find out
area_landing_gear = 0.125  # [m^2] # estimation assuming that LG is 0.5 m and folds in half
prop_diameter = 1  # TODO Find out
Seat_dimensions = [0.558, 1.01]  # m [width,height]
print("n_batteries", n_batteries)
print("volume batteries",volume_batteries)
# To pack: Crew, Batteries, propeller, engines, Landing gear


#### Landing gear sizing
xcg = np.arange(3.7, 4.0 + 0.1 , 0.1)  # [m]
ycg = 0.5  # [m]
Wto = 1220  # [kg]
Dtmg = 127*2  # [mm] assuming Dt = dt*2
Btmg = 127  # [mm]
volume_mw = np.pi * (Dtmg * 1E-3 * 0.5)**2 * (Btmg * 1E-3) # m^3
Dtng = 177.8  # [mm] assuming Dt = dt*2
Btng = 101.6  # [mm]
volume_nw = np.pi * (Dtng * 1E-3 * 0.5)**2 * (Btng * 1E-3) # m^3
ns = 1
x_ng = 1.29  # [m] Start of cockpit
bn = xcg[-1] - x_ng
l_ng = 0.8 # [m]
l_mg = l_ng # [m]
strut_length_nose = l_ng - (Dtng*1E-3 * 0.5) + (69.7E-3) # m
print("strut nose", strut_length_nose)
thick_fuselage = 7.686 # Position of thickest part of fuselage bottom used for ground clearance
ln = xcg - x_ng
print("Volume main gear", volume_mw)
print("volume nose gear", volume_nw)
# Longitudinal tip over
ycg_array = np.ones_like(xcg) * ycg
long_tip_over = ((ycg + l_ng) * np.tan(np.radians(15))) + xcg
long_tip_over_y = np.ones_like(long_tip_over) * -(l_ng)
x_mg = xcg + (np.tan(np.radians(15)) * (l_mg + ycg))
y_mg = np.ones_like(x_mg) * -l_mg
plt.plot([xcg, long_tip_over], [ycg_array, long_tip_over_y], label="lateral tip over")
for i in range(0, len(x_mg)):
    plt.plot(x_mg[i], y_mg[i], label=i, marker="x")

print("nose gear", x_ng, l_ng)
print("main gear", x_mg, l_mg)
# ground clearance
fuselage_left = thick_fuselage - x_mg
thick_part = np.ones_like(x_mg) * thick_fuselage
y_thick_part = l_mg + 0.576
ground_clearance = np.degrees(np.arctan(y_thick_part / fuselage_left))
print("ground clearance angles:", ground_clearance)
plt.plot([thick_part, x_mg], [np.zeros_like(y_mg), y_mg])
# lateral tip over
lm = x_mg - xcg
numer = ln + lm
denom = (ln ** 2 * (np.tan(np.radians(55))) ** 2 / (ycg + l_mg) ** 2) - 1
z_mg_min = numer / np.sqrt(denom)
print("z_mg_min", z_mg_min)
'''
d = x_mg - x_ng
alpha = np.degrees(np.arctan((ycg+l_mg) / d))
c = bn * np.sin(np.radians(alpha))
phi = np.degrees(np.arctan(7 * (ycg+l_mg) / c))
print("phi", phi)
'''

Pnw = (Wto * lm) / (lm + ln)
Pmw = (Wto * ln) / (2 * (lm + ln))
Loading_ratio = (Pnw / Wto) * 100
print("Pmw and Pnw are", Pmw, Pnw)
print("Loading ratio", Loading_ratio)
plt.legend()
plt.show()

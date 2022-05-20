import numpy as np
import matplotlib.pyplot as plt

n_crew = 2  # Requirements
weight_crew = 80 * n_crew
weight_batteries = 400  # [kg] Class 1
n_batteries = int(weight_batteries / 27.7) + 1  # ELEO series 35
volume_batteries = n_batteries * (0.864 * 0.303 * 0.08)  # [m^2]
weight_engines = 60  # [kg] Class 1
n_engines = 2
volume_1_engine = np.pi * (0.15 ** 2) * 0.30  # [m^2] Saluqi
diameter_addition = 0.0508  # [m] Roskam
weight_landing_gear = 1  # TODO Find out
area_landing_gear = 0.125  # [m^2] # estimation assuming that LG is 0.5 m and folds in half
prop_diameter = 1  # TODO Find out
Seat_dimensions = [0.558, 1.01]  # m [width,height]

# To pack: Crew, Batteries, propeller, engines, Landing gear


#### Landing gear sizing
xcg = np.arange(3.2, 3.4 + 0.1, 0.1)  # [m]
ycg = 0.5  # [m]
Wto = 1035  # [kg]
Dtmg = 431.8  # [mm]
Btmg = 152.4  # [mm]
Dtng = 12.5  # [mm]
Btng = 5  # [mm]


ns = 1
x_ng = 0.5  # [m] Start of cockpit
bn = xcg[-1] - x_ng
l_ng = 1.0
l_mg = l_ng
thick_fuselage = 5.5  # Position of thickest part of fuselage bottom used for ground clearance
ln = xcg-x_ng


# Longitudinal tip over
ycg_array = np.ones_like(xcg) * ycg
long_tip_over = ((ycg + l_ng) * np.tan(np.radians(15))) + xcg
long_tip_over_y = np.ones_like(long_tip_over) * -(l_ng)
x_mg = xcg + (np.tan(np.radians(15)) * (l_mg + ycg))
y_mg = np.ones_like(x_mg) * -l_mg
plt.plot([xcg, long_tip_over], [ycg_array, long_tip_over_y], label="lateral tip over")
for i in range(0, len(x_mg)):
    plt.plot(x_mg[i], y_mg[i], label=i, marker="x")


print("nose gear", x_ng,l_ng)
print("main gear", x_mg, l_mg)
# ground clearance
fuselage_left = thick_fuselage - x_mg
thick_part = np.ones_like(x_mg) * thick_fuselage
ground_clearance = np.degrees(np.arctan(l_mg / fuselage_left))
print("ground clearance angles:", ground_clearance)
plt.plot([thick_part, x_mg], [np.zeros_like(y_mg), y_mg])
#lateral tip over
lm = x_mg - xcg
numer = ln+lm
denom = (ln**2 * (np.tan(np.radians(55)))**2 / (ycg+l_mg)**2) -1
z_mg_min = numer / np.sqrt(denom)
print("z_mg_min", z_mg_min)
'''
d = x_mg - x_ng
alpha = np.degrees(np.arctan((ycg+l_mg) / d))
c = bn * np.sin(np.radians(alpha))
phi = np.degrees(np.arctan(7 * (ycg+l_mg) / c))
print("phi", phi)
'''
Pmw = 0.92 * Wto / 2
Pnw = 0.08 * Wto
Pnw = (Wto* lm)/(lm+ln)
Pmw = (Wto* ln)/(lm+ln)
Loading_ratio = Pnw/(Pnw+Pmw) *100
print("Pmw and Pnw are", Pmw, Pnw)
print("Loading ratio",Loading_ratio)
plt.legend()
plt.show()

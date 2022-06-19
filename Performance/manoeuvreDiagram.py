#%%
import numpy as np
import math as mt
import matplotlib.pyplot as plt

kts = 1/0.5144444
metersps = 1/kts

density = 1.225
MAC = 1.46
S = 12.3
W = 978 * 9.80065
CLmax = 1.6
CLalpha = 0.11 * (180 / np.pi)
a = CLalpha * 1.1
CNmax = CLmax * 1
n_max = 6
n_min = -4
wing_loading = W / S

mu_g = (2 * wing_loading) / (density * MAC * a * 9.80065)
k_g = (0.88 * mu_g) / (5.3 + mu_g)
U_de_cruise = 0.3048 * 50
U_de_dive = 0.3048 * 25

V = np.linspace(0, 160 * metersps, 1000)
V_cruise = 150 * 0.514444
V_dive = 1.25 * V_cruise
q = 0.5 * V ** 2  * density

n_pos = (q * CNmax * S) / W
n_neg = (q * -CNmax * S) / W

n_gust_cruise_pos = 1 + (
    (k_g * density * U_de_cruise * V_cruise * a) / (2 * wing_loading)
)
n_gust_cruise_neg = 1 - (
    (k_g * density * U_de_cruise * V_cruise * a) / (2 * wing_loading)
)

n_gust_dive_pos = 1 + ((k_g * density * U_de_dive * V_dive * a) / (2 * wing_loading))
n_gust_dive_neg = 1 - ((k_g * density * U_de_dive * V_dive * a) / (2 * wing_loading))

for i in range(len(V)):
    # if mt.isclose(n_pos[i], 1, rel_tol=0.01):
    #     V_stall = V[i]
    if mt.isclose(n_neg[i], n_min, rel_tol=0.01):
        V_prime = V[i]

V_stall = (W /(0.5*density*CLmax*S))**0.5
V_manu = V_stall * np.sqrt(n_max)

plt.xlim(0, V_dive*kts+10)
plt.ylim(-7, 9)

# Speed 
plt.vlines(V_stall*kts , ymax=1, ymin=0, color="red", label="Stall speed", linestyle="--")
# plt.vlines(
#     V_prime, ymax=0, ymin=n_min, color="black", linestyle="--"
# )
plt.vlines(
    V_cruise*kts ,
    ymax=n_max,
    ymin=n_min,
    color="fuchsia",
    label="Cruise speed",
    linestyle="--",
)

n_pos = [n_pos[i] for i in range(len(n_pos)) if n_pos[i] <=n_max]
n_neg = [n_neg[i] for i in range(len(n_neg)) if n_neg[i] >=n_min]

plt.vlines(V_dive*kts , ymax=n_max, ymin=n_min, color="Blue", label="Dive speed")

vplot = V *kts
# Manoeuvre Diagram
plt.hlines(n_max, xmin=V_manu*kts, xmax=V_dive*kts, color="purple", label="Manoeuvre")
plt.hlines(n_min, xmin=V_prime*kts, xmax=V_dive*kts, color="purple")
plt.plot(vplot[:len(n_pos)], n_pos, color="purple")
plt.plot(vplot[:len(n_neg)], n_neg, color="purple")

# Dust Diagram
plt.plot([0, V_cruise*kts], [1, n_gust_cruise_pos], color="orange", label="Gust")
plt.plot([0, V_cruise*kts], [1, n_gust_cruise_neg], color="orange")
plt.plot([0, V_dive*kts], [1, n_gust_dive_pos], color="orange", linestyle="--")
plt.plot([0, V_dive*kts], [1, n_gust_dive_neg], color="orange", linestyle="--")
plt.plot([V_cruise*kts, V_dive*kts], [n_gust_cruise_pos, n_gust_dive_pos], color="orange")
plt.plot([V_cruise*kts, V_dive*kts], [n_gust_cruise_neg, n_gust_dive_neg], color="orange")
plt.plot(
    [V_dive*kts, V_dive*kts], [n_gust_dive_pos, n_gust_dive_neg], color="orange", linestyle="--"
)

#Set the grid nicely
plt.grid(visible=True, which="major", color="#666666", linestyle="-")
plt.minorticks_on()
plt.grid(visible=True, which="minor", color="#999999", linestyle="-", alpha=0.2)
plt.yticks(np.arange(-6,9,1))

plt.legend(loc="upper left")
plt.xlabel("Operating Velocities [kts]")
plt.ylabel("G-load [-]")
plt.savefig("vndiagram2P")
plt.show()
#%%
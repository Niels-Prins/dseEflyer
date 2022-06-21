#%%
import numpy as np
import math as mt
import matplotlib.pyplot as plt

kts = 1/0.5144444
metersps = 1/kts

density = 1.225
MAC = 1.46
S = 12.3
W = 876 * 9.80065
CLmax = 1.6
CLalpha = 0.11 * (180 / np.pi)
a = CLalpha * 1.1
CNmax = CLmax * 1
n_max = 8
n_min = -6
wing_loading = W / S

CLmax_flap = 1.92

# Gust load parameters
mu_g = (2 * wing_loading) / (density * MAC * a * 9.80065)
k_g = (0.88 * mu_g) / (5.3 + mu_g)
U_de_cruise = 0.3048 * 50
U_de_dive = 0.3048 * 25

# Speed initializations
V = np.linspace(0, 160 * metersps, 1000)
V_cruise = 150 * 0.514444
V_dive = 1.25 * V_cruise
q = 0.5 * V ** 2  * density

n_posClean = (q * CNmax * S) / W
n_negClean = (q * -CNmax * S) / W

n_posFlap = (q * CLmax_flap * S) / W

n_gust_cruise_pos = 1 + (
    (k_g * density * U_de_cruise * V_cruise * a) / (2 * wing_loading)
)
n_gust_cruise_neg = 1 - (
    (k_g * density * U_de_cruise * V_cruise * a) / (2 * wing_loading)
)

n_gust_dive_pos = 1 + ((k_g * density * U_de_dive * V_dive * a) / (2 * wing_loading))
n_gust_dive_neg = 1 - ((k_g * density * U_de_dive * V_dive * a) / (2 * wing_loading))

for i in range(len(V)):
    # if mt.isclose(n_posClean[i], 1, rel_tol=0.01):
    #     V_stall = V[i]
    if mt.isclose(n_negClean[i], n_min, rel_tol=0.01):
        V_prime = V[i]
    if mt.isclose(n_posFlap[i], 2, rel_tol=0.01):
        V_sf = V[i]

V_stall = (W /(0.5*density*CLmax*S))**0.5
V_manu = V_stall * np.sqrt(n_max)
V_f = 1.8*V_sf

plt.xlim(0, V_dive*kts+10)
plt.ylim(-7, 9)

plt.vlines(V_stall*kts , ymax=1.01, ymin=0, color="brown", linestyle="--")
plt.vlines(
    V_cruise*kts ,
    ymax=n_max,
    ymin=n_min,
    color="brown",
    linestyle="--",
)
plt.vlines(V_dive*kts , ymax=n_max, ymin=n_min, color="blue")

n_posClean = [n_posClean[i] for i in range(len(n_posClean)) if n_posClean[i] <=n_max]
n_negClean = [n_negClean[i] for i in range(len(n_negClean)) if n_negClean[i] >=n_min]

n_posFlap = [n_posFlap[i] for i in range(len(n_posFlap)) if n_posFlap[i] <=2]



vplot = V *kts
# Manoeuvre Diagram
plt.hlines(n_max, xmin=V_manu*kts, xmax=V_dive*kts, color="blue", label="Manoeuvre Clean")
plt.hlines(n_min, xmin=V_prime*kts, xmax=V_dive*kts, color="blue")
plt.plot(vplot[:len(n_posClean)], n_posClean, color="blue")
plt.plot(vplot[:len(n_negClean)], n_negClean, color="blue")
plt.plot(vplot[:len(n_posFlap)], n_posFlap, color="red", label="Manoeuvre Flaps")
plt.hlines(2,vplot[:len(n_posFlap)][-1], V_f*kts, color="red",)
plt.vlines(V_f*kts,2,0,color="red")
# Dust Diagram
plt.plot([0, V_cruise*kts], [1, n_gust_cruise_pos], color="orange", label="Gust")
plt.plot([0, V_cruise*kts], [1, n_gust_cruise_neg], color="orange")
plt.plot([0, V_dive*kts], [1, n_gust_dive_pos], color="orange", linestyle="--")
plt.plot([0, V_dive*kts], [1, n_gust_dive_neg], color="orange", linestyle="--")
plt.plot([V_cruise*kts, V_dive*kts], [n_gust_cruise_pos, n_gust_dive_pos], color="orange")
plt.plot([V_cruise*kts, V_dive*kts], [n_gust_cruise_neg, n_gust_dive_neg], color="orange")

#Set the grid nicely
plt.grid(visible=True, which="major", color="#666666", linestyle="-")
plt.minorticks_on()
plt.grid(visible=True, which="minor", color="#999999", linestyle="-", alpha=0.2)
plt.yticks(np.arange(-6,9,1))

plt.text(V_stall*kts,0.1,"$V_s$")
plt.text(V_f*kts,0.1,"$V_f$")
plt.text(V_cruise*kts,0.1,"$V_c$")
plt.text(V_dive*kts,0.1,"$V_D$")

plt.legend(loc="upper left")
plt.xlabel("Operating Velocities [kts]")
plt.ylabel("G-load [-]")
plt.savefig("vndiagram")
plt.show()
#%%

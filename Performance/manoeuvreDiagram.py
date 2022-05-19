import numpy as np
import math as mt
import matplotlib.pyplot as plt

density = 1.225
MAC = 1.5
S = 13.1
W = 1035 * 9.80065
CLmax = 1.6
CLalpha = 0.11 * (180/np.pi)
a = CLalpha * 1.1
CNmax = CLmax * 1.1
n_max = 8
n_min = -6
wing_loading = W / S

mu_g = (2 * wing_loading) / density * MAC * a * 9.80065
k_g = (0.88 * mu_g)/(5.3 + mu_g)
U_de_cruise = 0.3048 * 50
U_de_dive = 0.3048 * 25

V = np.linspace(0, 160 * 0.51444, 1000)
V_cruise = 150 * 0.514444
V_dive = 1.25 * V_cruise
q = 0.5 * V ** 2 * S



n_pos = (q * CNmax) / wing_loading
n_neg = (q * -CNmax) / wing_loading

n_gust_cruise_pos = 1 + ((k_g * density * U_de_cruise * V_cruise * a)/(2*wing_loading))
n_gust_cruise_neg = 1 - ((k_g * density * U_de_cruise * V_cruise * a)/(2*wing_loading))

n_gust_dive_pos = 1 + ((k_g * density * U_de_dive * V_dive * a)/(2*wing_loading))
n_gust_dive_neg = 1 - ((k_g * density * U_de_dive * V_dive * a)/(2*wing_loading))

for i in range(len(V)):
    if mt.isclose(n_pos[i], 1, rel_tol=0.01):
        V_stall = V[i]
    if mt.isclose(n_neg[i], -6, rel_tol=0.01):
        V_prime = V[i]



V_manu = V_stall * np.sqrt(n_max)

plt.xlim(0, 200 * 0.51444)
plt.ylim(-7, 9)
plt.axhline(0, color='black')

plt.axvline(0, color='black')




plt.vlines(V_stall, ymax=1, ymin=0, color="black", label="stall speed", linestyle='--')
plt.vlines(V_prime, ymax=0, ymin=n_min, color="black", label="stall speed", linestyle='--')
plt.vlines(V_cruise, ymax=n_max, ymin=n_min, color="black", label="cruise speed", linestyle='--')
plt.vlines(V_dive, ymax=n_max, ymin=n_min, color="purple", label="dive speed")

plt.hlines(n_max, xmin=V_manu, xmax=V_dive, color='purple')
plt.hlines(n_min, xmin=V_prime, xmax=V_dive, color='purple')



plt.plot(V, n_pos, color='purple')
plt.plot(V, n_neg, color='purple')

print(n_gust_cruise_pos)
plt.plot([0, V_cruise],[1, n_gust_cruise_pos], color='orange')
plt.plot([0, V_cruise],[1, n_gust_cruise_neg], color='orange')
plt.plot([0, V_dive],[1, n_gust_dive_pos], color ='orange', linestyle='--')
plt.plot([0, V_dive],[1, n_gust_dive_neg], color='orange', linestyle='--')

plt.plot([V_cruise, V_dive], [n_gust_cruise_pos, n_gust_dive_pos], color='orange')
plt.plot([V_cruise, V_dive], [n_gust_cruise_neg, n_gust_dive_neg], color='orange')
plt.plot([V_dive, V_dive], [n_gust_dive_pos, n_gust_dive_neg], color='orange', linestyle='--')

plt.grid(b=True, which='major', color='#666666', linestyle='-')

# Show the minor grid lines with very faint and almost transparent grey lines
plt.minorticks_on()
plt.grid(b=True, which='minor', color='#999999', linestyle='-', alpha=0.2)


plt.show()
#%%
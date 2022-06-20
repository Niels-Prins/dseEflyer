"""
Wessel Kruidenier
Thijs van Lith
09/06/2022
Noise calculation and noise reduction due to ducted fan
"""
import numpy as np
import matplotlib.pyplot as plt

# Constants
D = 0.63
B = 13
n_p = 795
P_br = 260000
R = 287
gamma = 1.4
N = 1
V = 28.29
noise_perception = 14
duct_effect = 19.27
r_ft = 3
r_m = 0.3048 * r_ft

#ISA definition
def ISA(h):
    T0 = 288.15
    p0 = 101325.0
    g0 = 9.80665
    R = 287.0
    a1 = -0.0065
    T = T0 + a1 * (h - 0)
    p = p0 * (T / T0) ** (-g0 / (R * a1))
    rho1 = (p / (R * T))
    return [T, p, rho1]

#Calculation of the total perceived noise
c = np.sqrt(R * gamma * ISA(0)[0])
M_t = np.pi * n_p * D / (60 * c)
SPL = 83.4 + 15.3 * np.log10(P_br) - 20 * np.log10(D) \
      + 38.5 * M_t - 3 * (B - 2) + 10 * np.log10(N) - 20 * np.log10(r_m)
print("Noise: ", round(SPL - noise_perception - duct_effect, 2), "dB")

# Noise reduction calculation
fig = plt.figure()
ax = fig.add_subplot()
lst_x = [-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1]
lst_y = [10.8615, 10.3927, 12.4115, 14.1085, 17.5014, 19.2706, 20.5347, 21.0397, 19.7745, 18.5103, 16.8703]
plt.plot(lst_x, lst_y)
plt.title("Reduction in dB with a duct length of 2.22D")
plt.xlabel("Fore- and aft-position in terms of D")
plt.ylabel("Reduction in dB")
plt.vlines(0, 10, lst_y[5], colors="red", linestyles="dotted")
plt.hlines(19.27, -1, 0, colors="red", linestyles="dotted")
ax.text(-0.4, 10.25, "Upstream <-", fontsize=10)
ax.text(0.03, 10.25, "-> downstream", fontsize=10)
ax.text(-1.18, 19.1, "19.27", fontsize=10, color="red")
plt.xlim(-1, 1)
plt.ylim(10, 22)
plt.grid()
plt.show()

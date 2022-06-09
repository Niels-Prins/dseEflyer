"""
Wessel Kruidenier
09/06/2022
Propeller noise calculation with effect of duct
"""

import numpy as np
import pandas as pd

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
r_m = 0.3048*r_ft

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


c = np.sqrt(R * gamma * ISA(0)[0])
M_t = np.pi * n_p * D / (60 * c)
SPL = 83.4 + 15.3 * np.log10(P_br) - 20 * np.log10(D) \
      + 38.5 * M_t - 3 * (B - 2) + 10 * np.log10(N) - 20 * np.log10(r_m)
print("Noise: ", round(SPL - noise_perception - duct_effect, 2), "dB")

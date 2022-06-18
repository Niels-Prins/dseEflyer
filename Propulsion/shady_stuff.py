import numpy as np
import math as mt
import matplotlib.pyplot as plt

density = 1.225

S = 13.1
W = 1035 * 9.86005
CLmax = 1.6
n_max = 8
n_min = -6

wing_loading = W / S
V = np.linspace(0, 160, 1000) * 0.51444

V_cruise = 150 * 0.514444
V_dive = 1.25 * V_cruise
q = 0.5 * V ** 2 * S
n_pos = (q * CLmax) / wing_loading
n_neg = (q * -CLmax) / wing_loading
print(np.where(mt.isclose(n_pos, 1, rel_tol=0.2)))


plt.axhline(0, color='black')
plt.axhline(1, color='green', linestyle='--')
plt.axvline(0, color='black')
plt.xlim(0, 160)
plt.ylim(-7, 9)
plt.axhline(n_max, color='purple')
plt.axhline(n_min, color='cyan')
plt.plot(V, n_pos)
plt.plot(V, n_neg)
plt.grid(True)
plt.show()



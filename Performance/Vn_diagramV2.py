#%%
import matplotlib.pyplot as plt
import numpy as np  

metersps = 0.51444 # From kts to m/s
kts = 1/0.5144444# From m/2 to kts
feet = 1/ 0.3048 # from meter to feet 
lbs = 2.20462 # from kg to lbs

V_stall = 50 * metersps # m/s 

Clmax = 1.6
rho = 0.6601
S = 12.3 # m^2
mass1 = 888
mass2 = 978 
weight1 = mass1 * 9.81
weight2 = mass2 * 9.81

V_cruise = 150 * metersps
V_a = (8*weight1 /(0.5*rho*S *Clmax))**0.5

speeds = np.linspace(0, V_a, 1000)
gloads = 0.5 * rho * 
print(V_a*kts)

print(V_cruise*metersps)



# %%

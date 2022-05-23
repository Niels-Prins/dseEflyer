#%%
# Imports
import numpy as np
import matplotlib.pyplot as plt

"""Model using DAPCA IV method to estimate the Development and Procurement Cost of Aircraft. Taken from snorri"""

# Inputs
lbs = 2.20462 # [kg]

W_airframe = 1035 * lbs  # [lbs] Weight of structure
V_h = 150  # [kts]max airspeed
N = 20  # [-] Number of planned aircraft to be produced in 5 years
F_cert1 = 1  # Certification in CS-23
F_cf1 = 1  # 1.03 for a complex flap system,4 ¼1 if a simple flap system
composites = 0  # Fraction of composites of aircraft
F_comp1 = 1 + composites  #
F_press1 = 1  # Unpressuried aircraft


Q_m = 2  # production rate per month
F_cf2 = 1  # 1.02 for complex flap, 1 for simple flap
F_comp2 = 1 + composites
F_press2 = 1  # For unpressurised cabin
F_taper = 1  # tapered wing

F_cert3 = 1 # Certification
F_cf3 = 1 # Simple flap system
F_comp3 = 1 + 0.25*composites


#%%
# Formulas
# Engineering workhours
H_engr = (
    0.0396
    * (W_airframe ** 0.791)
    * (V_h ** 1.526)
    * (N ** 0.183)
    * F_cert1
    * F_comp1
    * F_press1
)

# Tooling workhours
H_tool = (
    1.0032
    * (W_airframe ** 0.764)
    * (V_h ** 0.899)
    * (N ** 0.178)
    * (Q_m ** 0.066)
    * F_cf2
    * F_comp1
    * F_press2
    * F_taper
)

# Manufacturing labour hours
H_mfg = 9.6613 * (W_airframe**0.74) * (V_h**0.543) * (N**0.524) * F_cert3 * F_cf3 * F_comp3



# %%

#%%
# Imports
import numpy as np
import matplotlib.pyplot as plt

"""Model using DAPCA IV method to estimate the Development and Procurement Cost of Aircraft. Taken from snorri"""
# EASTLAKE MODEL
# Inputs
lbs = 2.20462  # [kg]

W_airframe = 1035 * lbs  # [lbs] Weight of structure
V_h = 150  # [kts]max airspeed
N = 100  # [-] Number of planned aircraft to be produced in 5 years
F_cf1 = 1  # 1.03 for a complex flap system,4 ¼1 if a simple flap system
composites = 0.5  # Fraction of composites of aircraft
F_comp1 = 1 + composites  #

Q_m = 2  # production rate per month
F_taper = 1  # tapered wing

F_comp3 = 1 + 0.25 * composites

N_p = 1  # number of prototypes

F_cert6 = 5  # CS-23 certification

# hourly rates
R_eng = 92  # [$/hour] in 2012
R_tool = 61  # [$/hour] in 2012
R_man = 53  # [$/hour] in 2021

cpi = 1.021 * 1.015 * 1.0001 * 1.013 * 1.021 * 1.024 * 1.018 * 1.012 * 1.047

euro = 0.93

#%%
# Formulas
# Engineering workhours
H_engr = 0.0396 * (W_airframe ** 0.791) * (V_h ** 1.526) * (N ** 0.183) * F_comp1

# Tooling workhours
H_tool = (
    1.0032
    * (W_airframe ** 0.764)
    * (V_h ** 0.899)
    * (N ** 0.178)
    * (Q_m ** 0.066)
    * F_comp1
    * F_taper
)

# Manufacturing labour hours
H_mfg = 9.6613 * (W_airframe ** 0.74) * (V_h ** 0.543) * (N ** 0.524) * F_comp3
# Total cost of engineering
C_engr = H_engr * R_eng * cpi

# Total cost of development
C_dev = 0.06458 * (W_airframe ** 0.873) * (V_h ** 1.89) * (N_p ** 0.346) * cpi * F_comp1

# Total cost of Flight test Operations
C_ft = 0.009646 * (W_airframe ** 1.16) * (V_h ** 1.3718) * (N_p ** 1.281) * F_cert6

# Total cost of tooling
C_tool = H_tool * R_tool * cpi

# Total cost of manufacturing
C_mfg = H_mfg * R_man * cpi
# Total cost of quality control
C_qc = 0.13 * C_mfg * F_comp1
# Total cost of manufacturing
C_mat = 24.896 * (W_airframe ** 0.689) * (V_h ** 0.624) * (N ** 0.792) * cpi

batteries = (103572 / 1000) * 132 * euro  # Price of batteries

C_vsc = 10000 + 13000 + batteries
# Variable cost
C_var = (C_mfg + C_qc + C_mat) / N + C_vsc


# Fixed cost (cost to certify)
C_fix = (C_engr * 0 + C_dev * 0 + C_ft + C_tool) * euro

# Cost per unit
C_unit = (C_fix / N + C_var) * euro

print(
    f"Fixed Cost = €{round(C_fix)} \nVariable cost = €{round(C_var)} \nPrice per unit = €{round(C_unit)}"
)

# %%
# Operational cost

F_mf = 0.3 - 0.15
R_ap = 60
Q_flgt = 250
C_ap = F_mf * R_ap * Q_flgt
C_isnp = 500

C_year = C_ap + C_isnp

print(f"Yeary cost for maintenance and inspection = €{C_year}")
# %%

#%%
# Imports
import numpy as np
import matplotlib.pyplot as plt

"""Model using DAPCA IV method to estimate the Development and Procurement Cost of Aircraft. Taken from snorri"""
# EASTLAKE MODEL
# Inputs
lbs = 2.20462  # [kg]

W_airframe = 876 * lbs  # [lbs] Weight of structure
V_h = 150  # [kts]max airspeed
N = 75  # [-] Number of planned aircraft to be produced in 5 years
F_cf1 = 1  # 1.03 for a complex flap system,4 ¼1 if a simple flap system
composites = 1  # Fraction of composites of aircraft
F_comp1 = 1 + composites  #

Q_m = 5  # production rate per month
F_taper = 1  # tapered wing

F_comp3 = 1 + 0.25 * composites

N_p = 1  # number of prototypes

F_cert6 = 5  # CS-23 certification

# hourly rates
R_eng = 92  # [$/hour] in 2012
R_tool = 61  # [$/hour] in 2012
R_man = 53  # [$/hour] in 2012

cpi = 1.021 * 1.015 * 1.0001 * 1.013 * 1.021 * 1.024 * 1.018 * 1.012 * 1.047*1.08
euro = 0.93

#%%
# Formulas
# Engineering workhours
H_engr = 0.0396 * (W_airframe ** 0.791) * (V_h ** 1.526) * (N ** 0.183) * F_comp1 -4000

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
motor = 10000
avionics = 32500/euro
propeller = 9500
C_vsc = motor + avionics + batteries + propeller
# Variable cost
C_var = (C_mfg + C_qc + C_mat) / N + C_vsc


# Fixed cost (cost to certify)
C_fix = C_engr + C_dev + C_ft + C_tool

# Cost per unit
C_unit = C_fix / N + C_var
C_fixUnit = C_fix / N


print(
    f"Fixed Cost = €{round(C_fix * euro)} \nVariable cost = €{round(C_var * euro)} \nPrice per unit = €{round(C_unit * euro)}\n"
)
print(
    f"Engineering (\euro{round(R_eng*euro*cpi)}/hr) & {round(H_engr)} & {round(C_engr*euro)} & {round( C_engr*euro/N)} \\\ \hline\n Tooling  (\euro{round(R_tool*euro*cpi)}/hr)& {round(H_tool)} & {round(C_tool*euro)} & {round( C_tool*euro/N)}\\\ \hline\n Development Support & - & {round(C_dev*euro)} & {round( C_dev*euro/N)} \\\ \hline\n Flight Tests & - & {round(C_ft*euro)} & {round( C_engr*euro/N)}\\\ \hline \hline\n Total fixed cost & - & {round(C_fix*euro)} & {round(C_fixUnit*euro)} \\\ \hline \hline\n Manufacturing (\euro{round(R_man*euro*cpi)}/hr) & {round(H_mfg)} & {round(C_mfg*euro)} & {round(C_mfg*euro/N)} \\\ \hline\n Quality Control & - & {round(C_qc*euro)} & {round(C_qc*euro/N)} \\\ \hline\n Materials & - & {round(C_mat*euro)} & {round(C_mat*euro/N)}\\\ \hline\n Batteries & - &{round(batteries*euro*N)} &{round(batteries*euro)}\\\ \hline\n Motor & - & {round(motor*euro*N)} & {round(motor*euro)} \\\ \hline\n Avionics & - & {round(avionics * euro*N)} & {round(avionics*euro)} \\\ \hline\n Propeller & - & {round(propeller*euro*N)} & {round(propeller*euro)} \\\ \hline \hline\n Total Variable Cost & - & {round(C_var*euro*N)} & {round(C_var*euro)} \\\ \hline \hline \n Total & - & {round((C_var *N + C_fix)*euro)} & {round(C_unit*euro)} \\\ \hline"
)

# %%
# Operational cost

F_mf = 0.3 - 0.15 + 0.02
R_ap = 67*cpi*euro
Q_flgt = 250 *0.5
C_ap = F_mf * R_ap * Q_flgt
elec = 0.55 #€/kwh
capacity = 59.1

storage = 12*275*cpi

fuel = elec*capacity*0.85*250
inspection = 500*cpi*euro
C_year = C_ap +storage + inspection + fuel

print(f"Yeary cost for maintenance= €{C_ap} = {C_ap/Q_flgt}")
print(f"Inspection = {inspection} = {inspection/Q_flgt}")
print(f"Fuel = {fuel} = {fuel/Q_flgt}")
print(f"Storag = {storage} = {storage/Q_flgt}")
print(f"total = {C_year} = {C_year/Q_flgt}")


# %%

print("COMBUSTION ENGINE")
print(f"Yeary cost for maintenance= €{C_ap +12000 /1500} = {(C_ap)/Q_flgt + 12000 /1500}")
print(f"Inspection = {inspection} = {inspection/Q_flgt}")
print(f"Fuel = {elec*capacity*0.85*250} = {elec*capacity*0.85*250/Q_flgt}")
print(f"Storag = {storage} = {storage/Q_flgt}")
print(f"total = {C_year} = {C_year/Q_flgt}")
# %%

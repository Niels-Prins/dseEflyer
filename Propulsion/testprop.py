import numpy as np
import math as mt
import scipy.interpolate
import scipy as sc
import matplotlib.pyplot as plt
import pandas as pd

h = 0
v = 30

kappa = 1.4
cp = 1000
T0 = 288.15
p0 = 101325.0
g0 = 9.80665
R = 287.0
a1 = -0.0065
T = T0 + a1 * (h - 0)
p = p0 * (T / T0) ** (-g0 / (R * a1))
rho = (p / (R * T))


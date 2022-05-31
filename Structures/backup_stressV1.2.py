# Lennart Krieg
# DSE group 20: Aerobatic E-Flyer
# Wing load model, part of wing structural model

import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import scipy.interpolate
from scipy import integrate
from mpl_toolkits import mplot3d
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm
import matplotlib.cm as cm
import pandas as pd
from DSE_EFlyer_wing_properties import WingProperties
from DSE_EFlyer_wing_loads import WingLoading

#  Pycharm Project folder: C:\Users\lenna\PycharmProjects\DSE_22

DIRECTORY = 'C:/Users/lenna/Documents/DSE/Wing_structural_model'
INPUT_FOLDER = '/Input_Files/'
CL_dist_file = 'CL_dist_apha_clean.csv'
CL_alpha = ''
CD_dist_file = 'CD_dist_alpha_clean.csv'

# Airfoil files with coordinates in x/z fractions of chord.
# Before use these files will have to be preprocessed to have them start at the LE, go to TE and wrap around again
# to the LE. Text and other non-data will have to be removed and separation by spaces replaced by semicolon (;).
airfoil_file = 'NACA-0015_LE_71.dat'

# Wing geometry parameters
c_root = 2.18  # Root chord, [m]
c_tip = 0.98  # Tip chord, [m]
t_red_tip = 0.12/0.15  # Ratio of tip thickness (of % chord) over root thickness (of % chord), [-]
MAC = 1.58  # Mean aerodynamic chord, [m]
S = 14.536  # Total wing area, [m^2]

dihedral = 0  # Wing dihedral, [deg]
LE_sweep = 0  # Leading egde sweep, [deg]
half_span = 9.2/2  # Halfspan of the wing, [m]
AR = 5.8  # Aspect ratio, [-]

CL_0 = 0.05765
CL_10 = 0.82880

# Drag over alpha polynomial constants (can be taken from Excel)
# form of CD = a * alpha^2 + b * alpha + c
a = 0.0003
b = 0.0006
c = 0.0059

nr_span_sections = 10  # Number of cross-sections in which the wing is split up, [-]
ac_mass = 1063  # Aircraft manoeuvring mass, [kg]
g_load = 8  # Aircraft manoeuvring g load,[-]

rho = 1.225  # Air density, [kg/m^3]
V = 150  # Maneuvering speed, [kts]
nr_spars = 2  # Amount of spars, [-]
spar_locs = [0.25, 0.75]  # Location of spars as fraction of the chords, [-]. Input should be a list or an array.

# Thicknesses
t_le = 0  # Leading edge thickness, [mm]
t_top = 6  # Top skin thickness, [mm]
t_tr = 6  # Trailing edge thickness, [mm]
t_bot = 0  # Bottom skin thickness, [mm]
t_spar = [2.5, 2.5]  # Spar thicknesses, [mm]. Input must be list of same size as spar_locs list!
t_sp_flange = [3, 3]
w_sp_flange = [100, 100]

x_le = 0.25  # Fraction of chord where LE thickness becomes top or bottom skin thickness
x_te = 0.75  # Fraction of chord where top or bottom skin thickness become TE thickness

dens_spars = 1600  # Spar material density, [kg/m^3]
dens_skins = 1600  # Skin material density, [kg/m^3]


class WingStresses:
    def __init__(self, DIRECTORY, INPUT_FOLDER, airfoil_file, c_root, c_tip, MAC, dihedral, LE_sweep, half_span, AR,
                 nr_span_sections, nr_spars, spar_locs, t_le, t_top, t_tr, t_bot, t_spar, x_le, x_te, t_sp_flange,
                 w_sp_flange, t_red_tip, CL_dist_file, CD_dist_file, ac_mass, g_load, CL_0, CL_10,
                 S, rho, V, a, b, c, dens_spars, dens_skins):

        self.GeoProps = WingProperties(DIRECTORY, INPUT_FOLDER, airfoil_file, c_root, c_tip, MAC, dihedral,
                                       LE_sweep, half_span, AR, nr_span_sections, nr_spars, spar_locs, t_le,
                                       t_top, t_tr, t_bot, t_spar, x_le, x_te, t_sp_flange, w_sp_flange, t_red_tip,
                                       dens_spars, dens_skins)
        print(self.GeoProps.spar_coords)
        print('---')
        print(self.GeoProps.centroids)
        self.WingLoads = WingLoading(DIRECTORY, INPUT_FOLDER, CL_dist_file, CD_dist_file, c_root, c_tip, MAC, dihedral,
                                     LE_sweep, half_span, AR, nr_span_sections, ac_mass, g_load, CL_0, CL_10, S, rho, V,
                                     a, b, c)

        self.calc_max_min_bending_stresses()

    def bending_stress(self, y, x, z):
        siqma_y = (self.WingLoads.moment_dist_z(y)/self.GeoProps.Izz_dist(y)) * (x-self.GeoProps.centroid_dist(y)) + \
                  (self.WingLoads.moment_dist_x(y)/self.GeoProps.Ixx_dist(y)) * z
        return siqma_y

    def calc_max_min_bending_stresses(self):
        stresses_min = []
        stresses_max = []
        for span_idx, sec in enumerate(self.GeoProps.w_coords):
            span = self.GeoProps.sec_spans[span_idx]
            x_coords = sec[0]
            z_coords = sec[1]
            sec_stresses = []
            for idx_x, x in enumerate(x_coords):
                z = z_coords[idx_x]
                sigma = self.bending_stress(span, x, z)
                sec_stresses.append(sigma)
            sigma_min = min(sec_stresses)
            sigma_max = max(sec_stresses)
            stresses_min.append(sigma_min)
            stresses_max.append(sigma_max)

            norm = plt.Normalize(sigma_min, sigma_max)
            fig, ax = plt.subplots()
            plt.scatter(x_coords, z_coords, c=sec_stresses, cmap=plt.get_cmap('jet'))
            plt.colorbar(label='Bending stress [Pa]')
            ax.set_aspect('equal')
            plt.show()

        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_size_inches(13, 6)
        ax1.set_title('Compressive stress over wing halfspan.')
        ax1.set_xlabel('Wing halfspan, [m]')
        ax1.set_ylabel('Bending stress, [Nm]')
        ax2.set_title('Tensile stress over wing halfspan.')
        ax2.set_xlabel('Wing halfspan, [m]')
        ax2.set_ylabel('Bending stress, [Nm]')
        ax1.plot(self.GeoProps.sec_spans, stresses_min)
        ax2.plot(self.GeoProps.sec_spans, stresses_max)
        plt.tight_layout()
        plt.show()


WS = WingStresses(DIRECTORY, INPUT_FOLDER, airfoil_file, c_root, c_tip, MAC, dihedral, LE_sweep, half_span, AR,
                  nr_span_sections, nr_spars, spar_locs, t_le, t_top, t_tr, t_bot, t_spar, x_le, x_te, t_sp_flange,
                  w_sp_flange, t_red_tip, CL_dist_file, CD_dist_file, ac_mass, g_load, CL_0, CL_10,
                  S, rho, V, a, b, c, dens_spars, dens_skins)


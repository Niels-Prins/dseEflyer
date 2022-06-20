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
from DSE_EFlyer_CLT import LaminateTheory

#  Pycharm Project folder: C:\Users\lenna\PycharmProjects\DSE_22

DIRECTORY = 'C:/Users/lenna/Documents/DSE/Wing_structural_model'
INPUT_FOLDER = '/Input_Files/'
CL_dist_file = 'CL_dist_apha_clean_V2.csv'
CL_alpha = ''
CD_dist_file = 'CD_dist_alpha_clean_V2.csv'
CM_dist_file = 'CM_dist_alpha_clean.csv'
COP_dist_file = 'COP_dist_alpha_clean.csv'

# Airfoil files with coordinates in x/z fractions of chord.
# Before use these files will have to be preprocessed to have them start at the LE, go to TE and wrap around again
# to the LE. Text and other non-data will have to be removed and separation by spaces replaced by semicolon (;).
airfoil_file = 'NACA-0015_LE_71.dat'

# Wing geometry parameters
c_root = 2.01  # Root chord, [m]
c_tip = 0.90  # Tip chord, [m]
t_red_tip = 0.12/0.15  # Ratio of tip thickness (of % chord) over root thickness (of % chord), [-]
MAC = 1.46  # Mean aerodynamic chord, [m]
S = 12.3  # Total wing area, [m^2]

dihedral = 0  # Wing dihedral, [deg]
LE_sweep = 0  # Leading egde sweep, [deg]
half_span = 8.48/2  # Halfspan of the wing, [m]
AR = 5.8  # Aspect ratio, [-]

CL_0 = 0.05765
CL_10 = 0.82880

# Drag over alpha polynomial constants (can be taken from Excel)
# form of CD = a * alpha^2 + b * alpha + c
a = 0.0003
b = 0.0006
c = 0.0059

nr_span_sections = 100  # Number of cross-sections in which the wing is split up, [-]
ac_mass = 876  # Aircraft manoeuvring mass, [kg]
g_load = 8  # Aircraft manoeuvring g load,[-]

rho = 1.225  # Air density, [kg/m^3]
V = 155  # Maneuvering speed, [kts]
nr_spars = 2  # Amount of spars, [-]
spar_locs = [0.3, 0.75]  # Location of spars as fraction of the chords, [-]. Input should be a list or an array.
q_sec61_loc = 0.1  # Location of the LE vertical element in simplified shear calculations as fraction of chord, [-]
q_sec52_loc = spar_locs[0]
q_sec34_loc = spar_locs[1]

loc_ribs = [0, 0.3, 0.65, 1.1, 1.65, 2.7, half_span]  # Locations/span of the ribs, [m]
t_ribs = [1.92, 1.92, 1.92, 2.56, 1.92, 1.92, 1.92]  # Thickness of the ribs, [mm]

# Thicknesses
t_le = [1.92, 1.76, 1.76, 1.6, 1.28, 1.28]  # Leading edge thickness, [mm]
t_top = [2.24, 2.24, 1.92, 1.92, 1.6, 1.6]  # Top skin thickness, [mm]
t_tr = [0, 0, 0, 0, 0, 0]  # Trailing edge thickness, [mm]
t_bot = [2.24, 2.24, 1.92, 1.92, 1.6, 1.6]  # Bottom skin thickness, [mm]

t_fsp_web = [2.08, 2.08, 2.08, 2.08, 1.92, 1.44]
t_msp_web = [2.72, 2.72, 2.72, 2.56, 2.56, 1.92]
t_rsp_web = [2.56, 2.56, 2.4, 2.24, 1.92, 1.44]

t_spar = [3, 3]  # Spar thicknesses, [mm]. Input must be list of same size as spar_locs list!
t_sp_flange = [0, 0]
w_sp_flange = [0, 0]

x_le = 0.3  # Fraction of chord where LE thickness becomes top or bottom skin thickness
x_te = 0.75  # Fraction of chord where top or bottom skin thickness become TE thickness

# Number of stringers including the 'flanges' of the spars, so 4 stingers is 2 stringers on the skin and 1 at each spar
nr_str_top = 12
nr_str_bot = 12
nr_str_LE_top = 7
nr_str_LE_bot = 7

# The following lists contain the spans at which each stringer must end. Number of elements must equal number of
# stringers entered above.
str_top_end = [half_span, half_span, 2.7, half_span, 2.7, half_span, 2.7, half_span, 2.7, half_span, 2.7, half_span]
str_bot_end = [half_span, half_span, 2.7, half_span, 2.7, half_span, 2.7, half_span, 2.7, half_span, 2.7, half_span]
str_LE_top_end = [half_span, 2.7, half_span, 2.7, half_span, 2.7, half_span]
str_LE_bot_end = [half_span, 2.7, half_span, 2.7, half_span, 2.7, half_span]

# If set to true, the centroid, Ixx and Izz are calculated based on the simplified two-cell geometry. Value will be
# lower, but method for shear flow calculations will be consistent with the geometry. Only this method works with
# stringers at the moment
SMOA_simplified = True

dens_spars = 1480  # Spar material density, [kg/m^3]
dens_skins = 1480  # Skin material density, [kg/m^3]

safety_margin = 1.5


class WingStresses:
    def __init__(self, DIRECTORY, INPUT_FOLDER, airfoil_file, c_root, c_tip, MAC, dihedral, LE_sweep, half_span, AR,
                 nr_span_sections, nr_spars, spar_locs, t_le, t_top, t_tr, t_bot, t_spar, x_le, x_te, t_sp_flange,
                 w_sp_flange, t_red_tip, CL_dist_file, CD_dist_file, ac_mass, g_load, CL_0, CL_10,
                 S, rho, V, a, b, c, dens_spars, dens_skins, loc_ribs, q_sec61_loc, q_sec52_loc, q_sec34_loc, SFmargin,
                 t_msp_web, t_rsp_web, SMOA_simplified, nr_str_top, nr_str_bot, nr_str_LE_top, nr_str_LE_bot,
                 str_top_end, str_bot_end, str_LE_top_end, str_LE_bot_end, t_ribs, CM_dist_file, COP_dist_file, t_fsp_web):

        self.GeoProps = WingProperties(DIRECTORY, INPUT_FOLDER, airfoil_file, c_root, c_tip, MAC, dihedral,
                                       LE_sweep, half_span, AR, nr_span_sections, nr_spars, spar_locs, t_le,
                                       t_top, t_tr, t_bot, t_spar, x_le, x_te, t_sp_flange, w_sp_flange, t_red_tip,
                                       dens_spars, dens_skins, loc_ribs, q_sec61_loc, q_sec52_loc, q_sec34_loc,
                                       t_msp_web, t_rsp_web, SMOA_simplified, nr_str_top, nr_str_bot, nr_str_LE_top,
                                       nr_str_LE_bot, str_top_end, str_bot_end, str_LE_top_end, str_LE_bot_end, t_ribs, t_fsp_web)

        self.WingLoads = WingLoading(DIRECTORY, INPUT_FOLDER, CL_dist_file, CD_dist_file, c_root, c_tip, MAC, dihedral,
                                     LE_sweep, half_span, AR, nr_span_sections, ac_mass, g_load, CL_0, CL_10, S, rho, V,
                                     a, b, c, CM_dist_file, COP_dist_file)
        self.CLT = LaminateTheory()

        self.req_SF = SFmargin
        self.torques_y = []
        self.dvdy_dist = None
        self.v_z_dist = None
        self.thau_max_msp = None
        self.thau_max_rsp = None
        self.thau_max_fsp = None
        self.thau_crit_msp = None
        self.thau_crit_rsp = None
        self.thau_crit_fsp = None
        self.thau_max_msp_dist = None
        self.thau_max_rsp_dist = None
        self.thau_max_fsp_dist = None
        self.thau_crit_msp_dist = None
        self.thau_crit_rsp_dist = None
        self.thau_crit_fsp_dist = None

        self.str_top_stress = None
        self.top_sigma_crit = None
        self.str_LE_top_stress = None
        self.LE_top_sigma_crit = None

        # Shear web buckling constant
        asp_sb = [0, 1, 1.65, 2, 3, 3.75, 5, 50, 100]
        ks = [10, 10, 7, 6.2, 5.8, 5.8, 5.475, 5, 5]

        self.Ks_func = sc.interpolate.interp1d(asp_sb, ks, kind='linear')
        # Skin compressive buckling constant (simply supported edges)
        asp_cb = [0, 0.4, 0.6, 1, 1.4, 1.8, 100]
        kc = [8.4, 8.4, 5, 4.1, 4.5, 4.1, 4.1]
        self.Kc_func = sc.interpolate.interp1d(asp_cb, kc, kind='linear')

        self.est_wing_mass()
        self.calc_torsion_y()
        self.calc_shearstress_web()
        self.spar_web_shear_buckling()
        # self.calc_max_min_bending_stresses()
        self.plot_bending_stress(0)
        self.skin_compressive_stress()
        self.skin_compr_buckling()
        self.eval_dv_dy()
        self.eval_v_z()

    def bending_stress(self, y, x, z):
        siqma_y = (self.WingLoads.moment_dist_z(y)/self.GeoProps.Izz_dist(y)) * (x-self.GeoProps.centroid_dist(y)) + \
                  (self.WingLoads.moment_dist_x(y)/self.GeoProps.Ixx_dist(y)) * z
        return siqma_y

    def calc_torsion_y(self):
        tor = []
        for i, y in enumerate(self.GeoProps.sec_spans):
            sc = self.GeoProps.shear_centers[i]
            cp = (self.WingLoads.COP_dist(y)/100) * self.GeoProps.sec_chords[i]
            arm = cp - sc
            torque = self.WingLoads.shear_dist_z[i] * arm + self.WingLoads.aero_moment_dist[i]
            tor.append(torque)

        self.torques_y = np.array(tor)
        plt.plot(self.GeoProps.sec_spans, self.torques_y, c='orange')
        plt.xlabel('Wing halfspan, [m]')
        plt.ylabel('Torque, [Nm]')
        plt.title('Torque distribution over wing span')
        plt.minorticks_on()
        plt.grid(True, which='minor', color='#999999', linestyle='-', alpha=0.2)
        plt.grid(True, which='major', color='#666666', linestyle='-')
        plt.show()

    def calc_shearstress_web(self):
        # Main spar
        msp_shear_max_dist = []
        for i in range(self.GeoProps.nr_bsec):
            eval_z = np.linspace(self.GeoProps.xz5[i][1], self.GeoProps.xz2[i][1], num=21, endpoint=True)
            q_msp = []
            for z_idx, z in enumerate(eval_z):
                if z <= 0:
                    q_Sz = self.GeoProps.q_58fz(i, z) * self.WingLoads.shear_dist_z[i]
                    q_Sx = self.GeoProps.q_58fx(i, z) * self.WingLoads.shear_dist_x[i]
                    q_T_1 = self.GeoProps.tor_q_0[i][1] * self.torques_y[i]
                    q_T_0 = self.GeoProps.tor_q_0[i][0] * self.torques_y[i]
                    q_T = q_T_1 - q_T_0
                    # print('q_T_1', q_T_1)
                    # print('q_T_0', q_T_0)
                    q = q_Sz + q_Sx + q_T
                else:
                    q_Sz = self.GeoProps.q_82fz(i, z) * self.WingLoads.shear_dist_z[i]
                    q_Sx = self.GeoProps.q_82fx(i, z) * self.WingLoads.shear_dist_x[i]
                    q_T_1 = self.GeoProps.tor_q_0[i][1] * self.torques_y[i]
                    q_T_0 = self.GeoProps.tor_q_0[i][0] * self.torques_y[i]
                    q_T = q_T_1 - q_T_0
                    # print('q_T_1', q_T_1)
                    # print('q_T_0', q_T_0)
                    q = q_Sz + q_Sx + q_T
                q_msp.append(q)
            # q_msp = np.array(q_msp)
            q_max = np.abs(q_msp).max()
            thau_max = q_max/self.GeoProps.t58[i]
            # print('MSP max shear flow, [N/m]:', q_max)
            # print('MSP max shear stress, [MPa]:', thau_max/1000000)
            msp_shear_max_dist.append(thau_max)
        self.thau_max_msp = np.array(msp_shear_max_dist)
        self.thau_max_msp_dist = sc.interpolate.interp1d(self.GeoProps.sec_spans, self.thau_max_msp, kind='cubic')

        # Rear spar
        rsp_shear_max_dist = []
        for i in range(self.GeoProps.nr_bsec):
            eval_s = np.linspace(0, self.GeoProps.l34[i], num=21, endpoint=True)
            q_rsp = []
            for s in eval_s:
                q_Sz = self.GeoProps.q_34fz(i, s) * self.WingLoads.shear_dist_z[i]
                q_Sx = self.GeoProps.q_34fx(i, s) * self.WingLoads.shear_dist_x[i]
                q_T = self.GeoProps.tor_q_0[i][1] * self.torques_y[i]
                q_r = q_Sz + q_Sx + q_T
                q_rsp.append(q_r)
            # q_msp = np.array(q_msp)
            q_max_r = np.abs(q_rsp).max()
            thau_max_r = (q_max_r / self.GeoProps.t34[i])
            # print('RSP max shear flow, [N/m]:', q_max_r)
            # print('RSP max shear stress, [MPa]:', thau_max_r / 1000000)
            rsp_shear_max_dist.append(thau_max_r)
        self.thau_max_rsp = np.array(rsp_shear_max_dist)
        self.thau_max_rsp_dist = sc.interpolate.interp1d(self.GeoProps.sec_spans, self.thau_max_rsp, kind='cubic')

        # Front spar
        fsp_shear_max_dist = []
        for i in range(self.GeoProps.nr_bsec):
            eval_z = np.linspace(self.GeoProps.xz6[i][1], self.GeoProps.xz1[i][1], num=21, endpoint=True)
            q_fsp = []
            for z_idx, z in enumerate(eval_z):
                if z <= 0:
                    q_Sz = self.GeoProps.q_67fz(i, z) * self.WingLoads.shear_dist_z[i]
                    q_Sx = self.GeoProps.q_67fx(i, z) * self.WingLoads.shear_dist_x[i]
                    q_T_0 = self.GeoProps.tor_q_0[i][0] * self.torques_y[i]
                    q_T = q_T_0
                    # print('q_T_1', q_T_1)
                    # print('q_T_0', q_T_0)
                    q = q_Sz + q_Sx + q_T
                else:
                    q_Sz = self.GeoProps.q_71fz(i, z) * self.WingLoads.shear_dist_z[i]
                    q_Sx = self.GeoProps.q_71fx(i, z) * self.WingLoads.shear_dist_x[i]
                    q_T_0 = self.GeoProps.tor_q_0[i][0] * self.torques_y[i]
                    q_T = q_T_0
                    # print('q_T_1', q_T_1)
                    # print('q_T_0', q_T_0)
                    q = q_Sz + q_Sx + q_T
                q_fsp.append(q)
            # q_msp = np.array(q_msp)
            q_max = np.abs(q_fsp).max()
            thau_max = q_max / self.GeoProps.t67[i]
            # print('MSP max shear flow, [N/m]:', q_max)
            # print('MSP max shear stress, [MPa]:', thau_max/1000000)
            fsp_shear_max_dist.append(thau_max)
        self.thau_max_fsp = np.array(fsp_shear_max_dist)
        self.thau_max_fsp_dist = sc.interpolate.interp1d(self.GeoProps.sec_spans, self.thau_max_fsp, kind='cubic')

    def spar_web_shear_buckling(self):
        cap = 8
        # The ribs devide the wing spars and skins into smaller 'segments'
        # The critical buckling stresses within each of these segments are determined here
        thau_crit_msp = []
        thau_crit_rsp = []
        thau_crit_fsp = []
        for idx_spar, y_s in enumerate(self.GeoProps.ribs_loc[:-1]):
            y_e = self.GeoProps.ribs_loc[idx_spar + 1]
            dy = y_e - y_s
            rib_span_idxs = np.where(np.logical_and(self.GeoProps.sec_spans >= y_s, self.GeoProps.sec_spans <= y_e))
            # print(rib_span_idxs[0])
            # t_msp_s = self.GeoProps.xz2[rib_span_idxs[0]][1] * 2
            # aspect = dy/t_msp_s
            # Ks = self.Ks_func(aspect)
            # print('t_msp', t_msp_s)

            for sec_idx in rib_span_idxs[0]:
                # Main spar
                h_msp = self.GeoProps.xz2[sec_idx][1] * 2
                aspect_msp = dy / h_msp
                Ks_msp = self.Ks_func(aspect_msp)
                # print('Local aspect ratio:', aspect)
                # print('Local height:', h_msp_sec)
                # print('Local dy', dy)
                # print('Local Ks:', Ks)
                t_msp = self.GeoProps.t82[sec_idx]
                thau_cr_msp = ((np.pi**2 * Ks_msp * self.CLT.Exx)/(12 * (1 - self.CLT.vxy**2)))*(t_msp/h_msp)**2
                thau_crit_msp.append(thau_cr_msp)

                # Rear spar
                h_fsp = self.GeoProps.xz1[sec_idx][1] * 2
                aspect_fsp = dy / h_fsp
                Ks_fsp = self.Ks_func(aspect_fsp)
                t_fsp = self.GeoProps.t71[sec_idx]
                thau_cr_fsp = ((np.pi ** 2 * Ks_fsp * self.CLT.Exx) / (12 * (1 - self.CLT.vxy ** 2))) * (t_fsp / h_fsp) ** 2
                thau_crit_fsp.append(thau_cr_fsp)

                # Rear spar
                h_rsp = self.GeoProps.xz3[sec_idx][1] * 2
                aspect_rsp = dy / h_rsp
                Ks_rsp = self.Ks_func(aspect_rsp)
                # print('Local aspect ratio:', aspect)
                # print('Local height:', h_msp_sec)
                # print('Local dy', dy)
                # print('Local Ks:', Ks)
                t_rsp = self.GeoProps.t34[sec_idx]
                thau_cr_rsp = ((np.pi ** 2 * Ks_rsp * self.CLT.Exx) / (12 * (1 - self.CLT.vxy ** 2))) * (
                            t_rsp / h_rsp) ** 2
                thau_crit_rsp.append(thau_cr_rsp)

        self.thau_crit_msp = np.array(thau_crit_msp)
        self.thau_crit_rsp = np.array(thau_crit_rsp)
        self.thau_crit_fsp = np.array(thau_crit_fsp)

        # Main spar
        # plt.plot(self.GeoProps.sec_spans, self.thau_crit_msp / 1000000, c='red')
        # plt.plot(self.GeoProps.sec_spans, self.thau_max_msp / 1000000, c='blue')
        # for rib_y in self.GeoProps.ribs_loc:
        #     plt.axvline(x=rib_y, c='black', linestyle='--')
        # plt.xlabel('Wing halfspan, [m]')
        # plt.ylabel('MSP web shear stress, [MPa]')
        # plt.title('thau_max and thau_crit distribution over wing span')
        # plt.show()

        MSP_SF_capped = [max(min(thau, cap), 0) for thau in self.thau_crit_msp / self.thau_max_msp]

        # plt.plot(self.GeoProps.sec_spans, MSP_SF_capped, c='red')
        # plt.plot(self.GeoProps.sec_spans, np.full(len(self.GeoProps.sec_spans), self.req_SF), c='blue')
        # for rib_y in self.GeoProps.ribs_loc:
        #     plt.axvline(x=rib_y, c='black', linestyle='--')
        # plt.xlabel('Wing halfspan, [m]')
        # plt.ylabel('MSP web shear buckling SF, [-]')
        # plt.title('Safety margin distribution')
        # plt.show()

        # Rear spar
        # plt.plot(self.GeoProps.sec_spans, self.thau_crit_rsp / 1000000, c='red')
        # plt.plot(self.GeoProps.sec_spans, self.thau_max_rsp / 1000000, c='blue')
        # for rib_y in self.GeoProps.ribs_loc:
        #     plt.axvline(x=rib_y, c='black', linestyle='--')
        # plt.xlabel('Wing halfspan, [m]')
        # plt.ylabel('RSP web shear stress, [MPa]')
        # plt.title('thau_max and thau_crit distribution over wing span')
        # plt.show()

        RSP_SF_capped = [max(min(thau, cap), 0) for thau in self.thau_crit_rsp / self.thau_max_rsp]

        # plt.plot(self.GeoProps.sec_spans, RSP_SF_capped, c='red')
        # plt.plot(self.GeoProps.sec_spans, np.full(len(self.GeoProps.sec_spans), self.req_SF), c='blue')
        # for rib_y in self.GeoProps.ribs_loc:
        #     plt.axvline(x=rib_y, c='black', linestyle='--')
        # plt.xlabel('Wing halfspan, [m]')
        # plt.ylabel('RSP web shear buckling SF, [-]')
        # plt.title('Safety margin distribution')
        # plt.show()

        # Front spar
        # plt.plot(self.GeoProps.sec_spans, self.thau_crit_fsp / 1000000, c='red')
        # plt.plot(self.GeoProps.sec_spans, self.thau_max_fsp / 1000000, c='blue')
        # for rib_y in self.GeoProps.ribs_loc:
        #     plt.axvline(x=rib_y, c='black', linestyle='--')
        # plt.xlabel('Wing halfspan, [m]')
        # plt.ylabel('FSP web shear stress, [MPa]')
        # plt.title('thau_max and thau_crit distribution over wing span')
        # plt.show()

        FSP_SF_capped = [max(min(thau, cap), 0) for thau in self.thau_crit_fsp / self.thau_max_fsp]

        plt.figure(figsize=(8,6))
        plt.plot(self.GeoProps.sec_spans, MSP_SF_capped, c='red', label='MSP')
        plt.plot(self.GeoProps.sec_spans, RSP_SF_capped, c='purple', label='RSP')
        plt.plot(self.GeoProps.sec_spans, FSP_SF_capped, c='orange', label='FSP')
        plt.plot(self.GeoProps.sec_spans, np.full(len(self.GeoProps.sec_spans), self.req_SF), c='blue')
        for rib_y in self.GeoProps.ribs_loc:
            plt.axvline(x=rib_y, c='black', linestyle='--')
        plt.xlabel('Wing halfspan, [m]')
        plt.ylabel('SF, [-]')
        plt.title('Web buckling safety margin distribution')
        plt.minorticks_on()
        plt.grid(True, which='minor', color='#999999', linestyle='-', alpha=0.2)
        plt.grid(True, which='major', color='#666666', linestyle='-')
        plt.legend()
        plt.show()

    def skin_compressive_stress(self):
        str_top_stress = []
        for y_idx, y_sec in enumerate(self.GeoProps.sec_spans):
            yeet = []
            for cor_idx, cor in enumerate(self.GeoProps.str_top_coor[y_idx][:-1]):
                stress_buckpan = []
                cor_next = self.GeoProps.str_top_coor[y_idx][cor_idx + 1]
                eval_x = np.linspace(cor[0], cor_next[0], num=10, endpoint=True)
                eval_y = np.linspace(cor[1], cor_next[1], num=10, endpoint=True)
                for x_idx, x in enumerate(eval_x):
                    # print(x)
                    sigma = self.bending_stress(y_sec, x, eval_y[x_idx])
                    # print(sigma)
                    stress_buckpan.append(sigma)
                max_compr = min(stress_buckpan)
                yeet.append(max_compr)
            str_top_stress.append(yeet)

        self.str_top_stress = str_top_stress

        str_LE_top_stress = []
        for y_idx, y_sec in enumerate(self.GeoProps.sec_spans):
            LE_yeet = []
            for cor_idx, cor in enumerate(self.GeoProps.str_LE_top_coor[y_idx][:-1]):
                LE_stress_buckpan = []
                cor_next = self.GeoProps.str_LE_top_coor[y_idx][cor_idx + 1]
                eval_x = np.linspace(cor[0], cor_next[0], num=10, endpoint=True)
                eval_y = np.linspace(cor[1], cor_next[1], num=10, endpoint=True)
                for x_idx, x in enumerate(eval_x):
                    # print(x)
                    sigma = self.bending_stress(y_sec, x, eval_y[x_idx])
                    # print(sigma)
                    LE_stress_buckpan.append(sigma)
                max_compr = min(LE_stress_buckpan)
                LE_yeet.append(max_compr)
            str_LE_top_stress.append(LE_yeet)

        self.str_LE_top_stress = str_LE_top_stress

    def skin_compr_buckling(self):
        # The ribs devide the wing spars and skins into smaller 'segments'
        # The critical buckling stresses within each of these segments are determined here
        sigma_crit_top = []
        sigma_crit_LE_top = []

        for idx_spar, y_s in enumerate(self.GeoProps.ribs_loc[:-1]):
            y_e = self.GeoProps.ribs_loc[idx_spar + 1]
            dy = y_e - y_s
            rib_span_idxs = np.where(np.logical_and(self.GeoProps.sec_spans >= y_s, self.GeoProps.sec_spans <= y_e))

            for sec_idx in rib_span_idxs[0]:
                sigma_sec = []
                sigma_sec_LE = []
                for i in range(len(self.str_top_stress[sec_idx])):
                    w = np.sqrt((self.GeoProps.str_top_coor[sec_idx][i+1][0] - self.GeoProps.str_top_coor[sec_idx][i][0])**2 + (self.GeoProps.str_top_coor[sec_idx][i+1][1] - self.GeoProps.str_top_coor[sec_idx][i][1])**2)
                    aspect = dy/w
                    KC = self.Kc_func(aspect)
                    sigma_crit = (((np.pi**2) * KC * self.CLT.Exx)/(12*(1-self.CLT.vxy**2))) * (self.GeoProps.t23_act[sec_idx]/w)**2
                    sigma_sec.append(sigma_crit)

                for i in range(len(self.str_LE_top_stress[sec_idx])):
                    w_le = np.sqrt((self.GeoProps.str_LE_top_coor[sec_idx][i+1][0] - self.GeoProps.str_LE_top_coor[sec_idx][i][0])**2 + (self.GeoProps.str_LE_top_coor[sec_idx][i+1][1] - self.GeoProps.str_LE_top_coor[sec_idx][i][1])**2)
                    aspectle = dy/w_le
                    KC_le = self.Kc_func(aspectle)
                    sigma_crit_le = (((np.pi**2) * KC_le * self.CLT.Exx)/(12*(1-self.CLT.vxy**2))) * (self.GeoProps.t12_act[sec_idx]/w_le)**2
                    sigma_sec_LE.append(sigma_crit_le)

                sigma_crit_top.append(sigma_sec)
                sigma_crit_LE_top.append(sigma_sec_LE)

        self.top_sigma_crit = sigma_crit_top
        self.LE_top_sigma_crit = sigma_crit_LE_top

        max_compr_top = []
        compr_SF_top = []
        for idx, crit_sec in enumerate(self.top_sigma_crit):
            SF_arr = np.array(crit_sec)/np.abs(self.str_top_stress[idx])
            compr_SF_top.append(min(SF_arr))
            max_compr_top.append(max(np.abs(self.str_top_stress[idx])))

        # plt.plot(self.GeoProps.sec_spans, np.abs(self.str_top_stress) / 1000000, c='red')
        # plt.plot(self.GeoProps.sec_spans, self.top_sigma_crit / 1000000, c='blue')
        # for rib_y in self.GeoProps.ribs_loc:
        #     plt.axvline(x=rib_y, c='black', linestyle='--')
        # plt.xlabel('Wing halfspan, [m]')
        # plt.ylabel('Skin compressive stress, [MPa]')
        # plt.title('Skin compressive stress distribution over wing span')
        # plt.show()

        compr_SF_top_capped = [max(min(sigma, 10), 0) for sigma in compr_SF_top]

        # plt.plot(self.GeoProps.sec_spans, compr_SF_top_capped, c='red')
        # plt.plot(self.GeoProps.sec_spans, np.full(len(self.GeoProps.sec_spans), self.req_SF), c='blue')
        # for rib_y in self.GeoProps.ribs_loc:
        #     plt.axvline(x=rib_y, c='black', linestyle='--')
        # plt.xlabel('Wing halfspan, [m]')
        # plt.ylabel('Skin compressive buckling SF, [-]')
        # plt.title('Skin buckling safety margin distribution')
        # plt.show()

        max_compr_LE_top = []
        compr_SF_LE_top = []
        for idx, crit_sec in enumerate(self.LE_top_sigma_crit):
            LE_SF_arr = np.array(crit_sec) / np.abs(self.str_LE_top_stress[idx])
            compr_SF_LE_top.append(min(LE_SF_arr))
            max_compr_LE_top.append(max(np.abs(self.str_LE_top_stress[idx])))

        # plt.plot(self.GeoProps.sec_spans, np.abs(self.str_top_stress) / 1000000, c='red')
        # plt.plot(self.GeoProps.sec_spans, self.top_sigma_crit / 1000000, c='blue')
        # for rib_y in self.GeoProps.ribs_loc:
        #     plt.axvline(x=rib_y, c='black', linestyle='--')
        # plt.xlabel('Wing halfspan, [m]')
        # plt.ylabel('Skin compressive stress, [MPa]')
        # plt.title('Skin compressive stress distribution over wing span')
        # plt.show()

        compr_SF_LE_top_capped = [max(min(sigma, 10), 0) for sigma in compr_SF_LE_top]

        plt.figure(figsize=(8, 6))
        plt.plot(self.GeoProps.sec_spans, compr_SF_top_capped, c='red', label='Top skin')
        plt.plot(self.GeoProps.sec_spans, compr_SF_LE_top_capped, c='orange', label='LE skin')
        plt.plot(self.GeoProps.sec_spans, np.full(len(self.GeoProps.sec_spans), self.req_SF), c='blue')
        for rib_y in self.GeoProps.ribs_loc:
            plt.axvline(x=rib_y, c='black', linestyle='--')
        plt.xlabel('Wing halfspan, [m]')
        plt.ylabel('SF, [-]')
        plt.title('Skin buckling safety margin distribution')
        plt.legend()
        plt.minorticks_on()
        plt.grid(True, which='minor', color='#999999', linestyle='-', alpha=0.2)
        plt.grid(True, which='major', color='#666666', linestyle='-')
        plt.show()

    def plot_bending_stress(self, plot_idx, plot_3D=True):
        sigma_min = []
        sigma_max = []

        LE_LST = []
        MSP_LST = []
        RSP_LST = []
        LETOP_LST = []
        LEBOT_LST = []
        TOP_LST = []
        BOT_LST = []

        for span_idx in range(self.GeoProps.nr_bsec):
            y_span = self.GeoProps.sec_spans[span_idx]
            # LE spar
            LE_xarr = np.linspace(self.GeoProps.xz1[span_idx][0], self.GeoProps.xz6[span_idx][0], num=20, endpoint=True)
            LE_zarr = np.linspace(self.GeoProps.xz1[span_idx][1], self.GeoProps.xz6[span_idx][1], num=20, endpoint=True)

            LE_sigma = [self.bending_stress(y_span, LE_xarr[i], LE_zarr[i]) for i in range(len(LE_xarr))]
            LE_y = [y_span for i in range(len(LE_xarr))]
            LE_min = min(LE_sigma)
            LE_max = max(LE_sigma)

            LE_LST.append([LE_min, LE_max, LE_xarr, LE_zarr, LE_y, LE_sigma])

            # Main spar
            MSP_xarr = np.linspace(self.GeoProps.xz2[span_idx][0], self.GeoProps.xz5[span_idx][0], num=30, endpoint=True)
            MSP_zarr = np.linspace(self.GeoProps.xz2[span_idx][1], self.GeoProps.xz5[span_idx][1], num=30, endpoint=True)

            MSP_sigma = [self.bending_stress(y_span, MSP_xarr[i], MSP_zarr[i]) for i in range(len(MSP_xarr))]
            MSP_y = [y_span for i in range(len(MSP_xarr))]
            MSP_min = min(MSP_sigma)
            MSP_max = max(MSP_sigma)

            MSP_LST.append([MSP_min, MSP_max, MSP_xarr, MSP_zarr, MSP_y, MSP_sigma])

            # Main spar
            RSP_xarr = np.linspace(self.GeoProps.xz3[span_idx][0], self.GeoProps.xz4[span_idx][0], num=15, endpoint=True)
            RSP_zarr = np.linspace(self.GeoProps.xz3[span_idx][1], self.GeoProps.xz4[span_idx][1], num=15, endpoint=True)

            RSP_sigma = [self.bending_stress(y_span, RSP_xarr[i], RSP_zarr[i]) for i in range(len(RSP_xarr))]
            RSP_y = [y_span for i in range(len(RSP_xarr))]
            RSP_min = min(RSP_sigma)
            RSP_max = max(RSP_sigma)

            RSP_LST.append([RSP_min, RSP_max, RSP_xarr, RSP_zarr, RSP_y, RSP_sigma])

            # LE top
            LETOP_xarr = np.linspace(self.GeoProps.xz1[span_idx][0], self.GeoProps.xz2[span_idx][0], num=30, endpoint=True)
            LETOP_zarr = np.linspace(self.GeoProps.xz1[span_idx][1], self.GeoProps.xz2[span_idx][1], num=30, endpoint=True)

            LETOP_sigma = [self.bending_stress(y_span, LETOP_xarr[i], LETOP_zarr[i]) for i in range(len(LETOP_xarr))]
            LETOP_y = [y_span for i in range(len(LETOP_xarr))]
            LETOP_min = min(LETOP_sigma)
            LETOP_max = max(LETOP_sigma)

            LETOP_LST.append([LETOP_min, LETOP_max, LETOP_xarr, LETOP_zarr, LETOP_y, LETOP_sigma])

            # LE bot
            LEBOT_xarr = np.linspace(self.GeoProps.xz6[span_idx][0], self.GeoProps.xz5[span_idx][0], num=30, endpoint=True)
            LEBOT_zarr = np.linspace(self.GeoProps.xz6[span_idx][1], self.GeoProps.xz5[span_idx][1], num=30, endpoint=True)

            LEBOT_sigma = [self.bending_stress(y_span, LEBOT_xarr[i], LEBOT_zarr[i]) for i in range(len(LEBOT_xarr))]
            LEBOT_y = [y_span for i in range(len(LEBOT_xarr))]
            LEBOT_min = min(LEBOT_sigma)
            LEBOT_max = max(LEBOT_sigma)

            LEBOT_LST.append([LEBOT_min, LEBOT_max, LEBOT_xarr, LEBOT_zarr, LEBOT_y, LEBOT_sigma])

            # top
            TOP_xarr = np.linspace(self.GeoProps.xz2[span_idx][0], self.GeoProps.xz3[span_idx][0], num=60, endpoint=True)
            TOP_zarr = np.linspace(self.GeoProps.xz2[span_idx][1], self.GeoProps.xz3[span_idx][1], num=60, endpoint=True)

            TOP_sigma = [self.bending_stress(y_span, TOP_xarr[i], TOP_zarr[i]) for i in range(len(TOP_xarr))]
            TOP_y = [y_span for i in range(len(TOP_xarr))]
            TOP_min = min(TOP_sigma)
            TOP_max = max(TOP_sigma)

            TOP_LST.append([TOP_min, TOP_max, TOP_xarr, TOP_zarr, TOP_y, TOP_sigma])

            # bot
            BOT_xarr = np.linspace(self.GeoProps.xz5[span_idx][0], self.GeoProps.xz4[span_idx][0], num=60, endpoint=True)
            BOT_zarr = np.linspace(self.GeoProps.xz5[span_idx][1], self.GeoProps.xz4[span_idx][1], num=60, endpoint=True)

            BOT_sigma = [self.bending_stress(y_span, BOT_xarr[i], BOT_zarr[i]) for i in range(len(BOT_xarr))]
            BOT_y = [y_span for i in range(len(BOT_xarr))]
            BOT_min = min(BOT_sigma)
            BOT_max = max(BOT_sigma)

            BOT_LST.append([BOT_min, BOT_max, BOT_xarr, BOT_zarr, BOT_y, BOT_sigma])

            ALL_sigma = [*LE_sigma, *MSP_sigma, *RSP_sigma, *LETOP_sigma, *LEBOT_sigma, *TOP_sigma, *BOT_sigma]
            sig_min = min(ALL_sigma)
            sig_max = max(ALL_sigma)

            sigma_min.append(sig_min)
            sigma_max.append(sig_max)

            if span_idx == 0:
                wing_sig_max = sig_max
                wing_sig_min = sig_min

            if span_idx == plot_idx:
                ALL_xarr = [*LE_xarr, *MSP_xarr, *RSP_xarr, *LETOP_xarr, *LEBOT_xarr, *TOP_xarr, *BOT_xarr]
                ALL_zarr = [*LE_zarr, *MSP_zarr, *RSP_zarr, *LETOP_zarr, *LEBOT_zarr, *TOP_zarr, *BOT_zarr]
                wing_crossx = [x for x in self.GeoProps.w_coords[span_idx][0]]
                wing_crossz = [z for z in self.GeoProps.w_coords[span_idx][1]]

                norm = plt.Normalize(sig_min, sig_max)
                fig, ax = plt.subplots()
                plt.plot(wing_crossx, wing_crossz, c='black', linestyle='-.')
                plt.scatter(np.array(ALL_xarr), np.array(ALL_zarr), c=np.array(ALL_sigma)/1000000, cmap=plt.get_cmap('jet'))
                plt.colorbar(label='Bending stress [MPa]')
                ax.set_aspect('equal')
                plt.show()

        if plot_3D:
            ALL_y = []
            ALL_x = []
            ALL_z = []
            ALL_sigma = []
            for i_s in range(self.GeoProps.nr_bsec):
                for y in LE_LST[i_s][4]:
                    ALL_y.append(y)
                for y in MSP_LST[i_s][4]:
                    ALL_y.append(y)
                for y in RSP_LST[i_s][4]:
                    ALL_y.append(y)
                for y in LETOP_LST[i_s][4]:
                    ALL_y.append(y)
                for y in LEBOT_LST[i_s][4]:
                    ALL_y.append(y)
                for y in TOP_LST[i_s][4]:
                    ALL_y.append(y)
                for y in BOT_LST[i_s][4]:
                    ALL_y.append(y)

                for x in LE_LST[i_s][2]:
                    ALL_x.append(x)
                for x in MSP_LST[i_s][2]:
                    ALL_x.append(x)
                for x in RSP_LST[i_s][2]:
                    ALL_x.append(x)
                for x in LETOP_LST[i_s][2]:
                    ALL_x.append(x)
                for x in LEBOT_LST[i_s][2]:
                    ALL_x.append(x)
                for x in TOP_LST[i_s][2]:
                    ALL_x.append(x)
                for x in BOT_LST[i_s][2]:
                    ALL_x.append(x)

                for y in LE_LST[i_s][3]:
                    ALL_z.append(y)
                for y in MSP_LST[i_s][3]:
                    ALL_z.append(y)
                for y in RSP_LST[i_s][3]:
                    ALL_z.append(y)
                for y in LETOP_LST[i_s][3]:
                    ALL_z.append(y)
                for y in LEBOT_LST[i_s][3]:
                    ALL_z.append(y)
                for y in TOP_LST[i_s][3]:
                    ALL_z.append(y)
                for y in BOT_LST[i_s][3]:
                    ALL_z.append(y)

                for sig in LE_LST[i_s][5]:
                    ALL_sigma.append(sig)
                for sig in MSP_LST[i_s][5]:
                    ALL_sigma.append(sig)
                for sig in RSP_LST[i_s][5]:
                    ALL_sigma.append(sig)
                for sig in LETOP_LST[i_s][5]:
                    ALL_sigma.append(sig)
                for sig in LEBOT_LST[i_s][5]:
                    ALL_sigma.append(sig)
                for sig in TOP_LST[i_s][5]:
                    ALL_sigma.append(sig)
                for sig in BOT_LST[i_s][5]:
                    ALL_sigma.append(sig)

            # wing_crossx = [x for x in self.GeoProps.w_coords[i][0]]
            # wing_crossz = [z for z in self.GeoProps.w_coords[i][1]]
            fig = plt.figure()
            ax = plt.axes(projection='3d')
            ax.set_xlim3d([0, self.GeoProps.half_span*3])
            ax.set_ylim3d([0, self.GeoProps.half_span*3])
            ax.set_zlim3d([0, self.GeoProps.half_span*3])
            # make the panes transparent
            ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
            # make the grid lines transparent
            ax.xaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
            ax.yaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
            ax.zaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
            ax.set_axis_off()

            sc = ax.scatter(np.array(ALL_x), np.array(ALL_y), np.array(ALL_z), c=np.array(ALL_sigma)/1000000, cmap='jet')
            for yeet in range(self.GeoProps.nr_bsec):
                yarrr = np.repeat(self.GeoProps.sec_spans[yeet], self.GeoProps.nr_points)
                wing_crossx = [x for x in self.GeoProps.w_coords[yeet][0]]
                wing_crossz = [z for z in self.GeoProps.w_coords[yeet][1]]
                plt.plot(wing_crossx, yarrr, wing_crossz, c='black', linestyle='-.')
            plt.xlabel('Chord, [m]')
            plt.ylabel('Halfspan, [m]')
            fig.colorbar(sc, label='Bending stress [MPa]')
            # ax.set_aspect('equal')
            # for angle in range(0, 360):
            #     ax.view_init(30, angle)
            #     plt.draw()
            #     plt.pause(.1)
            plt.show()

        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_size_inches(13, 6)
        ax1.set_title('Max compressive stress over wing halfspan.')
        ax1.set_xlabel('Wing halfspan, [m]')
        ax1.set_ylabel('Bending stress, [MPa]')
        ax2.set_title('Max tensile stress over wing halfspan.')
        ax2.set_xlabel('Wing halfspan, [m]')
        ax2.set_ylabel('Tensile stress, [MPa]')
        ax1.plot(self.GeoProps.sec_spans, np.abs(sigma_min)/1000000)
        ax2.plot(self.GeoProps.sec_spans, np.array(sigma_max)/1000000)
        plt.tight_layout()
        plt.show()

    def stringer_column_buckling(self):
        pass

    def spar_column_buckling(self):
        pass

    def calc_max_min_bending_stresses(self, plot=True):
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

            if plot:
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

    def d2v_d2y(self, y):
        return -self.WingLoads.moment_dist_x(y)/(self.CLT.Exx * self.GeoProps.Ixx_dist(y))

    def dv_dy(self, y):
        dv_dy = sc.integrate.quad(self.d2v_d2y, 0, y)
        return dv_dy[0]

    def eval_dv_dy(self):
        y_arr = np.linspace(0, self.GeoProps.half_span, num=5, endpoint=True)  # Temporary smaller array to decrease computational time
        dvdy = []
        for y in y_arr:
            # print('Deflection slopes cumputed at span:', y)
            dvdy.append(self.dv_dy(y))

        self.dvdy_dist = sc.interpolate.interp1d(y_arr, np.array(dvdy), kind='cubic')
        plt.plot(y_arr, dvdy)
        plt.axis('scaled')
        plt.xlabel('Wing halfspan, [m]')
        plt.ylabel('Wing deflection slope, [-], ')
        plt.title('Wing upwards deflection under 8g steady load')
        plt.show()

    def v_z(self, y):
        v = sc.integrate.quad(self.dvdy_dist, 0, y)
        return v[0]

    def eval_v_z(self):
        y_arr = np.linspace(0, self.GeoProps.half_span, num=5, endpoint=True)  # Temporary smaller array to decrease computational time
        vz = []
        for y in y_arr:
            # print('Deflection slopes computed at span:', y)
            vz.append(self.v_z(y))

        self.v_z_dist = sc.interpolate.interp1d(y_arr, np.array(vz), kind='cubic')
        plt.plot(y_arr, vz)
        plt.axis('scaled')
        plt.xlabel('Wing halfspan, [m]')
        plt.ylabel('Wing deflection, [m], ')
        plt.title('Wing upwards deflection under 8g steady load')
        plt.show()

    def est_wing_mass(self):
        dy = self.GeoProps.half_span/(self.GeoProps.nr_bsec - 1)
        spars_mass = 0
        skins_mass = 0
        str_mass = 0
        ribs_mass = 0

        check_mass = 0

        for i in range(self.GeoProps.nr_bsec - 1):
            A_spars_i = 2 * self.GeoProps.A71[i] + 2 * self.GeoProps.A82[i] + self.GeoProps.A34[i]
            A_spars_ip1 = 2 * self.GeoProps.A71[i+1] + 2 * self.GeoProps.A82[i+1] + self.GeoProps.A34[i+1]
            A_spars_avg = (A_spars_i + A_spars_ip1)/2

            A_skins_i = self.GeoProps.A12[i] + self.GeoProps.A23[i] + self.GeoProps.A45[i] + self.GeoProps.A56[i]
            A_skins_ip1 = self.GeoProps.A12[i+1] + self.GeoProps.A23[i+1] + self.GeoProps.A45[i+1] + self.GeoProps.A56[i+1]
            A_skins_avg = (A_skins_i + A_skins_ip1)/2

            A_str_i = self.GeoProps.A_str * (len(self.GeoProps.str_top_coor[i]) + len(self.GeoProps.str_bot_coor[i]) + len(self.GeoProps.str_LE_top_coor[i]) + len(self.GeoProps.str_LE_bot_coor[i]))
            A_str_ip1 = self.GeoProps.A_str * (len(self.GeoProps.str_top_coor[i+1]) + len(self.GeoProps.str_bot_coor[i+1]) + len(self.GeoProps.str_LE_top_coor[i+1]) + len(self.GeoProps.str_LE_bot_coor[i+1]))
            A_str_avg = (A_str_i + A_str_ip1)/2

            spars_mass += A_spars_avg * dy * self.CLT.rho
            skins_mass += A_skins_avg * dy * self.CLT.rho
            str_mass += A_str_avg * dy * self.CLT.rho

            if i == 0 or i == (self.GeoProps.nr_bsec - 2):
                check_mass += (A_spars_avg + A_skins_avg + A_str_avg) * self.GeoProps.half_span * self.CLT.rho

        for rib_idx, y_rib in enumerate(self.GeoProps.ribs_loc):

            diff_arr = np.absolute(y_rib - self.GeoProps.sec_spans)
            index = diff_arr.argmin()
            A_rib = self.GeoProps.As_cell1[index] + self.GeoProps.As_cell2[index]
            print(A_rib)
            ribs_mass += A_rib * self.GeoProps.t_ribs[rib_idx] * self.CLT.rho * 1.25

        print('---------- Wing structure mass est. ----------')
        print('Spars mass, [kg]:', spars_mass)
        print('Skins mass, [kg]:', skins_mass)
        print('Stringers mass, [kg]:', str_mass)
        print('Ribs mass, [kg]:', ribs_mass)
        print('Total (half) wing structure mass, [kg]:', spars_mass + skins_mass + str_mass + ribs_mass)
        print('Total wing structure mass, [kg]:', (spars_mass + skins_mass + str_mass + ribs_mass) * 2)
        print('Check mass, [kg]:', check_mass/2)


WS = WingStresses(DIRECTORY, INPUT_FOLDER, airfoil_file, c_root, c_tip, MAC, dihedral, LE_sweep, half_span, AR,
                  nr_span_sections, nr_spars, spar_locs, t_le, t_top, t_tr, t_bot, t_spar, x_le, x_te, t_sp_flange,
                  w_sp_flange, t_red_tip, CL_dist_file, CD_dist_file, ac_mass, g_load, CL_0, CL_10,
                  S, rho, V, a, b, c, dens_spars, dens_skins, loc_ribs, q_sec61_loc, q_sec52_loc, q_sec34_loc,
                  safety_margin, t_msp_web, t_rsp_web, SMOA_simplified, nr_str_top, nr_str_bot, nr_str_LE_top,
                  nr_str_LE_bot, str_top_end, str_bot_end, str_LE_top_end, str_LE_bot_end, t_ribs,
                  CM_dist_file, COP_dist_file, t_fsp_web)


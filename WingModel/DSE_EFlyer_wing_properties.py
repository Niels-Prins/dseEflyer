# Lennart Krieg
# DSE group 20: Aerobatic E-Flyer
# Wing Discretization model, part of wing structural model

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
from scipy import integrate
from mpl_toolkits import mplot3d
import pandas as pd
import scipy as sc

#  Pycharm Project folder: C:\Users\lenna\PycharmProjects\DSE_22

DIRECTORY = 'C:/Users/lenna/Documents/DSE/Wing_structural_model'
INPUT_FOLDER = '/Input_Files/'

# Airfoil files with coordinates in x/z fractions of chord.
# Before use these files will have to be preprocessed to have them start at the LE, go to TE and wrap around again
# to the LE. Text and other non-data will have to be removed and separation by spaces replaced by semicolon (;).
airfoil_file = 'NACA-0015_LE_71.dat'

# Wing geometry parameters
c_root = 2.18  # Root chord, [m]
c_tip = 0.98  # Tip chord, [m]
t_red_tip = 0.12/0.15  # Ratio of tip thickness (of % chord) over root thickness (of % chord), [-]
MAC = 1.58  # Mean aerodynamic chord, [m]

dihedral = 0  # Wing dihedral, [deg]
LE_sweep = 0  # Leading egde sweep, [deg]
half_span = 9.2/2  # Halfspan of the wing, [m]
AR = 5.8  # Aspect ratio, [-]

nr_spars = 2  # Amount of spars, [-]
spar_locs = [0.3, 0.75]  # Location of spars as fraction of the chords, [-]. Input should be a list or an array.

nr_span_sections = 20  # Number of cross-sections in which the wing is split up, [-]

loc_ribs = [0, 0.3, 0.7, 1.2, 1.8, 3.0, half_span]  # Locations/span of the ribs, [m]

# # Thicknesses
# t_le = [2, 2, 2, 2, 2, 2]  # Leading edge thickness, [mm]
# t_top = [2, 2, 2, 2, 2, 2]  # Top skin thickness, [mm]
# t_tr = [0, 0, 0, 0, 0, 0]  # Trailing edge thickness, [mm]
# t_bot = [2, 2, 2, 2, 2, 2]  # Bottom skin thickness, [mm]
# t_spar = [3, 3]  # Spar thicknesses, [mm]. Input must be list of same size as spar_locs list!
# t_sp_flange = [3, 3]
# w_sp_flange = [100, 100]

t_le = 2  # Leading edge thickness, [mm]
t_top = 2  # Top skin thickness, [mm]
t_tr = 0  # Trailing edge thickness, [mm]
t_bot = 2  # Bottom skin thickness, [mm]
t_spar = [4, 4]  # Spar thicknesses, [mm]. Input must be list of same size as spar_locs list!
t_sp_flange = [3, 3]
w_sp_flange = [40, 40]

x_le = 0.3  # Fraction of chord where LE thickness becomes top or bottom skin thickness
x_te = 0.75  # Fraction of chord where top or bottom skin thickness become TE thickness

dens_spar = 2800  # Spar material density, [kg/m^3]
dens_skins = 2800  # Skin material density, [kg/m^3]

a_stfnrs = 112  # Stiffner/stringer area, [mm^2]

stf_top = [0.321, 0.393, 0.464, 0.536, 0.607, 0.678]
stf_bot = [0.321, 0.393, 0.464, 0.536, 0.607, 0.678]
stf_LEtop = [0.1786, 0.107]
stf_LEbot = [0.1786, 0.107]

stf_top_end = [half_span, half_span, half_span, half_span, half_span, half_span]
stf_bot_end = [half_span, half_span, half_span, half_span, half_span, half_span]
stf_LEtop_end = [half_span, half_span]
stf_LEbot_end = [half_span, half_span]


class WingProperties:
    def __init__(self, DIRECTORY, INPUT_FOLDER, airfoil_file, c_root, c_tip, MAC, dihedral, LE_sweep,
                 half_span, AR, span_sections, nr_spar, spar_loc, t_le, t_top, t_tr, t_bot, t_spar,
                 x_le, x_te, t_sp_flange, w_sp_flange, t_red_tip, dens_spars, dens_skins, loc_ribs,
                 q_sec61_loc, q_sec52_loc, q_sec34_loc, t_msp_web, t_rsp_web, SMOA_simplified,
                 nr_str_top, nr_str_bot, nr_str_LE_top, nr_str_LE_bot, str_top_end, str_bot_end,
                 str_LE_top_end, str_LE_bot_end, t_ribs, t_fsp_web):

        self.airfoil_name = airfoil_file[:-4]
        self.airfoil_input = pd.read_table(DIRECTORY + INPUT_FOLDER + airfoil_file, delimiter=';', header=None)
        self.airfoil_x = self.airfoil_input.values[:, 0]  # Airfoil x coordinates as fraction of chord
        self.airfoil_z = self.airfoil_input.values[:, 1]  # Airfoil z coordinates as fraction of chord
        self.nr_points = len(self.airfoil_x)
        self.airfoil_inter = None

        self.c_root = c_root  # Root chord, [m]
        self.c_tip = c_tip  # Tip chord, [m]
        self.MAC = MAC  # Mean aerodynamic chord, [m]

        self.dihedral = dihedral  # Wing dihedral, [deg]
        self.LE_sweep = LE_sweep  # Leading egde sweep, [deg]
        self.half_span = half_span  # Halfspan of the wing, [m]
        self.AR = AR  # Aspect ratio, [-]
        self.t_red_tip = t_red_tip  # Ratio of tip thickness (of % chord) over root thickness (of % chord), [-]

        self.nr_bsec = span_sections  # Number of cross-sections in which the wing is split up, [-]
        self.sec_spans = None  # Spans of discretized cross-sections. Declared here to be calculated later
        self.sec_chords = None  # chord at these discretized cross-sections. Declared here to be calculated later
        self.w_coords = None  # Discretized coordinate of wing cross-sections. Declared here to be calculated later
        self.ribs_loc = np.array(loc_ribs)
        self.t_ribs = np.array(t_ribs)/1000

        # Properties of simplified cross-section for shear flow calculations
        # 'Effective' thickness compensated with stringer area
        self.t82 = []
        self.t23 = []
        self.t34 = []
        self.t45 = []
        self.t58 = []
        self.t56 = []
        self.t67 = []
        self.t71 = []
        self.t12 = []

        # Actual thicknesses of the skins (no compensation for stringer areas)
        self.t23_act = []
        self.t45_act = []
        self.t56_act = []
        self.t12_act = []

        self.l82 = []
        self.l23 = []
        self.l34 = []
        self.l45 = []
        self.l58 = []
        self.l56 = []
        self.l67 = []
        self.l71 = []
        self.l12 = []

        self.A82 = []
        self.A23 = []
        self.A34 = []
        self.A45 = []
        self.A58 = []
        self.A56 = []
        self.A67 = []
        self.A71 = []
        self.A12 = []

        self.str_top_coor = []
        self.str_bot_coor = []
        self.str_LE_top_coor = []
        self.str_LE_bot_coor = []

        self.xz1 = []
        self.xz2 = []
        self.xz3 = []
        self.xz4 = []
        self.xz5 = []
        self.xz6 = []
        self.xz7 = []
        self.xz8 = []

        # Spars are defined based on start end coordinate, as the airfoil might be unsymmetrical.
        self.nr_spars = nr_spar  # Number of required spars, [-]
        self.spar_loc = spar_loc  # Fraction of chord at which spar(s) is/are required. Input given as list.
        self.spars_x = None
        self.spars_z = None
        self.spar_coords = None  # Coordinates of the start and end points of the spars. Declared to be calculated later
        self.spar_indices = None  # Indices of x coordinates at which spars are located.

        # Thicknesses
        self.t_le = np.array(t_le)/1000  # Leading edge thickness, [m]
        self.t_top = np.array(t_top)/1000  # Top skin thickness, [m]
        self.t_te = np.array(t_tr)/1000  # Trailing edge thickness, [m]
        self.t_bot = np.array(t_bot)/1000  # Bottom skin thickness, [m]
        self.t_spar = np.array(t_spar)/1000  # Spar thicknesses, [m]
        self.t_sp_flange = np.array(t_sp_flange)/1000  # Spars flange thickness, [m]
        self.w_sp_flange = np.array(w_sp_flange)/1000  # Spars flange width, [m]
        self.t_msp_web = np.array(t_msp_web)/1000
        self.t_rsp_web = np.array(t_rsp_web)/1000
        self.t_fsp_web = np.array(t_fsp_web)/1000
        self.t_distr = None  # Array with thicknesses

        self.nr_str_top = nr_str_top
        self.nr_str_bot = nr_str_bot
        self.nr_str_LE_top = nr_str_LE_top
        self.nr_str_LE_bot = nr_str_LE_bot
        self.str_top_end = str_top_end
        self.str_bot_end = str_bot_end
        self.str_LE_top_end = str_LE_top_end
        self.str_LE_bot_end = str_LE_bot_end

        self.x_le = x_le  # Fraction of chord where LE thickness becomes top or bottom skin thickness
        self.x_te = x_te  # Fraction of chord where top or bottom skin thickness become TE thickness

        self.spar_centroids = None
        self.spar_areas = None
        self.skin_centroids = None
        self.skin_areas = None
        self.seg_len = None

        t_str = 4 * 0.32/1000
        self.w_str = 15/1000
        self.h_str = 15/1000

        self.A_str = t_str*self.w_str*2 + t_str*(self.h_str-2*t_str)
        self.Ixx_str = 2 * (self.w_str*t_str**3)/12 + 2 * t_str * self.w_str * self.h_str**2 + (t_str * self.h_str**3)/12
        self.Izz_str = 2 * (t_str * self.w_str**3)/12 + 2 * t_str * self.w_str * (0.5 * self.w_str)**2 + (self.h_str*t_str**3)/12

        self.centroids = None
        self.centroid_dist = None
        self.Ixx = None
        self.Izz = None
        self.Ixx_dist = None
        self.Izz_dist = None
        self.A = None
        self.J = None
        self.zeta_x = None

        self.mass = 0
        self.q_0 = None
        self.shear_centers = []
        self.As_cell1 = []
        self.As_cell2 = []
        self.G_dtheta_dy = []
        self.tor_q_0 = []

        self.print_initials()
        self.add_spars()
        self.discretize_wing()
        self.shear_crosssections(q_sec61_loc, q_sec52_loc, q_sec34_loc)
        self.plot_wing_structure(q_sec61_loc, q_sec52_loc, q_sec34_loc)

        if SMOA_simplified:
            self.calc_SMOA_simp()
            print('Simplified geometry active')
        else:
            self.calc_centroid()
            self.calc_SMOA()

        self.interpolate_SMOAs()
        # self.est_mass(dens_spars, dens_skins)
        # self.plot_airfoil()
        self.calc_red_shear_flows()
        self.calc_shear_center()
        self.plot_crosssection(0)
        self.plot_wing()

    def print_initials(self):
        print('------------------------------------------------------------------------------')
        print('-------------------------- Wing geometry properties --------------------------')
        print('Imported airfoil geometry:', self.airfoil_name)
        print('Wingspan, [m]:', self.half_span * 2)
        print('Root chord, [m]:', self.c_root)
        print('Tip chord, [m]:', self.c_tip)
        # print('------------------------- Wing structural components -------------------------')
        # print('Nr. of spars:', self.nr_spars)
        # print('Spar locations, [% of c]:', self.spar_loc)
        # print('Spar web thicknesses, [mm]:', self.t_spar*1000)
        # print('Spar flange thicknesses, [mm]:', self.t_sp_flange*1000)
        # print('Spar flange widths, [mm]:', self.w_sp_flange * 1000)
        # print('Leading edge skin thickness, [mm]:', self.t_le*1000)
        # print('Trailing edge skin thickness, [mm]:', self.t_le * 1000)
        # print('Trailing edge skin thickness, [mm]:', self.t_le * 1000)
        # print('Top center skin thickness, [mm]:', self.t_top * 1000)
        # print('Top bottom skin thickness, [mm]:', self.t_bot * 1000)
        print('---------------------------- Wing discretization ----------------------------')
        print('Nr. of nodes:', self.nr_points)
        print('Nr of span sections:', self.nr_bsec)
        print('LE end, [% of c]:', self.x_le)
        print('TE start, [% of c]:', self.x_te)

    def scale_airfoil(self, chord, t_red):
        x_scaled = self.airfoil_x * chord
        z_scaled = self.airfoil_z * chord * t_red

        return x_scaled, z_scaled

    def scale_spar(self, chord, t_red):
        x_scaled = self.spars_x * chord
        z_scaled = self.spars_z * chord * t_red

        return x_scaled, z_scaled

    def add_spars(self):
        # Calculating the start and end z-coordinates of the spar(s) by determining closest airfoil coordinate and
        # interpolating between that coordinate and the previous/next one. This depends on whether the spar is in front
        # or behind the closest coordinate.

        spars_x_temp = []
        spars_z_temp = []
        self.airfoil_inter = sc.interpolate.interp1d(self.airfoil_x[:int((self.nr_points-1)/2)+1],
                                                     self.airfoil_z[:int((self.nr_points-1)/2)+1], kind='linear')

        for c_spar in self.spar_loc:
            z_spar = self.airfoil_inter(c_spar)
            spars_x_temp.append(c_spar)
            spars_z_temp.append(z_spar)

        self.spars_x = np.array(spars_x_temp)
        self.spars_z = np.array(spars_z_temp)

        self.spar_discretization_edit()

    def upper_lower_idx_relation(self, idx_upper):
        idx_low = self.nr_points - 1 - idx_upper
        return idx_low

    def spar_discretization_edit(self):
        # For each spar, the closest node is found and moved to the spar end locations.
        # This way, there is a clear distinction between LE, MID and TE sections
        for idx_sp, x_sp in enumerate(self.spars_x):
            z_sp = self.spars_z[idx_sp]
            # Finding nearest node at top half of the airfoil
            diff_arr = np.absolute(self.airfoil_x[:int((self.nr_points-1)/2)+1] - x_sp)
            index = diff_arr.argmin()

            self.airfoil_x[index] = x_sp
            self.airfoil_z[index] = z_sp
            # The nodes are symmetrically distributed around the x-axis, therefore the upper/lower wing indices are
            # relatable through the upper_lower_idx_relation function.
            self.airfoil_x[self.upper_lower_idx_relation(index)] = x_sp
            self.airfoil_z[self.upper_lower_idx_relation(index)] = -1 * z_sp

    def plot_airfoil(self):
        x = np.linspace(0, self.airfoil_x[:int((self.nr_points-1)/2)+1][-1], num=200, endpoint=True)
        plt.plot(x, self.airfoil_inter(x))
        plt.plot(x, -1 * self.airfoil_inter(x))
        for idx, sp_x in enumerate(self.spars_x):
            sp_z = self.spars_z[idx]
            plt.plot([sp_x, sp_x], [-1 * sp_z, sp_z], c='red')
            if idx < self.nr_spars - 1:
                plt.plot([sp_x, self.spars_x[idx+1]], [sp_z, self.spars_z[idx+1]], c='blue')
                plt.plot([sp_x, self.spars_x[idx + 1]], [-1 * sp_z, -1 * self.spars_z[idx + 1]], c='blue')

        plt.scatter(self.airfoil_x, self.airfoil_z, c='orange')
        plt.axis('scaled')
        plt.xlabel('Percentage of chord')
        plt.ylabel('Thickness as percentage of chord')
        plt.title('Wing cross-section at the root.')
        plt.show()

    def plot_crosssection(self, span_idx):
        x_arr = self.w_coords[span_idx][0]
        z_arr = self.w_coords[span_idx][1]

        # Plotting wing cross-section
        plt.plot(x_arr, z_arr, c='blue')

        # Plotting simplified structure cross-section
        plt.plot([self.xz1[span_idx][0], self.xz2[span_idx][0]], [self.xz1[span_idx][1], self.xz2[span_idx][1]], c='red')
        plt.plot([self.xz2[span_idx][0], self.xz3[span_idx][0]], [self.xz2[span_idx][1], self.xz3[span_idx][1]], c='red')
        plt.plot([self.xz3[span_idx][0], self.xz4[span_idx][0]], [self.xz3[span_idx][1], self.xz4[span_idx][1]], c='red')
        plt.plot([self.xz4[span_idx][0], self.xz5[span_idx][0]], [self.xz4[span_idx][1], self.xz5[span_idx][1]], c='red')
        plt.plot([self.xz5[span_idx][0], self.xz2[span_idx][0]], [self.xz5[span_idx][1], self.xz2[span_idx][1]], c='red')
        plt.plot([self.xz5[span_idx][0], self.xz6[span_idx][0]], [self.xz5[span_idx][1], self.xz6[span_idx][1]], c='red')
        plt.plot([self.xz1[span_idx][0], self.xz6[span_idx][0]], [self.xz1[span_idx][1], self.xz6[span_idx][1]], c='red')

        # Plotting stringers
        for cor in self.str_top_coor[span_idx]:
            plt.scatter(cor[0], cor[1], c='black')
        for cor in self.str_bot_coor[span_idx]:
            plt.scatter(cor[0], cor[1], c='black')
        for cor in self.str_LE_top_coor[span_idx]:
            plt.scatter(cor[0], cor[1], c='black')
        for cor in self.str_LE_bot_coor[span_idx]:
            plt.scatter(cor[0], cor[1], c='black')

        # Plot centroid
        plt.scatter(self.centroids[span_idx], 0, c='green', marker='+')

        # Plot shear center
        plt.scatter(self.shear_centers[span_idx], 0, c='magenta', marker='+')
        plt.minorticks_on()
        plt.grid(True, which='minor', color='#999999', linestyle='-', alpha=0.2)
        plt.grid(True, which='major', color='#666666', linestyle='-')
        plt.axis('scaled')
        plt.xlabel('Chord (x), [m]')
        plt.ylabel('Thickness (z), ')
        if span_idx == 0:
            plt.title('Wing cross-section at the root of the wing.')
        else:
            plt.title(f'Wing cross-section at {self.sec_spans[span_idx]} m from the root.')
        plt.show()

    def discretize_wing(self):
        self.sec_spans = np.linspace(0, self.half_span, num=self.nr_bsec, endpoint=True)  # Spans of the cross-sections
        chord_deltas = np.linspace(0, self.c_tip - self.c_root, num=self.nr_bsec,
                                   endpoint=True)  # Delta chord wrt, root at these cross-sections
        thick_reductions = np.linspace(1, self.t_red_tip, num=self.nr_bsec, endpoint=True)
        self.sec_chords = chord_deltas + self.c_root
        XZaf_temp = []
        XZsp_temp = []
        for idx, c in enumerate(self.sec_chords):
            x_af, z_af = self.scale_airfoil(c, thick_reductions[idx])
            XZaf_temp.append([x_af, z_af])

            x_sp, z_sp = self.scale_spar(c, thick_reductions[idx])
            XZsp_temp.append([x_sp, z_sp])

        self.w_coords = np.array(XZaf_temp)
        self.spar_coords = np.array(XZsp_temp)

    def calc_centroid(self):
        # To calculate the length of the linear segments/lines in which the airfoil is discretized.
        # The order of this list corresponds to the x and z arrays within w_coords, except for the last element
        # since the start and end points overlap.
        len_temp = []
        cen_temp = []
        t_temp = []
        for span_idx, sec in enumerate(self.w_coords):
            x = sec[0]
            z = sec[1]
            y = self.sec_spans[span_idx]
            diff_arr = np.absolute(self.ribs_loc[:-1] - y)
            idx_close = diff_arr.argmin()
            if y >= self.ribs_loc[idx_close]:
                t_idx = idx_close
            else:
                t_idx = idx_close - 1
            sec_len_temp = []
            sec_cen_temp = []
            sec_t_temp = []
            for idx, x_cor in enumerate(x):
                if idx < self.nr_points - 1:
                    x_next = x[idx + 1]
                    z_cor = z[idx]
                    z_next = z[idx + 1]
                    dx = x_next - x_cor
                    dz = z_next - z_cor
                    d = np.sqrt(dx**2 + dz**2)
                    c = (x_cor + dx / 2, z_cor + dz / 2)
                    sec_len_temp.append(d)
                    sec_cen_temp.append(c)

                    # Splitting up the wing cross-section into the different thicknesses.
                    # Storing these thicknesses in an array where for each line segment of each spanwise cross-section
                    # the thickness is stored at the same index as the w_coords array.
                    chord = self.sec_chords[span_idx]
                    if x_cor < chord*self.x_le:
                        t = self.t_le[t_idx]
                    elif x_cor > chord * self.x_te:
                        t = self.t_te[t_idx]
                    else:
                        if z_cor >= 0:
                            t = self.t_top[t_idx]
                        else:
                            t = self.t_bot[t_idx]

                    sec_t_temp.append(t)

            len_temp.append(sec_len_temp)
            cen_temp.append(sec_cen_temp)
            t_temp.append(sec_t_temp)

        self.seg_len = np.array(len_temp)
        self.skin_centroids = np.array(cen_temp)
        self.t_distr = np.array(t_temp)
        self.skin_areas = self.seg_len * self.t_distr

        # Calculating the centroid location of the wing structure.
        # Since this model is written for symmetrical airfoils, the centroid will lay on the symmetry axis (x-axis).
        # Spar centroids:
        spar_areas_temp = []
        spar_cens_temp = []

        for sec_spar in self.spar_coords:
            x = sec_spar[0]
            z = sec_spar[1]
            spar_a = []
            spar_c = []
            for idx, x_pos in enumerate(x):
                z_top = z[idx]
                # Area of spar web + flanges
                area = self.t_spar[idx] * 2 * z_top + 2 * self.t_sp_flange[idx] * self.w_sp_flange[idx]
                spar_a.append(area)
                spar_c.append(x_pos)
            spar_areas_temp.append(spar_a)
            spar_cens_temp.append(spar_c)

        self.spar_areas = np.array(spar_areas_temp)
        self.spar_centroids = np.array(spar_cens_temp)

        c_temp = []
        for idx, sec_area in enumerate(self.skin_areas):
            a_tot = np.sum(sec_area)
            ax_tot = 0
            for a_idx, a in enumerate(sec_area):
                ax_tot += a * self.skin_centroids[idx][a_idx][0]

            for idx_s, spar_area in enumerate(self.spar_areas[idx]):
                spar_cen = self.spar_centroids[idx][idx_s]
                a_tot += spar_area
                ax_tot += spar_area * spar_cen

            c = ax_tot / a_tot

            c_temp.append(c)

        self.centroids = np.array(c_temp)
        print('full xc', self.centroids[0])

    def calc_SMOA(self):
        Ixx_temp = []
        Izz_temp = []
        for idx, sec in enumerate(self.w_coords):
            x_arr = sec[0]
            z_arr = sec[1]
            Ixx = 0
            Izz = 0
            for idx_x, x in enumerate(x_arr[:int((self.nr_points-1)/2)+1]):
                if idx_x < int((self.nr_points-1)/2):
                    z = z_arr[idx_x]
                    x_next = x_arr[idx_x + 1]
                    z_next = z_arr[idx_x + 1]
                    dx = x_next - x
                    dz = z_next - z
                    beta = np.arctan(dz/dx)
                    l = np.sqrt(dx**2 + dz**2)
                    t = self.t_distr[idx][idx_x]
                    c = self.skin_centroids[idx][idx_x]
                    a = self.skin_areas[idx][idx_x]
                    Ixx += 2 * ((l**3 * t * np.sin(beta)**2)/12 + a * (c[1])**2)
                    Izz += 2 * ((l**3 * t * np.cos(beta)**2)/12 + a * (self.centroids[idx] - c[0])**2)

            for idx_s, spar_x in enumerate(self.spar_coords[idx][0]):
                spar_top_z = self.spar_coords[idx][1][idx_s]
                # for Ixx only spar flange Steiner term is taken, for Izz both about own centroid as well as Steiner
                Ixx += ((spar_top_z*2)**3 * self.t_spar[idx_s])/12 + 2 * self.t_sp_flange[idx_s] * self.w_sp_flange[idx_s] * spar_top_z**2
                Izz += (self.t_spar[idx_s]**3 * (spar_top_z*2))/12 + self.spar_areas[idx][idx_s] * (self.centroids[idx] - spar_x)**2
                Izz += 2 * (self.t_sp_flange[idx_s] * (self.w_sp_flange[idx_s]**3)/12 + self.t_sp_flange[idx_s] * self.w_sp_flange[idx_s] * (self.centroids[idx] - spar_x)**2)

            Ixx_temp.append(Ixx)
            Izz_temp.append(Izz)

        self.Ixx = np.array(Ixx_temp)
        self.Izz = np.array(Izz_temp)
        print('full Ixx', self.Ixx[0])
        print('full Izz', self.Izz[0])

    def interpolate_SMOAs(self):
        self.Ixx_dist = sc.interpolate.interp1d(self.sec_spans, self.Ixx, kind='cubic')
        self.Izz_dist = sc.interpolate.interp1d(self.sec_spans, self.Izz, kind='cubic')
        self.centroid_dist = sc.interpolate.interp1d(self.sec_spans, self.centroids, kind='linear')

    def est_mass(self, rho_spars, rho_skins):
        mass = 0
        # Iterating through the discretized span sections and estimate the mass based on the average of area between
        # the ith section and the i+1th.
        for idx, areas_skin_i in enumerate(self.skin_areas):
            if idx < self.nr_bsec - 1:
                a_skins_i = sum(areas_skin_i)
                a_spars_i = sum(self.spar_areas[idx])

                a_skins_ip1 = sum(self.skin_areas[idx+1])
                a_spars_ip1 = sum(self.spar_areas[idx+1])

                a_skin_avg = (a_skins_i+a_skins_ip1)
                a_spar_avg = (a_spars_i + a_spars_ip1)
                mass += (a_skin_avg * rho_skins + a_spar_avg * rho_spars) * self.half_span/(self.nr_bsec-1)

        self.mass = mass
        print('Estimated wing mass, [kg]:', round(self.mass, 2))

    # Unit shear flow distributions to calculate

    def shear_crosssections(self, x_c1, x_c2, x_c3):
        thick_red = np.linspace(1, self.t_red_tip, num=self.nr_bsec, endpoint=True)

        z_c1 = self.airfoil_inter(x_c1)
        z_c2 = self.airfoil_inter(x_c2)
        z_c3 = self.airfoil_inter(x_c3)

        for c_idx, c in enumerate(self.sec_chords):
            self.xz1.append([x_c1 * c, z_c1 * c * thick_red[c_idx]])
            self.xz2.append([x_c2 * c, z_c2 * c * thick_red[c_idx]])
            self.xz3.append([x_c3 * c, z_c3 * c * thick_red[c_idx]])

            self.xz4.append([x_c3 * c, -1 * z_c3 * c * thick_red[c_idx]])
            self.xz5.append([x_c2 * c, -1 * z_c2 * c * thick_red[c_idx]])
            self.xz6.append([x_c1 * c, -1 * z_c1 * c * thick_red[c_idx]])

            self.xz7.append([x_c1 * c, 0])
            self.xz8.append([x_c2 * c, 0])

            self.l82.append(np.sqrt((np.array(self.xz2[c_idx][0]) - np.array(self.xz8[c_idx][0])) ** 2 + (np.array(self.xz2[c_idx][1]) - np.array(self.xz8[c_idx][1])) ** 2))
            self.l23.append(np.sqrt((np.array(self.xz3[c_idx][0]) - np.array(self.xz2[c_idx][0])) ** 2 + (np.array(self.xz3[c_idx][1]) - np.array(self.xz2[c_idx][1])) ** 2))
            self.l34.append(np.sqrt((np.array(self.xz4[c_idx][0]) - np.array(self.xz3[c_idx][0])) ** 2 + (np.array(self.xz4[c_idx][1]) - np.array(self.xz3[c_idx][1])) ** 2))
            self.l45.append(np.sqrt((np.array(self.xz5[c_idx][0]) - np.array(self.xz4[c_idx][0])) ** 2 + (np.array(self.xz5[c_idx][1]) - np.array(self.xz4[c_idx][1])) ** 2))
            self.l58.append(np.sqrt((np.array(self.xz8[c_idx][0]) - np.array(self.xz5[c_idx][0])) ** 2 + (np.array(self.xz8[c_idx][1]) - np.array(self.xz5[c_idx][1])) ** 2))
            self.l56.append(np.sqrt((np.array(self.xz6[c_idx][0]) - np.array(self.xz5[c_idx][0])) ** 2 + (np.array(self.xz6[c_idx][1]) - np.array(self.xz5[c_idx][1])) ** 2))
            self.l67.append(np.sqrt((np.array(self.xz7[c_idx][0]) - np.array(self.xz6[c_idx][0])) ** 2 + (np.array(self.xz7[c_idx][1]) - np.array(self.xz6[c_idx][1])) ** 2))
            self.l71.append(np.sqrt((np.array(self.xz1[c_idx][0]) - np.array(self.xz7[c_idx][0])) ** 2 + (np.array(self.xz1[c_idx][1]) - np.array(self.xz7[c_idx][1])) ** 2))
            self.l12.append(np.sqrt((np.array(self.xz2[c_idx][0]) - np.array(self.xz1[c_idx][0])) ** 2 + (np.array(self.xz2[c_idx][1]) - np.array(self.xz1[c_idx][1])) ** 2))

            y = self.sec_spans[c_idx]
            diff_arr = np.absolute(self.ribs_loc[:-1] - y)
            idx_close = diff_arr.argmin()
            if y >= self.ribs_loc[idx_close]:
                t_idx = idx_close
            else:
                t_idx = idx_close - 1

            self.A82.append(self.t_msp_web[t_idx] * self.l82[c_idx])
            self.A23.append(self.t_top[t_idx] * self.l23[c_idx])
            self.A34.append(self.t_rsp_web[t_idx] * self.l34[c_idx])
            self.A45.append(self.t_bot[t_idx] * self.l45[c_idx])
            self.A58.append(self.t_msp_web[t_idx] * self.l58[c_idx])

            self.A56.append(self.t_le[t_idx] * self.l56[c_idx])
            self.A67.append(self.t_fsp_web[t_idx] * self.l67[c_idx])
            self.A71.append(self.t_fsp_web[t_idx] * self.l71[c_idx])
            self.A12.append(self.t_le[t_idx] * self.l12[c_idx])

            # Locating stringers on top skin
            dx_str_top = (self.xz3[c_idx][0] - self.xz2[c_idx][0]) / (self.nr_str_top - 1)
            dz_str_top = (self.xz3[c_idx][1] - self.xz2[c_idx][1]) / (self.nr_str_top - 1)

            sec_str_top_coor = []
            for i in range(self.nr_str_top):
                if self.sec_spans[c_idx] <= self.str_top_end[i]:
                    sec_str_top_coor.append([self.xz2[c_idx][0] + i * dx_str_top, self.xz2[c_idx][1] + i * dz_str_top])
            self.str_top_coor.append(sec_str_top_coor)

            # Locating stringers on bottom skin
            dx_str_bot = (self.xz4[c_idx][0] - self.xz5[c_idx][0]) / (self.nr_str_bot - 1)
            dz_str_bot = (self.xz4[c_idx][1] - self.xz5[c_idx][1]) / (self.nr_str_bot - 1)

            sec_str_bot_coor = []
            for i in range(self.nr_str_bot):
                if self.sec_spans[c_idx] <= self.str_bot_end[i]:
                    sec_str_bot_coor.append([self.xz5[c_idx][0] + i * dx_str_bot, self.xz5[c_idx][1] + i * dz_str_bot])
            self.str_bot_coor.append(sec_str_bot_coor)

            # Locating stringers on top LE skin
            dx_str_LE_top = (self.xz2[c_idx][0] - self.xz1[c_idx][0]) / (self.nr_str_LE_top - 1)
            dz_str_LE_top = (self.xz2[c_idx][1] - self.xz1[c_idx][1]) / (self.nr_str_LE_top - 1)

            sec_str_LE_top_coor = []
            for i in range(self.nr_str_LE_top):
                if self.sec_spans[c_idx] <= self.str_LE_top_end[i]:
                    sec_str_LE_top_coor.append([self.xz1[c_idx][0] + i * dx_str_LE_top, self.xz1[c_idx][1] + i * dz_str_LE_top])
            self.str_LE_top_coor.append(sec_str_LE_top_coor)

            # Locating stringers on top LE skin
            dx_str_LE_bot = (self.xz5[c_idx][0] - self.xz6[c_idx][0]) / (self.nr_str_LE_bot - 1)
            dz_str_LE_bot = (self.xz5[c_idx][1] - self.xz6[c_idx][1]) / (self.nr_str_LE_bot - 1)

            sec_str_LE_bot_coor = []
            for i in range(self.nr_str_LE_bot):
                if self.sec_spans[c_idx] <= self.str_LE_bot_end[i]:
                    sec_str_LE_bot_coor.append([self.xz6[c_idx][0] + i * dx_str_LE_bot, self.xz6[c_idx][1] + i * dz_str_LE_bot])
            self.str_LE_bot_coor.append(sec_str_LE_bot_coor)

            self.t82.append(self.t_msp_web[t_idx])
            self.t23.append(self.t_top[t_idx] + len(self.str_top_coor[c_idx]) * self.A_str/self.l23[c_idx])
            self.t34.append(self.t_rsp_web[t_idx])
            self.t45.append(self.t_bot[t_idx] + len(self.str_bot_coor[c_idx]) * self.A_str/self.l45[c_idx])
            self.t58.append(self.t_msp_web[t_idx])

            self.t56.append(self.t_le[t_idx] + len(self.str_LE_bot_coor[c_idx]) * self.A_str/self.l56[c_idx])
            self.t67.append(self.t_fsp_web[t_idx])
            self.t71.append(self.t_fsp_web[t_idx])
            self.t12.append(self.t_le[t_idx] + len(self.str_LE_top_coor[c_idx]) * self.A_str/self.l12[c_idx])

            self.t23_act.append(self.t_top[t_idx])
            self.t45_act.append(self.t_bot[t_idx])
            self.t56_act.append(self.t_le[t_idx])
            self.t12_act.append(self.t_le[t_idx])

            self.As_cell1.append((self.xz1[c_idx][1] + self.xz2[c_idx][1]) * (self.xz2[c_idx][0] - self.xz1[c_idx][0]))
            self.As_cell2.append((self.xz2[c_idx][1] + self.xz3[c_idx][1]) * (self.xz3[c_idx][0] - self.xz2[c_idx][0]))

        self.xz1 = np.array(self.xz1)
        self.xz2 = np.array(self.xz2)
        self.xz3 = np.array(self.xz3)

        self.xz4 = np.array(self.xz4)
        self.xz5 = np.array(self.xz5)
        self.xz6 = np.array(self.xz6)

        self.xz7 = np.array(self.xz7)
        self.xz8 = np.array(self.xz8)

    def calc_SMOA_simp(self):
        centroids_x = []
        IXX = []
        IZZ = []
        for i in range(self.nr_bsec):

            c_12 = ((self.xz1[i][0] + (self.xz2[i][0] - self.xz1[i][0]) / 2), (self.xz1[i][1] + (self.xz2[i][1] - self.xz1[i][1]) / 2))
            c_23 = ((self.xz2[i][0] + (self.xz3[i][0] - self.xz2[i][0]) / 2), (self.xz2[i][1] + (self.xz3[i][1] - self.xz2[i][1]) / 2))
            c_56 = ((self.xz6[i][0] + (self.xz5[i][0] - self.xz6[i][0]) / 2), (self.xz6[i][1] + (self.xz5[i][1] - self.xz6[i][1]) / 2))
            c_45 = ((self.xz5[i][0] + (self.xz4[i][0] - self.xz5[i][0]) / 2), (self.xz5[i][1] + (self.xz4[i][1] - self.xz5[i][1]) / 2))

            Ax_12 = self.A12[i] * c_12[0]
            Ax_23 = self.A23[i] * c_23[0]
            Ax_56 = self.A56[i] * c_56[0]
            Ax_45 = self.A45[i] * c_45[0]

            Ax_67 = self.A67[i] * self.xz7[i][0]
            Ax_71 = self.A71[i] * self.xz7[i][0]
            Ax_58 = self.A58[i] * self.xz8[i][0]
            Ax_82 = self.A82[i] * self.xz8[i][0]
            Ax_34 = self.A34[i] * self.xz3[i][0]

            Ax_str_top = 0
            for cor in self.str_top_coor[i]:
                Ax_str_top += self.A_str * cor[0]

            Ax_str_bot = 0
            for cor in self.str_bot_coor[i]:
                Ax_str_bot += self.A_str * cor[0]

            Ax_str_LE_top = 0
            for cor in self.str_LE_top_coor[i]:
                Ax_str_LE_top += self.A_str * cor[0]

            Ax_str_LE_bot = 0
            for cor in self.str_LE_bot_coor[i]:
                Ax_str_LE_bot += self.A_str * cor[0]

            # Calculating centroid location
            Ax = Ax_12 + Ax_23 + Ax_56 + Ax_45 + Ax_67 + Ax_71 + Ax_58 + Ax_82 + Ax_34 + \
                 Ax_str_top + Ax_str_bot + Ax_str_LE_top + Ax_str_LE_bot
            A = self.A12[i] + self.A23[i] + self.A56[i] + self.A45[i] + self.A67[i] + self.A71[i] + self.A58[i] + self.A82[i] + \
                self.A34[i] + len(self.str_top_coor[i]) * self.A_str + len(self.str_bot_coor[i]) * self.A_str + \
                len(self.str_LE_top_coor[i]) * self.A_str + len(self.str_LE_bot_coor[i]) * self.A_str
            c_x = Ax/A
            centroids_x.append(c_x)

            # Calculating SMOA's
            dx_12 = self.xz2[i][0] - self.xz1[i][0]
            dz_12 = self.xz2[i][1] - self.xz1[i][1]
            beta_12 = np.arctan(dz_12 / dx_12)

            Ixx_12 = (self.l12[i] ** 3 * self.t12_act[i] * np.sin(beta_12) ** 2) / 12 + self.A12[i] * c_12[1]**2
            Izz_12 = (self.l12[i] ** 3 * self.t12_act[i] * np.cos(beta_12) ** 2) / 12 + self.A12[i] * (c_x - c_12[0])**2

            dx_23 = self.xz3[i][0] - self.xz2[i][0]
            dz_23 = self.xz3[i][1] - self.xz2[i][1]
            beta_23 = np.arctan(dz_23 / dx_23)

            Ixx_23 = (self.l23[i] ** 3 * self.t23_act[i] * np.sin(beta_23) ** 2) / 12 + self.A23[i] * c_23[1]**2
            Izz_23 = (self.l23[i] ** 3 * self.t23_act[i] * np.cos(beta_23) ** 2) / 12 + self.A23[i] * (c_x - c_23[0])**2

            dx_56 = self.xz6[i][0] - self.xz5[i][0]
            dz_56 = self.xz6[i][1] - self.xz5[i][1]
            beta_56 = np.arctan(dz_56 / dx_56)

            Ixx_56 = (self.l56[i] ** 3 * self.t56_act[i] * np.sin(beta_56) ** 2) / 12 + self.A56[i] * c_56[1]**2
            Izz_56 = (self.l56[i] ** 3 * self.t56_act[i] * np.cos(beta_56) ** 2) / 12 + self.A56[i] * (c_x - c_56[0])**2

            dx_45 = self.xz5[i][0] - self.xz4[i][0]
            dz_45 = self.xz5[i][1] - self.xz4[i][1]
            beta_45 = np.arctan(dz_45 / dx_45)

            Ixx_45 = (self.l45[i] ** 3 * self.t45_act[i] * np.sin(beta_45) ** 2) / 12 + self.A45[i] * c_45[1] ** 2
            Izz_45 = (self.l45[i] ** 3 * self.t45_act[i] * np.cos(beta_45) ** 2) / 12 + self.A45[i] * (c_x - c_45[0]) ** 2

            Ixx_str_top = 0
            Izz_str_top = 0
            for cor in self.str_top_coor[i]:
                Ixx_str_top += self.Ixx_str + self.A_str * cor[1]**2
                Izz_str_top += self.Izz_str + self.A_str * (c_x - cor[0])**2

            Ixx_str_bot = 0
            Izz_str_bot = 0
            for cor in self.str_bot_coor[i]:
                Ixx_str_bot += self.Ixx_str + self.A_str * cor[1]**2
                Izz_str_bot += self.Izz_str + self.A_str * (c_x - cor[0])**2

            Ixx_str_LE_top = 0
            Izz_str_LE_top = 0
            for cor in self.str_LE_top_coor[i]:
                Ixx_str_LE_top += self.Ixx_str + self.A_str * cor[1]**2
                Izz_str_LE_top += self.Izz_str + self.A_str * (c_x - cor[0])**2

            Ixx_str_LE_bot = 0
            Izz_str_LE_bot = 0
            for cor in self.str_LE_bot_coor[i]:
                Ixx_str_LE_bot += self.Ixx_str + self.A_str * cor[1]**2
                Izz_str_LE_bot += self.Izz_str + self.A_str * (c_x - cor[0])**2

            Ixx_61 = (self.t67[i] * (self.l67[i] * 2)**3)/12
            Izz_61 = ((self.l67[i] * 2) * self.t67[i]**3)/12 + 2 * self.A67[i] * (c_x - self.xz7[i][0])**2

            Ixx_52 = (self.t58[i] * (self.l58[i] * 2) ** 3) / 12
            Izz_52 = ((self.l58[i] * 2) * self.t58[i] ** 3) / 12 + 2 * self.A58[i] * (c_x - self.xz8[i][0]) ** 2

            Ixx_34 = (self.t34[i] * (self.l34[i] * 2) ** 3) / 12
            Izz_34 = ((self.l34[i] * 2) * self.t34[i] ** 3) / 12 + 2 * self.A34[i] * (c_x - self.xz3[i][0]) ** 2

            ixx = Ixx_12 + Ixx_23 + Ixx_56 + Ixx_45 + Ixx_str_top + Ixx_str_bot + Ixx_str_LE_top + Ixx_str_LE_bot + \
                  Ixx_61 + Ixx_52 + Ixx_34
            izz = Izz_12 + Izz_23 + Izz_56 + Izz_45 + Izz_str_top + Izz_str_bot + Izz_str_LE_top + Izz_str_LE_bot + \
                  Izz_61 + Izz_52 + Izz_34

            IXX.append(ixx)
            IZZ.append(izz)

        self.centroids = np.array(centroids_x)
        self.Ixx = np.array(IXX)
        self.Izz = np.array(IZZ)

        print('full xc', self.centroids[0])
        print('full Ixx', self.Ixx[0])
        print('full Izz', self.Izz[0])

        fig, (ax1, ax2) = plt.subplots(1, 2)

        # fig.set_size_inches(13, 6)
        # fig.title('Ixx and Izz along the wing halfspan.')
        ax1.set_xlabel('Wing halfspan, [m]')
        ax1.set_ylabel(r'Ixx, [$m^4$]')
        ax2.set_xlabel('Wing halfspan, [m]')
        ax2.set_ylabel(r'Izz, [$m^4$]')
        ax1.plot(self.sec_spans, self.Ixx)
        ax2.plot(self.sec_spans, self.Izz)
        ax1.minorticks_on()
        ax1.grid(True, which='minor', color='#999999', linestyle='-', alpha=0.2)
        ax1.grid(True, which='major', color='#666666', linestyle='-')
        ax2.minorticks_on()
        ax2.grid(True, which='minor', color='#999999', linestyle='-', alpha=0.2)
        ax2.grid(True, which='major', color='#666666', linestyle='-')
        # plt.tight_layout()
        plt.show()

    def q_b82(self, span_idx, z):
        qb = -1 * (self.t82[span_idx]/(2 * self.Ixx[span_idx])) * z**2
        return qb

    def q_b82_int(self, span_idx, z):
        qb_inf = -1 * (self.t82[span_idx]/(6 * self.Ixx[span_idx])) * z**3
        return qb_inf

    def q_b23(self, span_idx, s):
        t = self.t23[span_idx]
        # print('span idx', span_idx)
        # print('s', s)
        # print('l82', self.l82[span_idx])
        # print('y2', self.xz2[span_idx][1])
        dy = self.xz3[span_idx][1] - self.xz2[span_idx][1]
        qb82_2 = self.q_b82(span_idx, self.l82[span_idx])
        qb12_2 = self.q_b12(span_idx, self.l12[span_idx])

        qb = -1 * (t/self.Ixx[span_idx]) * (self.xz2[span_idx][1] * s + dy/(2 * self.l23[span_idx])*s**2) + qb82_2 + qb12_2
        return qb

    def q_b23_int(self, span_idx, s):
        t = self.t23[span_idx]
        dy = self.xz3[span_idx][1] - self.xz2[span_idx][1]
        qb82_2 = self.q_b82(span_idx, self.l82[span_idx])
        qb12_2 = self.q_b12(span_idx, self.l12[span_idx])

        qb_int = -1 * (t/self.Ixx[span_idx]) * (0.5 * self.xz2[span_idx][1] * s**2 + dy/(6 * self.l23[span_idx])*s**3) + qb82_2 * s + qb12_2 * s
        return qb_int

    def q_b34(self, span_idx, s):
        t = self.t34[span_idx]
        qb23_3 = self.q_b23(span_idx, self.l23[span_idx])
        qb = -1 * (t/self.Ixx[span_idx]) * (self.xz3[span_idx][1] * s - 0.5 * s**2) + qb23_3
        return qb

    def q_b34_int(self, span_idx, s):
        t = self.t34[span_idx]
        qb23_3 = self.q_b23(span_idx, self.l23[span_idx])
        qb_int = -1 * (t/self.Ixx[span_idx]) * (0.5 * self.xz3[span_idx][1] * s**2 - (s**3)/6) + qb23_3 * s
        return qb_int

    def q_b45(self, span_idx, s):
        t = self.t45[span_idx]
        dy = self.xz5[span_idx][1] - self.xz4[span_idx][1]
        qb34_4 = self.q_b34(span_idx, self.l34[span_idx])
        qb = -1 * (t/self.Ixx[span_idx]) * (self.xz4[span_idx][1] * s + dy/(2 * self.l45[span_idx]) * s**2) + qb34_4
        return qb

    def q_b45_int(self, span_idx, s):
        t = self.t45[span_idx]
        dy = self.xz5[span_idx][1] - self.xz4[span_idx][1]
        qb34_4 = self.q_b34(span_idx, self.l34[span_idx])
        qb_int = -1 * (t/self.Ixx[span_idx]) * (0.5 * self.xz4[span_idx][1] * s**2 + dy/(6 * self.l45[span_idx]) * s**3) + qb34_4 * s
        return qb_int

    def q_b58(self, span_idx, z):
        # t = self.t58[span_idx]
        # qb = (t/(2 * self.Ixx[span_idx])) * (self.l58[span_idx] - s)**2
        qb = -1 * (self.t58[span_idx] / (2 * self.Ixx[span_idx])) * z**2
        return qb

    def q_b58_int(self, span_idx, z):
        qb_int = (self.t58[span_idx] / (6 * self.Ixx[span_idx])) * z**3
        return qb_int

    def q_b56(self, span_idx, s):
        t = self.t56[span_idx]
        dy = self.xz6[span_idx][1] - self.xz5[span_idx][1]
        qb45_5 = self.q_b45(span_idx, self.l45[span_idx])
        qb58_5 = self.q_b58(span_idx, self.xz5[span_idx][1])
        qb = -t/self.Ixx[span_idx] * (self.xz5[span_idx][1] * s + dy/(2 * self.l56[span_idx]) * s**2) + qb45_5 - qb58_5
        return qb

    def q_b56_int(self, span_idx, s):
        t = self.t56[span_idx]
        dy = self.xz6[span_idx][1] - self.xz5[span_idx][1]
        qb45_5 = self.q_b45(span_idx, self.l45[span_idx])
        qb58_5 = self.q_b58(span_idx, self.xz5[span_idx][1])
        qb_int = -t/self.Ixx[span_idx] * (0.5 * self.xz5[span_idx][1] * s**2 + dy/(6 * self.l56[span_idx]) * s**3) + qb45_5 * s - qb58_5 * s
        return qb_int

    def q_b67(self, span_idx, s):
        t = self.t67[span_idx]
        qb56_6 = self.q_b56(span_idx, self.l56[span_idx])
        qb = -1 * t/self.Ixx[span_idx] * (self.xz6[span_idx][1] * s + 0.5 * s**2) + qb56_6
        return qb

    def q_b67_int(self, span_idx, s):
        t = self.t67[span_idx]
        qb56_6 = self.q_b56(span_idx, self.l56[span_idx])
        qb_int = -1 * t/self.Ixx[span_idx] * (0.5 * self.xz6[span_idx][1] * s**2 + (s**3)/6) + qb56_6 * s
        return qb_int

    def q_b71(self, span_idx, s):
        t = self.t71[span_idx]
        qb = -t/(2 * self.Ixx[span_idx]) * s**2
        return qb

    def q_b71_int(self, span_idx, s):
        t = self.t71[span_idx]
        qb_int = -t/(6 * self.Ixx[span_idx]) * s**3
        return qb_int

    def q_b12(self, span_idx, s):
        t = self.t12[span_idx]
        dy = self.xz2[span_idx][1] - self.xz1[span_idx][1]
        qb71_1 = self.q_b71(span_idx, self.l71[span_idx])
        qb = -t/self.Ixx[span_idx] * (self.xz1[span_idx][1] * s + dy/(2 * self.l12[span_idx]) * s**2) + qb71_1
        return qb

    def q_b12_int(self, span_idx, s):
        t = self.t12[span_idx]
        dy = self.xz2[span_idx][1] - self.xz1[span_idx][1]
        qb71_1 = self.q_b71(span_idx, self.l71[span_idx])
        qb_int = -t/self.Ixx[span_idx] * (0.5 * self.xz1[span_idx][1] * s**2 + dy/(6 * self.l12[span_idx]) * s**3) + qb71_1 * s
        return qb_int

    def calc_red_shear_flows(self):
        # Calculating redundant shear flow based on rate of twist. Since vertical shear force is assumed to be
        # acting through the shear center, no twist occurs and therefore for both cells d_theta/d_y = 0
        # Solving the following matrix equation for redundant shear flows q_01 and q_02
        # | A11 A12 |  | q_01 | = | -B1 |
        # | A21 A22 |  | q_02 |   | -B2 |

        self.q_0 = []
        self.tor_q_0 = []
        self.G_dtheta_dy = []

        for span_idx in range(len(self.sec_spans)):
            A11 = self.l56[span_idx]/self.t56[span_idx] + self.l67[span_idx]/self.t67[span_idx] + \
                  self.l71[span_idx]/self.t71[span_idx] + self.l12[span_idx]/self.t12[span_idx] + \
                  self.l82[span_idx]/self.t82[span_idx] + self.l58[span_idx]/self.t58[span_idx]
            A12 = -1 * self.l82[span_idx]/self.t82[span_idx] - self.l58[span_idx]/self.t58[span_idx]
            A21 = -1 * self.l58[span_idx]/self.t58[span_idx] - self.l82[span_idx]/self.t82[span_idx]
            A22 = self.l23[span_idx]/self.t23[span_idx] + self.l34[span_idx]/self.t34[span_idx] + \
                  self.l45[span_idx]/self.t45[span_idx] + self.l58[span_idx]/self.t58[span_idx] + \
                  self.l82[span_idx]/self.t82[span_idx]

            B1 = self.q_b56_int(span_idx, self.l56[span_idx])/self.t56[span_idx] + \
                 self.q_b67_int(span_idx, self.l67[span_idx])/self.t67[span_idx] + \
                 self.q_b71_int(span_idx, self.l71[span_idx])/self.t71[span_idx] + \
                 self.q_b12_int(span_idx, self.l12[span_idx])/self.t12[span_idx] - \
                 self.q_b82_int(span_idx, self.xz2[span_idx][1])/self.t82[span_idx] - \
                 self.q_b58_int(span_idx, self.xz5[span_idx][1])/self.t58[span_idx]

            B2 = self.q_b23_int(span_idx, self.l23[span_idx])/self.t23[span_idx] + \
                 self.q_b34_int(span_idx, self.l34[span_idx])/self.t34[span_idx] + \
                 self.q_b45_int(span_idx, self.l45[span_idx])/self.t45[span_idx] + \
                 self.q_b82_int(span_idx, self.xz2[span_idx][1])/self.t82[span_idx] + \
                 self.q_b58_int(span_idx, self.xz5[span_idx][1])/self.t58[span_idx]

            A = np.array([[A11, A12], [A21, A22]])
            B = np.array([-B1, -B2])

            sol = np.linalg.solve(A, B)
            self.q_0.append([sol[0], sol[1]])

            # The same A matrix entries are used for calculating the cross-section unit shear flow due to
            # a 1 Nm torsion. Additionally torsional stiffness of the section is derived from this.
            A1 = self.As_cell1[span_idx]
            A2 = self.As_cell2[span_idx]
            A = np.array([[A1, A2, 0],
                          [A11/(2*A1), A12/(2*A1), -1],
                          [A21/(2*A2), A22/(2*A2), -1]])

            B = np.array([1, 0, 0])

            sol_tor = np.linalg.solve(A, B)
            self.G_dtheta_dy.append(sol_tor[2])
            self.tor_q_0.append([sol_tor[0], sol_tor[1]])
            # print('GJ:', 26900000000 * 1/sol_tor[2])

    # Shear flow distributions due to unit shear force Sz

    def q_82fz(self, span_idx, z):
        return self.q_b82(span_idx, z) - self.q_0[span_idx][0] + self.q_0[span_idx][1]

    def q_23fz(self, span_idx, s):
        return self.q_b23(span_idx, s) + self.q_0[span_idx][1]

    def q_34fz(self, span_idx, s):
        return self.q_b34(span_idx, s) + self.q_0[span_idx][1]

    def q_45fz(self, span_idx, s):
        return self.q_b45(span_idx, s) + self.q_0[span_idx][1]

    def q_58fz(self, span_idx, z):
        return self.q_b58(span_idx, z) - self.q_0[span_idx][0] + self.q_0[span_idx][1]

    def q_56fz(self, span_idx, s):
        return self.q_b56(span_idx, s) + self.q_0[span_idx][0]

    def q_67fz(self, span_idx, s):
        return self.q_b67(span_idx, s) + self.q_0[span_idx][0]

    def q_71fz(self, span_idx, s):
        return self.q_b71(span_idx, s) + self.q_0[span_idx][0]

    def q_12fz(self, span_idx, s):
        return self.q_b12(span_idx, s) + self.q_0[span_idx][0]

    def calc_shear_center(self):
        # Internal moment due to shear flows is evaluated around point 5, therefore the shear flows in
        # sections 45, 56, 58 and 82 can be ignored as they act through point 5.
        for i in range(len(self.sec_spans)):
            r_67 = self.xz5[i][0] - self.xz6[i][0]
            r_71 = r_67
            r_12 = np.abs(np.cross(self.xz2[i]-self.xz1[i], self.xz1[i]-self.xz5[i]))/np.linalg.norm(self.xz2[i]-self.xz1[i])
            r_23 = np.abs(np.cross(self.xz3[i]-self.xz2[i], self.xz2[i]-self.xz5[i]))/np.linalg.norm(self.xz3[i]-self.xz2[i])
            r_34 = self.xz4[i][0] - self.xz5[i][0]

            M67 = r_67 * sc.integrate.quad(lambda s: self.q_67fz(i, s), 0, self.l67[i])[0]
            M71 = r_71 * sc.integrate.quad(lambda s: self.q_71fz(i, s), 0, self.l71[i])[0]
            M12 = r_12 * sc.integrate.quad(lambda s: self.q_12fz(i, s), 0, self.l12[i])[0]
            M23 = r_23 * sc.integrate.quad(lambda s: self.q_23fz(i, s), 0, self.l23[i])[0]
            M34 = r_34 * sc.integrate.quad(lambda s: self.q_34fz(i, s), 0, self.l34[i])[0]

            zeta_5 = -(M67 + M71 + M12 + M23 + M34)

            shear_center = self.xz5[i][0] + zeta_5
            # print('shear center, [m]:', shear_center)
            self.shear_centers.append(shear_center)

    # Shear flow distributions due to unit shear
    # No redundant shear flow calculation present due to symmetry assumption and Sx acting through the symmetry axis

    def q_82fx(self, span_idx, z):
        qb = -1 * (self.t82[span_idx] / self.Izz[span_idx]) * (self.xz8[span_idx][0] * z - self.centroids[span_idx] * z)
        return qb

    def q_23fx(self, span_idx, s):
        t = self.t23[span_idx]
        dx = self.xz3[span_idx][0] - self.xz2[span_idx][0]
        q82_2 = self.q_82fx(span_idx, self.xz2[span_idx][1])
        q12_2 = self.q_12fx(span_idx, self.l12[span_idx])
        qb = -1 * (t / self.Izz[span_idx]) * (self.xz2[span_idx][0] * s + dx / (2 * self.l23[span_idx]) * s**2 - self.centroids[span_idx] * s) + q82_2 + q12_2
        return qb

    def q_34fx(self, span_idx, s):
        t = self.t34[span_idx]
        q23_3 = self.q_23fx(span_idx, self.l23[span_idx])
        qb = -1 * (t/self.Izz[span_idx]) * (self.xz3[span_idx][0] * s - self.centroids[span_idx] * s) + q23_3
        return qb

    def q_45fx(self, span_idx, s):
        t = self.t45[span_idx]
        q34_4 = self.q_34fx(span_idx, self.l34[span_idx])
        dx = self.xz5[span_idx][0] - self.xz4[span_idx][0]
        qb = -1 * (t/self.Izz[span_idx]) * (self.xz4[span_idx][0] * s + (dx/(2 * self.l45[span_idx])) * s**2 - self.centroids[span_idx] * s) + q34_4
        return qb

    def q_58fx(self, span_idx, z):
        qb = -1 * (self.t58[span_idx] / self.Izz[span_idx]) * (self.xz5[span_idx][0] * z - self.centroids[span_idx] * z)
        return qb

    def q_56fx(self, span_idx, s):
        t = self.t56[span_idx]
        dx = self.xz6[span_idx][0] - self.xz5[span_idx][0]
        q45_5 = self.q_45fx(span_idx, self.l45[span_idx])
        q58_5 = self.q_58fx(span_idx, self.xz5[span_idx][1])
        qb = -1 * (t/self.Izz[span_idx]) * (self.xz5[span_idx][0] * s + dx/(2 * self.l56[span_idx]) * s**2 - self.centroids[span_idx] * s) + q45_5 - q58_5
        return qb

    def q_67fx(self, span_idx, s):
        t = self.t67[span_idx]
        q56_6 = self.q_56fx(span_idx, self.l56[span_idx])
        qb = -1 * (t/self.Izz[span_idx]) * (self.xz6[span_idx][0] * s - self.centroids[span_idx] * s) + q56_6
        return qb

    def q_71fx(self, span_idx, s):
        t = self.t71[span_idx]
        qb = -1 * (t/self.Izz[span_idx]) * (self.xz7[span_idx][0] * s - self.centroids[span_idx] * s)
        return qb

    def q_12fx(self, span_idx, s):
        t = self.t12[span_idx]
        q71_1 = self.q_71fx(span_idx, self.l71[span_idx])
        dx = self.xz2[span_idx][0] - self.xz1[span_idx][0]
        qb = -1 * (t/self.Izz[span_idx]) * (self.xz1[span_idx][0] * s + dx/(2 * self.l12[span_idx]) * s**2 - self.centroids[span_idx] * s) + q71_1
        return qb

    def plot_wing(self):
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.set_xlim3d([0, self.half_span*2])
        ax.set_ylim3d([0, self.half_span*2])
        ax.set_zlim3d([0, self.half_span*2])

        for idx, sec in enumerate(self.w_coords):
            x = sec[0]
            z = sec[1]
            y = np.repeat(self.sec_spans[idx], self.nr_points)
            ax.plot(x, y, z, c='blue')

        for idx, sec_spar in enumerate(self.spar_coords):
            x = sec_spar[0]
            z = sec_spar[1]
            y = np.repeat(self.sec_spans[idx], self.nr_spars)
            for idx, x_spar in enumerate(x):
                z_spar = z[idx]
                y_spar = y[idx]
                ax.plot([x_spar, x_spar], [y_spar, y_spar], [z_spar, -1 * z_spar])

        for idx, sec_spar in enumerate(self.spar_coords):
            x = sec_spar[0]
            z = sec_spar[1]
            y = np.repeat(self.sec_spans[idx], self.nr_spars)

            if idx < self.nr_bsec-1:
                x_next = self.spar_coords[idx+1][0]
                z_next = self.spar_coords[idx + 1][1]
                y_next = np.repeat(self.sec_spans[idx+1], self.nr_spars)
                for idx, x_spar in enumerate(x):
                    z_spar = z[idx]
                    y_spar = y[idx]
                    ax.plot([x_spar, x_next[idx]], [y_spar, y_next[idx]], [z_spar, z_next[idx]], c='r')
                    ax.plot([x_spar, x_next[idx]], [y_spar, y_next[idx]], [-1 * z_spar, -1 * z_next[idx]], c='r')

            for idx, x_spar in enumerate(x):
                z_spar = z[idx]
                y_spar = y[idx]
                ax.plot([x_spar, x_spar], [y_spar, y_spar], [z_spar, -1 * z_spar], c='r')

        plt.title('3D view of wing discretization')
        plt.xlabel('Chord, [m]')
        plt.ylabel('Halfspan, [m]')
        plt.show()

    def plot_SMOA_dist(self):
        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_size_inches(13, 6)
        ax1.set_title('Ixx over span')
        ax1.set_xlabel('Wing halfspan, [m]')
        ax1.set_ylabel('Second moment of area, [m^4]')
        ax2.set_title('Izz over span')
        ax2.set_xlabel('Wing halfspan, [m]')
        ax2.set_ylabel('Second moment of area, [m^4]')
        ax1.plot(self.sec_spans, self.Ixx)
        ax2.plot(self.sec_spans, self.Izz)
        plt.show()

    def check_SMOA(self):
        # Function to verify that the centroid location. For this, put LE and TE thicknesses to 0. Use two spars
        # sufficiently distanced with a thickness of two.
        x_spar1 = self.spar_coords[0][0][0]
        x_spar2 = self.spar_coords[0][0][1]
        z_spar1 = self.spar_coords[0][1][0]
        z_spar2 = self.spar_coords[0][1][1]

        dx = x_spar2 - x_spar1
        dz = z_spar2 - z_spar1

        beta = np.arctan(dz / dx)
        l = np.sqrt(dx**2 + dz**2)
        t = 2/1000

        centroid_top_x = x_spar1 + dx / 2
        centroid_top_z = z_spar1 + dz / 2

        a_spar1 = 2 * z_spar1 * 2/1000
        a_spar2 = 2 * z_spar2 * 2/1000
        a_top = l * 2 / 1000
        a_tot = a_spar1 + a_spar2 + 2 * a_top

        centr = (a_spar1 * x_spar1 + a_spar2 * x_spar2 + 2 * centroid_top_x * a_top)/a_tot

        print('-----------------------------------------------------------------------------')
        print('Model calculated centroid location, [m]:', self.centroids[0])
        print('Verification estimate centroid location, [m]:', centr)
        print('Percentage difference wrt ver. est. [%]:', self.centroids[0]/centr*100 - 100)

        Ixx = 2 * ((l ** 3 * t * np.sin(beta) ** 2) / 12 + a_top * centroid_top_z**2) + (t*(2*z_spar1)**3)/12 + (t*(2*z_spar2)**3)/12
        Izz = 2 * ((l ** 3 * t * np.cos(beta) ** 2) / 12 + a_top * (centr - centroid_top_x) ** 2) + a_spar1 * (centr - x_spar1) ** 2 + a_spar2 * (centr - x_spar2) ** 2

        print('-----------------------------------------------------------------------------')
        print('Model calculated Ixx, [m]:', self.Ixx[0])
        print('Verification estimate Ixx, [m]:', Ixx)
        print('Percentage difference wrt ver. est. [%]:', self.Ixx[0] / Ixx * 100 - 100)

        print('-----------------------------------------------------------------------------')
        print('Model calculated Izz, [m]:', self.Izz[0])
        print('Verification estimate Izz, [m]:', Izz)
        print('Percentage difference wrt ver. est. [%]:', self.Izz[0] / Izz * 100 - 100)

    def plot_wing_structure(self, x_c1, x_c2, x_c3):
        opc = 0.5
        opc_skin = 0.4

        color_le = ['magenta', 'grey', 'grey', 'yellow', 'green', 'green']  # Leading edge thickness, [mm]
        color_top = ['red', 'red', 'magenta', 'magenta', 'yellow', 'yellow']

        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.set_xlim3d([0, self.half_span * 2])
        ax.set_ylim3d([0, self.half_span * 2])
        ax.set_zlim3d([0, self.half_span * 2])
        ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
        # make the grid lines transparent
        ax.xaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
        ax.yaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
        ax.zaxis._axinfo["grid"]['color'] = (1, 1, 1, 0)
        ax.set_axis_off()

        thick_red = np.linspace(1, self.t_red_tip, num=self.nr_bsec, endpoint=True)

        t_red_func = sc.interpolate.interp1d(self.sec_spans, thick_red, kind='linear')
        chord_func = sc.interpolate.interp1d(self.sec_spans, self.sec_chords, kind='linear')

        z_c1 = self.airfoil_inter(x_c1)
        z_c2 = self.airfoil_inter(x_c2)
        z_c3 = self.airfoil_inter(x_c3)

        x_r_1 = x_c1 * self.c_root
        x_r_2 = x_c2 * self.c_root
        x_r_3 = x_c3 * self.c_root

        z_r_1 = z_c1 * self.c_root
        z_r_2 = z_c2 * self.c_root
        z_r_3 = z_c3 * self.c_root

        plt.scatter(0, 0, 0, color='red', label='2.24 mm')
        plt.scatter(0, 0, 0, color='magenta', label='1.92 mm')
        plt.scatter(0, 0, 0, color='grey', label='1.76 mm')
        plt.scatter(0, 0, 0, color='yellow', label='1.60 mm')
        plt.scatter(0, 0, 0, color='green', label='1.28 mm')

        dx_str_top_r = (x_r_3 - x_r_2) / (self.nr_str_top - 1)
        dz_str_top_r = (z_r_3 - z_r_2) / (self.nr_str_top - 1)
        beta_top = np.arctan(dz_str_top_r / dx_str_top_r)

        dx_str_LE_top_r = (x_r_2 - x_r_1) / (self.nr_str_LE_top - 1)
        dz_str_LE_top_r = (z_r_2 - z_r_1) / (self.nr_str_LE_top - 1)
        beta_LE_top = np.arctan(dz_str_LE_top_r / dx_str_LE_top_r)

        sec_str_top_coor_r = []
        str_p1r = []
        str_p3r = []
        str_p4r = []
        for i in range(self.nr_str_top):
            x_str = x_r_2 + i * dx_str_top_r
            z_str = z_r_2 + i * dz_str_top_r
            x_str3 = x_str + self.h_str * np.sin(beta_top)
            z_str3 = z_str - self.h_str * np.cos(beta_top)
            sec_str_top_coor_r.append([x_str, z_str])
            str_p1r.append([x_str - self.w_str * np.cos(beta_top), z_str - self.w_str * np.sin(beta_top)])
            str_p3r.append([x_str3, z_str3])
            str_p4r.append([x_str3 + self.w_str * np.cos(beta_top), z_str3 + self.w_str * np.sin(beta_top)])

        sec_str_LE_top_coor_r = []
        str_LE_p1r = []
        str_LE_p3r = []
        str_LE_p4r = []
        for i in range(self.nr_str_LE_top):
            x_str = x_r_1 + i * dx_str_LE_top_r
            z_str = z_r_1 + i * dz_str_LE_top_r
            x_str3 = x_str + self.h_str * np.sin(beta_top)
            z_str3 = z_str - self.h_str * np.cos(beta_top)
            sec_str_LE_top_coor_r.append([x_str, z_str])
            str_LE_p1r.append([x_str - self.w_str * np.cos(beta_top), z_str + self.w_str * np.sin(beta_top)])
            str_LE_p3r.append([x_str3, z_str3])
            str_LE_p4r.append([x_str3 + self.w_str * np.cos(beta_top), z_str3 + self.w_str * np.sin(beta_top)])

            # sec_str_LE_top_coor_r.append([x_r_1 + i * dx_str_LE_top_r, z_r_1 + i * dz_str_LE_top_r])

        for idx, y_rib in enumerate(self.ribs_loc):
            # rib_chords.append(chord_func(y_rib))
            c = chord_func(y_rib)
            red = t_red_func(y_rib)
            x, z = self.scale_airfoil(c, red)
            y = np.full(len(x), y_rib)
            plt.plot(x, y, z, c='blue', linewidth=1)

            z1 = z_c1 * c * red
            z2 = z_c2 * c * red
            z3 = z_c3 * c * red
            x1 = x_c1 * c
            x2 = x_c2 * c
            x3 = x_c3 * c

            if idx > 0:
                y_prev = self.ribs_loc[idx - 1]
                c_prev = chord_func(y_prev)
                red_prev = t_red_func(y_prev)
                z1_prev = z_c1 * c_prev * red_prev
                z2_prev = z_c2 * c_prev * red_prev
                z3_prev = z_c3 * c_prev * red_prev
                x1_prev = x_c1 * c_prev
                x2_prev = x_c2 * c_prev
                x3_prev = x_c3 * c_prev

                ax.plot_trisurf(np.array([x1, x2, x1_prev]), np.array([y_rib, y_rib, y_prev]), np.array([z1, z2, z1_prev]), color=color_le[idx - 1], alpha=opc_skin)
                ax.plot_trisurf(np.array([x2, x1_prev, x2_prev]), np.array([y_rib, y_prev, y_prev]), np.array([z2, z1_prev, z2_prev]), color=color_le[idx - 1], alpha=opc_skin)
                ax.plot_trisurf(np.array([x1, x2, x1_prev]), np.array([y_rib, y_rib, y_prev]), np.array([-z1, -z2, -z1_prev]), color=color_le[idx - 1], alpha=opc_skin)
                ax.plot_trisurf(np.array([x2, x1_prev, x2_prev]), np.array([y_rib, y_prev, y_prev]), np.array([-z2, -z1_prev, -z2_prev]), color=color_le[idx - 1], alpha=opc_skin)

                ax.plot_trisurf(np.array([x2, x3, x2_prev]), np.array([y_rib, y_rib, y_prev]), np.array([z2, z3, z2_prev]), color=color_top[idx - 1], alpha=opc_skin)
                ax.plot_trisurf(np.array([x3, x2_prev, x3_prev]), np.array([y_rib, y_prev, y_prev]), np.array([z3, z2_prev, z3_prev]), color=color_top[idx - 1], alpha=opc_skin)
                ax.plot_trisurf(np.array([x2, x3, x2_prev]), np.array([y_rib, y_rib, y_prev]), np.array([-z2, -z3, -z2_prev]), color=color_top[idx - 1], alpha=opc_skin)
                ax.plot_trisurf(np.array([x3, x2_prev, x3_prev]), np.array([y_rib, y_prev, y_prev]), np.array([-z3, -z2_prev, -z3_prev]), color=color_top[idx - 1], alpha=opc_skin)

                ax.plot_trisurf(np.array([x2, x2/0.999999, x2_prev]), np.array([y_rib, y_rib, y_prev]), np.array([z2, -z2, z2_prev]), color='blue', alpha=opc_skin)
                ax.plot_trisurf(np.array([x2, x2_prev, x2_prev/0.999999]), np.array([y_rib, y_prev, y_prev]), np.array([-z2, z2_prev, -z2_prev]), color='blue', alpha=opc_skin)

                ax.plot_trisurf(np.array([x3, x3 / 0.999999, x3_prev]), np.array([y_rib, y_rib, y_prev]), np.array([z3, -z3, z3_prev]), color='blue', alpha=opc_skin)
                ax.plot_trisurf(np.array([x3, x3_prev, x3_prev / 0.999999]), np.array([y_rib, y_prev, y_prev]), np.array([-z3, z3_prev, -z3_prev]), color='blue', alpha=opc_skin)

                ax.plot_trisurf(np.array([x1, x1 / 0.999999, x1_prev]), np.array([y_rib, y_rib, y_prev]), np.array([z1, -z1, z1_prev]), color='blue', alpha=opc_skin)
                ax.plot_trisurf(np.array([x1, x1_prev, x1_prev / 0.999999]), np.array([y_rib, y_prev, y_prev]), np.array([-z1, z1_prev, -z1_prev]), color='blue', alpha=opc_skin)
                #ax.plot_trisurf(np.array([x2, x3, x2_prev]), np.array([y_rib, y_rib, y_prev]), np.array([-z2, -z3, -z2_prev]), color=color_le[idx - 1], alpha=opc_skin)
                #ax.plot_trisurf(np.array([x3, x2_prev, x3_prev]), np.array([y_rib, y_prev, y_prev]), np.array([-z3, -z2_prev, -z3_prev]), color=color_le[idx - 1], alpha=opc_skin)

                plt.plot([x1, x1_prev], [y_rib, y_prev], [z1, z1_prev], c='black', linewidth=3)
                plt.plot([x2, x2_prev], [y_rib, y_prev], [z2, z2_prev], c='black', linewidth=3)
                plt.plot([x3, x3_prev], [y_rib, y_prev], [z3, z3_prev], c='black', linewidth=3)
                plt.plot([x1, x1_prev], [y_rib, y_prev], [-z1, -z1_prev], c='black', linewidth=3)
                plt.plot([x2, x2_prev], [y_rib, y_prev], [-z2, -z2_prev], c='black', linewidth=3)
                plt.plot([x3, x3_prev], [y_rib, y_prev], [-z3, -z3_prev], c='black', linewidth=3)

            plt.plot([x1, x1], [y_rib, y_rib], [z1, -z1], c='black')
            plt.plot([x2, x2], [y_rib, y_rib], [z2, -z2], c='black')
            plt.plot([x3, x3], [y_rib, y_rib], [z3, -z3], c='black')

            dx_str_top = (x3 - x2) / (self.nr_str_top - 1)
            dz_str_top = (z3 - z2) / (self.nr_str_top - 1)
            beta_top = np.arctan(dz_str_top / dx_str_top)

            dx_str_LE_top = (x2 - x1) / (self.nr_str_LE_top - 1)
            dz_str_LE_top = (z2 - z1) / (self.nr_str_LE_top - 1)
            beta_LE_top = np.arctan(dz_str_LE_top / dx_str_LE_top)

            for i in range(self.nr_str_top):
                if y_rib == self.str_top_end[i]:
                    x_str_2 = x2 + i * dx_str_top
                    z_str_2 = z2 + i * dz_str_top
                    x_str_1 = x_str_2 - self.w_str * np.cos(beta_top)
                    z_str_1 = z_str_2 - self.w_str * np.sin(beta_top)
                    x_str_3 = x_str_2 + self.h_str * np.sin(beta_top)
                    z_str_3 = z_str_2 - self.h_str * np.cos(beta_top)
                    x_str_4 = x_str_2 + self.w_str * np.cos(beta_top)
                    z_str_4 = z_str_3 + self.w_str * np.sin(beta_top)

                    x_str_1r = str_p1r[i][0]
                    z_str_1r = str_p1r[i][1]
                    x_str_2r = sec_str_top_coor_r[i][0]
                    z_str_2r = sec_str_top_coor_r[i][1]
                    x_str_3r = str_p3r[i][0]
                    z_str_3r = str_p3r[i][1]
                    x_str_4r = str_p4r[i][0]
                    z_str_4r = str_p4r[i][1]

                    if not i==0 and not i== self.nr_str_top - 1:
                        # ax.plot_trisurf(np.array([x_str_1r, x_str_1, x_str_2r]), np.array([0, y_rib, 0]), np.array([-z_str_1r, -z_str_1, -z_str_2r]), color='black', alpha=opc)
                        # ax.plot_trisurf(np.array([x_str_1, x_str_2, x_str_2r]), np.array([y_rib, y_rib, 0]), np.array([-z_str_1, -z_str_2, -z_str_1r]), color='black', alpha=opc)
                        # ax.plot_trisurf(np.array([x_str_2r, x_str_2, x_str_3r]), np.array([0, y_rib, 0]), np.array([-z_str_2r, -z_str_2, -z_str_3r]), color='black', alpha=opc)
                        # ax.plot_trisurf(np.array([x_str_2, x_str_3, x_str_3r]), np.array([y_rib, y_rib, 0]), np.array([-z_str_2, -z_str_3, -z_str_3r]), color='black', alpha=opc)
                        # ax.plot_trisurf(np.array([x_str_3r, x_str_3, x_str_4r]), np.array([0, y_rib, 0]), np.array([-z_str_3r, -z_str_3, -z_str_4r]), color='black', alpha=opc)
                        # ax.plot_trisurf(np.array([x_str_3, x_str_4, x_str_4r]), np.array([y_rib, y_rib, 0]), np.array([-z_str_3, -z_str_4, -z_str_4r]), color='black', alpha=opc)

                        ax.plot_trisurf(np.array([x_str_1r, x_str_1, x_str_2r]), np.array([0, y_rib, 0]), np.array([z_str_1r, z_str_1, z_str_2r]), color='black', alpha=opc)
                        ax.plot_trisurf(np.array([x_str_1, x_str_2, x_str_2r]), np.array([y_rib, y_rib, 0]), np.array([z_str_1, z_str_2, z_str_1r]), color='black', alpha=opc)
                        ax.plot_trisurf(np.array([x_str_2r, x_str_2, x_str_3r]), np.array([0, y_rib, 0]), np.array([z_str_2r, z_str_2, z_str_3r]), color='black', alpha=opc)
                        ax.plot_trisurf(np.array([x_str_2, x_str_3, x_str_3r]), np.array([y_rib, y_rib, 0]), np.array([z_str_2, z_str_3, z_str_3r]), color='black', alpha=opc)
                        ax.plot_trisurf(np.array([x_str_3r, x_str_3, x_str_4r]), np.array([0, y_rib, 0]), np.array([z_str_3r, z_str_3, z_str_4r]), color='black', alpha=opc)
                        ax.plot_trisurf(np.array([x_str_3, x_str_4, x_str_4r]), np.array([y_rib, y_rib, 0]), np.array([z_str_3, z_str_4, z_str_4r]), color='black', alpha=opc)

                    plt.plot([x_str_1r, x_str_1], [0, y_rib], [z_str_1r, z_str_1], c='black', linewidth=1)
                    plt.plot([x_str_2r, x_str_2], [0, y_rib], [z_str_2r, z_str_2], c='black', linewidth=1)
                    plt.plot([x_str_3r, x_str_3], [0, y_rib], [z_str_3r, z_str_3], c='black', linewidth=1)
                    plt.plot([x_str_4r, x_str_4], [0, y_rib], [z_str_4r, z_str_4], c='black', linewidth=1)

                    # plt.plot([x_str_1r, x_str_1], [0, y_rib], [-z_str_1r, -z_str_1], c='black', linewidth=1)
                    plt.plot([x_str_2r, x_str_2], [0, y_rib], [-z_str_2r, -z_str_2], c='black', linewidth=1)
                    # plt.plot([x_str_3r, x_str_3], [0, y_rib], [-z_str_3r, -z_str_3], c='black', linewidth=1)
                    # plt.plot([x_str_4r, x_str_4], [0, y_rib], [-z_str_4r, -z_str_4], c='black', linewidth=1)

            for i in range(self.nr_str_LE_top):
                if y_rib == self.str_LE_top_end[i]:
                    x_str_LE2 = x1 + i * dx_str_LE_top
                    z_str_LE2 = z1 + i * dz_str_LE_top
                    x_str_LE1 = x_str_LE2 - self.w_str * np.cos(beta_LE_top)
                    z_str_LE1 = z_str_LE2 - self.w_str * np.sin(beta_LE_top)
                    x_str_LE3 = x_str_LE2 + self.h_str * np.sin(beta_LE_top)
                    z_str_LE3 = z_str_LE2 - self.h_str * np.cos(beta_LE_top)
                    x_str_LE4 = x_str_LE3 + self.w_str * np.cos(beta_LE_top)
                    z_str_LE4 = z_str_LE3 + self.w_str * np.sin(beta_LE_top)

                    x_str_LE1r = str_LE_p1r[i][0]
                    z_str_LE1r = str_LE_p1r[i][1]
                    x_str_LE2r = sec_str_LE_top_coor_r[i][0]
                    z_str_LE2r = sec_str_LE_top_coor_r[i][1]
                    x_str_LE3r = str_LE_p3r[i][0]
                    z_str_LE3r = str_LE_p3r[i][1]
                    x_str_LE4r = str_LE_p4r[i][0]
                    z_str_LE4r = str_LE_p4r[i][1]

                    plt.plot([x_str_LE2r, x_str_LE2], [0, y_rib], [z_str_LE2r, z_str_LE2], c='black')
                    plt.plot([x_str_LE2r, x_str_LE2], [0, y_rib], [-z_str_LE2r, -z_str_LE2], c='black')

                    if not i==0 and not i== self.nr_str_LE_top - 1:
                        ax.plot_trisurf(np.array([x_str_LE1r, x_str_LE1, x_str_LE2r]), np.array([0, y_rib, 0]), np.array([z_str_LE1r, z_str_LE1, z_str_LE2r]), color='black', alpha=opc)
                        ax.plot_trisurf(np.array([x_str_LE1, x_str_LE2, x_str_LE2r]), np.array([y_rib, y_rib, 0]), np.array([z_str_LE1, z_str_LE2, z_str_LE1r]), color='black', alpha=opc)
                        ax.plot_trisurf(np.array([x_str_LE2r, x_str_LE2, x_str_LE3r]), np.array([0, y_rib, 0]), np.array([z_str_LE2r, z_str_LE2, z_str_LE3r]), color='black', alpha=opc)
                        ax.plot_trisurf(np.array([x_str_LE2, x_str_LE3, x_str_LE3r]), np.array([y_rib, y_rib, 0]), np.array([z_str_LE2, z_str_LE3, z_str_LE3r]), color='black', alpha=opc)
                        ax.plot_trisurf(np.array([x_str_LE3r, x_str_LE3, x_str_LE4r]), np.array([0, y_rib, 0]), np.array([z_str_LE3r, z_str_LE3, z_str_LE4r]), color='black', alpha=opc)
                        ax.plot_trisurf(np.array([x_str_LE3, x_str_LE4, x_str_LE4r]), np.array([y_rib, y_rib, 0]), np.array([z_str_LE3, z_str_LE4, z_str_LE4r]), color='black', alpha=opc)

                        # ax.plot_trisurf(np.array([x_str_LE1r, x_str_LE1, x_str_LE2r]), np.array([0, y_rib, 0]), np.array([-z_str_LE1r, -z_str_LE1, -z_str_LE2r]), color='black', alpha=opc)
                        # ax.plot_trisurf(np.array([x_str_LE1, x_str_LE2, x_str_LE2r]), np.array([y_rib, y_rib, 0]), np.array([-z_str_LE1, -z_str_LE2, -z_str_LE1r]), color='black', alpha=opc)
                        # ax.plot_trisurf(np.array([x_str_LE2r, x_str_LE2, x_str_LE3r]), np.array([0, y_rib, 0]), np.array([-z_str_LE2r, -z_str_LE2, -z_str_LE3r]), color='black', alpha=opc)
                        # ax.plot_trisurf(np.array([x_str_LE2, x_str_LE3, x_str_LE3r]), np.array([y_rib, y_rib, 0]), np.array([-z_str_LE2, -z_str_LE3, -z_str_LE3r]), color='black', alpha=opc)
                        # ax.plot_trisurf(np.array([x_str_LE3r, x_str_LE3, x_str_LE4r]), np.array([0, y_rib, 0]), np.array([-z_str_LE3r, -z_str_LE3, -z_str_LE4r]), color='black', alpha=opc)
                        # ax.plot_trisurf(np.array([x_str_LE3, x_str_LE4, x_str_LE4r]), np.array([y_rib, y_rib, 0]), np.array([-z_str_LE3, -z_str_LE4, -z_str_LE4r]), color='black', alpha=opc)

                        plt.plot([x_str_LE1r, x_str_LE1], [0, y_rib], [z_str_LE1r, z_str_LE1], c='black', linewidth=1)
                        plt.plot([x_str_LE2r, x_str_LE2], [0, y_rib], [z_str_LE2r, z_str_LE2], c='black', linewidth=1)
                        plt.plot([x_str_LE3r, x_str_LE3], [0, y_rib], [z_str_LE3r, z_str_LE3], c='black', linewidth=1)
                        plt.plot([x_str_LE4r, x_str_LE4], [0, y_rib], [z_str_LE4r, z_str_LE4], c='black', linewidth=1)

                        # plt.plot([x_str_LE1r, x_str_LE1], [0, y_rib], [-z_str_LE1r, -z_str_LE1], c='black', linewidth=1)
                        plt.plot([x_str_LE2r, x_str_LE2], [0, y_rib], [-z_str_LE2r, -z_str_LE2], c='black', linewidth=1)
                        # plt.plot([x_str_LE3r, x_str_LE3], [0, y_rib], [-z_str_LE3r, -z_str_LE3], c='black', linewidth=1)
                        # plt.plot([x_str_LE4r, x_str_LE4], [0, y_rib], [-z_str_LE4r, -z_str_LE4], c='black', linewidth=1)

        plt.legend()
        plt.show()

# WP = WingProperties(DIRECTORY, INPUT_FOLDER, airfoil_file, c_root, c_tip, MAC, dihedral,
#                     LE_sweep, half_span, AR, nr_span_sections, nr_spars, spar_locs, t_le,
#                     t_top, t_tr, t_bot, t_spar, x_le, x_te, t_sp_flange, w_sp_flange, t_red_tip)


# WP.plot_airfoil()
# WP.plot_wing()
# WP.plot_SMOA_dist()
# WP.check_SMOA()





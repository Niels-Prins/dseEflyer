# Lennart Krieg
# DSE group 20: Aerobatic E-Flyer
# Wing Discretization model, part of wing structural model

import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate
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

# Thicknesses
t_le = 1  # Leading edge thickness, [mm]
t_top = 2  # Top skin thickness, [mm]
t_tr = 1  # Trailing edge thickness, [mm]
t_bot = 2  # Bottom skin thickness, [mm]
t_spar = [4, 4]  # Spar thicknesses, [mm]. Input must be list of same size as spar_locs list!
t_sp_flange = [3, 3]
w_sp_flange = [40, 40]

x_le = 0.3  # Fraction of chord where LE thickness becomes top or bottom skin thickness
x_te = 0.75  # Fraction of chord where top or bottom skin thickness become TE thickness

dens_spar = 2800  # Spar material density, [kg/m^3]
dens_skins = 2800  # Skin material density, [kg/m^3]

class WingProperties:
    def __init__(self, DIRECTORY, INPUT_FOLDER, airfoil_file, c_root, c_tip, MAC, dihedral, LE_sweep,
                 half_span, AR, span_sections, nr_spar, spar_loc, t_le, t_top, t_tr, t_bot, t_spar,
                 x_le, x_te, t_sp_flange, w_sp_flange, t_red_tip, dens_spars, dens_skins):

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

        # Spars are defined based on start end coordinate, as the airfoil might be unsymmetrical.
        self.nr_spars = nr_spar  # Number of required spars, [-]
        self.spar_loc = spar_loc  # Fraction of chord at which spar(s) is/are required. Input given as list.
        self.spars_x = None
        self.spars_z = None
        self.spar_coords = None  # Coordinates of the start and end points of the spars. Declared to be calculated later
        self.spar_indices = None  # Indices of x coordinates at which spars are located.

        # Thicknesses
        self.t_le = t_le/1000  # Leading edge thickness, [m]
        self.t_top = t_top/1000  # Top skin thickness, [m]
        self.t_te = t_tr/1000  # Trailing edge thickness, [m]
        self.t_bot = t_bot/1000  # Bottom skin thickness, [m]
        self.t_spar = np.array(t_spar)/1000  # Spar thicknesses, [m]
        self.t_sp_flange = np.array(t_sp_flange)/1000  # Spars flange thickness, [m]
        self.w_sp_flange = np.array(w_sp_flange)/1000  # Spars flange width, [m]
        self.t_distr = None  # Array with

        self.x_le = x_le  # Fraction of chord where LE thickness becomes top or bottom skin thickness
        self.x_te = x_te  # Fraction of chord where top or bottom skin thickness become TE thickness

        self.spar_centroids = None
        self.spar_areas = None
        self.skin_centroids = None
        self.skin_areas = None
        self.seg_len = None

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

        self.print_initials()
        self.add_spars()
        self.discretize_wing()
        self.calc_centroid()
        self.calc_SMOA()
        self.interpolate_SMOAs()
        self.est_mass(dens_spars, dens_skins)
        self.plot_airfoil()
        self.plot_wing()

    def print_initials(self):
        print('------------------------------------------------------------------------------')
        print('-------------------------- Wing geometry properties --------------------------')
        print('Imported airfoil geometry:', self.airfoil_name)
        print('Wingspan, [m]:', self.half_span * 2)
        print('Root chord, [m]:', self.c_root)
        print('Tip chord, [m]:', self.c_tip)
        print('------------------------- Wing structural components -------------------------')
        print('Nr. of spars:', self.nr_spars)
        print('Spar locations, [% of c]:', self.spar_loc)
        print('Spar web thicknesses, [mm]:', self.t_spar*1000)
        print('Spar flange thicknesses, [mm]:', self.t_sp_flange*1000)
        print('Spar flange widths, [mm]:', self.w_sp_flange * 1000)
        print('Leading edge skin thickness, [mm]:', self.t_le*1000)
        print('Trailing edge skin thickness, [mm]:', self.t_le * 1000)
        print('Trailing edge skin thickness, [mm]:', self.t_le * 1000)
        print('Top center skin thickness, [mm]:', self.t_top * 1000)
        print('Top bottom skin thickness, [mm]:', self.t_bot * 1000)
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
                        t = self.t_le
                    elif x_cor > chord * self.x_te:
                        t = self.t_te
                    else:
                        if z_cor >= 0:
                            t = self.t_top
                        else:
                            t = self.t_bot

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
        self.centroid_dist = sc.interpolate.interp1d(self.sec_spans, self.centroids, kind='linear')

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

    def interpolate_SMOAs(self):
        self.Ixx_dist = sc.interpolate.interp1d(self.sec_spans, self.Ixx, kind='cubic')
        self.Izz_dist = sc.interpolate.interp1d(self.sec_spans, self.Izz, kind='cubic')

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

    def plot_wing(self):
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        ax.set_xlim3d([0, self.half_span])
        ax.set_ylim3d([0, self.half_span])
        ax.set_zlim3d([0, self.half_span])

        for idx, sec in enumerate(self.w_coords):
            x = sec[0]
            z = sec[1]
            y = np.repeat(self.sec_spans[idx], self.nr_points)
            ax.plot(x, y, z)

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


# WP = WingProperties(DIRECTORY, INPUT_FOLDER, airfoil_file, c_root, c_tip, MAC, dihedral,
#                     LE_sweep, half_span, AR, nr_span_sections, nr_spars, spar_locs, t_le,
#                     t_top, t_tr, t_bot, t_spar, x_le, x_te, t_sp_flange, w_sp_flange, t_red_tip)


# WP.plot_airfoil()
# WP.plot_wing()
# WP.plot_SMOA_dist()
# WP.check_SMOA()





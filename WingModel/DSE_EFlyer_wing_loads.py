# Lennart Krieg
# DSE group 20: Aerobatic E-Flyer
# Wing load model, part of wing structural model

import numpy as np
import matplotlib.pyplot as plt
import scipy as sc
import scipy.interpolate
from scipy import integrate
from mpl_toolkits import mplot3d
import pandas as pd

#  Pycharm Project folder: C:\Users\lenna\PycharmProjects\DSE_22

DIRECTORY = 'C:/Users/lenna/Documents/DSE/Wing_structural_model'
INPUT_FOLDER = '/Input_Files/'

# Airfoil files with coordinates in x/z fractions of chord.
# Before use these files will have to be preprocessed to have them start at the LE, go to TE and wrap around again
# to the LE. Text and other non-data will have to be removed and separation by spaces replaced by semicolon (;).
CL_dist_file = 'CL_dist_apha_clean_V2.csv'
CL_alpha = ''
CD_dist_file = 'CD_dist_alpha_clean_V2.csv'
CM_dist_file = 'CM_dist_alpha_clean.csv'

# Wing geometry parameters
c_root = 2.01  # Root chord, [m]
c_tip = 0.98  # Tip chord, [m]
MAC = 1.58  # Mean aerodynamic chord, [m]
S = 14.536  # Total wing area, [m^2]

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

nr_span_sections = 20  # Number of cross-sections in which the wing is split up, [-]
ac_mass = 1063  # Aircraft manoeuvring mass, [kg]
g_load = 8  # Aircraft manoeuvring g load,[-]

rho = 1.225  # Air density, [kg/m^3]
V = 150  # Maneuvering speed, [kts]


class WingLoading:
    def __init__(self, DIRECTORY, INPUT_FOLDER, CL_dist_file, CD_dist_file, c_root, c_tip, MAC, dihedral, LE_sweep,
                 half_span, AR, span_sections, ac_mass, g_load, CL_0, CL_10, S, rho, V, a, b, c,
                 CM_dist_file, COP_dist_file):

        self.CL_distributions = pd.read_csv(DIRECTORY+INPUT_FOLDER+CL_dist_file, delimiter=';')
        self.CD_distributions = pd.read_csv(DIRECTORY + INPUT_FOLDER + CD_dist_file, delimiter=';')
        self.CM_distributions = pd.read_csv(DIRECTORY + INPUT_FOLDER + CM_dist_file, delimiter=';')
        self.COP_distributions = pd.read_csv(DIRECTORY + INPUT_FOLDER + COP_dist_file, delimiter=';')

        self.c_root = c_root  # Root chord, [m]
        self.c_tip = c_tip  # Tip chord, [m]
        self.MAC = MAC  # Mean aerodynamic chord, [m]

        self.dihedral = dihedral  # Wing dihedral, [deg]
        self.LE_sweep = LE_sweep  # Leading egde sweep, [deg]
        self.half_span = half_span  # Halfspan of the wing, [m]
        self.AR = AR  # Aspect ratio, [-]

        self.nr_bsec = span_sections  # Number of cross-sections in which the wing is split up, [-]
        self.sec_spans = None  # Spans of discretized cross-sections. Declared here to be calculated later
        self.sec_chords = None  # chord at these discretized cross-sections. Declared here to be calculated later

        self.ac_mass = ac_mass
        self.g_load = g_load
        self.V = V * 0.515  # Flight speed, [m/s]
        self.rho = rho  # Air density, [kg/m^3]
        self.S = S  # Total wing area (both wing halves), [m^2]

        self.CL_dist = None
        self.CL_0 = CL_0
        self.CL_10 = CL_10
        self.CL_alpha = (CL_10 - CL_0)/10
        self.CL_req = ac_mass * 9.81 * g_load / (0.5 * rho * S * (V*0.515)**2)
        self.alpha_req = (self.CL_req - self.CL_0)/self.CL_alpha
        print('Required alpha at 150 kt', self.alpha_req)
        self.a = a
        self.b = b
        self.c = c

        self.CD_dist = None
        self.CM_dist = None
        self.COP_dist = None

        self.shear_dist_z = None
        self.shear_dist_x = None

        self.shear_func_z = None
        self.shear_func_x = None

        self.aero_moment_dist = None
        self.aero_moment_func = None

        self.moment_dist_x = None
        self.moment_dist_z = None
        # Temporary quarter chord center of pressure assumption until proper data is imported
        self.cps = [0.25 for i in range(self.nr_bsec)]

        self.discretize_wing()

        if self.alpha_req > 17.84:
            self.extrapolate_CL_dist(self.alpha_req)
        else:
            self.select_CL_dist(self.alpha_req)

        self.plot_decomposed_forces()
        self.plot_shear_forces()
        self.plot_bending_moments()
        self.plot_COP_dist()

    def discretize_wing(self):
        self.sec_spans = np.linspace(0, self.half_span, num=self.nr_bsec,
                                     endpoint=True)  # Spans of the cross-sections
        chord_deltas = np.linspace(0, self.c_tip - self.c_root, num=self.nr_bsec,
                                   endpoint=True)  # Delta chord wrt, root at these cross-sections
        self.sec_chords = chord_deltas + self.c_root

    def select_CL_dist(self, alpha):
        print('----- Plotting CL distribution over halfspan -----')
        span_CL = self.CL_distributions['Span']
        CL_locals = self.CL_distributions[str(round(alpha, 1))]

        self.CL_dist = sc.interpolate.interp1d(span_CL, CL_locals, kind='cubic')

        span_CD = self.CD_distributions['Span']
        CD_locals = self.CD_distributions[str(round(alpha, 1))]

        self.CD_dist = sc.interpolate.interp1d(span_CD, CD_locals, kind='cubic')

        span_CM = self.CM_distributions['Span']
        CM_locals = self.CM_distributions[str(round(alpha, 1))]

        self.CM_dist = sc.interpolate.interp1d(span_CM, CM_locals, kind='cubic')

        span_COP = self.COP_distributions['Span']
        COP_locals = self.COP_distributions[str(round(alpha, 1))]

        self.COP_dist = sc.interpolate.interp1d(span_COP, COP_locals, kind='cubic')

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
        # fig.set_size_inches(13, 6)
        ax1.set_title('CL distribution over halfspan.')
        ax1.set_xlabel('Wing halfspan, [m]')
        ax1.set_ylabel('Lift coefficient, [-]')
        ax2.set_title('CD distribution over halfspan.')
        ax2.set_xlabel('Wing halfspan, [m]')
        ax2.set_ylabel('Drag coefficient, [-]')
        ax3.set_title('CM distribution over halfspan.')
        ax3.set_xlabel('Wing halfspan, [m]')
        ax3.set_ylabel('Moment coefficient, [-]')
        ax1.plot(self.sec_spans, self.CL_dist(self.sec_spans))
        ax2.plot(self.sec_spans, self.CD_dist(self.sec_spans))
        ax3.plot(self.sec_spans, self.CM_dist(self.sec_spans))
        plt.show()

        # plt.plot(self.sec_spans, self.CL_dist(self.sec_spans))
        # plt.title('CL distribution over halfspan.')
        # plt.xlabel('Wing halfspan [m]')
        # plt.ylabel('Lift coefficient, CL [-]')
        # plt.show()

    def extrapolate_CL_dist(self, alpha):
        print('----- Plotting CL distribution over halfspan -----')
        span = self.CL_distributions['Span']
        CL_dist_0 = self.CL_distributions[str(0)].to_numpy()
        CL_dist_10 = self.CL_distributions[str(10)].to_numpy()
        CL_req_arr = CL_dist_0 + ((self.CL_req - self.CL_0)/(self.CL_10 - self.CL_0))*(CL_dist_10 - CL_dist_0)
        self.CL_dist = sc.interpolate.interp1d(span, CL_req_arr, kind='cubic')

        span_CD = self.CD_distributions['Span']
        CD_locals_0 = self.CD_distributions[str(0)]

        CD_scale = (self.a * alpha ** 2 + self.b * alpha + self.c) / self.c
        CD_locals_alpha = CD_locals_0 * CD_scale

        self.CD_dist = sc.interpolate.interp1d(span_CD, CD_locals_alpha, kind='cubic')

        span_CM = self.CM_distributions['Span']
        CM_locals = self.CM_distributions[str(round(alpha, 1))]

        self.CM_dist = sc.interpolate.interp1d(span_CM, CM_locals, kind='cubic')

        span_COP = self.COP_distributions['Span']
        COP_locals = self.COP_distributions[str(round(alpha, 1))]

        self.COP_dist = sc.interpolate.interp1d(span_COP, COP_locals, kind='cubic')

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
        # fig.set_size_inches(13, 6)
        ax1.set_title('CL distribution over halfspan.')
        ax1.set_xlabel('Wing halfspan, [m]')
        ax1.set_ylabel('Lift coefficient, [-]')
        ax2.set_title('CD distribution over halfspan.')
        ax2.set_xlabel('Wing halfspan, [m]')
        ax2.set_ylabel('Drag coefficient, [-]')
        ax3.set_title('CM distribution over halfspan.')
        ax3.set_xlabel('Wing halfspan, [m]')
        ax3.set_ylabel('Moment coefficient, [-]')
        ax1.plot(self.sec_spans, self.CL_dist(self.sec_spans))
        ax2.plot(self.sec_spans, self.CD_dist(self.sec_spans))
        ax3.plot(self.sec_spans, self.CM_dist(self.sec_spans))
        plt.show()

        # plt.plot(self.sec_spans, self.CL_dist(self.sec_spans))
        # plt.title('CL distribution over halfspan.')
        # plt.xlabel('Wing halfspan [m]')
        # plt.ylabel('Lift coefficient, CL [-]')
        # plt.show()

    def battery_load(self, y):
        y_bat_end = 0.3  # [m]
        bat_mass_tot = 165.6  # [kg]
        if isinstance(y, float) or isinstance(y, int):
            if y <= y_bat_end:
                return -(bat_mass_tot * 0.5 * 9.81) / y_bat_end
            else:
                return 0
        else:
            F_arr = []
            for y_val in y:
                if y_val <= y_bat_end:
                    F_arr.append(-(bat_mass_tot * 0.5 * 9.81) / y_bat_end)
                else:
                    F_arr.append(0)
            return np.array(F_arr)

    def normal_force(self, y):
        chords = sc.interpolate.interp1d(self.sec_spans, self.sec_chords, kind='cubic')
        L = self.CL_dist(y) * np.cos(np.deg2rad(self.alpha_req)) * chords(y) * 0.5 * self.rho * self.V**2
        D = self.CD_dist(y) * np.sin(np.deg2rad(self.alpha_req)) * chords(y) * 0.5 * self.rho * self.V**2
        F_bat = self.battery_load(y) * np.cos(np.deg2rad(self.alpha_req))
        return L + D + F_bat

    def tangential_force(self, y):
        chords = sc.interpolate.interp1d(self.sec_spans, self.sec_chords, kind='cubic')
        L = self.CL_dist(y) * np.sin(np.deg2rad(self.alpha_req)) * chords(y) * 0.5 * self.rho * self.V ** 2
        D = self.CD_dist(y) * np.cos(np.deg2rad(self.alpha_req)) * chords(y) * 0.5 * self.rho * self.V ** 2
        F_bat = self.battery_load(y) * np.sin(np.deg2rad(self.alpha_req))
        return L + D + F_bat

    def aero_moment(self, y):
        chords = sc.interpolate.interp1d(self.sec_spans, self.sec_chords, kind='cubic')
        M = self.CM_dist(y) * chords(y) * 0.5 * self.rho * self.V ** 2
        return M

    def shear_force_z(self, y_span):
        S_z = sc.integrate.quad(self.normal_force, y_span, self.half_span)
        return S_z[0]

    def shear_force_x(self, y_span):
        S_x = sc.integrate.quad(self.tangential_force, y_span, self.half_span)
        return S_x[0]

    def aero_torsion(self, y_span):
        M_y = sc.integrate.quad(self.aero_moment, y_span, self.half_span)
        return M_y[0]

    def bending_moment_x(self, y_span):
        # Bending moment around x-axis
        # M_x = sc.integrate.quad(self.shear_force_z, y_span, self.half_span)
        M_x = sc.integrate.quad(self.shear_func_z, y_span, self.half_span)
        return -M_x[0]

    def bending_moment_z(self, y_span):
        # Bending moment around z-axis
        # M_z = sc.integrate.quad(self.shear_force_x, y_span, self.half_span)
        M_z = sc.integrate.quad(self.shear_func_x, y_span, self.half_span)
        return M_z[0]

    def plot_COP_dist(self):
        y_arr = np.linspace(0, self.half_span, num=self.nr_bsec, endpoint=True)
        COP = self.COP_dist(y_arr)/100

        plt.plot(y_arr, COP)
        plt.xlabel('Wing halfspan, [m]')
        plt.ylabel('COP location, [% of chord]')
        plt.title('Center of pressure locations along the wing halfspan')
        plt.show()

    def plot_decomposed_forces(self):
        print('----- Plotting distribution of lift per unit wingspan -----')
        y_arr = np.linspace(0, self.half_span, num=self.nr_bsec, endpoint=True)
        Fz = self.normal_force(y_arr)
        Fx = self.tangential_force(y_arr)
        My = self.aero_moment(y_arr)

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
        # fig.set_size_inches(13, 6)
        ax1.set_title('Distribution of normal force per unit wingspan.')
        ax1.set_xlabel('Wing halfspan, [m]')
        ax1.set_ylabel('Force per unit wingspan [N/m]')
        ax2.set_title('Distribution of tangential force per unit wingspan.')
        ax2.set_xlabel('Wing halfspan, [m]')
        ax2.set_ylabel('Force per unit wingspan [N/m]')
        ax3.set_title('Distribution of aero. moment per unit wingspan.')
        ax3.set_xlabel('Wing halfspan, [m]')
        ax3.set_ylabel('Aero. moment per unit wingspan [N]')
        ax1.plot(y_arr, Fz)
        ax2.plot(y_arr, Fx)
        ax3.plot(y_arr, My)
        plt.tight_layout()
        plt.show()

        # plt.plot(y_arr, Fz)
        # plt.title('Distribution of normal force per unit wingspan.')
        # plt.xlabel('Wing halfspan [m]')
        # plt.ylabel('Lift per unit wingspan [N/m]')
        # plt.show()

    def plot_shear_forces(self):
        print('----- Plotting internal shear distribution in z direction -----')
        y_arr = np.linspace(0, self.half_span, num=self.nr_bsec, endpoint=True)
        # y_arr = np.linspace(0, self.half_span, num=5, endpoint=True)

        shear_dist_z = []
        shear_dist_x = []
        aero_moment_y = []
        for y in y_arr:
            shear_dist_z.append(self.shear_force_z(y))
            shear_dist_x.append(self.shear_force_x(y))
            aero_moment_y.append((self.aero_torsion(y)))

        self.shear_dist_z = np.array(shear_dist_z)
        self.shear_dist_x = np.array(shear_dist_x)
        self.aero_moment_dist = np.array(aero_moment_y)

        self.shear_func_z = sc.interpolate.interp1d(y_arr, self.shear_dist_z, kind='cubic')
        self.shear_func_x = sc.interpolate.interp1d(y_arr, self.shear_dist_x, kind='cubic')
        self.aero_moment_func = sc.interpolate.interp1d(y_arr, self.aero_moment_dist, kind='cubic')

        fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
        # fig.set_size_inches(13, 6)
        ax1.set_title('Distributions of Sz over wing halfspan.')
        ax1.set_xlabel('Wing halfspan, [m]')
        ax1.set_ylabel('Sz [N]')
        ax2.set_title('Distributions of Sx over wing halfspan.')
        ax2.set_xlabel('Wing halfspan, [m]')
        ax2.set_ylabel('Sx [N]')
        ax3.set_title('Distributions of M_aero over wing halfspan.')
        ax3.set_xlabel('Wing halfspan, [m]')
        ax3.set_ylabel('M_aero [N]')
        ax1.plot(y_arr, shear_dist_z)
        ax2.plot(y_arr, shear_dist_x)
        ax3.plot(y_arr, aero_moment_y)
        plt.tight_layout()
        plt.show()

        # plt.plot(y_arr, shear_dist, c='r')
        # plt.title('Internal shear distribution in z direction.')
        # plt.xlabel('Wing halfspan [m]')
        # plt.ylabel('Shear force, [N]')
        # plt.show()

    def plot_bending_moments(self):
        print('----- Plotting bending moment distribution along x-axis -----')
        # y_arr = np.linspace(0, self.half_span, num=self.nr_bsec, endpoint=True)
        y_arr = np.linspace(0, self.half_span, num=self.nr_bsec, endpoint=True)  # Temporary smaller array to decrease computational time
        moment_dist_x = []
        moment_dist_z = []
        for y in y_arr:
            print('Moments are being computed at span:', y)
            moment_dist_x.append(self.bending_moment_x(y))
            moment_dist_z.append(self.bending_moment_z(y))

        self.moment_dist_x = sc.interpolate.interp1d(y_arr, np.array(moment_dist_x), kind='cubic')
        self.moment_dist_z = sc.interpolate.interp1d(y_arr, np.array(moment_dist_z), kind='cubic')

        fig, (ax1, ax2) = plt.subplots(1, 2)
        fig.set_size_inches(13, 6)
        ax1.set_title('Distributions of Mx over wing halfspan.')
        ax1.set_xlabel('Wing halfspan, [m]')
        ax1.set_ylabel('Bending moment Mx, [Nm]')
        ax2.set_title('Distributions of Mz over wing halfspan.')
        ax2.set_xlabel('Wing halfspan, [m]')
        ax2.set_ylabel('Bending moment Mz, [Nm]')
        ax1.plot(y_arr, moment_dist_x)
        ax2.plot(y_arr, moment_dist_z)
        plt.tight_layout()
        plt.show()

        # plt.plot(y_arr, moment_dist_x, c='r')
        # plt.title('Bending moment distribution in z direction.')
        # plt.xlabel('Wing halfspan [m]')
        # plt.ylabel('Bending moment, [Nm]')
        # plt.show()


# WL = WingLoading(DIRECTORY, INPUT_FOLDER, CL_dist_file, CD_dist_file, c_root, c_tip, MAC, dihedral,
#                  LE_sweep, half_span, AR, nr_span_sections, ac_mass, g_load, CL_0, CL_10, S, rho, V,
#                  a, b, c)

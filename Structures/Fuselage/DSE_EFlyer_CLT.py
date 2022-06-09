# Lennart Krieg
# DSE group 20: Aerobatic E-Flyer
# Classical laminate theory module as part of wing structural model

import numpy as np


class LaminateTheory:
    def __init__(self):
        self.E11 = 56.9E9  # [Pa]
        self.E22 = 56.9E9  # [Pa]
        self.v12 = 0.033  # [-]
        self.v21 = (self.E22/self.E11) * self.v12  # [-]
        self.G12 = 5.5E9  # [Pa]
        self.t_ply = 0.08/1000  # [mm]
        self.layup = [0, 45, 45, 0]  # [deg]

        self.tens_max_11 = 1259E6  # [Pa]
        self.tens_max_22 = 1259E6  # [Pa]
        self.compr_max_11 = 755E6  # [Pa]
        self.compr_max_22 = 755E6  # [Pa]
        self.shear_max12 = 102E6  # [Pa]

        self.C11 = self.E11/(1 - self.v12 * self.v21)
        self.C12 = self.v12 * self.E22/(1 - self.v12 * self.v21)
        self.C22 = self.E22/(1 - self.v12 * self.v21)
        self.C66 = self.G12

        self.Q = None
        self.t_lam = 0
        self.h = []
        self.t_half = 0

        self.A = None
        self.B = None
        self.D = None

        self.Exx = 0
        self.Eyy = 0
        self.Gxy = 0
        self.vxy = 0
        self.vyx = 0

        self.tens_max_xx = 0
        self.tens_max_yy = 0
        self.compr_max_xx = 0
        self.compr_max_yy = 0

        self.geometric_props()
        self.Q_matrices()
        self.ABD()
        self.in_plane_properties()
        self.max_stresses()
        self.print_properties()

    def geometric_props(self):
        self.t_lam = self.t_ply * len(self.layup)
        self.t_half = self.t_lam/2

        for i in range(len(self.layup)):
            ply_top = (i + 1) * self.t_ply - self.t_half
            ply_bot = i * self.t_ply - self.t_half
            self.h.append([ply_top, ply_bot])

    def Q_matrices(self):
        Qs = []
        for theta in self.layup:

            if theta != 0:

                m = np.cos(theta * np.pi / 180)
                n = np.sin(theta * np.pi / 180)

                Q11 = (self.C11 * (m ** 4)) + ((2 * (self.C12 + (2 * self.C66))) * (m ** 2) * (n ** 2)) + (self.C22 * (n ** 4))
                Q12 = ((self.C11 + self.C22 - (4 * self.C66)) * (m ** 2) * (n ** 2)) + (self.C12 * ((m ** 4) + (n ** 4)))
                Q22 = (self.C11 * (n ** 4)) + ((2 * (self.C12 + (2 * self.C66))) * (m ** 2) * (n ** 2)) + (self.C22 * (m ** 4))
                Q16 = ((self.C11 - self.C12 - (2 * self.C66)) * n * (m ** 3)) + ((self.C12 - self.C22 + (2 * self.C66)) * (n ** 3) * m)
                Q26 = ((self.C11 - self.C12 - (2 * self.C66)) * m * (n ** 3)) + ((self.C12 - self.C22 + (2 * self.C66)) * (m ** 3) * n)
                Q66 = ((self.C11 + self.C22 - (2 * self.C12) - (2 * self.C66)) * (n ** 2) * (m ** 2)) + (self.C66 * ((n ** 4) + (m ** 4)))

            else:
                Q11 = self.C11
                Q12 = self.C12
                Q22 = self.C22
                Q16 = 0
                Q26 = 0
                Q66 = self.C66

            Q = np.array([[Q11, Q12, Q16],
                          [Q12, Q22, Q26],
                          [Q16, Q26, Q66]], dtype='object')
            Qs.append(Q)
        self.Q = np.array(Qs)

    def ABD(self):

        A11 = 0
        A12 = 0
        A22 = 0
        A16 = 0
        A26 = 0
        A66 = 0

        B11 = 0
        B12 = 0
        B22 = 0
        B16 = 0
        B26 = 0
        B66 = 0

        D11 = 0
        D12 = 0
        D22 = 0
        D16 = 0
        D26 = 0
        D66 = 0

        for idx, z in enumerate(self.h):
            A11 += self.Q[idx][0][0] * (z[0] - z[1])
            A12 += self.Q[idx][0][1] * (z[0] - z[1])
            A22 += self.Q[idx][1][1] * (z[0] - z[1])
            A16 += self.Q[idx][0][2] * (z[0] - z[1])
            A26 += self.Q[idx][1][2] * (z[0] - z[1])
            A66 += self.Q[idx][2][2] * (z[0] - z[1])

            B11 += self.Q[idx][0][0] * (z[0]**2 - z[1]**2)/2
            B12 += self.Q[idx][0][1] * (z[0]**2 - z[1]**2)/2
            B22 += self.Q[idx][1][1] * (z[0]**2 - z[1]**2)/2
            B16 += self.Q[idx][0][2] * (z[0]**2 - z[1]**2)/2
            B26 += self.Q[idx][1][2] * (z[0]**2 - z[1]**2)/2
            B66 += self.Q[idx][2][2] * (z[0]**2 - z[1]**2)/2

            D11 += self.Q[idx][0][0] * (z[0]**3 - z[1]**3)/3
            D12 += self.Q[idx][0][1] * (z[0]**3 - z[1]**3)/3
            D22 += self.Q[idx][1][1] * (z[0]**3 - z[1]**3)/3
            D16 += self.Q[idx][0][2] * (z[0]**3 - z[1]**3)/3
            D26 += self.Q[idx][1][2] * (z[0]**3 - z[1]**3)/3
            D66 += self.Q[idx][2][2] * (z[0]**3 - z[1]**3)/3

        self.A = np.array([[A11, A12, A16],
                           [A12, A22, A26],
                           [A16, A26, A66]])

        self.B = np.array([[B11, B12, B16],
                           [B12, B22, B26],
                           [B16, B26, B66]])

        self.D = np.array([[D11, D12, D16],
                           [D12, D22, D26],
                           [D16, D26, D66]])

    def in_plane_properties(self):
        self.Exx = (self.A[0][0] * self.A[1][1] - self.A[0][1]**2) / (self.A[1][1] * self.t_lam)
        self.Eyy = (self.A[0][0] * self.A[1][1] - self.A[0][1]**2) / (self.A[0][0] * self.t_lam)
        self.Gxy = self.A[2][2] / self.t_lam
        self.vxy = self.A[0][1] / self.A[1][1]
        self.vyx = self.A[0][1] / self.A[0][0]

    def max_stresses(self):
        theta_largest = max(self.layup)

        if self.tens_max_11 <= self.tens_max_22:
            self.tens_max_xx = self.tens_max_11 * np.cos(theta_largest * np.pi / 180)
            self.tens_max_yy = self.tens_max_xx
        else:
            self.tens_max_xx = self.tens_max_22 * np.cos(theta_largest * np.pi / 180)
            self.tens_max_yy = self.tens_max_xx

        if self.compr_max_11 <= self.compr_max_22:
            self.compr_max_xx = self.compr_max_11 * np.cos(theta_largest * np.pi / 180)
            self.compr_max_yy = self.compr_max_xx
        else:
            self.compr_max_xx = self.compr_max_22 * np.cos(theta_largest * np.pi / 180)
            self.compr_max_yy = self.compr_max_xx

    def print_properties(self):
        print('---------- Laminate Properties ----------')
        print('Number of plies:', len(self.layup))
        print('ply thickness, [mm]:', self.t_ply*1000)
        print('Laminate thickness, [mm]:', self.t_lam * 1000)
        print('Orientations, [deg]:', self.layup)
        print('Exx, [GPa]:', self.Exx / 1000000000)
        print('Eyy, [GPa]:', self.Eyy / 1000000000)
        print('Gxy, [GPa]:', self.Gxy / 1000000000)
        print('Poissant vxy, [-]:', self.vxy)
        print('Poissant vyx, [-]:', self.vyx)
        print('Max tensile stress, [MPa]:', self.tens_max_xx/1000000)
        print('Max compressive stress, [MPa]:', self.compr_max_xx/1000000)
        print('Max shear stress, [MPa]:', self.shear_max12/1000000)


CLT = LaminateTheory()

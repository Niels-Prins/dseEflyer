import unittest

import numpy as np

from Stability.Code.aircraft import Aircraft

test_case = Aircraft(path='Aircraft/Cessna', example=True)
test_step = 0.1


class TestFuselage(unittest.TestCase):

    def test_initialization(self, margin=0.01):
        length_true = 6.900
        length_calc = round(test_case.fuselage.length, 3)

        diameter_true = 1.449
        diameter_calc = round(test_case.fuselage.diameter, 3)

        nose_front = 0.800
        nose_width_sides = 2 * 1.520
        nose_height_sides = 2 * 2.000

        area_nose_true = nose_front + nose_width_sides + nose_height_sides
        area_nose_calc = round(test_case.fuselage.area_nose, 3)

        main_width_sides = 2 * 1.870
        main_height_sides = 2 * 2.550

        area_main_true = main_width_sides + main_height_sides
        area_main_calc = round(test_case.fuselage.area_main, 3)

        tail_rear = 0.200
        tail_width_sides = 2 * 2.700
        tail_height_sides = 2 * 3.600

        area_tail_true = tail_rear + tail_width_sides + tail_height_sides
        area_tail_calc = round(test_case.fuselage.area_tail, 3)

        area_true = area_nose_true + area_main_true + area_tail_true
        area_calc = round(test_case.fuselage.area, 3)

        apparent_mass_true = 0.855
        apparent_mass_calc = round(test_case.fuselage.apparent_mass, 3)

        net_area_true = 14.263
        net_area_calc = test_case.wing.area_net

        cross_area_true = 1.649
        cross_area_calc = round(test_case.fuselage.area_cross, 3)

        self.assertAlmostEqual(length_true, length_calc, delta=margin * length_true)
        self.assertAlmostEqual(diameter_true, diameter_calc, delta=margin * diameter_true)
        self.assertAlmostEqual(area_nose_true, area_nose_calc, delta=margin * area_nose_true)
        self.assertAlmostEqual(area_main_true, area_main_calc, delta=margin * area_main_true)
        self.assertAlmostEqual(area_tail_true, area_tail_calc, delta=margin * area_tail_true)
        self.assertAlmostEqual(area_true, area_calc, delta=margin * area_true)
        self.assertAlmostEqual(apparent_mass_true, apparent_mass_calc, delta=margin * apparent_mass_true)
        self.assertAlmostEqual(net_area_true, net_area_calc, delta=margin * net_area_true)
        self.assertAlmostEqual(cross_area_true, cross_area_calc, delta=margin * cross_area_true)

    def test_washes(self, margin=0.01):
        d_x = abs(test_case.wing.X_ac_cg - test_case.aircraft[1].X_ac_cg)
        d_z = abs(test_case.wing.Z_ac_cg - test_case.aircraft[1].Z_ac_cg)

        K_ar = (1 / test_case.wing.ar) - (1 / (1 + (test_case.wing.ar ** 1.7)))
        K_taper = (10 - (3 * test_case.wing.taper)) / 7
        K_height = (1 - (d_z / test_case.wing.span)) / ((2 * d_x) / test_case.wing.span) ** (1 / 3)

        downwash_true = round(4.44 * (K_ar * K_taper * K_height) ** 1.19, 3)
        downwash_calc = round(test_case.aircraft[1].downwash, 3)

        wing_height = - 0.75
        fuselage_max_height = 1.50

        sidewash_true = round(((0.724 + ((3.06 / 2) * (test_case.aircraft[2].area / test_case.wing.area)))
                              + (0.4 * (wing_height / fuselage_max_height)) + (0.009 * test_case.wing.ar)), 3)
        sidewash_calc = round(test_case.aircraft[2].sidewash, 3)

        self.assertAlmostEqual(downwash_true, downwash_calc, delta=margin * downwash_true)
        self.assertAlmostEqual(sidewash_true, sidewash_calc, delta=margin * sidewash_true)

    def test_wing_corrections(self, margin=0.01):
        C_L_alpha = 4.5863
        C_L_alpha_cor = 5.1129

        K_n = 0.0431
        K_wf = 1.0830
        K_fw = 0.13734

        C_L_correction_true = round(C_L_alpha * (((K_n + K_wf + K_fw) * 0.88234) - 1), 3)
        C_L_correction_calc = round(test_case.C_L_correction, 3)

        C_D_skin = 0.009005
        C_D_base = 0.000420

        C_D_correction_true = round(C_D_skin + C_D_base, 3)
        C_D_correction_calc = round(test_case.C_D_correction, 3)

        C_L_nose = 0.1744
        C_L_wf = 4.3825
        C_L_fw = 0.5558

        X_ac_nose = -0.2243
        X_ac_wf = 0.2840
        X_ac_fw = 0.25

        X_ac_root = ((((C_L_nose * X_ac_nose) + (C_L_wf * X_ac_wf) + (C_L_fw * X_ac_fw)) / C_L_alpha_cor)
                     * test_case.wing.chord_root)

        C_M_ac_correction_true = -0.024
        C_M_ac_correction_calc = round(test_case.C_M_correction, 3)

        X_ac_correction_true = np.round((- 1.6 - X_ac_root) + 2.091, 3)
        X_ac_correction_calc = round(test_case.ac_correction, 3)

        C_Y_beta_true = -0.121
        C_Y_beta_calc = round(test_case.side_force_correction, 3)

        self.assertAlmostEqual(C_L_correction_true, C_L_correction_calc, delta=margin * C_L_correction_true)
        self.assertAlmostEqual(C_D_correction_true, C_D_correction_calc, delta=margin * C_D_correction_true)
        self.assertAlmostEqual(C_M_ac_correction_true, C_M_ac_correction_calc, delta=margin * C_M_ac_correction_true)
        self.assertAlmostEqual(X_ac_correction_true, X_ac_correction_calc, delta=margin * X_ac_correction_true)
        self.assertAlmostEqual(C_Y_beta_true, C_Y_beta_calc, delta=margin * C_Y_beta_true)

    def test_velocity_corrections(self):
        velocity = np.array([[8], [1], [0]])

        # Kelvin tests.
        width, height = 2, 4

        position_left = [- width, 0]
        position_right = [width, 0]
        position_down = [0, - height]
        position_up = [0, height]

        _, Y_dot_left, Z_dot_left = (test_case.fuselage.velocity_correction
                                     (velocity, position_left, height=height, width=width))
        _, Y_dot_right, Z_dot_right = (test_case.fuselage.velocity_correction
                                       (velocity, position_right, height=height, width=width))
        _, Y_dot_down, Z_dot_down = (test_case.fuselage.velocity_correction
                                     (velocity, position_down, height=height, width=width))
        _, Y_dot_up, Z_dot_up = (test_case.fuselage.velocity_correction
                                 (velocity, position_up, height=height, width=width))

        self.assertTrue(np.round(Y_dot_left, 3) == 0, np.round(Z_dot_left, 3) == 0)
        self.assertTrue(np.round(Y_dot_right, 3) == 0, np.round(Z_dot_right, 3) == 0)
        self.assertTrue(np.round(Y_dot_down, 3) == velocity[1, 0], np.round(Z_dot_down, 3) == 0)
        self.assertTrue(np.round(Y_dot_up, 3) == velocity[1, 0], np.round(Z_dot_up, 3) == 0)

        # Rankine tests.
        width, height = 4, 2

        position_left = [- width, 0]
        position_right = [width, 0]
        position_down = [0, - height]
        position_up = [0, height]

        _, Y_dot_left, Z_dot_left = (test_case.fuselage.velocity_correction
                                     (velocity, position_left, height=height, width=width))
        _, Y_dot_right, Z_dot_right = (test_case.fuselage.velocity_correction
                                       (velocity, position_right, height=height, width=width))
        _, Y_dot_down, Z_dot_down = (test_case.fuselage.velocity_correction
                                     (velocity, position_down, height=height, width=width))
        _, Y_dot_up, Z_dot_up = (test_case.fuselage.velocity_correction
                                 (velocity, position_up, height=height, width=width))

        self.assertTrue(np.round(Y_dot_left, 3) == 0, np.round(Z_dot_left, 3) == 0)
        self.assertTrue(np.round(Y_dot_right, 3) == 0, np.round(Z_dot_right, 3) == 0)
        self.assertTrue(np.round(Y_dot_down, 3) == velocity[1, 0], np.round(Z_dot_down, 3) == 0)
        self.assertTrue(np.round(Y_dot_up, 3) == velocity[1, 0], np.round(Z_dot_up, 3) == 0)

        # Doublet tests.
        width, height = 4, 4

        position_left = [- width, 0]
        position_right = [width, 0]
        position_down = [0, - height]
        position_up = [0, height]

        _, Y_dot_left, Z_dot_left = (test_case.fuselage.velocity_correction
                                     (velocity, position_left, height=height, width=width))
        _, Y_dot_right, Z_dot_right = (test_case.fuselage.velocity_correction
                                       (velocity, position_right, height=height, width=width))
        _, Y_dot_down, Z_dot_down = (test_case.fuselage.velocity_correction
                                     (velocity, position_down, height=height, width=width))
        _, Y_dot_up, Z_dot_up = (test_case.fuselage.velocity_correction
                                 (velocity, position_up, height=height, width=width))

        self.assertTrue(np.round(Y_dot_left, 3) == 0, np.round(Z_dot_left, 3) == 0)
        self.assertTrue(np.round(Y_dot_right, 3) == 0, np.round(Z_dot_right, 3) == 0)
        self.assertTrue(np.round(Y_dot_down, 3) == velocity[1, 0], np.round(Z_dot_down, 3) == 0)
        self.assertTrue(np.round(Y_dot_up, 3) == velocity[1, 0], np.round(Z_dot_up, 3) == 0)


if __name__ == '__main__':
    unittest.main()

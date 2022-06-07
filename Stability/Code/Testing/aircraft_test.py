import numpy as np

import unittest

from Stability.Code.aircraft import Aircraft

test_case = Aircraft(path='Aircraft/Cessna', example=True)
test_step = 0.1

gravity = 9.80665
mass = 1199.980
weight = mass * gravity
density = 1.060
area = 16.165
velocity = 60.000
mac = 1.494
span = 10.922


class TestAircraft(unittest.TestCase):

    def test_initialization(self, margin=0.01):
        m_c_true = 46.957   # calculated by hand using table 4-5 Flight dynamics.
        m_c_calc = round(test_case.m_c, 3)

        m_b_true = 6.423    # calculated by hand using table 4-8 Flight dynamics.
        m_b_calc = round(test_case.m_b, 3)

        K_XX_true = 0.009   # calculated by hand using table 4-8 Flight dynamics.
        K_XX_calc = round(test_case.K_XX, 3)

        K_YY_true = 0.681   # calculated by hand using table 4-5 Flight dynamics.
        K_YY_calc = round(test_case.K_YY, 3)

        K_ZZ_true = 0.019   # calculated by hand using table 4-8 Flight dynamics.
        K_ZZ_calc = round(test_case.K_ZZ, 3)

        K_XZ_true = 0.000   # calculated by hand using table 4-8 Flight dynamics.
        K_XZ_calc = round(test_case.K_XZ, 3)

        density_true = 1.060
        density_calc = round(test_case.rho, 3)

        magnitude_true = 60.000  # calculated by hand using Pythagoras Theorem.
        magnitude_calc = round(test_case.magnitude, 3)

        self.assertAlmostEqual(m_c_true, m_c_calc, delta=margin * m_c_true)
        self.assertAlmostEqual(m_b_true, m_b_calc, delta=margin * m_c_true)

        self.assertAlmostEqual(K_XX_true, K_XX_calc, delta=margin * K_XX_true)
        self.assertAlmostEqual(K_YY_true, K_YY_calc, delta=margin * K_YY_true)
        self.assertAlmostEqual(K_ZZ_true, K_ZZ_calc, delta=margin * K_ZZ_true)
        self.assertAlmostEqual(K_XZ_true, K_XZ_calc, delta=margin * K_XZ_true)

        self.assertAlmostEqual(density_true, density_calc, delta=margin * density_true)
        self.assertAlmostEqual(magnitude_true, magnitude_calc, delta=margin * magnitude_true)

    def test_summation(self, margin=0.01):
        # Outputs from each surface has been printed, summarized by Numpy and compared to the reference output.
        output_wing = np.array([[-22.447, -22.447, 0.000],
                                [189.690, -189.69, -2963.631],
                                [-6271.444, -6271.444, 0.000]])

        output_horizontal = np.array([[-32.329, 0.000, -729.573],
                                     [0.000, 0.000, 2088.987],
                                     [475.163, 0.000, -49.639]])

        output_vertical = np.array([[-13.795, 0.000, 0.000],
                                    [0.000, 0.000, 13.549],
                                    [0.000, 0.000, -0.177]])

        summation = output_wing + output_horizontal + output_vertical

        forces = (summation[:, 0] + summation[:, 1]).reshape(-1, 1)
        moments = summation[:, 2].reshape(-1, 1)

        summation_true = np.hstack((forces, moments))
        summation_calc = test_case.outputs_reference

        self.assertTrue(np.allclose(summation_true, summation_calc, rtol=margin))

    def test_symmetric(self, margin=0.01):
        symmetric_derivatives_calc = test_case.derivatives_symmetric

        output_reference = np.array([[-91.018, -729.573],
                                     [0.000, -861.096],
                                     [-12067.725, -49.816]])

        C_X_0, C_Z_0, C_M_0 = 0, (- 2 * weight) / (density * velocity ** 2 * area), 0

        initial_conditions_true = np.array([C_X_0, C_Z_0, C_M_0])
        initial_conditions_calc = symmetric_derivatives_calc[0]

        output_X_dot = np.array([[-92.535, -733.766],
                                 [0.000, -846.539],
                                 [-12094.187, -49.942]]) - output_reference

        X_X_dot, Z_X_dot, M_Y_X_dot = output_X_dot[0, 0], output_X_dot[2, 0], output_X_dot[1, 1]

        C_X_X_dot = (2 * X_X_dot) / (density * velocity * area * test_step)
        C_Z_X_dot = (2 * Z_X_dot) / (density * velocity * area * test_step)
        C_M_Y_X_dot = (2 * M_Y_X_dot) / (density * velocity * area * mac * test_step)

        X_dot_derivatives_true = np.round(np.array([C_X_X_dot, C_Z_X_dot, C_M_Y_X_dot]), 3)
        X_dot_derivatives_calc = np.round(symmetric_derivatives_calc[1, :], 3)

        output_Z_dot = np.array([[-65.201, -692.901],
                                 [0.000, -976.338],
                                 [-12356.75, -50.593]]) - output_reference

        X_Z_dot, Z_Z_dot, M_Y_Z_dot = output_Z_dot[0, 0], output_Z_dot[2, 0], output_Z_dot[1, 1]

        C_X_Z_dot = (2 * X_Z_dot) / (density * velocity * area * test_step)
        C_Z_Z_dot = (2 * Z_Z_dot) / (density * velocity * area * test_step)
        C_M_Y_Z_dot = (2 * M_Y_Z_dot) / (density * velocity * area * mac * test_step)

        Z_dot_derivatives_true = np.round(np.array([C_X_Z_dot, C_Z_Z_dot, C_M_Y_Z_dot]), 3)
        Z_dot_derivatives_calc = np.round(symmetric_derivatives_calc[2, :], 3)

        output_Z_dot_dot = np.array([[-91.008, -728.102],
                                     [0.000, -865.287],
                                     [-12068.683, -49.801]]) - output_reference

        X_Z_dot_dot, Z_Z_dot_dot, M_Y_Z_dot_dot = output_Z_dot_dot[0, 0], output_Z_dot_dot[2, 0], output_Z_dot_dot[1, 1]

        C_X_Z_dot_dot = (2 * X_Z_dot_dot) / (density * area * mac * test_step)
        C_Z_Z_dot_dot = (2 * Z_Z_dot_dot) / (density * area * mac * test_step)
        C_M_Y_Z_dot_dot = (2 * M_Y_Z_dot_dot) / (density * area * mac ** 2 * test_step)

        Z_dot_dot_derivatives_true = np.round(np.array([C_X_Z_dot_dot, C_Z_Z_dot_dot, C_M_Y_Z_dot_dot]), 3)
        Z_dot_dot_derivatives_calc = np.round(symmetric_derivatives_calc[3, :], 3)

        output_pitch_dot = np.array([[-93.893, -569.292],
                                     [0.000, -1316.107],
                                     [-12160.026, -52.469]]) - output_reference

        X_pitch_dot, Z_pitch_dot, M_Y_pitch_dot = output_pitch_dot[0, 0], output_pitch_dot[2, 0], output_pitch_dot[1, 1]

        C_X_pitch_dot = (2 * X_pitch_dot) / (density * velocity * area * mac * test_step)
        C_Z_pitch_dot = (2 * Z_pitch_dot) / (density * velocity * area * mac * test_step)
        C_M_Y_pitch_dot = (2 * M_Y_pitch_dot) / (density * velocity * area * mac ** 2 * test_step)

        pitch_dot_derivatives_true = np.round(np.array([C_X_pitch_dot, C_Z_pitch_dot, C_M_Y_pitch_dot]), 3)
        pitch_dot_derivatives_calc = np.round(symmetric_derivatives_calc[4, :], 3)

        output_elevator = np.array([[-62.929, 568.176],
                                    [0.000, -4771.282],
                                    [-12912.935, -6.688]]) - output_reference

        X_elevator, Z_elevator, M_Y_elevator = output_elevator[0, 0], output_elevator[2, 0], output_elevator[1, 1]

        C_X_elevator = (2 * X_elevator) / (density * velocity ** 2 * area * test_step)
        C_Z_elevator = (2 * Z_elevator) / (density * velocity ** 2 * area * test_step)
        C_M_Y_elevator = (2 * M_Y_elevator) / (density * velocity ** 2 * area * mac * test_step)

        elevator_derivatives_true = np.round(np.array([C_X_elevator, C_Z_elevator, C_M_Y_elevator]), 3)
        elevator_derivatives_calc = np.round(symmetric_derivatives_calc[5, :], 3)

        self.assertTrue(np.allclose(initial_conditions_true, initial_conditions_calc, rtol=margin))
        self.assertTrue(np.allclose(X_dot_derivatives_true, X_dot_derivatives_calc, rtol=margin))
        self.assertTrue(np.allclose(Z_dot_derivatives_true, Z_dot_derivatives_calc, rtol=margin))
        self.assertTrue(np.allclose(Z_dot_dot_derivatives_true, Z_dot_dot_derivatives_calc, rtol=margin))
        self.assertTrue(np.allclose(pitch_dot_derivatives_true, pitch_dot_derivatives_calc, rtol=margin))
        self.assertTrue(np.allclose(elevator_derivatives_true, elevator_derivatives_calc, rtol=margin))

    def test_asymmetric(self, margin=0.01):
        asymmetric_derivatives_calc = test_case.derivatives_asymmetric

        output_reference = np.array([[-91.018, -729.573],
                                     [0.000, -861.096],
                                     [-12067.725, -49.816]])

        output_beta = np.array([[-91.012, -866.419],
                                [-14.146, -861.091],
                                [-12067.755, -2.878]]) - output_reference

        Y_beta, M_X_beta, M_Z_beta = output_beta[1, 0], output_beta[0, 1], output_beta[2, 1]

        C_Y_beta = (2 * Y_beta) / (density * velocity * area * test_step)
        C_M_X_beta = (2 * M_X_beta) / (density * velocity * area * span * test_step)
        C_M_Z_beta = (2 * M_Z_beta) / (density * velocity * area * span * test_step)

        beta_derivatives_true = np.round(np.array([C_Y_beta + test_case.side_force_correction,
                                                   C_M_X_beta, C_M_Z_beta]), 3)
        beta_derivatives_calc = np.round(asymmetric_derivatives_calc[0, :], 3)

        output_beta_dot = np.array([[-91.022, -728.598],
                                    [0.993, -861.092],
                                    [-12067.725, -54.430]]) - output_reference

        Y_beta_dot, M_X_beta_dot, M_Z_beta_dot = output_beta_dot[1, 0], output_beta_dot[0, 1], output_beta_dot[2, 1]

        C_Y_beta_dot = (2 * Y_beta_dot) / (density * area * span * test_step)
        C_M_X_beta_dot = (2 * M_X_beta_dot) / (density * area * span ** 2 * test_step)
        C_M_Z_beta_dot = (2 * M_Z_beta_dot) / (density * area * span ** 2 * test_step)

        beta_dot_derivatives_true = np.round(np.array([C_Y_beta_dot, C_M_X_beta_dot, C_M_Z_beta_dot]), 3)
        beta_dot_derivatives_calc = np.round(asymmetric_derivatives_calc[1, :], 3)

        output_roll_dot = np.array([[-86.636, -2584.65],
                                    [-29.5, -703.819],
                                    [-12031.035, -182.249]]) - output_reference

        Y_roll_dot, M_X_roll_dot, M_Z_roll_dot = output_roll_dot[1, 0], output_roll_dot[0, 1], output_roll_dot[2, 1]

        C_Y_roll_dot = (4 * Y_roll_dot) / (density * velocity * area * span * test_step)
        C_M_X_roll_dot = (4 * M_X_roll_dot) / (density * velocity * area * span ** 2 * test_step)
        C_M_Z_roll_dot = (4 * M_Z_roll_dot) / (density * velocity * area * span ** 2 * test_step)

        roll_dot_derivatives_true = np.round(np.array([C_Y_roll_dot, C_M_X_roll_dot, C_M_Z_roll_dot]), 3)
        roll_dot_derivatives_calc = np.round(asymmetric_derivatives_calc[2, :], 3)

        output_yaw_dot = np.array([[-91.018, -481.673],
                                   [61.299, -842.768],
                                   [-12063.641, -333.485]]) - output_reference

        Y_yaw_dot, M_X_yaw_dot, M_Z_yaw_dot = output_yaw_dot[1, 0], output_yaw_dot[0, 1], output_yaw_dot[2, 1]

        C_Y_yaw_dot = (4 * Y_yaw_dot) / (density * velocity * area * span * test_step)
        C_M_X_yaw_dot = (4 * M_X_yaw_dot) / (density * velocity * area * span ** 2 * test_step)
        C_M_Z_yaw_dot = (4 * M_Z_yaw_dot) / (density * velocity * area * span ** 2 * test_step)

        yaw_dot_derivatives_true = np.round(np.array([C_Y_yaw_dot, C_M_X_yaw_dot, C_M_Z_yaw_dot]), 3)
        yaw_dot_derivatives_calc = np.round(asymmetric_derivatives_calc[3, :], 3)

        output_ailerons = np.array([[71.446, 9145.643],
                                    [114.707, -1004.866],
                                    [-12042.754, 910.703]]) - output_reference

        Y_ailerons, M_X_ailerons, M_Z_ailerons = output_ailerons[1, 0], output_ailerons[0, 1], output_ailerons[2, 1]

        C_Y_ailerons = (2 * Y_ailerons) / (density * velocity ** 2 * area * test_step)
        C_M_X_ailerons = (2 * M_X_ailerons) / (density * velocity ** 2 * area * span * test_step)
        C_M_Z_ailerons = (2 * M_Z_ailerons) / (density * velocity ** 2 * area * span * test_step)

        ailerons_derivatives_true = np.round(np.array([C_Y_ailerons, C_M_X_ailerons, C_M_Z_ailerons]), 3)
        ailerons_derivatives_calc = np.round(asymmetric_derivatives_calc[4, :], 3)

        output_rudder = np.array([[-52.579, 40.097],
                                  [783.634, -898.85],
                                  [-12067.725, -3889.281]]) - output_reference

        Y_rudder, M_X_rudder, M_Z_rudder = output_rudder[1, 0], output_rudder[0, 1], output_rudder[2, 1]

        C_Y_rudder = (2 * Y_rudder) / (density * velocity ** 2 * area * test_step)
        C_M_X_rudder = (2 * M_X_rudder) / (density * velocity ** 2 * area * span * test_step)
        C_M_Z_rudder = (2 * M_Z_rudder) / (density * velocity ** 2 * area * span * test_step)

        rudder_derivatives_true = np.round(np.array([C_Y_rudder, C_M_X_rudder, C_M_Z_rudder]), 3)
        rudder_derivatives_calc = np.round(asymmetric_derivatives_calc[5, :], 3)

        self.assertTrue(np.allclose(beta_derivatives_true, beta_derivatives_calc, rtol=margin))
        self.assertTrue(np.allclose(beta_dot_derivatives_true, beta_dot_derivatives_calc, rtol=margin))
        self.assertTrue(np.allclose(roll_dot_derivatives_true, roll_dot_derivatives_calc, rtol=margin))
        self.assertTrue(np.allclose(yaw_dot_derivatives_true, yaw_dot_derivatives_calc, rtol=margin))
        self.assertTrue(np.allclose(ailerons_derivatives_true, ailerons_derivatives_calc, rtol=margin))
        self.assertTrue(np.allclose(rudder_derivatives_true, rudder_derivatives_calc, rtol=margin))


if __name__ == '__main__':
    unittest.main()

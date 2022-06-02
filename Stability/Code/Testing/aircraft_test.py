import numpy as np

import unittest

from Stability.Code.aircraft import Aircraft

test_case = Aircraft(path='Aircraft/Cessna', example=True)


class TestAircraft(unittest.TestCase):

    def test_initialization(self, margin=0.01):
        m_c_true = 46.957   # calculated by hand using table 4-5 Flight dynamics.
        m_c_calc = test_case.m_c

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

        density_true = 1.058    # ...
        density_calc = round(test_case.rho, 3)

        magnitude_true = 60.075  # calculated by hand using Pythagoras Theorem
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
        # Outputs from each surface has been printed, summarized and compared to the reference output.
        output_wing = np.array([[-186.745, -186.745, 0.000],
                                [195.645, -195.645, -3404.751],
                                [-6468.334, -6468.334, 0.000]])
        output_horizontal = np.array([[-33.143, 0.000, -694.572],
                                      [0.000, 0.000, 1963.187],
                                      [452.367, 0.000, -50.888]])
        output_vertical = np.array([[-13.807, 0.000, 0.000],
                                    [0.000, 0.000, 5.277],
                                    [0.000, 0.000, -0.177]])

        summation = output_wing + output_horizontal + output_vertical

        forces = (summation[:, 0] + summation[:, 1]).reshape(-1, 1)
        moments = summation[:, 2].reshape(-1, 1)

        summation_true = np.hstack((forces, moments))
        summation_calc = test_case.outputs_reference

        self.assertTrue(np.allclose(summation_true, summation_calc, rtol=margin))

    def test_symmetric(self, margin=0.01):
        symmetric_derivatives_calc = test_case.derivatives_symmetric

        output_reference = test_case.outputs_reference

        gravity = 9.80665
        mass = 1200.000
        weight = mass * gravity
        density = 1.060
        area = 16.165
        velocity = 60.000
        mac = 1.494

        C_X_0, C_Z_0, C_M_0 = 0, (- 2 * weight) / (density * velocity ** 2 * area), 0

        initial_conditions_true = np.array([C_X_0, C_Z_0, C_M_0])  # C_Z_0 is equal to - C_L.
        initial_conditions_calc = symmetric_derivatives_calc[0]

        output_X_dot = np.array([[-447.842, -736.389],
                                 [0.000, -1455.15],
                                 [-12756.617, -52.41]]) - output_reference

        # print(np.round(output_X_dot))

        X, Z, M_Y = output_X_dot[0, 0], output_X_dot[2, 0], output_X_dot[1, 1]

        print(X, Z, M_Y)

        C_X_u = (2 * X) / (density * velocity * area)
        C_Z_u = (2 * Z) / (density * velocity * area)
        C_M_u = (2 * M_Y) / (density * velocity * area * mac)

        x_dot_derivatives_true = np.round(np.array([C_X_u, C_Z_u, C_M_u]), 3)
        x_dot_derivatives_calc = np.round(symmetric_derivatives_calc[1, :], 3)

        # z_dot_derivatives_true = np.array([None, None, None])
        # z_dot_derivatives_calc = None
        #
        # z_dot_dot_derivatives_true = np.array([None, None, None])
        # z_dot_dot_derivatives_calc = None
        #
        # pitch_dot_derivatives_true = np.array([None, None, None])
        # pitch_dot_derivatives_calc = None
        #
        # elevator_effectiveness_true = np.array([None, None, None])
        # elevator_effectiveness_calc = None
        #
        # self.assertTrue(np.allclose(initial_conditions_true, initial_conditions_calc, rtol=margin))
        # self.assertTrue(np.allclose(x_dot_derivatives_true, x_dot_derivatives_calc, rtol=margin))
        # self.assertTrue(np.allclose(z_dot_derivatives_true, z_dot_derivatives_calc, rtol=margin))
        # self.assertTrue(np.allclose(z_dot_dot_derivatives_true, z_dot_dot_derivatives_calc, rtol=margin))
        # self.assertTrue(np.allclose(pitch_dot_derivatives_true, pitch_dot_derivatives_calc, rtol=margin))
        # self.assertTrue(np.allclose(elevator_effectiveness_true, elevator_effectiveness_calc, rtol=margin))

    # def test_asymmetric(self, margin=0.01):
    #     # Per derivative, get outputs per derivative, make them dimensionless and crosscheck.
    #
    #     beta_derivatives_true = np.array([None, None, None])
    #     beta_derivatives_calc = None
    #
    #     beta_dot_derivatives_true = np.array([None, None, None])
    #     beta_dot_derivatives_calc = None
    #
    #     roll_dot_derivatives_true = np.array([None, None, None])
    #     roll_dot_derivatives_calc = None
    #
    #     yaw_dot_derivatives_true = np.array([None, None, None])
    #     yaw_dot_derivatives_calc = None
    #
    #     ailerons_effectiveness_true = np.array([None, None, None])
    #     ailerons_effectiveness_calc = None
    #
    #     rudder_effectiveness_true = np.array([None, None, None])
    #     rudder_effectiveness_calc = None
    #
    #     self.assertTrue(np.allclose(beta_derivatives_true, beta_derivatives_calc, rtol=margin))
    #     self.assertTrue(np.allclose(beta_dot_derivatives_true, beta_dot_derivatives_calc, rtol=margin))
    #     self.assertTrue(np.allclose(roll_dot_derivatives_true, roll_dot_derivatives_calc, rtol=margin))
    #     self.assertTrue(np.allclose(yaw_dot_derivatives_true, yaw_dot_derivatives_calc, rtol=margin))
    #     self.assertTrue(np.allclose(ailerons_effectiveness_true, ailerons_effectiveness_calc, rtol=margin))
    #     self.assertTrue(np.allclose(rudder_effectiveness_true, rudder_effectiveness_calc, rtol=margin))


if __name__ == '__main__':
    unittest.main()

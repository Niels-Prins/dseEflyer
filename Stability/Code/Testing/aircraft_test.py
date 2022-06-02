import numpy as np

import unittest

from Stability.Code.aircraft import Aircraft

test_case = Aircraft(path='Aircraft/Cessna', example=True)


class TestAircraft(unittest.TestCase):

    def test_initialization(self, margin=0.01):
        m_c_true = None
        m_c_calc = test_case.m_c

        m_b_true = None
        m_b_calc = test_case.m_b

        K_XX_true = None
        K_XX_calc = test_case.K_XX

        K_YY_true = None
        K_YY_calc = test_case.K_YY

        K_ZZ_true = None
        K_ZZ_calc = test_case.K_ZZ

        K_XZ_true = 0.000
        K_XZ_calc = test_case.K_XZ

        density_true = None
        density_calc = test_case.rho

        magnitude_true = None
        magnitude_calc = test_case.magnitude

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

        gravity = 9.80665
        mass = 1200.000
        weight = mass * gravity
        density = 1.058
        area = 16.165
        velocity = 60.000

        initial_conditions_true = np.array([0.000, -0.382, 0.000])  # C_Z_0 is equal to - C_L.
        initial_conditions_calc = symmetric_derivatives_calc[0]

        x_dot_derivatives_true = np.array([None, None, None])
        x_dot_derivatives_calc = None

        z_dot_derivatives_true = np.array([None, None, None])
        z_dot_derivatives_calc = None

        z_dot_dot_derivatives_true = np.array([None, None, None])
        z_dot_dot_derivatives_calc = None

        pitch_dot_derivatives_true = np.array([None, None, None])
        pitch_dot_derivatives_calc = None

        elevator_effectiveness_true = np.array([None, None, None])
        elevator_effectiveness_calc = None

        self.assertTrue(np.allclose(initial_conditions_true, initial_conditions_calc, rtol=margin))
        self.assertTrue(np.allclose(x_dot_derivatives_true, x_dot_derivatives_calc, rtol=margin))
        self.assertTrue(np.allclose(z_dot_derivatives_true, z_dot_derivatives_calc, rtol=margin))
        self.assertTrue(np.allclose(z_dot_dot_derivatives_true, z_dot_dot_derivatives_calc, rtol=margin))
        self.assertTrue(np.allclose(pitch_dot_derivatives_true, pitch_dot_derivatives_calc, rtol=margin))
        self.assertTrue(np.allclose(elevator_effectiveness_true, elevator_effectiveness_calc, rtol=margin))

    def test_asymmetric(self, margin=0.01):
        # Per derivative, get outputs per derivative, make them dimensionless and crosscheck.

        beta_derivatives_true = np.array([None, None, None])
        beta_derivatives_calc = None

        beta_dot_derivatives_true = np.array([None, None, None])
        beta_dot_derivatives_calc = None

        roll_dot_derivatives_true = np.array([None, None, None])
        roll_dot_derivatives_calc = None

        yaw_dot_derivatives_true = np.array([None, None, None])
        yaw_dot_derivatives_calc = None

        ailerons_effectiveness_true = np.array([None, None, None])
        ailerons_effectiveness_calc = None

        rudder_effectiveness_true = np.array([None, None, None])
        rudder_effectiveness_calc = None

        self.assertTrue(np.allclose(beta_derivatives_true, beta_derivatives_calc, rtol=margin))
        self.assertTrue(np.allclose(beta_dot_derivatives_true, beta_dot_derivatives_calc, rtol=margin))
        self.assertTrue(np.allclose(roll_dot_derivatives_true, roll_dot_derivatives_calc, rtol=margin))
        self.assertTrue(np.allclose(yaw_dot_derivatives_true, yaw_dot_derivatives_calc, rtol=margin))
        self.assertTrue(np.allclose(ailerons_effectiveness_true, ailerons_effectiveness_calc, rtol=margin))
        self.assertTrue(np.allclose(rudder_effectiveness_true, rudder_effectiveness_calc, rtol=margin))


if __name__ == '__main__':
    unittest.main()

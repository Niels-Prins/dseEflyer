import numpy as np

import unittest

from Stability.Code.aircraft import Aircraft


class TestAircraft(unittest.TestCase):

    def test_initialization(self, margin=0.01):
        m_c_true = None
        m_c_calc = None

        m_b_true = None
        m_b_calc = None

        K_XX_true = None
        K_XX_calc = None

        K_YY_true = None
        K_YY_calc = None

        K_ZZ_true = None
        K_ZZ_calc = None

        K_XZ_true = None
        K_XZ_calc = None

        density_true = None
        density_calc = None

        magnitude_true = None
        magnitude_calc = None

        self.assertAlmostEqual(m_c_true, m_c_calc, delta=margin * m_c_true)
        self.assertAlmostEqual(m_b_true, m_b_calc, delta=margin * m_c_true)

        self.assertAlmostEqual(K_XX_true, K_XX_calc, delta=margin * K_XX_true)
        self.assertAlmostEqual(K_YY_true, K_YY_calc, delta=margin * K_YY_true)
        self.assertAlmostEqual(K_ZZ_true, K_ZZ_calc, delta=margin * K_ZZ_true)
        self.assertAlmostEqual(K_XZ_true, K_XZ_calc, delta=margin * K_XZ_true)

        self.assertAlmostEqual(density_true, density_calc, delta=margin * density_true)
        self.assertAlmostEqual(magnitude_true, magnitude_calc, delta=margin * magnitude_true)

    def test_symmetric(self, margin=0.01):
        # Per derivative, get outputs per derivative, make them dimensionless and crosscheck.

        initial_conditions_true = np.array([None, None, None])
        initial_conditions_calc = None

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

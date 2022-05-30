import unittest
import numpy as np

from Code.surface import AerodynamicSurface


class TestAerodynamicSurface(unittest.TestCase):
    """
    All formulas used during testing have an identifier and can be found in the formulas file.
    """
    def load_surface(self, symmetric=True, vertical=False, main_body=None, simplify=False, controls=0):
        if not simplify:
            inputs = np.array([[0.349, 0.070, 15.215, 0.100],
                               [0.175, 0.035, 13.108, 0.200],
                               [0.698, 0.035, -1.300, 0.300]])
        else:
            inputs = np.array([[0.000, 0.000, 20.000, 0.000],
                               [0.000, 0.000, 0.000, 0.000],
                               [0.000, 0.000, 0.000, 0.000]])

        input_rho = 1.225
        input_magnitude = np.linalg.norm(inputs[:, 2])

        if symmetric:
            surface = AerodynamicSurface('Aircraft/Symmetric',
                                         symmetric=symmetric, vertical=vertical, downwash_body=main_body)
        else:
            if vertical:
                surface = AerodynamicSurface('Aircraft/Vertical',
                                             symmetric=symmetric, vertical=vertical, downwash_body=main_body)
            else:
                surface = AerodynamicSurface('Aircraft/Asymmetric',
                                             symmetric=symmetric, vertical=vertical, downwash_body=main_body)

        surface.calculate_outputs(inputs, input_magnitude, input_rho, controls=controls)

        return surface

    def test_initialization(self, margin=0.01):

        def test_symmetric():
            surface = self.load_surface()

            area_true = 0.180  # calculated by hand using formula [1].
            area_calc = round(surface.area, 3)

            mac_true = 0.156  # calculated by hand using formula [2].
            mac_calc = round(surface.mac, 3)

            X_ac_true = -0.044  # calculated by hand using formula [3.1].
            X_ac_calc = round(surface.X_ac_rel, 3)

            Y_ac_true = -0.267  # calculated by hand using formula [3.2].
            Y_ac_calc = round(surface.Y_ac_rel, 3)

            Z_ac_true = -0.023  # calculated by hand using formula [3.3].
            Z_ac_calc = round(surface.Z_ac_rel, 3)

            alpha_true = 0.048  # calculated by averaging the first column in the reference file.
            alpha_calc = round(np.mean(surface.alpha), 3)

            C_L_alpha_true = 0.082  # calculated by averaging the third column gradient in the reference file.
            C_L_alpha_calc = round(2 * np.mean(np.gradient(surface.C_L)), 3)

            C_D_alpha_true = 0.001  # calculated by averaging the sixth column gradient in the reference file.
            C_D_alpha_calc = round(2 * np.mean(np.gradient(surface.C_D)), 3)

            C_M_true = 0.052  # calculated by averaging the ninth column in the reference file.
            C_M_calc = round(surface.C_M, 3)

            downwash_true = 0.000  # must be zero as no main body is selected.
            downwash_calc = round(surface.downwash, 3)

            # Geometry tests.
            self.assertAlmostEqual(area_true, area_calc, delta=margin * area_true)
            self.assertAlmostEqual(mac_true, mac_calc, delta=margin * mac_true)
            self.assertAlmostEqual(X_ac_true, X_ac_calc, delta=margin * X_ac_true)
            self.assertAlmostEqual(Y_ac_true, Y_ac_calc, delta=margin * Y_ac_true)
            self.assertAlmostEqual(Z_ac_true, Z_ac_calc, delta=margin * Z_ac_true)

            # Aerodynamic tests.
            self.assertAlmostEqual(alpha_true, alpha_calc, delta=margin * alpha_true)
            self.assertAlmostEqual(C_L_alpha_true, C_L_alpha_calc, delta=margin * C_L_alpha_true)
            self.assertAlmostEqual(C_D_alpha_true, C_D_alpha_calc, delta=margin * C_D_alpha_true)
            self.assertAlmostEqual(C_M_true, C_M_calc, delta=margin * C_M_true)
            self.assertAlmostEqual(downwash_true, downwash_calc, delta=margin * downwash_true)

        def test_asymmetric():
            surface = self.load_surface(symmetric=False)

            area_true = 0.090  # calculated by hand using formula [1].
            area_calc = round(surface.area, 3)

            mac_true = 0.156  # calculated by hand using formula [2].
            mac_calc = round(surface.mac, 3)

            X_ac_true = -0.044  # calculated by hand using formula [3.1].
            X_ac_calc = round(surface.X_ac_rel, 3)

            Y_ac_true = -0.267  # calculated by hand using formula [3.2].
            Y_ac_calc = round(surface.Y_ac_rel, 3)

            Z_ac_true = -0.023  # calculated by hand using formula [3.3].
            Z_ac_calc = round(surface.Z_ac_rel, 3)

            # Tests.
            self.assertAlmostEqual(area_true, area_calc, delta=margin * area_true)
            self.assertAlmostEqual(mac_true, mac_calc, delta=margin * mac_true)
            self.assertAlmostEqual(X_ac_true, X_ac_calc, delta=margin * X_ac_true)
            self.assertAlmostEqual(Y_ac_true, Y_ac_calc, delta=margin * Y_ac_true)
            self.assertAlmostEqual(Z_ac_true, Z_ac_calc, delta=margin * Z_ac_true)

        def test_transformed():
            surface = self.load_surface(symmetric=False, vertical=True)

            X_ac_true = -0.184  # calculated by hand using formula [3.1] and [4.1].
            X_ac_calc = round(surface.X_ac_rel, 3)

            Y_ac_true = -0.317  # calculated by hand using formula [3.2] and [4.1].
            Y_ac_calc = round(surface.Y_ac_rel, 3)

            Z_ac_true = 0.277  # calculated by hand using formula [3.3] and [4.1].
            Z_ac_calc = round(surface.Z_ac_rel, 3)

            # Tests.
            self.assertAlmostEqual(X_ac_true, X_ac_calc, delta=margin * X_ac_true)
            self.assertAlmostEqual(Y_ac_true, Y_ac_calc, delta=margin * Y_ac_true)
            self.assertAlmostEqual(Z_ac_true, Z_ac_calc, delta=margin * Z_ac_true)

        test_symmetric()
        test_asymmetric()
        test_transformed()

    def test_transformation_input(self, margin=0.01):

        def test_untransformed():
            surface = self.load_surface()

            rho_true = 1.225
            rho_calc = round(surface.rho, 3)

            magnitude_true = 20.125
            magnitude_calc = round(surface.magnitude, 3)

            angles_true = np.array([[0.349], [0.175], [0.698]])  # unmodified input angles.
            angles_calc = np.round(surface.angles, 3)

            rates_true = np.array([[0.070], [0.035], [0.035]])  # unmodified input rates.
            rates_calc = np.round(surface.rates, 3)

            velocity_true = np.array([[20.001], [1.006], [1.992]])  # calculated by hand using formula [4].
            velocity_calc = np.round(surface.velocity, 3)

            acceleration_true = np.array([[0.150], [0.197], [0.281]])  # calculated by hand using formula [4].
            acceleration_calc = np.round(surface.acceleration, 3)

            # Tests.
            self.assertAlmostEqual(rho_true, rho_calc, delta=margin * rho_true)
            self.assertAlmostEqual(magnitude_true, magnitude_calc, delta=margin * magnitude_true)

            self.assertTrue(np.allclose(angles_true, angles_calc, rtol=margin))
            self.assertTrue(np.allclose(rates_true, rates_calc, rtol=margin))
            self.assertTrue(np.allclose(velocity_true, velocity_calc, rtol=margin))
            self.assertTrue(np.allclose(acceleration_true, acceleration_calc, rtol=margin))

        def test_transformed():
            surface = self.load_surface(vertical=True)

            angles_true = np.array([[0.349], [0.698], [-0.175]])  # calculated by hand using formula [4.1].
            angles_calc = np.round(surface.angles, 3)

            rates_true = np.array([[0.070], [0.035], [-0.035]])  # calculated by hand using formula [4.1].
            rates_calc = np.round(surface.rates, 3)

            velocity_true = np.array([[20.076], [1.195], [-0.720]])  # calculated by hand using formula [4].
            velocity_calc = np.round(surface.velocity, 3)

            acceleration_true = np.array([[0.164], [0.252], [-0.223]])  # calculated by hand using formula [4].
            acceleration_calc = np.round(surface.acceleration, 3)

            # Tests.
            self.assertTrue(np.allclose(angles_true, angles_calc, rtol=margin))
            self.assertTrue(np.allclose(rates_true, rates_calc, rtol=margin))
            self.assertTrue(np.allclose(velocity_true, velocity_calc, rtol=margin))
            self.assertTrue(np.allclose(acceleration_true, acceleration_calc, rtol=margin))

        test_untransformed()
        test_transformed()

    def test_calculate_local_velocity(self, margin=0.01):
        surface = self.load_surface()
        surface_variables = surface.calculate_local_velocity

        local_velocity_true = np.array([[19.980], [0.832], [2.124]])  # calculated by hand using formula [6].
        local_velocity_calc = surface_variables[0]

        local_magnitude_true = 20.110  # calculated by hand.
        local_magnitude_calc = round(surface_variables[1], 3)

        local_alpha_true = 0.106  # calculated by hand.
        local_alpha_calc = round(surface_variables[2], 3)

        local_alpha_dot_true = 0.013  # calculated by hand.
        local_alpha_dot_calc = np.round(surface_variables[3], 3)

        local_beta_true = 0.042  # calculated by hand.
        local_beta_calc = np.round(surface_variables[4], 3)

        local_beta_dot_true = 0.009  # calculated by hand.
        local_beta_dot_calc = np.round(surface_variables[5], 3)

        # Tests.
        self.assertTrue(np.allclose(local_velocity_true, local_velocity_calc, rtol=margin))
        self.assertAlmostEqual(local_magnitude_true, local_magnitude_calc, delta=margin * local_magnitude_true)
        self.assertAlmostEqual(local_alpha_true, local_alpha_calc, delta=margin * local_alpha_true)
        self.assertAlmostEqual(local_alpha_dot_true, local_alpha_dot_calc, delta=margin * local_alpha_dot_true)
        self.assertAlmostEqual(local_beta_true, local_beta_calc, delta=margin * local_beta_true)
        self.assertAlmostEqual(local_beta_dot_true, local_beta_dot_calc, delta=margin * local_beta_dot_true)

    def test_calculate_local_forces(self, margin=0.02):
        surface = self.load_surface()

        # Forces have been re-calculated by hand using the same steps as in the program.
        F_X_true, F_Y_true, F_Z_true = [0.814, -1.009, -11.567]
        F_X_calc, F_Y_calc, F_Z_calc = np.round(surface.calculate_local_forces, 3)

        # Tests.
        self.assertAlmostEqual(F_X_true, F_X_calc, delta=margin * abs(F_X_true))
        self.assertAlmostEqual(F_Y_true, F_Y_calc, delta=margin * abs(F_Y_true))
        self.assertAlmostEqual(F_Z_true, F_Z_calc, delta=margin * abs(F_Z_true))

    def test_calculate_forces(self, margin=0.02):

        def test_symmetric():
            surface = self.load_surface(simplify=True)

            # Forces have been calculated by simplifying the inputs and applying the lift formula.
            F_X_true, F_Y_true, F_Z_true = -0.322, 0, 2.770
            F_X_calc, F_Y_calc, F_Z_calc = surface.outputs[:, 0] + surface.outputs[:, 1]

            # Tests.
            self.assertAlmostEqual(F_X_true, F_X_calc, delta=margin * abs(F_X_true))
            self.assertAlmostEqual(F_Y_true, F_Y_calc, delta=margin)
            self.assertAlmostEqual(F_Z_true, F_Z_calc, delta=margin * abs(F_Z_true))

        def test_asymmetric():
            surface = self.load_surface(symmetric=False, simplify=True)

            # Forces have been calculated by simplifying the inputs and applying the lift formula.
            F_X_true, F_Y_true, F_Z_true = -0.161, -0.121, 1.385
            F_X_calc, F_Y_calc, F_Z_calc = surface.outputs[:, 0] + surface.outputs[:, 1]

            # Tests.
            self.assertAlmostEqual(F_X_true, F_X_calc, delta=margin * abs(F_X_true))
            self.assertAlmostEqual(F_Y_true, F_Y_calc, delta=margin * abs(F_Y_true))
            self.assertAlmostEqual(F_Z_true, F_Z_calc, delta=margin * abs(F_Z_true))

        test_symmetric()
        test_asymmetric()

    def test_calculate_moments(self, margin=0.02):
        surface = self.load_surface(symmetric=False, simplify=True)

        # Moments have been calculated by simplifying the inputs and applying torques.
        M_Y_due_F_X = 0.0037
        M_Z_due_F_X = - 0.0430

        M_X_due_F_Y = -0.0028
        M_Z_due_F_Y = 0.0053

        M_X_due_F_Z = -0.3700
        M_Y_due_F_Z = 0.0609

        M_Y_ac = 0.1790

        M_X_true = M_X_due_F_Y + M_X_due_F_Z
        M_Y_true = M_Y_due_F_X + M_Y_due_F_Z + M_Y_ac
        M_Z_true = M_Z_due_F_X + M_Z_due_F_Y

        M_X_calc, M_Y_calc, M_Z_calc = surface.outputs[:, 2]

        # Tests.
        self.assertAlmostEqual(M_X_true, M_X_calc, delta=margin * abs(M_X_true))
        self.assertAlmostEqual(M_Y_true, M_Y_calc, delta=margin * abs(M_Y_true))
        self.assertAlmostEqual(M_Z_true, M_Z_calc, delta=margin * abs(M_Z_true))

    def test_outputs(self, margin=0.02):

        def test_untransformed():
            surface = self.load_surface(simplify=True)

            # Same method as previous tests to obtain forces and moments.
            outputs_true = np.array([[-0.161, -0.161, 0.000],
                                     [-0.121, 0.121, 0.485],
                                     [1.385, 1.385, 0.000]])
            outputs_calc = surface.outputs

            # Test.
            self.assertTrue(np.allclose(outputs_true, outputs_calc, rtol=margin))

        def test_transformed():
            surface = self.load_surface(symmetric=False, vertical=True, simplify=True)

            # Same method as previous tests to obtain forces and moments.
            outputs_true = np.array([[-0.161, 0.000, -0.405],
                                     [-1.385, 0.000, 0.029],
                                     [-0.121, 0.000, 0.388]])
            outputs_calc = surface.outputs

            # Test.
            self.assertTrue(np.allclose(outputs_true, outputs_calc, rtol=margin))

        test_untransformed()
        test_transformed()


if __name__ == '__main__':
    unittest.main()

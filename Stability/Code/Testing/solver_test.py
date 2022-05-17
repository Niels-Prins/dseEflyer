import unittest
import numpy as np
import sympy as sym

from Code.solver import AircraftStability


class TestAircraftStability(unittest.TestCase):
    """
    All formulas used during testing have an identifier and can be found in the formulas file.
    """
    def test_state_space(self, margin=0.01):

        def test_symmetric():
            # Manual calculated values.
            P = - np.array([[6.93354, 0, 0, 0],
                            [0, 6.98181, 0, 0],
                            [0, 0, 0.03376, 0],
                            [0, 0.12490, 0, 6.79487]])

            Q = - np.array([[-0.2199, 0.4653, -1.1360, 0],
                            [-2.2720, -5.1600, 0, 201.5400],
                            [0, 0, 0, 1],
                            [0, -0.4300, 0, -7.0400]])

            R = - np.array([[0, 0],
                            [-0.6238, 0],
                            [0, 0],
                            [-1.5530, 0]])

            # Calculate matrix A and B.
            A_true = np.round(np.dot(np.linalg.inv(P), Q), 4)
            B_true = np.round(np.dot(np.linalg.inv(P), R), 4)

            A_calc = np.round(test_case.A_symmetric, 4)
            B_calc = np.round(test_case.B_symmetric, 4)

            # Tests.
            self.assertTrue(np.allclose(A_calc, A_true, rtol=margin))
            self.assertTrue(np.allclose(B_calc, B_true, rtol=margin))

        def test_asymmetric():
            # Manual calculated values.
            P = - np.array([[6.91419, 0, 0, 0],
                            [0, 0.11152, 0, 0],
                            [0, 0, 0.16594, -0.02766],
                            [0, 0, -0.02766, 0.51165]])

            Q = - np.array([[-0.9896, 1.1360, -0.0870, -61.5700],
                            [0, 0, 1, 0],
                            [-0.0772, 0, -0.3444, 0.2800],
                            [0.1638, 0, -0.0108, -0.1930]])

            R = - np.array([[0, 0.3037],
                            [0, 0],
                            [-0.2349, 0.0286],
                            [0.0286, -0.1261]])

            # Calculate matrix A and B.
            A_true = np.round(np.dot(np.linalg.inv(P), Q), 4)
            B_true = np.round(np.dot(np.linalg.inv(P), R), 4)

            A_calc = np.round(test_case.A_asymmetric, 4)
            B_calc = np.round(test_case.B_asymmetric, 4)

            # Tests.
            self.assertTrue(np.allclose(A_calc, A_true, rtol=margin))
            self.assertTrue(np.allclose(B_calc, B_true, rtol=margin))

        test_case = AircraftStability(example=True)
        test_symmetric()
        test_asymmetric()

    def test_eigenvalues(self, margin=0.01):

        def test_symmetric():
            # Forces and moments coefficients.
            C_X_0, C_Z_0, C_M_0 = symmetric_coefficients[0, :]
            C_X_u, C_Z_u, C_M_u = symmetric_coefficients[1, :]
            C_X_alpha, C_Z_alpha, C_M_alpha = symmetric_coefficients[2, :]
            C_X_alpha_dot, C_Z_alpha_dot, C_M_alpha_dot = symmetric_coefficients[3, :]
            C_X_q, C_Z_q, C_M_q = symmetric_coefficients[4, :]

            # Aerodynamic & geometric coefficients.
            velocity, m_c, mac, K_YY = symmetric_misc_coefficients

            # Symmetric equations of motion matrix.
            row_1 = [C_X_u - (2 * m_c * eigenvalue), C_X_alpha, C_Z_0, 0]
            row_2 = [C_Z_u, C_Z_alpha + (eigenvalue * (C_Z_alpha_dot - (2 * m_c))), - C_X_0, C_Z_q + (2 * m_c)]
            row_3 = [0, 0, - eigenvalue, 1]
            row_4 = [C_M_u, C_M_alpha + (eigenvalue * C_M_alpha_dot), 0, C_M_q - (2 * m_c * K_YY * eigenvalue)]

            eom = [row_1, row_2, row_3, row_4]
            eom_sym = sym.Matrix(eom)

            # Calculating true eigenvalues.
            eigenvalues_true = sym.solve(eom_sym.det(), eigenvalue)
            eigenvalues_true = [complex(i * (velocity / mac)) for i in eigenvalues_true]
            eigenvalues_true = np.sort(np.round(eigenvalues_true, 4))

            # Calculated eigenvalues by program.
            eigenvalues_calc = test_case.eigenvalues_symmetric
            eigenvalues_calc = np.sort(np.round(eigenvalues_calc, 4))

            # Test.
            self.assertTrue(np.allclose(eigenvalues_calc, eigenvalues_true, rtol=margin))
            print('Symmetric test completed!')

        def test_asymmetric():
            # Forces and moments coefficients.
            C_Y_beta, C_l_beta, C_n_beta = asymmetric_coefficients[0, :]
            C_Y_beta_dot, C_l_beta_dot, C_n_beta_dot = asymmetric_coefficients[1, :]
            C_Y_p, C_l_p, C_n_p = asymmetric_coefficients[2, :]
            C_Y_r, C_l_r, C_n_r = asymmetric_coefficients[3, :]

            # Aerodynamic & geometric coefficients.
            velocity, C_L, m_b, span, K_XX, K_ZZ, K_XZ = asymmetric_misc_coefficients

            # Asymmetric equations of motion.
            row_1 = [C_Y_beta + (C_Y_beta_dot - 2 * m_b) * eigenvalue, C_L, C_Y_p, C_Y_r - 4 * m_b]
            row_2 = [0, - eigenvalue / 2, 1, 0]
            row_3 = [C_l_beta, 0, C_l_p - (4 * m_b * K_XX * eigenvalue), C_l_r + (4 * m_b * K_XZ * eigenvalue)]
            row_4 = [C_n_beta + (C_n_beta_dot * eigenvalue), 0, C_n_p + (4 * m_b * K_XZ * eigenvalue),
                     C_n_r - (4 * m_b * K_ZZ * eigenvalue)]

            eom = [row_1, row_2, row_3, row_4]
            eom_sym = sym.Matrix(eom)

            # Calculating true eigenvalues.
            eigenvalues_true = sym.solve(eom_sym.det(), eigenvalue)
            eigenvalues_true = [complex(i * (velocity / span)) for i in eigenvalues_true]
            eigenvalues_true = np.sort(np.round(eigenvalues_true, 4))

            # Calculated eigenvalues by program.
            eigenvalues_calc = test_case.eigenvalues_asymmetric
            eigenvalues_calc = np.sort(np.round(eigenvalues_calc, 4))

            # Test.
            self.assertTrue(np.allclose(eigenvalues_calc, eigenvalues_true, rtol=margin))
            print('Asymmetric test completed!')

        test_case = AircraftStability(example=True)
        eigenvalue = sym.Symbol('eigenvalue')

        symmetric_coefficients = test_case.symmetric_coefficients
        symmetric_misc_coefficients = test_case.symmetric_misc_coefficients

        asymmetric_coefficients = test_case.asymmetric_coefficients
        asymmetric_misc_coefficients = test_case.asymmetric_misc_coefficients

        test_symmetric()
        test_asymmetric()


if __name__ == '__main__':
    unittest.main()

import matplotlib.pyplot as plt
import numpy as np
import control as ctrl

import os

from Stability.Code.aircraft import Aircraft


class AircraftStability():

    def __init__(self, path, fuselage=True, example=False, coefficients=None):
        self.example = example

        # Obtaining aircraft.
        if not self.example:
            self.aircraft = Aircraft(path, fuselage=fuselage)
        else:
            (self.symmetric_coefficients, self.symmetric_misc_coefficients,
             self.asymmetric_coefficients, self.asymmetric_misc_coefficients) = coefficients

        # Equations of motion in state-space format, calculated by state-space functions.
        self.A_symmetric = None
        self.B_symmetric = None
        self.C_symmetric = None
        self.D_symmetric = None

        self.A_asymmetric = None
        self.B_asymmetric = None
        self.C_asymmetric = None
        self.D_asymmetric = None

        self.state_space_symmetric_eom()
        self.state_space_asymmetric_eom()

        # Eigenvalues of the equations of motion, calculated by the eigenvalue function.
        self.eigenvalues_symmetric = None
        self.eigenvalues_asymmetric = None

        self.period_symmetric = None
        self.period_asymmetric = None

        self.half_amplitude_symmetric = None
        self.half_amplitude_asymmetric = None

        self.eigenvalues()

        if not self.example:
            if not os.path.isdir(f'Results/{self.aircraft.name}/Eigenvalues'):
                os.makedirs(f'Results/{self.aircraft.name}/Eigenvalues')

            if not os.path.isdir(f'Results/{self.aircraft.name}/Responses'):
                os.makedirs(f'Results/{self.aircraft.name}/Responses')

    def state_space_symmetric_eom(self):
        if not self.example:
            # Forces and moments coefficients.
            C_X_0, C_Z_0, C_M_0 = self.aircraft.derivatives_symmetric[0, :]
            C_X_u, C_Z_u, C_M_u = self.aircraft.derivatives_symmetric[1, :]
            C_X_alpha, C_Z_alpha, C_M_alpha = self.aircraft.derivatives_symmetric[2, :]
            C_X_alpha_dot, C_Z_alpha_dot, C_M_alpha_dot = self.aircraft.derivatives_symmetric[3, :]
            C_X_q, C_Z_q, C_M_q = self.aircraft.derivatives_symmetric[4, :]
            C_X_elevator, C_Z_elevator, C_M_elevator = self.aircraft.derivatives_symmetric[5, :]

            # Aerodynamic coefficients.
            velocity = self.aircraft.velocity

            # Geometry coefficients.
            m_c, mac, K_YY = self.aircraft.m_c, self.aircraft.wing.mac, self.aircraft.K_YY

        else:
            # Symmetric coefficients.
            C_X_0, C_Z_0, C_M_0 = self.symmetric_coefficients[0, :]
            C_X_u, C_Z_u, C_M_u = self.symmetric_coefficients[1, :]
            C_X_alpha, C_Z_alpha, C_M_alpha = self.symmetric_coefficients[2, :]
            C_X_alpha_dot, C_Z_alpha_dot, C_M_alpha_dot = self.symmetric_coefficients[3, :]
            C_X_q, C_Z_q, C_M_q = self.symmetric_coefficients[4, :]
            C_X_elevator, C_Z_elevator, C_M_elevator = self.symmetric_coefficients[5, :]

            # Miscellaneous coefficients.
            velocity, m_c, mac, K_YY = self.symmetric_misc_coefficients

        # Setting up the matrices to transform to the state-space format.
        P = - (mac / velocity) * np.array([[2 * m_c, 0, 0, 0],
                                           [0, ((2 * m_c) - C_Z_alpha_dot), 0, 0],
                                           [0, 0, 1, 0],
                                           [0, - C_M_alpha_dot, 0, 2 * m_c * K_YY]])

        Q = - np.array([[C_X_u, C_X_alpha, C_Z_0, 0],
                        [C_Z_u, C_Z_alpha, - C_X_0, C_Z_q + (2 * m_c)],
                        [0, 0, 0, 1],
                        [C_M_u, C_M_alpha, 0, C_M_q]])

        R = - np.array([[C_X_elevator, 0],
                        [C_Z_elevator, 0],
                        [0, 0],
                        [C_M_elevator, 0]])

        # Assigning state-space matrices to class attributes.
        self.A_symmetric = np.dot(np.linalg.inv(P), Q)
        self.B_symmetric = np.dot(np.linalg.inv(P), R)

        self.C_symmetric = np.array([[1, 0, 0, 0],   # dimensionless x-velocity
                                     [0, 1, 0, 0],   # alpha in radians
                                     [0, 0, 1, 0],   # theta in radians
                                     [0, 0, 0, 1]])  # dimensionless pitch rate

        self.D_symmetric = np.array([[0, 0],
                                     [0, 0],
                                     [0, 0],
                                     [0, 0]])

    def state_space_asymmetric_eom(self):
        if not self.example:
            # Forces and moments coefficients.
            C_Y_beta, C_l_beta, C_n_beta = self.aircraft.derivatives_asymmetric[0, :]
            C_Y_beta_dot, C_l_beta_dot, C_n_beta_dot = self.aircraft.derivatives_asymmetric[1, :]
            C_Y_p, C_l_p, C_n_p = self.aircraft.derivatives_asymmetric[2, :]
            C_Y_r, C_l_r, C_n_r = self.aircraft.derivatives_asymmetric[3, :]
            C_Y_ailerons, C_l_ailerons, C_n_ailerons = self.aircraft.derivatives_asymmetric[4, :]
            C_Y_rudder, C_l_rudder, C_n_rudder = self.aircraft.derivatives_asymmetric[5, :]

            # Aerodynamic coefficients.
            velocity, C_L, = self.aircraft.velocity, - self.aircraft.derivatives_symmetric[0, 1]

            # Geometry coefficients.
            m_b, span, K_XX, K_ZZ, K_XZ = (self.aircraft.m_b, self.aircraft.wing.span,
                                           self.aircraft.K_XX, self.aircraft.K_ZZ, self.aircraft.K_XZ)

        else:
            # Asymmetric coefficients.
            C_Y_beta, C_l_beta, C_n_beta = self.asymmetric_coefficients[0, :]
            C_Y_beta_dot, C_l_beta_dot, C_n_beta_dot = self.asymmetric_coefficients[1, :]
            C_Y_p, C_l_p, C_n_p = self.asymmetric_coefficients[2, :]
            C_Y_r, C_l_r, C_n_r = self.asymmetric_coefficients[3, :]
            C_Y_ailerons, C_l_ailerons, C_n_ailerons = self.asymmetric_coefficients[4, :]
            C_Y_rudder, C_l_rudder, C_n_rudder = self.asymmetric_coefficients[5, :]

            # Miscellaneous coefficients.
            velocity, C_L, m_b, span, K_XX, K_ZZ, K_XZ = self.asymmetric_misc_coefficients

        # Setting up the matrices to transform to the state-space format.
        P = - (span / velocity) * np.array([[(2 * m_b) - C_Y_beta_dot, 0, 0, 0],
                                            [0, 1 / 2, 0, 0],
                                            [0, 0, 4 * m_b * K_XX, -4 * m_b * K_XZ],
                                            [- C_n_beta_dot, 0, -4 * m_b * K_XZ, 4 * m_b * K_ZZ]])

        Q = - np.array([[C_Y_beta, C_L, C_Y_p, C_Y_r - (4 * m_b)],
                        [0, 0, 1, 0],
                        [C_l_beta, 0, C_l_p, C_l_r],
                        [C_n_beta, 0, C_n_p, C_n_r]])

        R = - np.array([[C_Y_ailerons, C_Y_rudder],
                        [0, 0],
                        [C_l_ailerons, C_l_rudder],
                        [C_n_ailerons, C_n_rudder]])

        # Assigning state-space matrices to class attributes.
        self.A_asymmetric = np.dot(np.linalg.inv(P), Q)
        self.B_asymmetric = np.dot(np.linalg.inv(P), R)

        self.C_asymmetric = np.array([[1, 0, 0, 0],   # beta in radians
                                      [0, 1, 0, 0],   # roll in radians
                                      [0, 0, 1, 0],   # dimensionless roll rate
                                      [0, 0, 0, 1]])  # dimensionless yaw rate

        self.D_asymmetric = np.array([[0, 0],
                                      [0, 0],
                                      [0, 0],
                                      [0, 0]])

    def eigenvalues(self):
        # Calculating and assigning eigenvalues to class attributes.
        self.eigenvalues_symmetric = np.linalg.eigvals(self.A_symmetric)
        self.eigenvalues_asymmetric = np.linalg.eigvals(self.A_asymmetric)

        # Extracting unique eigenvalues.
        eigenvalues_symmetric = np.array([[i.real, abs(i.imag)] for i in self.eigenvalues_symmetric])
        eigenvalues_symmetric = np.unique(eigenvalues_symmetric, axis=0)

        eigenvalues_asymmetric = np.array([[i.real, abs(i.imag)] for i in self.eigenvalues_asymmetric])
        eigenvalues_asymmetric = np.unique(eigenvalues_asymmetric, axis=0)

        # Calculating period from unique eigenvalues and assign infinite if eigenvalue only has real part.
        self.period_symmetric = [2 * np.pi / i if i != 0
                                 else np.inf for i in eigenvalues_symmetric[:, 1]]
        self.period_asymmetric = [2 * np.pi / i if i != 0
                                  else np.inf for i in eigenvalues_asymmetric[:, 1]]

        # Calculating half amplitude time from unique eigenvalues.
        self.half_amplitude_symmetric = - np.log(2) / eigenvalues_symmetric[:, 0]
        self.half_amplitude_asymmetric = - np.log(2) / eigenvalues_asymmetric[:, 0]

    def responses(self, deflections, deflections_time, start=0, stop=200, step=0.01, symmetric=True):
        # Relevant aerodynamic & geometric coefficients.
        velocity, mac, span = self.aircraft.velocity, self.aircraft.wing.mac, self.aircraft.wing.span

        # Setting up the time array.
        time = np.arange(start, stop + step, step)

        # Setting up the input array.
        inputs_on = np.ones((len(deflections), round(deflections_time / step)))
        inputs_off = np.zeros((len(deflections), len(time) - len(inputs_on[0])))

        inputs = np.hstack((inputs_on, inputs_off))
        inputs = [inputs[i] * deflections[i] for i in range(len(deflections))]

        # Creating the control system.
        if symmetric:
            system = ctrl.ss(self.A_symmetric, self.B_symmetric, self.C_symmetric, self.D_symmetric)
            time, response = ctrl.forced_response(system, time, inputs)

            # Give responses their relevant dimension and starting value.
            response[0] = (response[0] * velocity) + velocity
            response[1] = response[1] * (180 / np.pi)
            response[2] = response[2] * (180 / np.pi)
            response[3] = response[3] * (180 / np.pi) * (velocity / mac)

        else:
            system = ctrl.ss(self.A_asymmetric, self.B_asymmetric, self.C_asymmetric, self.D_asymmetric)
            time, response = ctrl.forced_response(system, time, inputs)

            # Give responses their relevant dimension and starting value.
            response[0] = response[0] * (180 / np.pi)
            response[1] = response[1] * (180 / np.pi)
            response[2] = response[2] * (180 / np.pi) * ((2 * velocity) / span)
            response[3] = response[3] * (180 / np.pi) * ((2 * velocity) / span)

        # Give inputs_reference their relevant dimension.
        inputs[0] = inputs[0] * (180 / np.pi)
        inputs[1] = inputs[1] * (180 / np.pi)

        return time, inputs, response

    def plot_eigenvalues(self, symmetric=True, show=False):
        plt.close()

        if symmetric:
            part_real = [i.real for i in self.eigenvalues_symmetric]
            part_imaginary = [i.imag for i in self.eigenvalues_symmetric]
        else:
            part_real = [i.real for i in self.eigenvalues_asymmetric]
            part_imaginary = [i.imag for i in self.eigenvalues_asymmetric]

        plt.scatter(part_real, part_imaginary)
        plt.xlabel('Real [-]')
        plt.ylabel('Imaginary [-]')

        plt.axhline(0, color='black')
        plt.axvline(0, color='black')

        if symmetric:
            plt.savefig(f'Results/{self.aircraft.name}/Eigenvalues/symmetric')
        else:
            plt.savefig(f'Results/{self.aircraft.name}/Eigenvalues/asymmetric')

        if show:
            plt.show()

    def plot_responses(self, deflections, deflections_time, stop=30, symmetric=True, show=False):
        time, inputs, response = self.responses(deflections, deflections_time, stop=stop, symmetric=symmetric)

        plt.close()
        figure, axis = plt.subplots(2, 3, tight_layout=True)

        if symmetric:
            axis[0, 0].set(xlabel='Time [s]', ylabel='Velocity [m/s]', title='Velocity')
            axis[0, 1].set(xlabel='Time [s]', ylabel='Alpha [deg]', title='Alpha')
            axis[0, 2].set(xlabel='Time [s]', ylabel='Elevator [deg]', title='Elevator')
            axis[1, 0].set(xlabel='Time [s]', ylabel='Theta [deg]', title='Theta')
            axis[1, 1].set(xlabel='Time [s]', ylabel='Pitch rate [deg/s]', title='Pitch Rate')
            axis[1, 2].set(xlabel='Time [s]', ylabel='Trim [deg]', title='Trim')

        else:
            axis[0, 0].set(xlabel='Time [s]', ylabel='Beta [deg]', title='Beta')
            axis[0, 1].set(xlabel='Time [s]', ylabel='Roll [deg]', title='Roll')
            axis[0, 2].set(xlabel='Time [s]', ylabel='Ailerons [deg]', title='Ailerons')
            axis[1, 0].set(xlabel='Time [s]', ylabel='Roll rate [deg/s]', title='Roll Rate')
            axis[1, 1].set(xlabel='Time [s]', ylabel='Yaw rate [deg/s]', title='Yaw Rate')
            axis[1, 2].set(xlabel='Time [s]', ylabel='Rudder [deg]', title='Rudder')

        axis[0, 0].plot(time, response[0])
        axis[0, 1].plot(time, response[1])
        axis[0, 2].plot(time, inputs[0])
        axis[1, 0].plot(time, response[2])
        axis[1, 1].plot(time, response[3])
        axis[1, 2].plot(time, inputs[1])

        if symmetric:
            plt.savefig(f'Results/{self.aircraft.name}/Responses/symmetric')
        else:
            plt.savefig(f'Results/{self.aircraft.name}/Responses/asymmetric')

        if show:
            plt.show()


def main(path, interval=100, fuselage=True):
    solver = AircraftStability(path=path, fuselage=fuselage)

    solver.plot_eigenvalues(show=False)
    solver.plot_responses(deflections=[0.1, 0], deflections_time=2, stop=interval, show=True)

    solver.plot_eigenvalues(symmetric=False, show=False)
    solver.plot_responses(deflections=[0.1, 0], deflections_time=0.5, stop=interval / 10, symmetric=False, show=True)

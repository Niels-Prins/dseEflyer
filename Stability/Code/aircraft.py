import numpy as np
import pandas as pd

import os


class Aircraft:
    """
    For the usage, assumptions and limitations of this class consult the readme file.
    """
    def __init__(self, path, aircraft):
        # Creating the initial aircraft provided by user.
        self.aircraft = aircraft
        self.wing = aircraft[0]
        self.name = path.split('/')[-1]

        # Obtain general aircraft information provided by user.
        with open(f'{path}/aircraft.txt') as file:
            aircraft_data = np.genfromtxt(file, dtype=str)

        self.g = 9.80665
        self.mass = float(aircraft_data[0, 1])
        self.weight = self.mass * self.g

        self.I_XX = float(aircraft_data[1, 1])
        self.I_YY = float(aircraft_data[2, 1])
        self.I_ZZ = float(aircraft_data[3, 1])
        self.I_XZ = float(aircraft_data[4, 1])

        self.rho = float(aircraft_data[5, 1])
        self.velocity = float(aircraft_data[6, 1])

        # Dimensionless mass coefficients.
        self.m_c = self.mass / (self.rho * self.wing.area * self.wing.mac)
        self.m_b = self.mass / (self.rho * self.wing.area * self.wing.span)

        self.K_XX = self.I_XX / (self.mass * (self.wing.span ** 2))
        self.K_YY = self.I_YY / (self.mass * (self.wing.mac ** 2))
        self.K_ZZ = self.I_ZZ / (self.mass * (self.wing.span ** 2))
        self.K_XZ = self.I_XZ / (self.mass * (self.wing.span ** 2))

        # Inputs and outputs at equilibrium condition used for reference.
        self.inputs_reference = np.array([[0, 0, self.velocity, 0],
                                          [0, 0, 0, 0],
                                          [0, 0, self.velocity / 20, 0]])
        self.outputs_reference = self.calculate_outputs(self.inputs_reference, self.velocity)

        if not os.path.isdir(f'Results/{self.name}/Coefficients'):
            os.makedirs(f'Results/{self.name}/Coefficients')

        # Calculate stability derivatives.
        self.derivatives_symmetric = np.zeros((6, 3))
        self.derivatives_asymmetric = np.zeros((6, 3))

        self.calculate_symmetric_derivatives()
        self.calculate_asymmetric_derivatives()

    def calculate_outputs(self, inputs, magnitude, motion=None, controls=None):
        outputs_aircraft = np.zeros((3, 2))

        # Ensures controls remain uncoupled during different motions.
        for surface in self.aircraft:
            if surface.motion == motion:
                outputs_surface = surface.calculate_outputs(inputs, magnitude, self.rho, controls)
            else:
                outputs_surface = surface.calculate_outputs(inputs, magnitude, self.rho, controls=0)

            # Sum force outputs and add to combined body forces.
            outputs_aircraft[0, 0] += outputs_surface[0, 0] + outputs_surface[0, 1]
            outputs_aircraft[1, 0] += outputs_surface[1, 0] + outputs_surface[1, 1]
            outputs_aircraft[2, 0] += outputs_surface[2, 0] + outputs_surface[2, 1]

            # Add moment outputs to combined body moments.
            outputs_aircraft[0, 1] += outputs_surface[0, 2]
            outputs_aircraft[1, 1] += outputs_surface[1, 2]
            outputs_aircraft[2, 1] += outputs_surface[2, 2]

        return outputs_aircraft

    def calculate_derivatives(self, derivative, entry, inputs_change, c_inputs,
                              symmetric=False, controls=False, motion=None):
        inputs = self.inputs_reference.copy()

        if not controls:
            inputs[entry] += inputs_change
            controls = 0
        else:
            controls = inputs_change

        magnitude = np.linalg.norm(inputs[:, 2])

        outputs = self.calculate_outputs(inputs, magnitude, motion, controls)
        outputs_change = (outputs - self.outputs_reference) / inputs_change

        if symmetric:
            c_outputs_1 = (1 / 2) * self.rho * magnitude ** 2 * self.wing.area
            c_outputs_2 = (1 / 2) * self.rho * magnitude ** 2 * self.wing.area
            c_outputs_3 = (1 / 2) * self.rho * magnitude ** 2 * self.wing.area * self.wing.mac

            self.derivatives_symmetric[derivative, 0] = outputs_change[0, 0] / (c_inputs * c_outputs_1)
            self.derivatives_symmetric[derivative, 1] = outputs_change[2, 0] / (c_inputs * c_outputs_2)
            self.derivatives_symmetric[derivative, 2] = outputs_change[1, 1] / (c_inputs * c_outputs_3)

        else:
            c_outputs_1 = (1 / 2) * self.rho * magnitude ** 2 * self.wing.area
            c_outputs_2 = (1 / 2) * self.rho * magnitude ** 2 * self.wing.area * self.wing.span
            c_outputs_3 = (1 / 2) * self.rho * magnitude ** 2 * self.wing.area * self.wing.span

            self.derivatives_asymmetric[derivative, 0] = outputs_change[1, 0] / (c_inputs * c_outputs_1)
            self.derivatives_asymmetric[derivative, 1] = outputs_change[0, 1] / (c_inputs * c_outputs_2)
            self.derivatives_asymmetric[derivative, 2] = outputs_change[2, 1] / (c_inputs * c_outputs_3)

    def calculate_symmetric_derivatives(self):
        # Initial conditions (C_X_0, C_Z_0, C_M_Y_0).
        self.derivatives_symmetric[0, 0] = 0
        self.derivatives_symmetric[0, 1] = - (2 * self.weight) / (self.rho * (self.velocity ** 2) * self.wing.area)
        self.derivatives_symmetric[0, 2] = 0

        # Velocity derivatives (C_X_u, C_Z_u, C_M_Y_u).
        X_dot, c_X_dot = 1, 1 / self.velocity
        self.calculate_derivatives(1, (0, 2), X_dot, c_X_dot, symmetric=True)

        # Alpha derivatives (C_X_alpha, C_Z_alpha, C_M_Y_alpha).
        Z_dot, c_Z_dot = 1, 1 / self.velocity
        self.calculate_derivatives(2, (2, 2), Z_dot, c_Z_dot, symmetric=True)

        # Alpha dot derivatives (C_X_alpha_dot, C_Z_alpha_dot, C_M_Y_alpha_dot).
        Z_dot_dot, c_Z_dot_dot = 1, self.wing.mac / (self.velocity ** 2)
        self.calculate_derivatives(3, (2, 3), Z_dot_dot, c_Z_dot_dot, symmetric=True)

        # Pitch dot derivatives (C_X_pitch_dot, C_Z_pitch_dot, C_M_Y_pitch_dot).
        pitch_dot, c_pitch_dot = 0.1, self.wing.mac / self.velocity
        self.calculate_derivatives(4, (1, 1), pitch_dot, c_pitch_dot, symmetric=True)

        # Elevator effectiveness (C_X_elevator, C_Z_elevator, C_M_Y_elevator).
        elevator, c_elevator = 0.1, 1
        self.calculate_derivatives(5, None, elevator, c_elevator, controls=True, symmetric=True, motion='Y')

        row_labels = np.array(['0', 'X dot', 'alpha', 'alpha dot', 'pitch dot', 'elevator'])
        column_labels = np.array(['X', 'Z', 'M_Y'])

        dataframe = pd.DataFrame(np.round(self.derivatives_symmetric, 3), columns=column_labels, index=row_labels)
        dataframe.to_csv(f'Results/{self.name}/Coefficients/symmetric.csv')

        print(f'\nSymmetric derivatives: \n \n{dataframe}')

    def calculate_asymmetric_derivatives(self):
        # Beta derivatives (C_Y_beta, C_M_X_beta, C_M_Z_beta).
        Y_dot, c_Y_dot = 1, 1 / self.velocity
        self.calculate_derivatives(0, (1, 2), Y_dot, c_Y_dot)

        # Beta dot derivatives (C_Y_beta_dot, C_M_X_beta_dot, C_M_Z_beta_dot).
        Y_dot_dot, c_Y_dot_dot = 1, self.wing.span / (self.velocity ** 2)
        self.calculate_derivatives(1, (1, 3), Y_dot_dot, c_Y_dot_dot)

        # Roll dot derivatives (C_Y_roll_dot, C_M_X_roll_dot, C_M_Z_roll_dot).
        roll_dot, c_roll_dot = 0.1, self.wing.span / (2 * self.velocity)
        self.calculate_derivatives(2, (0, 1), roll_dot, c_roll_dot)

        # Yaw dot derivatives (C_Y_yaw_dot, C_M_X_yaw_dot, C_M_Z_yaw_dot).
        yaw_dot, c_yaw_dot = 0.1, self.wing.span / (2 * self.velocity)
        self.calculate_derivatives(3, (2, 1), yaw_dot, c_yaw_dot)

        # Ailerons effectiveness (C_Y_ailerons, C_M_X_ailerons, C_M_Z_ailerons).
        ailerons, c_ailerons = 0.1, 1
        self.calculate_derivatives(4, None, ailerons, c_ailerons, controls=True, motion='X')

        # Rudder effectiveness (C_Y_rudder, C_M_X_rudder, C_M_Z_rudder).
        rudder, c_rudder = 0.1, 1
        self.calculate_derivatives(5, None, rudder, c_rudder, controls=True, motion='Z')

        row_labels = np.array(['beta', 'beta dot', 'roll dot', 'yaw dot', 'ailerons', 'rudder'])
        column_labels = np.array(['Y', 'M_X', 'M_Z'])

        dataframe = pd.DataFrame(np.round(self.derivatives_asymmetric, 3), columns=column_labels, index=row_labels)
        dataframe.to_csv(f'Results/{self.name}/Coefficients/asymmetric.csv')

        print(f'\nAsymmetric derivatives: \n \n{dataframe}')

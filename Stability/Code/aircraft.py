import numpy as np
import pandas as pd

from scipy.optimize import fsolve

import os

from Stability.Code.surface import AerodynamicSurface


class Aircraft:

    def __init__(self, path, fuselage=True, example=False):
        self.example = example

        # Obtain general aircraft information provided by user.
        with open(f'{path}/aircraft.txt') as file:
            aircraft_data = np.genfromtxt(file, dtype=str)

        # Obtain flight conditions provided by user.
        with open(f'{path}/conditions.txt') as file:
            conditions_data = np.genfromtxt(file, dtype=str)

        self.g = 9.80665
        self.mass = float(aircraft_data[0, 1])
        self.weight = self.mass * self.g

        self.X_cg = float(aircraft_data[1, 1])
        self.Y_cg = float(aircraft_data[2, 1])
        self.Z_cg = float(aircraft_data[3, 1])

        self.I_XX = float(aircraft_data[4, 1])
        self.I_YY = float(aircraft_data[5, 1])
        self.I_ZZ = float(aircraft_data[6, 1])
        self.I_XZ = float(aircraft_data[7, 1])

        self.velocity = float(conditions_data[0, 1])
        self.altitude = float(conditions_data[1, 1])

        temperature = 15.04 + 273.1 - (0.00649 * self.altitude)
        pressure = 101.29 * (temperature / 288.08) ** 5.256

        self.rho = pressure / (0.2869 * temperature)
        self.magnitude = self.velocity

        surfaces = np.array([name for name in os.listdir(path) if os.path.isdir(f'{path}/{name}')])
        cg = [self.X_cg, self.Y_cg, self.Z_cg]

        self.wing = AerodynamicSurface(f'{path}/Wing', cg, symmetric=True, vertical=False)
        self.aircraft = [self.wing]
        self.name = path.split('/')[-1]

        self.C_L = (2 * self.weight) / (self.rho * self.velocity ** 2 * self.wing.area)
        self.alpha = np.interp(self.C_L, self.wing.C_L, self.wing.alpha)

        if fuselage:
            self.fuselage = Fuselage(path, self.wing)

            (self.C_L_correction, self.C_D_correction,
             self.C_M_correction, self.ac_correction, self.side_force_correction) = self.fuselage.wing_corrections()

            self.wing.C_L += self.C_L_correction * self.wing.alpha
            self.wing.C_D += self.C_D_correction
            self.wing.C_M += self.C_M_correction
            self.wing.C_Y_beta = self.side_force_correction

            self.wing.X_ac_nose -= self.ac_correction
            self.wing.X_ac_cg += self.ac_correction

            self.wing.interference_body = self.fuselage

        for surface in surfaces:
            if surface.lower() != 'wing':
                with open(f'{path}/{surface}/values.txt') as file:
                    surface_data = np.genfromtxt(file, dtype=str)

                symmetric = True if surface_data[0, 0] == 'True' else False
                vertical = True if surface_data[1, 1] == 'True' else False
                downwash = True if surface_data[2, 1] == 'True' else False
                sidewash = True if surface_data[3, 1] == 'True' else False

                if downwash:
                    new_surface = AerodynamicSurface(f'{path}/{surface}', cg,
                                                     symmetric=symmetric, vertical=vertical, downwash=self.wing)
                elif sidewash and fuselage:
                    new_surface = AerodynamicSurface(f'{path}/{surface}', cg,
                                                     symmetric=symmetric, vertical=vertical,
                                                     sidewash=[self.wing, self.fuselage])
                else:
                    new_surface = AerodynamicSurface(f'{path}/{surface}', cg,
                                                     symmetric=symmetric, vertical=vertical)
                self.aircraft.append(new_surface)

        # Dimensionless mass coefficients.
        self.m_c = self.mass / (self.rho * self.wing.area * self.wing.mac)
        self.m_b = self.mass / (self.rho * self.wing.area * self.wing.span)

        self.K_XX = self.I_XX / (self.mass * (self.wing.span ** 2))
        self.K_YY = self.I_YY / (self.mass * (self.wing.mac ** 2))
        self.K_ZZ = self.I_ZZ / (self.mass * (self.wing.span ** 2))
        self.K_XZ = self.I_XZ / (self.mass * (self.wing.span ** 2))

        # Inputs and outputs at equilibrium condition used for reference.
        self.inputs_reference = np.array([[0, 0, self.velocity * np.cos(self.alpha), 0],
                                          [0, 0, 0, 0],
                                          [0, 0, self.velocity * np.sin(self.alpha), 0]])
        self.outputs_reference = self.calculate_outputs(self.inputs_reference, self.velocity)

        if not self.example:
            if not os.path.isdir(f'Results/{self.name}/Coefficients'):
                os.makedirs(f'Results/{self.name}/Coefficients')

        # Calculate stability derivatives.
        self.derivatives_symmetric = np.zeros((6, 3))
        self.derivatives_asymmetric = np.zeros((6, 3))

        self.calculate_symmetric_derivatives()
        self.calculate_asymmetric_derivatives()
        self.calculate_rates()

    def calculate_outputs(self, inputs, motion=None, controls=None):
        outputs_aircraft = np.zeros((3, 2))

        # Ensures controls remain uncoupled during different motions.
        for surface in self.aircraft:
            if surface.motion == motion:
                outputs_surface = surface.calculate_outputs(inputs, self.magnitude, self.rho, controls)
            else:
                outputs_surface = surface.calculate_outputs(inputs, self.magnitude, self.rho, controls=0)

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

        outputs = self.calculate_outputs(inputs, motion, controls)
        outputs_change = (outputs - self.outputs_reference) / inputs_change

        if symmetric:
            c_outputs_1 = (self.rho / 2) * self.magnitude ** 2 * self.wing.area
            c_outputs_2 = (self.rho / 2) * self.magnitude ** 2 * self.wing.area
            c_outputs_3 = (self.rho / 2) * self.magnitude ** 2 * self.wing.area * self.wing.mac

            self.derivatives_symmetric[derivative, 0] = outputs_change[0, 0] / (c_inputs * c_outputs_1)
            self.derivatives_symmetric[derivative, 1] = outputs_change[2, 0] / (c_inputs * c_outputs_2)
            self.derivatives_symmetric[derivative, 2] = outputs_change[1, 1] / (c_inputs * c_outputs_3)

        else:
            c_outputs_1 = (self.rho / 2) * self.magnitude ** 2 * self.wing.area
            c_outputs_2 = (self.rho / 2) * self.magnitude ** 2 * self.wing.area * self.wing.span
            c_outputs_3 = (self.rho / 2) * self.magnitude ** 2 * self.wing.area * self.wing.span

            self.derivatives_asymmetric[derivative, 0] = outputs_change[1, 0] / (c_inputs * c_outputs_1)
            self.derivatives_asymmetric[derivative, 1] = outputs_change[0, 1] / (c_inputs * c_outputs_2)
            self.derivatives_asymmetric[derivative, 2] = outputs_change[2, 1] / (c_inputs * c_outputs_3)

    def calculate_symmetric_derivatives(self):
        # Initial conditions (C_X_0, C_Z_0, C_M_Y_0).
        self.derivatives_symmetric[0, 0] = 0
        self.derivatives_symmetric[0, 1] = - (2 * self.weight) / (self.rho * (self.velocity ** 2) * self.wing.area)
        self.derivatives_symmetric[0, 2] = 0

        # Velocity derivatives (C_X_u, C_Z_u, C_M_Y_u).
        X_dot, c_X_dot = 0.1, 1 / self.velocity
        self.calculate_derivatives(1, (0, 2), X_dot, c_X_dot, symmetric=True)

        # Alpha derivatives (C_X_alpha, C_Z_alpha, C_M_Y_alpha).
        Z_dot, c_Z_dot = 0.1, 1 / self.velocity
        self.calculate_derivatives(2, (2, 2), Z_dot, c_Z_dot, symmetric=True)

        # Alpha dot derivatives (C_X_alpha_dot, C_Z_alpha_dot, C_M_Y_alpha_dot).
        Z_dot_dot, c_Z_dot_dot = 0.1, self.wing.mac / (self.velocity ** 2)
        self.calculate_derivatives(3, (2, 3), Z_dot_dot, c_Z_dot_dot, symmetric=True)

        # Pitch dot derivatives (C_X_pitch_dot, C_Z_pitch_dot, C_M_Y_pitch_dot).
        pitch_dot, c_pitch_dot = 0.1, self.wing.mac / self.velocity
        self.calculate_derivatives(4, (1, 1), pitch_dot, c_pitch_dot, symmetric=True)

        # Elevator effectiveness (C_X_elevator, C_Z_elevator, C_M_Y_elevator).
        elevator, c_elevator = 0.1, 1
        self.calculate_derivatives(5, None, elevator, c_elevator, controls=True, symmetric=True, motion='Y')

        if not self.example:
            row_labels = np.array(['0', 'X dot', 'alpha', 'alpha dot', 'pitch dot', 'elevator'])
            column_labels = np.array(['X', 'Z', 'M_Y'])

            dataframe = pd.DataFrame(np.round(self.derivatives_symmetric, 3), columns=column_labels, index=row_labels)
            dataframe.to_csv(f'Results/{self.name}/Coefficients/symmetric.csv')

            print(f'\nSymmetric derivatives: \n \n{dataframe}')

    def calculate_asymmetric_derivatives(self):
        # Beta derivatives (C_Y_beta, C_M_X_beta, C_M_Z_beta).
        Y_dot, c_Y_dot = 0.1, 1 / self.velocity
        self.calculate_derivatives(0, (1, 2), Y_dot, c_Y_dot)

        # Fuselage contribution to C_Y_beta.
        if self.wing.C_Y_beta is not None:
            self.derivatives_asymmetric[0, 0] += self.wing.C_Y_beta

        # Beta dot derivatives (C_Y_beta_dot, C_M_X_beta_dot, C_M_Z_beta_dot).
        Y_dot_dot, c_Y_dot_dot = 0.1, self.wing.span / (self.velocity ** 2)
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

        if not self.example:
            row_labels = np.array(['beta', 'beta dot', 'roll dot', 'yaw dot', 'ailerons', 'rudder'])
            column_labels = np.array(['Y', 'M_X', 'M_Z'])

            dataframe = pd.DataFrame(np.round(self.derivatives_asymmetric, 3), columns=column_labels, index=row_labels)
            dataframe.to_csv(f'Results/{self.name}/Coefficients/asymmetric.csv')

            print(f'\nAsymmetric derivatives: \n \n{dataframe}')

    def calculate_rates(self, deflection=30):
        roll = (- 2 * (self.derivatives_asymmetric[4, 1] / self.derivatives_asymmetric[2, 1])
                * (self.velocity / self.wing.span) * deflection)
        pitch = ((self.derivatives_symmetric[5, 2] / self.derivatives_symmetric[4, 2])
                 * (self.velocity / self.wing.mac) * deflection)
        yaw = (2 * (self.derivatives_asymmetric[5, 2] / self.derivatives_asymmetric[3, 2])
               * (self.velocity / self.wing.span) * deflection)

        if not self.example:
            data = np.round(np.array([[roll], [pitch], [yaw]]))
            dataframe = pd.DataFrame(data, columns=['Rates [deg/s]'], index=['Roll', 'Pitch', 'Yaw'])

            print(f'\n{dataframe}')


class Fuselage:

    def __init__(self, path, wing):
        # Obtain fuselage information provided by user.
        with open(f'{path}/fuselage.txt') as file:
            data_fuselage = np.genfromtxt(file, dtype=str)

        self.nose_length = float(data_fuselage[0, 1])
        self.nose_height_start = float(data_fuselage[1, 1])
        self.nose_height_end = float(data_fuselage[2, 1])
        self.nose_height_slope = ((self.nose_height_end - self.nose_height_start) / self.nose_length)

        self.nose_width_start = float(data_fuselage[3, 1])
        self.nose_width_end = float(data_fuselage[4, 1])
        self.nose_width_slope = (self.nose_width_end - self.nose_width_start) / self.nose_length
        self.nose_sweep = float(data_fuselage[5, 1])

        self.main_length = float(data_fuselage[6, 1])
        self.main_height = self.nose_height_end
        self.main_width = self.nose_width_end

        self.tail_length = float(data_fuselage[7, 1])
        self.tail_height_start = self.main_height
        self.tail_height_end = float(data_fuselage[8, 1])
        self.tail_height_slope = ((self.tail_height_end - self.tail_height_start) / self.tail_length)

        self.tail_width_start = self.main_width
        self.tail_width_end = float(data_fuselage[9, 1])
        self.tail_width_slope = (self.tail_width_end - self.tail_width_start) / self.tail_length
        self.tail_sweep = float(data_fuselage[10, 1])

        self.length = self.nose_length + self.main_length + self.tail_length
        self.diameter = np.sqrt((4 * self.main_height * self.main_width) / np.pi)

        self.area_nose = ((self.nose_height_start * self.nose_width_start)
                          + self.nose_length * (self.nose_width_start + self.nose_width_end)
                          + self.nose_length * (self.nose_height_start + self.nose_height_end))
        self.area_main = 2 * self.main_length * (self.main_height + self.main_width)
        self.area_tail = ((self.tail_height_end * self.tail_width_end)
                          + self.tail_length * (self.tail_width_start + self.tail_width_end)
                          + self.tail_length * (self.tail_height_start + self.tail_height_end))

        self.area = self.area_nose + self.area_main + self.area_tail
        self.area_cross = (np.pi / 4) * self.diameter ** 2

        self.apparent_mass = 1 - ((10 * self.main_width) / (11 * self.length))

        # Fuselage attributes extended.
        self.wing = wing
        self.wing.area_net = self.wing.area - (self.main_width * self.wing.chord_root)

    def wing_corrections(self):

        def lift():
            C_L_alpha_nose = 2 * self.apparent_mass * ((self.main_width * self.main_height) / self.wing.area)

            K_nose = (C_L_alpha_nose / self.wing.C_L_alpha) * (self.wing.area / self.wing.area_net)

            K_wing = (0.1714 * (self.main_width / self.wing.span) ** 2
                      + 0.8326 * (self.main_width / self.wing.span) + 0.9974)

            K_fuselage = (0.7810 * (self.main_width / self.wing.span) ** 2
                          + 1.1976 * (self.main_width / self.wing.span) + 0.0088)

            C_L_alpha_wf = (self.wing.C_L_alpha * (K_nose + K_wing + K_fuselage)
                            * (self.wing.area_net / self.wing.area))

            return C_L_alpha_wf, C_L_alpha_nose, K_wing, K_fuselage

        def drag():
            interference_factor = 1.050
            friction_coefficient = 0.003

            tail_diameter = np.sqrt((4 * self.tail_width_end * self.tail_height_end) / np.pi)

            C_D_skin = (interference_factor * friction_coefficient
                        * (1 + (60 / ((self.length / self.diameter) ** 3))
                           + (0.0025 * (self.length / self.diameter)))
                        * (self.area / self.wing.area))

            C_D_base = ((0.029 * ((tail_diameter / self.diameter) ** 3))
                        / np.sqrt((C_D_skin * (self.wing.area / self.area_cross)))
                        * (self.area_cross / self.wing.area))

            return C_D_skin + C_D_base

        def aerodynamic_moment():
            nose_width = (self.nose_width_end + self.nose_width_start) / 2
            tail_width = (self.tail_width_end + self.tail_width_start) / 2

            C_L_signs = np.sign(self.wing.C_L)
            alpha = self.wing.alpha[np.where(C_L_signs > 0)[0][0]] * (180 / np.pi)

            c_m_nose = self.nose_length * (alpha + self.nose_sweep) * nose_width ** 2
            c_m_tail = self.tail_length * (alpha + self.tail_sweep) * tail_width ** 2

            c_m_ac_f = (self.apparent_mass / (36.5 * self.wing.area * self.wing.mac)) * (c_m_nose + c_m_tail)

            return c_m_ac_f

        def aerodynamic_center():
            C_L_alpha_wf, C_L_alpha_nose, K_wing, K_fuselage = lift()

            # All measured from the exposed root chord!
            nose_first = (2 * self.nose_height_slope * self.nose_width_slope
                            * self.nose_length ** 2 * ((abs(self.wing.X_le) / 2) - (self.nose_length / 3)))
            nose_second = (self.nose_length * (abs(self.wing.X_le) - (self.nose_length / 2))
                           * ((self.nose_height_start * self.nose_width_slope)
                              + (self.nose_width_start * self.nose_height_slope)))

            X_ac_nose = - (nose_first + nose_second) / (self.wing.chord_root * self.area_cross)

            # Should have sweep correction implemented!
            X_ac_wing = abs(self.wing.X_ac_nose - self.wing.X_le) / self.wing.chord_root
            C_L_alpha_wing_fuselage = self.wing.C_L_alpha * K_wing * (self.wing.area_net / self.wing.area)

            X_ac_fuselage = ((((5 * (self.wing.span - self.main_width)) / (21 * self.wing.chord_root))
                             * np.sqrt(self.main_width / self.wing.span) * np.tan(self.wing.sweep)) + 0.25)

            C_L_alpha_fuselage_wing = self.wing.C_L_alpha * K_fuselage * (self.wing.area_net / self.wing.area)

            X_ac_wf = ((((X_ac_nose * C_L_alpha_nose)
                         + (X_ac_wing * C_L_alpha_wing_fuselage)
                         + (X_ac_fuselage * C_L_alpha_fuselage_wing)) / C_L_alpha_wf) * self.wing.chord_root)

            return (self.wing.X_le - X_ac_wf) - self.wing.X_ac_nose

        def side_force():
            x = (self.tail_length / 2) - (self.length / 10)
            x_area = (self.main_height + (self.tail_height_slope * x)) * (self.main_width + (self.tail_width_slope * x))

            if self.wing.position == 0:
                interference_factor = 1.50
            elif self.wing.position == 1:
                interference_factor = 1.00
            else:
                interference_factor = 1.85

            return - interference_factor * (x_area / self.wing.area)

        C_L_alpha_correction = lift()[0] - self.wing.C_L_alpha
        C_D_correction = drag()
        C_M_correction = aerodynamic_moment()
        ac_correction = aerodynamic_center()

        side_force_correction = side_force()

        return C_L_alpha_correction, C_D_correction, C_M_correction, ac_correction, side_force_correction

    def velocity_correction(self, velocity, position, height=None, width=None):

        def kelvin():

            def equation_kelvin(ratio):
                return ((velocity[1, 0] * height) +
                        (((height ** 2 * ratio ** 2 + width ** 2) * velocity[1, 0])
                         / (2 * ratio * height)) * np.log((1 - ratio) / (1 + ratio)))

            # Use Scipy to solve for location and calculate required strength.
            location = fsolve(equation_kelvin, 0.5)[0] * height
            strength = (((location ** 2 + width ** 2) * np.pi * velocity[1, 0]) / location)

            # Determine the velocities in Y and Z directions.
            u = velocity[1, 0] + ((strength / (2 * np.pi)) * (((y - location) / (x ** 2 + (y - location) ** 2))
                                                              - ((y + location) / (x ** 2 + (y + location) ** 2))))
            v = (strength / (2 * np.pi)) * ((x / (x ** 2 + (y + location) ** 2)) - (x / (x ** 2 + (y - location) ** 2)))

            return u, v

        def rankine():

            def equation_rankine(ratio):
                return ((height ** 2 * ratio ** 2)
                        + ((height ** 2 * ratio) / (np.arctan(ratio))) - width ** 2)

            # Use Scipy to solve for location and calculate required strength.
            location = fsolve(equation_rankine, 1.0)[0] * height
            strength = (velocity[1, 0] * height * np.pi) / np.arctan(location / height)

            # Determine the velocities in Y and Z directions.
            u = velocity[1, 0] + ((strength / (2 * np.pi))) * (((x + location) / ((x + location) ** 2 + y ** 2))
                                                               - ((x - location) / ((x - location) ** 2 + y ** 2)))
            v = ((strength / (2 * np.pi))) * ((y / ((x + location) ** 2 + y ** 2))
                                              - (y / ((x - location) ** 2 + y ** 2)))

            return u, v

        def doublet():
            # Doublet does not require location solution and strength can be calculated directly.
            strength = 2 * np.pi * velocity[1, 0] * width ** 2

            # Determine the velocities in Y and Z directions.
            u = velocity[1, 0] + ((strength * (y ** 2 - x ** 2)) / (2 * np.pi * (x ** 2 + y ** 2) ** 2))
            v = - ((strength * x * y) / (np.pi * (x ** 2 + y ** 2) ** 2))

            return u, v

        # Else statement is only followed for verification, otherwise always if statement.
        if height is None and width is None:
            height = self.main_height / 2
            width = self.main_width / 2
            correction = True

            wing_locations = [- height, 0, height]
            y = wing_locations[self.wing.position]
            x = position

            if abs(y) > 0:
                if self.main_height / self.main_width > 1:
                    Y_dot, Z_dot = kelvin()

                elif self.main_height / self.main_width < 1:
                    Y_dot, Z_dot = rankine()

                else:
                    Y_dot, Z_dot = doublet()
            else:
                Y_dot, Z_dot = 0, 0

        else:
            x, y = position
            correction = False

            if height / width > 1:
                Y_dot, Z_dot = kelvin()

            elif height / width < 1:
                Y_dot, Z_dot = rankine()

            else:
                Y_dot, Z_dot = doublet()

        # Scaling in case solution is larger than original components.
        if abs(Y_dot) > abs(velocity[1, 0]):
            Y_dot = (Y_dot / abs(Y_dot)) * velocity[1, 0]
        if abs(Z_dot) > abs(velocity[1, 0]):
            Z_dot = (Z_dot / abs(Z_dot)) * velocity[1, 0]

        # Empirical scaling based on validation.
        if correction:
            Y_dot /= 2
            Z_dot /= 2

        # Create new array with changed velocities (verification array is mutable).
        velocity_x = velocity[0, 0]
        velocity_y = Y_dot
        velocity_z = velocity[2, 0] - Z_dot

        return np.array([[velocity_x], [velocity_y], [velocity_z]])

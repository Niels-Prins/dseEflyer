import numpy as np


class AerodynamicSurface:
    """
    For the usage, assumptions and limitations of this class consult the readme file.
    """
    def __init__(self, path, symmetric=True, vertical=False, downwash_body=None, sidewash_body=None):

        # Load polars data provided by user.
        with open(f'{path}/polars.txt') as file:
            self.data_polars = np.genfromtxt(file, skip_header=True)

        # Load values data provided by user.
        with open(f'{path}/values.txt') as file:
            self.data_values = np.genfromtxt(file, dtype=str)

        # Attributes loaded from the values file.
        self.chord_root = float(self.data_values[0, 1])
        self.span = float(self.data_values[1, 1])
        self.sweep = float(self.data_values[2, 1]) * (np.pi / 180)
        self.taper = float(self.data_values[3, 1])
        self.dihedral = float(self.data_values[4, 1]) * (np.pi / 180)

        self.X_le = float(self.data_values[5, 1])
        self.Y_le = float(self.data_values[6, 1])
        self.Z_le = float(self.data_values[7, 1])

        self.start_control = float(self.data_values[8, 1])
        self.end_control = float(self.data_values[9, 1])
        self.length_control = self.end_control - self.start_control

        self.motion = self.data_values[10, 1]
        self.deflection = True if self.data_values[11, 1] == 'True' else False

        self.alpha_controls = float(self.data_values[12, 1])
        self.moment_controls = float(self.data_values[13, 1]) * (180 / np.pi)

        # Attributes loaded from the polars file.
        self.alpha = self.data_polars[:, 0] * (np.pi / 180)
        self.C_L = self.data_polars[:, 2]
        self.C_D = self.data_polars[:, 5]
        self.C_M = np.mean(self.data_polars[:, 8])

        # Attributes from class input.
        self.symmetric = symmetric
        self.vertical = vertical

        self.downwash_body = downwash_body
        self.sidewash_body = sidewash_body

        # Attributes calculated by the initialization function.
        self.area = None
        self.mac = None
        self.downwash = None

        self.X_ac = None
        self.Y_ac = None
        self.Z_ac = None

        self.initialize()

        # Attributes in the body frame calculated by the transformation function.
        self.angles = None
        self.rates = None
        self.velocity = None
        self.acceleration = None

        # Attributes related to dynamic pressure.
        self.magnitude = None
        self.rho = None

        # Attributes calculated and transformed by the forces and moments functions.
        self.outputs = np.zeros((3, 3))

        # Code below only used for verification tests.
        self.calculate_local_velocity = None
        self.calculate_local_forces = None

    def initialize(self):

        def area():
            self.area = (self.span / 2) * self.chord_root * (1 + self.taper)

        def positioning():
            # Local Y-coordinate of aerodynamic center.
            if self.symmetric:
                Y_0 = - (self.span * (1 + (2 * self.taper))) / (6 * (1 + self.taper))
            else:
                Y_0 = - (self.span * (1 + (2 * self.taper))) / (3 * (1 + self.taper))

            # Initial local X and Z-coordinates of aerodynamic center.
            X_0 = - (self.chord_root / 2)
            Z_0 = 0

            # Mean aerodynamic chord located at Y-coordinate of aerodynamic center.
            self.mac = self.chord(Y_0)

            # Shift due to geometry and position of aerodynamic center along the airfoil.
            X_shift_due_ac = self.mac / 4
            X_shift_due_sweep = - abs(Y_0) * np.sin(self.sweep)
            Z_shift_due_dihedral = - abs(Y_0) * self.dihedral

            # Summing effects with the initial coordinates.
            X_final = X_0 + X_shift_due_ac + X_shift_due_sweep
            Y_final = Y_0
            Z_final = Z_0 + Z_shift_due_dihedral

            # Transform position in case of a vertical surface which is not in the body reference frame.
            if self.vertical:
                self.X_le, self.Y_le, self.Z_le = self.transformation(np.array([[self.X_le], [self.Y_le], [self.Z_le]]),
                                                                      transform_vertical=True)[:, 0]

            # Sum with the leading edge to get the aerodynamic center with respect to the center of gravity.
            self.X_ac = self.X_le + X_final
            self.Y_ac = self.Y_le + Y_final
            self.Z_ac = self.Z_le + Z_final

        def downwash():
            self.downwash = 0

            # Downwash only applies if behind downwash body.
            if self.downwash_body is not None and self.X_ac < self.downwash_body.X_ac:
                C_L_alpha = ((self.downwash_body.C_L[-1] - self.downwash_body.C_L[0]) /
                             (self.downwash_body.alpha[-1] - self.downwash_body.alpha[0]))
                ar = self.downwash_body.span / self.downwash_body.mac

                # Empirical downwash estimation.
                nominator = (7 * C_L_alpha)
                denominator = (4 * np.pi * ar * (1 + abs(self.Z_ac - self.downwash_body.Z_ac)) *
                               (self.downwash_body.taper * abs(self.X_ac - self.downwash_body.X_ac)) ** (1 / 4))

                self.downwash = nominator / denominator

        area()
        positioning()
        downwash()

    def chord(self, coordinate):
        if self.symmetric:
            return self.chord_root - self.chord_root * (1 - self.taper) * ((2 * abs(coordinate)) / self.span)
        else:
            return self.chord_root + self.chord_root * (1 - self.taper) * (coordinate / self.span)

    def transformation(self, inputs, transform_vertical=False, earth_to_body=True):

        def transformation_vertical(rotation=1):
            return np.array([[1, 0, 0],
                             [0, 0, rotation],
                             [0, - rotation, 0]])

        def transformation_roll(roll):
            return np.array([[1, 0, 0],
                             [0, np.cos(roll), np.sin(roll)],
                             [0, - np.sin(roll), np.cos(roll)]])

        def transformation_pitch(pitch):
            return np.array([[np.cos(pitch), 0, - np.sin(pitch)],
                             [0, 1, 0],
                             [np.sin(pitch), 0, np.cos(pitch)]])

        def transformation_yaw(yaw):
            return np.array([[np.cos(yaw), np.sin(yaw), 0],
                             [- np.sin(yaw), np.cos(yaw), 0],
                             [0, 0, 1]])

        # Check if vertical transformation applies due to surface rotation.
        if transform_vertical:
            if earth_to_body:
                transformation_matrix = transformation_vertical()
            else:
                transformation_matrix = transformation_vertical(rotation=-1)

        # Check if transformation from the Earth to the body reference frame applies (order: roll, pitch, yaw).
        elif earth_to_body:
            transformation_matrix = np.dot(transformation_roll(self.angles[0, 0]),
                                           transformation_pitch(self.angles[1, 0]))
            transformation_matrix = np.dot(transformation_matrix, transformation_yaw(self.angles[2, 0]))

        # Check if transformation from the body to the Earth reference frame applies (order: yaw, pitch, roll).
        else:
            transformation_matrix = np.dot(transformation_yaw(- self.angles[2, 0]),
                                           transformation_pitch(- self.angles[1, 0]))
            transformation_matrix = np.dot(transformation_matrix, transformation_roll(- self.angles[0, 0]))

        return np.dot(transformation_matrix, inputs)

    def calculate_forces(self, controller):

        def roll_dot_effect(coordinate):
            return coordinate * self.rates[0, 0]

        def pitch_dot_effect():
            return - self.X_ac * self.rates[1, 0]

        def yaw_dot_effect(coordinate):
            return - coordinate * self.rates[2, 0]

        def calculate_local_velocity(coordinate):
            # Local effects due to turn rates.
            X_dot_due_yaw_dot = yaw_dot_effect(coordinate)
            Z_dot_due_roll_dot = roll_dot_effect(coordinate)
            Z_dot_due_pitch_dot = pitch_dot_effect()

            # Local effects due to dihedral.
            Y_dot_due_dihedral = self.velocity[2, 0] * self.dihedral * - (coordinate / abs(coordinate))
            Z_dot_due_dihedral = self.velocity[1, 0] * self.dihedral * (coordinate / abs(coordinate))

            X_dot = self.velocity[0, 0] + X_dot_due_yaw_dot
            Y_dot = self.velocity[1, 0] + Y_dot_due_dihedral
            Z_dot = self.velocity[2, 0] + Z_dot_due_roll_dot + Z_dot_due_pitch_dot + Z_dot_due_dihedral

            # Local velocity and local angles.
            local_velocity = np.array([[X_dot], [Y_dot], [Z_dot]])
            local_magnitude = np.linalg.norm(local_velocity)

            local_alpha = np.arctan(local_velocity[2, 0] / local_velocity[0, 0])
            local_alpha_dot = (np.arctan((local_velocity[2, 0] + self.acceleration[2, 0]) /
                                         (local_velocity[0, 0] + self.acceleration[0, 0])) - local_alpha)

            local_beta = np.arctan(local_velocity[1, 0] / local_velocity[0, 0])
            local_beta_dot = (np.arctan((local_velocity[1, 0] + self.acceleration[1, 0]) /
                                        (local_velocity[0, 0] + self.acceleration[0, 0])) - local_beta)

            return local_velocity, local_magnitude, local_alpha, local_alpha_dot, local_beta, local_beta_dot

        def calculate_local_forces(coordinate):
            # Local velocities.
            local_velocity, local_magnitude, local_alpha, local_alpha_dot, local_beta, local_beta_dot = \
                calculate_local_velocity(coordinate)

            # Local geometry.
            local_chord = self.chord(coordinate)

            # Local downwash.
            if self.downwash_body is not None:
                local_delay = local_alpha_dot * ((self.downwash_body.X_ac - self.X_ac) / local_velocity[0, 0])
                local_downwash = (local_alpha - local_delay) * self.downwash
            else:
                local_downwash = 0

            # Local sidewash.
            # TODO: implement sidewash.

            # Sidewash correction to Y-component velocity.
            # TODO: implement velocity tilt near fuselage.

            # Correct lift and drag coefficients for control input.
            if self.start_control <= abs(coordinate) <= self.end_control:
                if not self.deflection and coordinate > 0:
                    local_alpha -= controller * self.alpha_controls
                else:
                    local_alpha += controller * self.alpha_controls

            # Local lift and drag coefficients.
            local_C_L = np.interp(local_alpha - local_downwash, self.alpha, self.C_L)
            local_C_D = np.interp(local_alpha - local_downwash, self.alpha, self.C_D)

            # Local lift and drag forces.
            local_L = (self.rho / 2) * local_magnitude ** 2 * local_C_L * local_chord * \
                      (1 + ((coordinate / abs(coordinate)) * local_beta * np.sin(2 * self.sweep)))
            local_D = (self.rho / 2) * local_magnitude ** 2 * local_C_D * local_chord * \
                      (1 + ((coordinate / abs(coordinate)) * local_beta * np.sin(2 * self.sweep)))

            # Local normal and tangential forces.
            local_F_T = np.sin(local_alpha) * local_L - np.cos(local_alpha) * local_D
            local_F_N = np.cos(local_alpha) * local_L + np.sin(local_alpha) * local_D

            # Local forces in the body reference frame.
            F_X = local_F_T
            F_Y = local_F_N * self.dihedral * - (coordinate / abs(coordinate))
            F_Z = - local_F_N

            return np.array([F_X, F_Y, F_Z])

        # Riemann sum of local forces.
        step_size = 0.001

        if self.symmetric:
            steps = int((self.span - self.span % step_size) / (2 * step_size)) + 1

            # Separate into left and right part to facilitate moment calculation later.
            for i in range(-steps, 0):
                step = (i * step_size) + (step_size / 2)
                self.outputs[:, 0] += calculate_local_forces(step) * step_size

            for i in range(0, steps):
                step = (i * step_size) + (step_size / 2)
                self.outputs[:, 1] += calculate_local_forces(step) * step_size

        else:
            steps = int((self.span - self.span % step_size) / step_size) + 1

            for i in range(-steps, 0):
                step = (i * step_size) + (step_size / 2)
                self.outputs[:, 0] += calculate_local_forces(step) * step_size

        # Code below only used for verification tests.
        self.calculate_local_velocity = calculate_local_velocity(coordinate=self.span / 2)
        self.calculate_local_forces = calculate_local_forces(coordinate=self.span / 2)

    def calculate_moments(self, controller):
        # Moments due to forces.
        M_Y_due_F_X = self.Z_ac * (self.outputs[0, 0] + self.outputs[0, 1])
        M_Z_due_F_X = self.Y_ac * (self.outputs[0, 1] - self.outputs[0, 0])

        M_X_due_F_Y = self.Z_ac * (self.outputs[1, 0] + self.outputs[1, 1]) * - 1
        M_Z_due_F_Y = self.X_ac * (self.outputs[1, 0] + self.outputs[1, 1])

        M_X_due_F_Z = self.Y_ac * (self.outputs[2, 0] - self.outputs[2, 1])
        M_Y_due_F_Z = self.X_ac * (self.outputs[2, 0] + self.outputs[2, 1]) * - 1

        # Correct aerodynamic center moment coefficient for control input.
        if self.deflection:
            if self.symmetric:
                C_M_Y_ac = self.C_M + (((2 * self.length_control) / self.span) * self.moment_controls * controller)
            else:
                C_M_Y_ac = self.C_M + ((self.length_control / self.span) * self.moment_controls * controller)
        else:
            C_M_Y_ac = self.C_M

        # Aerodynamic center moment.
        M_Y_due_ac = (C_M_Y_ac / 2) * self.rho * self.magnitude ** 2 * self.area * self.mac

        # Sum all moments and assign to output attribute.
        self.outputs[0, 2] = M_X_due_F_Y + M_X_due_F_Z
        self.outputs[1, 2] = M_Y_due_F_X + M_Y_due_F_Z + M_Y_due_ac
        self.outputs[2, 2] = M_Z_due_F_X + M_Z_due_F_Y

    def calculate_outputs(self, inputs, magnitude, rho, controls):
        # Reset outputs and assign inputs to class attributes.
        self.outputs = np.zeros((3, 3))
        self.magnitude = magnitude
        self.rho = rho

        # Check if vertical transformation from the body frame applies.
        if self.vertical:
            inputs = self.transformation(inputs, transform_vertical=True)

        # Assign transformed inputs to class attributes.
        self.angles = inputs[:, 0].reshape(-1, 1)
        self.rates = inputs[:, 1].reshape(-1, 1)

        self.velocity = self.transformation(inputs[:, 2].reshape(-1, 1))
        self.acceleration = self.transformation(inputs[:, 3].reshape(-1, 1))

        # Calculate outputs.
        self.calculate_forces(controls)
        self.calculate_moments(controls)

        # Check if vertical transformation to the body frame applies.
        if self.vertical:
            self.outputs = self.transformation(self.outputs, transform_vertical=True, earth_to_body=False)

        return self.outputs

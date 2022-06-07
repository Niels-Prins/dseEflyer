import numpy as np


class AerodynamicSurface:

    def __init__(self, path, cg, symmetric=True, vertical=False,
                 downwash=None, sidewash=None, interference=None, step=0.01):
        # Load polars data provided by user.
        with open(f'{path}/polars.txt') as file:
            self.data_polars = np.genfromtxt(file, skip_header=True)

        # Load values data provided by user.
        with open(f'{path}/values.txt') as file:
            self.data_values = np.genfromtxt(file, dtype=str)

        # Attributes loaded from the values file.
        self.position = int(self.data_values[4, 1])
        self.chord_root = float(self.data_values[5, 1])
        self.span = float(self.data_values[6, 1])
        self.sweep = float(self.data_values[7, 1]) * (np.pi / 180)
        self.taper = float(self.data_values[8, 1])
        self.dihedral = float(self.data_values[9, 1]) * (np.pi / 180)
        self.incidence = float(self.data_values[10, 1]) * (np.pi / 180)

        self.X_le = float(self.data_values[11, 1])
        self.Y_le = float(self.data_values[12, 1])
        self.Z_le = float(self.data_values[13, 1])

        self.start_control = float(self.data_values[14, 1])
        self.end_control = float(self.data_values[15, 1])
        self.length_control = self.end_control - self.start_control

        self.motion = self.data_values[16, 1]
        self.deflection = True if self.data_values[17, 1] == 'True' else False

        self.alpha_controls = float(self.data_values[18, 1])
        self.moment_controls = float(self.data_values[19, 1]) * (180 / np.pi)

        self.area = (self.span / 2) * self.chord_root * (1 + self.taper)
        self.ar = (self.span ** 2) / self.area

        # Attributes loaded from the polars file.
        self.alpha = self.data_polars[:, 0] * (np.pi / 180)
        self.C_L = self.data_polars[:, 2]
        self.C_D = self.data_polars[:, 5]
        self.C_M = np.mean(self.data_polars[:, 8])

        alpha_step = (self.alpha[1] - self.alpha[0]) * (180 / np.pi)

        self.C_L_alpha = np.mean(np.gradient(self.C_L)) * (180 / (alpha_step * np.pi))
        self.C_D_alpha = np.mean(np.gradient(self.C_D)) * (180 / (alpha_step * np.pi))

        # Attributes from class input.
        self.symmetric = symmetric
        self.vertical = vertical
        self.X_cg, self.Y_cg, self.Z_cg = cg

        self.downwash_body = downwash
        self.interference_body = interference

        if sidewash is not None:
            self.sidewash_wing, self.sidewash_fuselage = sidewash
        else:
            self.sidewash_wing, self.sidewash_fuselage = None, None

        self.step = step

        # Attributes calculated by the initialization function.
        self.mac = None
        self.downwash = None
        self.sidewash = None

        self.X_ac_cg = None
        self.Y_ac_cg = None
        self.Z_ac_cg = None

        self.X_ac_nose = None
        self.Y_ac_nose = None
        self.Z_ac_nose = None

        self.initialize()

        # Attributes in the body frame calculated by the transformation function.
        self.angles = None
        self.rates = None
        self.velocity = None
        self.acceleration = None

        # Attributes related to dynamic pressure.
        self.magnitude = None
        self.rho = None

        # Attributes used for fuselage corrections.
        self.C_Y_beta = None

        # Attributes calculated and transformed by the forces and moments functions.
        self.outputs = np.zeros((3, 3))

        # Code below only used for verification tests.
        self.calculate_local_velocity = None
        self.calculate_local_forces = None

    def initialize(self):

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

            # Transform position in case of a vertical aircraft which is not in the body reference frame.
            if self.vertical:
                self.X_le, self.Y_le, self.Z_le = self.transformation(np.array([[self.X_le], [self.Y_le], [self.Z_le]]),
                                                                      transform_vertical=True)[:, 0]
                self.X_cg, self.Y_cg, self.Z_cg = self.transformation(np.array([[self.X_cg], [self.Y_cg], [self.Z_cg]]),
                                                                      transform_vertical=True)[:, 0]

            # Aerodynamic center with respect to the center of gravity.
            self.X_ac_cg = self.X_le - self.X_cg + X_final
            self.Y_ac_cg = self.Y_le - self.Y_cg + Y_final
            self.Z_ac_cg = self.Z_le - self.Z_cg + Z_final

            # Aerodynamic center with respect to the aircraft nose.
            self.X_ac_nose = self.X_le + X_final
            self.Y_ac_nose = self.Y_le + Y_final
            self.Z_ac_nose = self.Z_le + Z_final

        def downwash():
            self.downwash = 0

            # Downwash only applies if behind downwash body.
            if self.downwash_body is not None and self.X_ac_cg < self.downwash_body.X_ac_cg:
                d_x = abs(self.X_ac_cg - self.downwash_body.X_ac_cg)
                d_z = abs(self.Z_ac_cg - self.downwash_body.Z_ac_cg)

                K_a = (1 / self.downwash_body.ar) - (1 / (1 + self.downwash_body.ar ** 1.7))
                K_lambda = (10 - (3 * self.downwash_body.taper)) / 7
                K_h = (1 - (d_z / self.downwash_body.span)) / (((2 * d_x) / self.downwash_body.span) ** (1 / 3))

                # Empirical relation for downwash based on previous parameters.
                self.downwash = 4.44 * (K_a * K_lambda * K_h * np.sqrt(np.cos(self.downwash_body.sweep))) ** 1.19

        def sidewash():
            self.sidewash = 0

            # Sidewash only applies if behind sidewash wing and within the fuselage width.
            if self.sidewash_wing is not None and self.X_ac_cg < self.sidewash_wing.X_ac_cg \
                    and abs(self.Z_ac_cg) < self.sidewash_fuselage.main_width:
                d_z_options = [self.sidewash_fuselage.main_height / 2, 0, -self.sidewash_fuselage.main_height / 2]
                d_z = d_z_options[self.sidewash_wing.position]

                # Empirical relation for sidewash based on previous parameters.
                self.sidewash = (0.724 + ((3.06 / (1 + np.cos(self.sidewash_wing.sweep)))
                                 * (self.area / self.sidewash_wing.area))
                                 + ((0.4 * d_z) / self.sidewash_fuselage.main_height)
                                 + 0.009 * self.sidewash_wing.ar)

        positioning()
        downwash()
        sidewash()

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

        # Check if vertical transformation applies due to aircraft rotation.
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
            return - self.X_ac_cg * self.rates[1, 0]

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

            # Local interference from fuselage (wing only).
            if self.interference_body is not None:
                local_velocity = self.interference_body.velocity_correction(local_velocity, coordinate)

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

            # Local downwash (horizontal tail only).
            if self.downwash_body is not None and not self.vertical:
                local_delay = local_alpha_dot * ((self.downwash_body.X_ac_cg - self.X_ac_cg) / local_velocity[0, 0])
                local_downwash = (local_alpha - local_delay) * self.downwash

            else:
                local_downwash = 0

            # Local sidewash (vertical tail only).
            if self.sidewash_wing is not None and self.vertical:
                local_delay = local_alpha_dot * ((self.sidewash_wing.X_ac_cg - self.X_ac_cg) / local_velocity[0, 0])
                local_alpha = (local_alpha - local_delay) * self.sidewash

            local_alpha += self.incidence

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
            local_L = ((self.rho / 2) * local_magnitude ** 2 * local_C_L * local_chord
                       * (1 + ((coordinate / abs(coordinate)) * local_beta * np.sin(2 * self.sweep))))
            local_D = ((self.rho / 2) * local_magnitude ** 2 * local_C_D * local_chord
                       * (1 + ((coordinate / abs(coordinate)) * local_beta * np.sin(2 * self.sweep))))

            # Local normal and tangential forces.
            local_F_T = np.sin(local_alpha) * local_L - np.cos(local_alpha) * local_D
            local_F_N = np.cos(local_alpha) * local_L + np.sin(local_alpha) * local_D

            # Local forces in the body reference frame.
            F_X = local_F_T
            F_Y = local_F_N * self.dihedral * - (coordinate / abs(coordinate))
            F_Z = - local_F_N

            return np.array([F_X, F_Y, F_Z])

        # Riemann sum of local forces.
        step_size = self.step

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
        M_Y_due_F_X = self.Z_ac_cg * (self.outputs[0, 0] + self.outputs[0, 1])
        M_Z_due_F_X = self.Y_ac_cg * (self.outputs[0, 0] - self.outputs[0, 1]) * - 1

        M_X_due_F_Y = self.Z_ac_cg * (self.outputs[1, 0] + self.outputs[1, 1]) * - 1
        M_Z_due_F_Y = self.X_ac_cg * (self.outputs[1, 0] + self.outputs[1, 1])

        M_X_due_F_Z = self.Y_ac_cg * (self.outputs[2, 0] - self.outputs[2, 1])
        M_Y_due_F_Z = self.X_ac_cg * (self.outputs[2, 0] + self.outputs[2, 1]) * - 1

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

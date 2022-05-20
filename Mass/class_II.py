import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from class_I import class_I

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


class MassMethods:

    def __init__(self):
        # Class I attributes.
        self.mass_empty_1, self.mass_takeoff_1, self.mass_battery, self.mass_motor, self.mass_occupants = class_I()

        # Load factor attributes.
        self.safety_margin = 1.5
        self.load_factor = 8
        self.load_factor_ultimate = self.load_factor * self.load_factor

        # Wing attributes.
        self.wing_aspect_ratio = 5.80
        self.wing_area = 13.10
        self.wing_span = 8.72
        self.wing_sweep_half = -7.45 * (np.pi / 180)
        self.wing_sweep_quarter = -3.74 * (np.pi / 180)
        self.wing_taper = 0.45
        self.wing_t_to_c = 0.15
        self.wing_t_max = 0.31
        self.wing_chord_root = 2.07
        self.wing_MAC = ((2 / 3) * self.wing_chord_root *
                         ((1 + self.wing_taper + self.wing_taper ** 2) / (1 + self.wing_taper)))
        self.wing_C_L_alpha = 6.56
        self.wing_X_ac = 0.25

        # Horizontal tail attributes.
        self.h_tail_aspect_ratio = 5.15  # 5.00
        self.h_tail_area = 2.79  # 2.00
        self.h_tail_span = 3.79  # 3.80
        self.h_tail_sweep_half = 0
        self.h_tail_sweep_quarter = 1.78 * (np.pi / 180)
        self.h_tail_taper = 0.73  # 0.70
        self.h_tail_t_max = 0.12
        self.h_tail_C_L_alpha = None

        # Vertical tail attributes.
        self.v_tail_aspect_ratio = 1.65
        self.v_tail_area = 1.31
        self.v_tail_span = 1.47
        self.v_tail_sweep_quarter = 27 * (np.pi / 180)
        self.v_tail_taper = 0.45
        self.v_tail_t_max = 0.12
        self.v_tail_C_L_alpha = None

        # Fuselage attributes.
        self.fuselage_area = 21.00
        self.fuselage_length = 6.50
        self.fuselage_height = 1.58
        self.fuselage_width = 0.80
        self.fuselage_radius = max(self.fuselage_height, self.fuselage_width) / 2

        # Gear attributes.
        self.gear_length = 0.50
        self.gear_load_factor = 5.5

        # Empirical attributes.
        self.h_tail_arm = 0.40 * self.wing_span

        # Conversion factors.
        self.kg_to_pounds = 2.2046
        self.meters_to_feet = 3.2808
        self.pascal_to_psf = 0.0209
        self.ms_to_knots = 1 / 0.514444444

        # Flight condition attributes.
        self.velocity = 77
        self.density = 1.225
        self.pressure = (self.density / 2) * self.velocity ** 2

        # Mass fractions.
        self.mass_bat_1 = (self.mass_battery * 2) / 4
        self.mass_bat_2 = (self.mass_battery * 2) / 4

        self.mass_occupant_1 = self.mass_occupants / 2
        self.mass_occupant_2 = self.mass_occupants / 2

        # Fuselage group moment arms as fractions of fuselage length.
        self.arm_h_tail = 0.95 * self.fuselage_length
        self.arm_v_tail = 0.95 * self.fuselage_length
        self.arm_fuselage = 0.35 * self.fuselage_length
        self.arm_bat_1 = 0.10 * self.fuselage_length
        self.arm_bat_2 = 0.60 * self.fuselage_length
        self.arm_motor = 0.70 * self.fuselage_length
        self.arm_control = 0.60 * self.fuselage_length
        self.arm_electric = 0.30 * self.fuselage_length
        self.arm_misc = 0.70 * self.fuselage_length
        self.arm_occupant_1 = 0.30 * self.fuselage_length
        self.arm_occupant_2 = 0.50 * self.fuselage_length

        # Wing group moment arms as fraction of MAC.
        self.arm_wing = 0.40 * self.wing_MAC
        self.arm_gear = 0.50 * self.wing_MAC
        self.arm_EOM = 0.35 * self.wing_MAC

        # To be calculated attributes.
        self.mass_empty_2 = None
        self.mass_takeoff_2 = None

        self.wing_X_LE_correction = 0.2
        self.wing_X_LE = None
        self.aircraft_X_CG = None
        self.aircraft_C_L_alpha = None

    def cessna(self):
        mass_wing = ((self.load_factor_ultimate * self.mass_takeoff_1 * self.kg_to_pounds) ** 0.397
                     * (self.wing_area * self.meters_to_feet ** 2) ** 0.360
                     * self.wing_aspect_ratio ** 1.712
                     * (0.04674 / self.kg_to_pounds))

        mass_h_tail = ((self.mass_takeoff_1 * self.kg_to_pounds) ** 0.887
                       * (self.h_tail_area * self.meters_to_feet ** 2) ** 0.101
                       * ((self.h_tail_aspect_ratio ** 0.138) / ((self.h_tail_t_max * self.meters_to_feet) ** 0.223))
                       * (0.0554 / self.kg_to_pounds))

        mass_v_tail = ((self.mass_takeoff_1 * self.kg_to_pounds) ** 0.567
                       * (self.v_tail_area * self.meters_to_feet ** 2) ** 0.149
                       * ((self.v_tail_aspect_ratio ** 0.482) / ((self.v_tail_t_max * self.meters_to_feet) ** 0.747))
                       * (1 / (np.cos(self.v_tail_sweep_quarter) ** 0.882))
                       * (0.1077 / self.kg_to_pounds))

        mass_fuselage = ((self.mass_takeoff_1 * self.kg_to_pounds) ** 0.692
                         * 2 ** 0.374
                         * (self.fuselage_length * self.meters_to_feet) ** 0.590
                         * (0.04682 / self.kg_to_pounds))

        mass_gear = ((0.013 * (self.mass_takeoff_1 * self.kg_to_pounds)
                     + 0.146 * (self.mass_takeoff_1 * self.kg_to_pounds) ** 0.417 * self.gear_load_factor ** 0.950 * (
                                 self.gear_length * self.meters_to_feet) ** 0.183
                     + 6.2 + 0.0013 * (self.mass_takeoff_1 * self.kg_to_pounds) + 0.000143 * (
                                 self.mass_takeoff_1 * self.kg_to_pounds) ** 0.749 * self.gear_load_factor * (
                                 self.gear_length * self.meters_to_feet) ** 0.788
                     + 0.014 * (self.mass_takeoff_1 * self.kg_to_pounds))) * (1 / self.kg_to_pounds)

        mass_control = 0.0168 * self.mass_takeoff_1
        mass_electric = 0.0268 * self.mass_takeoff_1
        mass_misc = (0.911 * (self.mass_takeoff_1 * self.kg_to_pounds) ** 0.489) / self.kg_to_pounds

        return np.round(np.array([mass_wing, mass_h_tail, mass_v_tail, mass_fuselage, mass_gear, mass_control,
                                  mass_electric, mass_misc, self.mass_battery, self.mass_motor, self.mass_occupants]))

    def raymer(self):
        mass_wing = (0.036 * (self.meters_to_feet ** 2 * self.wing_area) ** 0.758
                     * (self.wing_aspect_ratio / ((np.cos(self.wing_sweep_quarter)) ** 2)) ** 0.6
                     * (self.pressure * self.pascal_to_psf) ** 0.006
                     * self.wing_taper ** 0.04
                     * ((100 * self.wing_t_to_c) / (np.cos(self.wing_sweep_quarter))) ** -0.3
                     * (self.load_factor_ultimate * self.mass_takeoff_1 * self.kg_to_pounds) ** 0.49
                     * 1 / self.kg_to_pounds)

        mass_h_tail = (0.016 * (self.load_factor_ultimate * self.mass_takeoff_1 * self.kg_to_pounds) ** 0.414
                       * (self.pressure * self.pascal_to_psf) ** 0.168
                       * (self.h_tail_area * self.meters_to_feet ** 2) ** 0.896
                       * ((100 * self.wing_t_to_c) / (np.cos(self.wing_sweep_quarter))) ** -0.12
                       * (self.h_tail_aspect_ratio / ((np.cos(self.h_tail_sweep_quarter)) ** 2)) ** 0.043
                       * self.h_tail_taper ** -0.02
                       * 1 / self.kg_to_pounds)

        mass_v_tail = (0.073 * (self.load_factor_ultimate * self.mass_takeoff_1 * self.kg_to_pounds) ** 0.376
                       * (self.pressure * self.pascal_to_psf) ** 0.112
                       * (self.v_tail_area * self.meters_to_feet ** 2) ** 0.873
                       * ((100 * self.wing_t_to_c) / (np.cos(self.v_tail_sweep_quarter))) ** -0.49
                       * (self.wing_aspect_ratio / (np.cos(self.v_tail_sweep_quarter)) ** 2) ** 0.357
                       * self.v_tail_taper ** 0.039
                       * 1 / self.kg_to_pounds)

        mass_fuselage = (0.052 * (self.fuselage_area * self.meters_to_feet ** 2) ** 1.086
                         * (self.load_factor_ultimate * self.mass_takeoff_1 * self.kg_to_pounds) ** 0.177
                         * (self.h_tail_arm * self.meters_to_feet) ** -0.051
                         * (self.fuselage_length / self.fuselage_width) ** -0.072
                         * (self.pressure * self.pascal_to_psf) ** 0.241
                         * (1 / self.kg_to_pounds))

        mass_gear_main = (0.095 * (self.gear_load_factor * self.mass_takeoff_1 * self.kg_to_pounds) ** 0.768
                          * (self.gear_length * self.meters_to_feet) ** 0.409
                          * 1 / self.kg_to_pounds)

        mass_gear_nose = (0.125 * (self.gear_load_factor * self.mass_takeoff_1 * self.kg_to_pounds) ** 0.566
                          * (self.gear_length * self.meters_to_feet) ** 0.845
                          * 1 / self.kg_to_pounds)

        mass_gear = mass_gear_main + mass_gear_nose

        mass_control = (0.053 * (self.fuselage_length * self.meters_to_feet) ** 1.536
                        * (self.wing_span * self.meters_to_feet) ** 0.371
                        * (self.load_factor_ultimate * self.mass_takeoff_1 * self.kg_to_pounds * 10 ** -4) ** 0.8
                        * 1 / self.kg_to_pounds)

        mass_electric = 0

        mass_misc = (((0.0582 * self.mass_takeoff_1 * self.kg_to_pounds) - 65)
                     * 1 / self.kg_to_pounds)

        return np.round(np.array([mass_wing, mass_h_tail, mass_v_tail, mass_fuselage, mass_gear, mass_control,
                                  mass_electric, mass_misc, self.mass_battery, self.mass_motor, self.mass_occupants]))

    def torenbeek(self):
        mass_wing = (self.load_factor_ultimate ** 0.55 * ((self.wing_span * self.wing_area * self.meters_to_feet ** 3) /
                                                          (self.wing_t_max * self.meters_to_feet * self.mass_takeoff_1 *
                                                           self.kg_to_pounds * np.cos(self.wing_sweep_half))) ** 0.30
                     * ((self.wing_span * self.meters_to_feet) / np.cos(self.wing_sweep_half)) ** 0.75
                     * (1 + np.sqrt((6.3 * np.cos(self.wing_sweep_half)) / (self.wing_span * self.meters_to_feet)))
                     * (0.00125 / self.kg_to_pounds) * (self.mass_takeoff_1 * self.kg_to_pounds))

        mass_h_tail = ((((self.h_tail_area + self.v_tail_area) * self.meters_to_feet ** 2) ** 2) ** 0.75
                       * (0.02 / self.kg_to_pounds) * self.load_factor_ultimate ** 0.75)

        mass_v_tail = mass_h_tail

        mass_fuselage = 0

        mass_gear_main = ((33 + 0.04 * (self.mass_takeoff_1 * self.kg_to_pounds) ** 0.75
                           + 0.021 * (self.mass_takeoff_1 * self.kg_to_pounds))
                          * (1 / self.kg_to_pounds))

        mass_gear_nose = ((12 + 0.06 * (self.mass_takeoff_1 * self.kg_to_pounds) ** 0.75)
                          * (1 / self.kg_to_pounds))

        mass_gear = mass_gear_main + mass_gear_nose
        mass_control = (0.33 / self.kg_to_pounds) * (self.mass_takeoff_1 * self.kg_to_pounds) ** (2 / 3)
        mass_electric = (0.0078 * (self.mass_takeoff_1 * self.kg_to_pounds) ** 1.2) / self.kg_to_pounds
        mass_misc = 0

        return np.round(np.array([mass_wing, mass_h_tail, mass_v_tail, mass_fuselage, mass_gear, mass_control,
                                  mass_electric, mass_misc, self.mass_battery, self.mass_motor, self.mass_occupants]))

    def usaf(self):
        mass_wing = (96.948 * (((self.load_factor_ultimate * self.mass_takeoff_1 * self.kg_to_pounds) / (10 ** 5)) ** 0.65
                               * (self.wing_aspect_ratio / (np.cos(self.wing_sweep_quarter))) ** 0.57
                               * ((self.wing_area * self.meters_to_feet ** 2) / 100) ** 0.61
                               * ((1 + self.wing_taper) / (2 * self.wing_t_to_c)) ** 0.36
                               * np.sqrt(1 + ((self.velocity * self.ms_to_knots) / 500))) * 0.993
                     * (1 / self.kg_to_pounds))

        mass_h_tail = (127 * (((self.load_factor_ultimate * self.mass_takeoff_1 * self.kg_to_pounds) / (10 ** 5)) ** 0.87
                           * ((self.h_tail_area * self.meters_to_feet ** 2) / 100) ** 1.2 * 0.289
                           * (self.h_tail_arm * self.meters_to_feet / 10) ** 0.483
                           * np.sqrt((self.h_tail_span * self.meters_to_feet) / (self.h_tail_t_max * self.meters_to_feet))) ** 0.458
                    * (1 / self.kg_to_pounds))

        mass_v_tail = (
                    98.5 * (((self.load_factor_ultimate * self.mass_takeoff_1 * self.kg_to_pounds) / (10 ** 5)) ** 0.87
                            * ((self.v_tail_area * self.meters_to_feet ** 2) / 100) ** 1.2 * 0.289
                            * np.sqrt((self.v_tail_span * self.meters_to_feet) / (self.v_tail_t_max * self.meters_to_feet))) ** 0.458
                    * 1 / self.kg_to_pounds)

        mass_fuselage = (200 * (((self.load_factor_ultimate * self.mass_takeoff_1 * self.kg_to_pounds) / (10 ** 5)) ** 0.286
                           * ((self.fuselage_length * self.meters_to_feet) / 10) ** 0.857
                           * (((self.fuselage_width + self.fuselage_height) * self.meters_to_feet) / 10)
                           * ((self.velocity * self.ms_to_knots) / 100) ** 0.338) ** 1.1
                    * 1 / self.kg_to_pounds)

        mass_gear = (0.054 * (self.gear_load_factor * self.mass_takeoff_1 * self.kg_to_pounds) ** 0.684
                     * (self.gear_length * self.meters_to_feet) ** 0.501
                     * 1 / self.kg_to_pounds)

        mass_control = (1.066 * (self.mass_takeoff_1 * self.kg_to_pounds) ** 0.626
                        * 1 / self.kg_to_pounds)

        mass_electric = 0

        mass_misc = (34.5 * 2 * (self.pressure * self.pascal_to_psf) ** 0.25
                     * (1 / self.kg_to_pounds))

        return np.round(np.array([mass_wing, mass_h_tail, mass_v_tail, mass_fuselage, mass_gear, mass_control,
                                  mass_electric, mass_misc, self.mass_battery, self.mass_motor, self.mass_occupants]))

    def combine(self):
        data_cessna = self.cessna()
        data_raymer = self.raymer()
        data_torenbeek = self.torenbeek()
        data_usaf = self.usaf()

        data = np.vstack((data_cessna, data_raymer, data_torenbeek, data_usaf))
        data_totals = np.round(np.array([np.sum(data[i]) for i in range(len(data))])).reshape(-1, 1)
        data = np.hstack((data, data_totals))

        # Correction factors.
        data[:, 0] = np.round(data[:, 0] * 0.50)
        data[:, 3] = np.round(data[:, 3] * 0.50)
        data[:, 4] = np.round(data[:, 4] * 0.60)
        data[:, 5] = np.round(data[:, 5] * 0.60)
        data[:, 7] = np.round(data[:, 7] * 0.50)

        average = []

        for i in range(len(data[0]) - 1):
            average_column = np.sum(data[:, i]) / len(np.where(data[:, i] > 0.1 * np.sum(data[:, i]))[0])
            average.append(round(average_column))

        average = np.array(average)
        relative = np.round((average / np.sum(average) * 100))

        data_average = np.vstack((average, relative))

        return np.sum(data_average[0]), data_average, data

    def positioning(self):
        # Fuselage group.
        moment_bat = (self.mass_bat_1 * self.arm_bat_1) + (self.mass_bat_2 * self.arm_bat_2)
        moment_motor = self.mass_motor * self.arm_motor

        moment_h_tail = self.mass_takeoff_2[0, 1] * self.arm_h_tail
        moment_v_tail = self.mass_takeoff_2[0, 2] * self.arm_v_tail
        moment_fuselage = self.mass_takeoff_2[0, 3] * self.arm_fuselage
        moment_control = self.mass_takeoff_2[0, 5] * self.arm_control
        moment_electric = self.mass_takeoff_2[0, 6] * self.arm_electric
        moment_misc = self.mass_takeoff_2[0, 7] * self.arm_misc

        moment_fuselage = (moment_bat + moment_motor + moment_h_tail + moment_v_tail
                           + moment_fuselage + moment_control + moment_electric + moment_misc)
        mass_fuselage = (np.sum(self.mass_takeoff_2[0])
                         - self.mass_takeoff_2[0, 0] - self.mass_takeoff_2[0, 4] - self.mass_takeoff_2[0, -1])
        arm_fuselage = moment_fuselage / mass_fuselage

        # Wing group.
        moment_wing = self.mass_takeoff_2[0, 0] * self.arm_wing
        moment_gear = self.mass_takeoff_2[0, 4] * self.arm_gear

        mass_wing = self.mass_takeoff_2[0, 0] + self.mass_takeoff_2[0, 4]
        arm_wing = (moment_wing + moment_gear) / mass_wing

        mass_eom = mass_wing + mass_fuselage
        arm_eom = 0.35 * self.wing_MAC

        self.wing_X_LE = arm_fuselage - arm_eom + (mass_wing / mass_fuselage) * ((arm_wing - arm_eom) * self.wing_MAC) \
                         + self.wing_X_LE_correction

        CG_0 = np.round(((((self.wing_X_LE + arm_wing) * mass_wing) + (arm_fuselage * mass_fuselage))
                         / (mass_wing + mass_fuselage)), 2)
        CG_0_mass = mass_eom
        CG_0_MAC = np.round(((CG_0 - self.wing_X_LE) / self.wing_MAC) * 100)

        CG_1_front = np.round((((mass_eom * CG_0) + self.arm_occupant_1 * self.mass_occupant_1)
                               / (mass_eom + self.mass_occupant_1)), 2)
        CG_1_front_mass = mass_eom + self.mass_occupant_1
        CG_1_front_MAC = np.round(((CG_1_front - self.wing_X_LE) / self.wing_MAC) * 100)

        CG_1_rear = np.round((((mass_eom * CG_0) + self.arm_occupant_2 * self.mass_occupant_2)
                              / (mass_eom + self.mass_occupant_2)), 2)
        CG_1_rear_mass = mass_eom + self.mass_occupant_2
        CG_1_rear_MAC = np.round(((CG_1_rear - self.wing_X_LE) / self.wing_MAC) * 100)

        CG_2 = np.round((((mass_eom * CG_0) + (self.arm_occupant_1 * self.mass_occupant_1) +
                          (self.arm_occupant_2 * self.mass_occupant_2)) / (mass_eom + self.mass_occupants)), 2)
        CG_2_mass = mass_eom + self.mass_occupants
        CG_2_MAC = np.round(((CG_2 - self.wing_X_LE) / self.wing_MAC) * 100)

        data_CG = np.array([[CG_0_mass, CG_0, CG_0_MAC],
                            [CG_1_front_mass, CG_1_front, CG_1_front_MAC],
                            [CG_1_rear_mass, CG_1_rear, CG_1_rear_MAC],
                            [CG_2_mass, CG_2, CG_2_MAC]])

        row_labels = np.array(['0 occupants', 'Front occupant', 'Rear occupant', '2 occupants'])
        column_labels = np.array(['Mass [kg]', 'CG location from nose [m]', 'CG location [% of MAC]'])

        dataframe = pd.DataFrame(data_CG, columns=column_labels, index=row_labels)

        self.aircraft_X_CG = data_CG

        print(dataframe)
        print()
        print(f'Wing leading edge position: {round(self.wing_X_LE, 2)} [m]')

    def scissors(self, ratio=0.15):
        wing_area_net = self.wing_area - (self.fuselage_width * self.wing_chord_root)

        beta = np.sqrt(1 - 0.232 ** 2)

        self.h_tail_C_L_alpha = ((2 * np.pi * self.h_tail_aspect_ratio)
                                 / (2 + np.sqrt(4 + (self.h_tail_aspect_ratio * beta / 0.95) ** 2
                                    * (1 + (np.tan(self.h_tail_sweep_half) ** 2) / beta ** 2))))

        self.aircraft_C_L_alpha = (self.wing_C_L_alpha * (1 + 2.15 * (self.fuselage_width / self.wing_span))
                                   * (wing_area_net / self.wing_area)
                                   + (np.pi / 2) * ((self.fuselage_width ** 2) / self.wing_area))

        X_ac_fuselage_1 = ((-1.8 * self.fuselage_width * self.fuselage_height * self.wing_X_LE)
                           / (self.aircraft_C_L_alpha * self.wing_area * self.wing_MAC))
        X_ac_fuselage_2 = (((0.273 * self.fuselage_width * self.wing_area * (self.wing_span - self.fuselage_width))
                            * np.tan(self.wing_sweep_quarter))
                           / ((1 + self.wing_taper) * self.wing_span * self.wing_MAC ** 2
                              * (self.wing_span + 2.15 * self.fuselage_width)))

        X_ac = self.wing_X_ac + X_ac_fuselage_1 + X_ac_fuselage_2
        X_cg = np.arange(0, 1, 0.01)

        C_m_ac = 0.0
        C_L_aircraft = 1.6
        C_L_h = -0.35 * self.h_tail_aspect_ratio ** (1 / 3)

        X_tail = self.arm_h_tail

        h_tail_arm = (X_tail - (X_ac * self.wing_MAC) - self.wing_X_LE)
        h_tail_speed_ratio = 0.85

        downwash = ((7 * self.wing_C_L_alpha)
                    / (4 * np.pi * self.h_tail_aspect_ratio * (self.wing_taper * h_tail_arm) ** 0.25))

        denominator = ((self.wing_C_L_alpha / self.h_tail_C_L_alpha)
                       * (1 - downwash) * (h_tail_arm / self.wing_MAC) * h_tail_speed_ratio)

        area_ratio_stability = (X_cg / denominator) - ((X_ac - 0.05) / denominator)
        area_ratio_control = (X_cg - X_ac + (C_m_ac / C_L_aircraft)) * ((C_L_aircraft * self.wing_MAC) /
                                                                        (C_L_h * h_tail_arm * h_tail_speed_ratio))

        min_cg = np.min(self.aircraft_X_CG[:, 2] / 100)
        max_cg = np.max(self.aircraft_X_CG[:, 2] / 100)

        plt.plot(X_cg, area_ratio_stability, label='Stability', c='blue')
        plt.plot(X_cg, area_ratio_stability - 0.05, label='Stab. without margin', c='green')
        plt.plot(X_cg, area_ratio_control, label='Controllability', c='red')
        plt.plot([min_cg, max_cg], [ratio, ratio], label='CG range', c='black')
        plt.plot()
        plt.xlabel('Xcg [% MAC]')
        plt.ylabel('Sh/S [-]')
        plt.xlim(0, 1)
        plt.ylim(0, 0.5)
        plt.legend()
        plt.show()

    def main(self, iterations=10):
        self.mass_takeoff_1, self.mass_takeoff_2, methods_data = self.combine()

        for i in range(iterations):
            self.mass_takeoff_1, self.mass_takeoff_2, methods_data = self.combine()

        row_labels = np.array(['Cessna', 'Raymer', 'Torenbeek', 'USAF'])
        column_labels = np.array(['Wing', 'H-tail', 'V-tail', 'Fuselage', 'Gear', 'Control',
                                  'Electric', 'Misc', 'Batteries', 'Motor', 'Occupants', 'Totals'])

        dataframe = pd.DataFrame(methods_data, columns=column_labels, index=row_labels)

        print()
        print(dataframe)
        print()

        row_labels = np.array(['Average [kg]', 'Average [%]'])
        column_labels = np.array(['Wing', 'H-tail', 'V-tail', 'Fuselage', 'Gear', 'Control',
                                  'Electric', 'Misc', 'Batteries', 'Motor', 'Occupants'])

        dataframe = pd.DataFrame(self.mass_takeoff_2, columns=column_labels, index=row_labels)

        print(dataframe)
        print()

        self.positioning()


if __name__ == '__main__':
    class_II = MassMethods()
    class_II.main()
    class_II.scissors()

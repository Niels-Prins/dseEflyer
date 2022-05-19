import numpy as np
import pandas as pd

from class_I import class_I

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


class MassMethods:

    def __init__(self):
        # Load factor attributes.
        self.safety_margin = 1.5
        self.load_factor = 4
        self.load_factor_ultimate = self.load_factor * self.load_factor

        # Class I attributes.
        # self.mass_empty_1, self.mass_takeoff_1, self.mass_battery, self.mass_motor, self.mass_occupants = class_I()

        self.mass_empty_1, self.mass_takeoff_1, self.mass_battery, self.mass_motor, self.mass_occupants = 1029, 1633, 308, 60, 180

        # Wing attributes.
        self.wing_aspect_ratio = 5.8
        self.wing_area = 13.1
        self.wing_span = 8.72
        self.wing_sweep_half = -7.45 * (np.pi / 180)
        self.wing_sweep_quarter = -3.74 * (np.pi / 180)
        self.wing_taper = 0.45
        self.wing_t_to_c = 0.15
        self.wing_t_max = 0.31
        self.wing_chord_root = 1.58
        self.wing_MAC = ((2 / 3) * self.wing_chord_root *
                         ((1 + self.wing_taper + self.wing_taper ** 2) / (1 + self.wing_taper)))
        self.wing_lift_alpha = 6.56

        # Horizontal tail attributes.
        self.h_tail_aspect_ratio = 5.56
        self.h_tail_area = 2.80
        self.h_tail_span = 4.00
        self.h_tail_sweep_quarter = 0
        self.h_tail_taper = 0.65
        self.h_tail_t_max = 0.12

        # Vertical tail attributes.
        self.v_tail_aspect_ratio = 2.78
        self.v_tail_area = 1.40
        self.v_tail_span = 2.00
        self.v_tail_sweep_quarter = 10 * (np.pi / 180)
        self.v_tail_taper = 0.65
        self.v_tail_t_max = 0.12

        # Fuselage attributes.
        self.fuselage_area = 21
        self.fuselage_length = 7.92
        self.fuselage_height = 1.3
        self.fuselage_width = 1.3
        self.fuselage_radius = max(self.fuselage_height, self.fuselage_width) / 2

        # Gear attributes.
        self.gear_length = 0.5
        self.gear_load_factor = 4

        # Empirical attributes.
        self.h_tail_arm = 0.40 * self.wing_span

        # Conversion factors.
        self.kg_to_pounds = 2.2046
        self.meters_to_feet = 3.2808
        self.pascal_to_psf = 0.0209

        # Flight condition attributes.
        self.velocity = 77
        self.density = 1.225
        self.pressure = (self.density / 2) * self.velocity ** 2

        # Mass fractions.
        self.mass_bat_1 = self.mass_battery / 3
        self.mass_bat_2 = self.mass_battery / 3
        self.mass_bat_3 = self.mass_battery / 3

        self.mass_occupant_1 = self.mass_occupants / 2
        self.mass_occupant_2 = self.mass_occupants / 2

        # Fuselage group moment arms as fractions of fuselage length.
        self.arm_h_tail = 1.00 * self.fuselage_length
        self.arm_v_tail = 1.00 * self.fuselage_length
        self.arm_fuselage = 0.35 * self.fuselage_length
        self.arm_bat_1 = 0.20 * self.fuselage_length
        self.arm_bat_2 = 0.40 * self.fuselage_length
        self.arm_bat_3 = 0.60 * self.fuselage_length
        self.arm_motor = 0.75 * self.fuselage_length
        self.arm_control = 0.60 * self.fuselage_length
        self.arm_electric = 0.30 * self.fuselage_length
        self.arm_misc = 0.75 * self.fuselage_length
        self.arm_occupant_1 = 0.30 * self.fuselage_length
        self.arm_occupant_2 = 0.50 * self.fuselage_length

        # Wing group moment arms as fraction of MAC.
        self.arm_wing = 0.40 * self.wing_MAC
        self.arm_gear = 0.50 * self.wing_MAC
        self.arm_EOM = 0.35 * self.wing_MAC

        # Stability attributes.

        self.wing_area_net = self.wing_area - (self.fuselage_width * self.wing_chord_root)

        self.X_ac_wing = 0.25
        self.X_ac_fuselage_1 = None
        self.X_ac_fuselage_2 = None
        self.X_ac = self.X_ac_wing + self.X_ac_fuselage_1 + self.X_ac_fuselage_2
        self.X_ac_h_tail = None

        self.h_tail_arm = self.X_ac_h_tail - self.X_ac
        self.h_tail_speed_ratio = 0.85

        # To be calculated attributes.
        self.mass_empty_2 = None
        self.mass_takeoff_2 = None

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
                       * (self.v_tail_area * self.meters_to_feet ** 2) ** 0.125
                       * ((self.v_tail_aspect_ratio ** 0.482) / ((self.v_tail_t_max * self.meters_to_feet) ** 0.747))
                       * (1 / (np.cos(self.v_tail_sweep_quarter) ** 0.882))
                       * (0.1077 / self.kg_to_pounds))

        mass_fuselage = ((self.mass_takeoff_1 * self.kg_to_pounds) ** 0.692
                         * (self.fuselage_radius * self.meters_to_feet) ** 0.374
                         * (self.fuselage_length * self.meters_to_feet) ** 0.590
                         * (0.04682 / self.kg_to_pounds))

        mass_gear = (0.013 * (self.mass_takeoff_1 * self.kg_to_pounds)
                     + 0.362 * (self.mass_takeoff_1 * self.kg_to_pounds) ** 0.417 * self.gear_load_factor ** 0.950 * (
                                 self.gear_length * self.meters_to_feet) ** 0.183
                     + 6.2 + 0.0013 * (self.mass_takeoff_1 * self.kg_to_pounds) + 0.007157 * (
                                 self.mass_takeoff_1 * self.kg_to_pounds) ** 0.749 * self.gear_load_factor * (
                                 self.gear_length * self.meters_to_feet) ** 0.788
                     + 0.014 * (self.mass_takeoff_1 * self.kg_to_pounds)) / self.kg_to_pounds

        mass_control = 0.0168 * self.mass_takeoff_1
        mass_electric = 0.0268 * self.mass_takeoff_1
        mass_misc = (0.0911 * (self.mass_takeoff_1 * self.kg_to_pounds) ** 0.489) / self.kg_to_pounds

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
                     * (self.wing_span * self.meters_to_feet) / np.cos(self.wing_sweep_half) ** 0.75
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
        mass_control = (0.33 / self.kg_to_pounds) * (self.mass_takeoff_1 * self.kg_to_pounds) ** 0.667
        mass_electric = (0.0078 * (self.mass_takeoff_1 * self.kg_to_pounds) ** 1.2) / self.kg_to_pounds
        mass_misc = 0

        return np.round(np.array([mass_wing, mass_h_tail, mass_v_tail, mass_fuselage, mass_gear, mass_control,
                                  mass_electric, mass_misc, self.mass_battery, self.mass_motor, self.mass_occupants]))

    def usaf(self):
        mass_wing = (96.948 * (((self.load_factor_ultimate * self.mass_takeoff_1 * self.kg_to_pounds) / 10 ** 5) ** 0.65
                               * (self.wing_aspect_ratio / (np.cos(self.wing_sweep_quarter)) ** 2) ** 0.57
                               * ((self.wing_area * self.meters_to_feet ** 2) / 100) ** 0.61
                               * ((1 + self.wing_taper) / (2 * self.wing_t_to_c)) ** 0.36
                               * np.sqrt(1 + (self.velocity * self.meters_to_feet) / 500)) * 0.993
                     * (1 / self.kg_to_pounds))

        mass_h_tail = (
                    127 * (((self.load_factor_ultimate * self.mass_takeoff_1 * self.kg_to_pounds) / (10 ** 5)) ** 0.87
                           * ((self.h_tail_area * self.meters_to_feet ** 2) / 100) ** 1.2 * 0.289
                           * ((self.h_tail_arm * self.meters_to_feet / 10)) ** 0.483
                           * np.sqrt(
                        (self.h_tail_span * self.meters_to_feet) / (self.h_tail_t_max * self.meters_to_feet))) ** 0.458
                    * (1 / self.kg_to_pounds))

        mass_v_tail = (
                    98.5 * (((self.load_factor_ultimate * self.mass_takeoff_1 * self.kg_to_pounds) / (10 ** 5)) ** 0.87
                            * ((self.v_tail_area * self.meters_to_feet ** 2) / 100) ** 1.2 * 0.289
                            * np.sqrt(
                        (self.v_tail_span * self.meters_to_feet) / (self.v_tail_t_max * self.meters_to_feet))) ** 0.458
                    * 1 / self.kg_to_pounds)

        mass_fuselage = (
                    200 * (((self.load_factor_ultimate * self.mass_takeoff_1 * self.kg_to_pounds) / 10 ** 5) ** 0.286
                           * ((self.fuselage_length * self.meters_to_feet) / 10) ** 0.857
                           * (((self.fuselage_width + self.fuselage_height) * self.meters_to_feet) / 10)
                           * ((self.velocity * self.meters_to_feet) / 100) ** 0.338) ** 1.1
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
        average = []

        for i in range(len(data[0])):
            average_column = np.sum(data[:, i]) / len(np.where(data[:, i] > 0.1 * np.sum(data[:, i]))[0])
            average.append(round(average_column))

        average = np.array(average)
        relative = np.round((average / np.sum(average) * 100))

        data_average = np.vstack((average, relative))

        return np.sum(data_average), data_average, data

    def positioning(self):
        # Fuselage group.
        moment_bat = ((self.mass_bat_1 * self.arm_bat_1)
                      + (self.mass_bat_2 * self.arm_bat_2)
                      + (self.mass_bat_3 * self.arm_bat_3))
        moment_occupants = (self.mass_occupant_1 * self.arm_occupant_1) + (self.mass_occupant_2 * self.arm_occupant_2)
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

        arm_wing_LE = arm_fuselage - arm_eom + (mass_wing / mass_fuselage) * ((arm_wing - arm_eom) * self.wing_MAC)

        CG_0 = np.round(((((arm_wing_LE + arm_wing) * mass_wing) + (arm_fuselage * mass_fuselage))
                         / (mass_wing + mass_fuselage)), 2)
        CG_0_mass = mass_eom
        CG_0_MAC = np.round(((CG_0 - arm_wing_LE) / self.wing_MAC) * 100, 2)

        CG_1_front = np.round((((mass_eom * CG_0) + self.arm_occupant_1 * self.mass_occupant_1)
                               / (mass_eom + self.mass_occupant_1)), 2)
        CG_1_front_mass = mass_eom + self.mass_occupant_1
        CG_1_front_MAC = np.round(((CG_1_front - arm_wing_LE) / self.wing_MAC) * 100, 2)

        CG_1_rear = np.round((((mass_eom * CG_0) + self.arm_occupant_2 * self.mass_occupant_2)
                              / (mass_eom + self.mass_occupant_2)), 2)
        CG_1_rear_mass = mass_eom + self.mass_occupant_2
        CG_1_rear_MAC = np.round(((CG_1_rear - arm_wing_LE) / self.wing_MAC) * 100, 2)

        CG_2 = np.round((((mass_eom * CG_0) + (self.arm_occupant_1 * self.mass_occupant_1) +
                          (self.arm_occupant_2 * self.mass_occupant_2)) / (mass_eom + self.mass_occupants)), 2)
        CG_2_mass = mass_eom + self.mass_occupants
        CG_2_MAC = np.round(((CG_2 - arm_wing_LE) / self.wing_MAC) * 100, 2)

        data_CG = np.array([[CG_0_mass, CG_0, CG_0_MAC],
                            [CG_1_front_mass, CG_1_front, CG_1_front_MAC],
                            [CG_1_rear_mass, CG_1_rear, CG_1_rear_MAC],
                            [CG_2_mass, CG_2, CG_2_MAC]])

        row_labels = np.array(['0 occupants', 'Front occupant', 'Rear occupant', '2 occupants'])
        column_labels = np.array(['Mass [kg]', 'CG location from nose [m]', 'CG location [% of MAC]'])

        dataframe = pd.DataFrame(data_CG, columns=column_labels, index=row_labels)

        print(dataframe)

    def scissors(self):
        h_tail_C_L_alpha = (2 * np.pi * self.h_tail_aspect_ratio) / (2 + np.sqrt(4 + (self.h_tail_aspect_ratio / 0.95) ** 2))

    def main(self, iteration=0.02):

        while self.mass_takeoff_1 / self.combine()[0] > (1 + iteration):
            self.mass_takeoff_1, self.mass_takeoff_2, methods_data = self.combine()

        row_labels = np.array(['Cessna', 'Raymer', 'Torenbeek', 'USAF'])
        column_labels = np.array(['Wing', 'H-tail', 'V-tail', 'Fuselage', 'Gear', 'Control',
                                  'Electric', 'Misc', 'Batteries', 'Motor', 'Occupants'])

        dataframe = pd.DataFrame(methods_data, columns=column_labels, index=row_labels)

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

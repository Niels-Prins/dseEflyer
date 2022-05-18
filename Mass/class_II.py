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
        # self.mass_empty, self.mass_takeoff, self.mass_battery, self.mass_motor = class_I()

        self.mass_empty, self.mass_takeoff, self.mass_battery, self.mass_motor = 1029, 1633, 308, 60

        # Wing attributes.
        self.wing_aspect_ratio = 10.12
        self.wing_area = 13.46
        self.wing_span = 11.67
        self.wing_sweep_half = 0
        self.wing_sweep_quarter = 0
        self.wing_taper = 0.5
        self.wing_t_to_c = 0.15
        self.wing_t_max = 0.24

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

    def cessna(self):
        mass_wing = ((self.load_factor_ultimate * self.mass_takeoff * self.kg_to_pounds) ** 0.397
                     * (self.wing_area * self.meters_to_feet ** 2) ** 0.360
                     * self.wing_aspect_ratio ** 1.712
                     * (0.04674 / self.kg_to_pounds))

        mass_h_tail = ((self.mass_takeoff * self.kg_to_pounds) ** 0.887
                       * (self.h_tail_area * self.meters_to_feet ** 2) ** 0.101
                       * ((self.h_tail_aspect_ratio ** 0.138) / ((self.h_tail_t_max * self.meters_to_feet) ** 0.223))
                       * (0.0183 / self.kg_to_pounds))

        mass_v_tail = ((self.mass_takeoff * self.kg_to_pounds) ** 0.567
                       * (self.v_tail_area * self.meters_to_feet ** 2) ** 0.125
                       * ((self.v_tail_aspect_ratio ** 0.482) / ((self.v_tail_t_max * self.meters_to_feet) ** 0.747))
                       * (1 / (np.cos(self.v_tail_sweep_quarter) ** 0.882))
                       * (0.0026 / self.kg_to_pounds))

        mass_fuselage = ((self.mass_takeoff * self.kg_to_pounds) ** 0.692
                         * (self.fuselage_radius * self.meters_to_feet) ** 0.374
                         * (self.fuselage_length * self.meters_to_feet) ** 0.590
                         * (0.0468 / self.kg_to_pounds))

        mass_gear = ((6.2 + 0.0143 * self.mass_takeoff * self.kg_to_pounds
                      + 0.3620 * (self.mass_takeoff * self.kg_to_pounds) ** 0.417
                      * self.gear_load_factor ** 0.950 * (self.gear_length * self.meters_to_feet) ** 0.183
                      + 0.0072 * (self.mass_takeoff * self.kg_to_pounds) ** 0.749
                      * self.load_factor_ultimate * (self.gear_length * self.meters_to_feet) ** 0.788)
                     * (1 / self.kg_to_pounds))

        mass_control = 0.0168 * self.mass_takeoff
        mass_electric = 0.0268 * self.mass_takeoff
        mass_misc = (0.0911 * (self.mass_takeoff * self.kg_to_pounds) ** 0.489) / self.kg_to_pounds

        return np.round(np.array([mass_wing, mass_h_tail, mass_v_tail, mass_fuselage, mass_gear,
                                  mass_control, mass_electric, mass_misc, self.mass_battery, self.mass_motor]))

    def raymer(self):
        mass_wing = (0.036 * (self.meters_to_feet ** 2 * self.wing_area) ** 0.758
                     * (self.wing_aspect_ratio / ((np.cos(self.wing_sweep_quarter)) ** 2)) ** 0.6
                     * (self.pressure * self.pascal_to_psf) ** 0.006
                     * self.wing_taper ** 0.04
                     * ((100 * self.wing_t_to_c) / (np.cos(self.wing_sweep_quarter))) ** -0.3
                     * (self.load_factor_ultimate * self.mass_takeoff * self.kg_to_pounds) ** 0.49
                     * 1 / self.kg_to_pounds)

        mass_h_tail = (0.016 * (self.load_factor_ultimate * self.mass_takeoff * self.kg_to_pounds) ** 0.414
                       * (self.pressure * self.pascal_to_psf) ** 0.168
                       * (self.h_tail_area * self.meters_to_feet ** 2) ** 0.896
                       * ((100 * self.wing_t_to_c) / (np.cos(self.wing_sweep_quarter))) ** -0.12
                       * (self.h_tail_aspect_ratio / ((np.cos(self.h_tail_sweep_quarter)) ** 2)) ** 0.043
                       * self.h_tail_taper ** -0.02
                       * 1 / self.kg_to_pounds)

        mass_v_tail = (0.073 * (self.load_factor_ultimate * self.mass_takeoff * self.kg_to_pounds) ** 0.376
                       * (self.pressure * self.pascal_to_psf) ** 0.112
                       * (self.v_tail_area * self.meters_to_feet ** 2) ** 0.873
                       * ((100 * self.wing_t_to_c) / (np.cos(self.v_tail_sweep_quarter))) ** -0.49
                       * (self.wing_aspect_ratio / (np.cos(self.v_tail_sweep_quarter)) ** 2) ** 0.357
                       * self.v_tail_taper ** 0.039
                       * 1 / self.kg_to_pounds)

        mass_fuselage = (0.052 * (self.fuselage_area * self.meters_to_feet ** 2) ** 1.086
                         * (self.load_factor_ultimate * self.mass_takeoff * self.kg_to_pounds) ** 0.177
                         * (self.h_tail_arm * self.meters_to_feet) ** -0.051
                         * (self.fuselage_length / self.fuselage_width) ** -0.072
                         * (self.pressure * self.pascal_to_psf) ** 0.241
                         * (1 / self.kg_to_pounds))

        mass_gear_main = (0.095 * (self.gear_load_factor * self.mass_takeoff * self.kg_to_pounds) ** 0.768
                          * (self.gear_length * self.meters_to_feet) ** 0.409
                          * 1 / self.kg_to_pounds)

        mass_gear_nose = (0.125 * (self.gear_load_factor * self.mass_takeoff * self.kg_to_pounds) ** 0.566
                          * (self.gear_length * self.meters_to_feet) ** 0.845
                          * 1 / self.kg_to_pounds)

        mass_gear = mass_gear_main + mass_gear_nose

        mass_control = (0.053 * (self.fuselage_length * self.meters_to_feet) ** 1.536
                        * (self.wing_span * self.meters_to_feet) ** 0.371
                        * (self.load_factor_ultimate * self.mass_takeoff * self.kg_to_pounds * 10 ** -4) ** 0.8
                        * 1 / self.kg_to_pounds)

        mass_electric = 0

        mass_misc = (((0.0582 * self.mass_takeoff * self.kg_to_pounds) - 65)
                     * 1 / self.kg_to_pounds)

        return np.round(np.array([mass_wing, mass_h_tail, mass_v_tail, mass_fuselage, mass_gear,
                                  mass_control, mass_electric, mass_misc, self.mass_battery, self.mass_motor]))

    def torenbeek(self):
        mass_wing = (self.load_factor_ultimate ** 0.55 * ((self.wing_span * self.wing_area * self.meters_to_feet ** 3) /
                                                          (self.wing_t_max * self.meters_to_feet * self.mass_takeoff *
                                                           self.kg_to_pounds * np.cos(self.wing_sweep_half))) ** 0.30
                     * (self.wing_span * self.meters_to_feet) / np.cos(self.wing_sweep_half) ** 0.75
                     * (1 + np.sqrt((6.3 * np.cos(self.wing_sweep_half)) / (self.wing_span * self.meters_to_feet)))
                     * (0.00125 / self.kg_to_pounds) * (self.mass_takeoff * self.kg_to_pounds))

        mass_h_tail = ((((self.h_tail_area + self.v_tail_area) * self.meters_to_feet ** 2) ** 2) ** 0.75
                       * (0.02 / self.kg_to_pounds) * self.load_factor_ultimate ** 0.75)

        mass_v_tail = mass_h_tail

        mass_fuselage = 0

        mass_gear_main = ((33 + 0.04 * (self.mass_takeoff * self.kg_to_pounds) ** 0.75
                           + 0.021 * (self.mass_takeoff * self.kg_to_pounds))
                          * (1 / self.kg_to_pounds))

        mass_gear_nose = ((12 + 0.06 * (self.mass_takeoff * self.kg_to_pounds) ** 0.75)
                          * (1 / self.kg_to_pounds))

        mass_gear = mass_gear_main + mass_gear_nose
        mass_control = (0.23 / self.kg_to_pounds) * (self.mass_takeoff * self.kg_to_pounds) ** 0.667
        mass_electric = (0.0078 * (self.mass_takeoff * self.kg_to_pounds) ** 1.2) / self.kg_to_pounds
        mass_misc = 0

        return np.round(np.array([mass_wing, mass_h_tail, mass_v_tail, mass_fuselage, mass_gear,
                                  mass_control, mass_electric, mass_misc, self.mass_battery, self.mass_motor]))

    def usaf(self):
        mass_wing = (96.948 * (((self.load_factor_ultimate * self.mass_takeoff * self.kg_to_pounds) / 10 ** 5) ** 0.65
                     * (self.wing_aspect_ratio / (np.cos(self.wing_sweep_quarter)) ** 2) ** 0.57
                     * ((self.wing_area * self.meters_to_feet ** 2) / 100) ** 0.61
                     * ((1 + self.wing_taper) / (2 * self.wing_t_to_c)) ** 0.36
                     * np.sqrt(1 + (self.velocity * self.meters_to_feet) / 500)) * 0.993
                     * (1 / self.kg_to_pounds))

        mass_h_tail = (71.927 * (((self.load_factor_ultimate * self.mass_takeoff * self.kg_to_pounds) / (10 ** 5)) ** 0.87
                       * ((self.h_tail_area * self.meters_to_feet ** 2) / 100) ** 1.2
                       * ((self.h_tail_arm * self.meters_to_feet / 10)) ** 0.483
                       * np.sqrt((self.h_tail_span * self.meters_to_feet) / (self.h_tail_t_max * self.meters_to_feet)))
                       * (0.458 / self.kg_to_pounds))

        mass_v_tail = (55.786 * (((self.load_factor_ultimate * self.mass_takeoff * self.kg_to_pounds) / (10 ** 5)) ** 0.87
                       * ((self.v_tail_area * self.meters_to_feet ** 2) / 100) ** 1.2
                       * np.sqrt((self.v_tail_span * self.meters_to_feet) / (self.v_tail_t_max * self.meters_to_feet))) ** 0.458
                       * 1 / self.kg_to_pounds)

        mass_fuselage = (200 * (((self.load_factor_ultimate * self.mass_takeoff * self.kg_to_pounds) / 10 ** 5) ** 0.286
                         * ((self.fuselage_length * self.meters_to_feet) / 10) ** 0.857
                         * (((self.fuselage_width + self.fuselage_height) * self.meters_to_feet) / 10)
                         * ((self.velocity * self.meters_to_feet) / 100) ** 0.338) ** 1.1
                         * 1 / self.kg_to_pounds)

        mass_gear = (0.054 * (self.gear_load_factor * self.mass_takeoff * self.kg_to_pounds) ** 0.684
                     * (self.gear_length * self.meters_to_feet) ** 0.501
                     * 1 / self.kg_to_pounds)

        mass_control = (1.066 * (self.mass_takeoff * self.kg_to_pounds) ** 0.626
                        * 1 / self.kg_to_pounds)

        mass_electric = 0

        mass_misc = (34.5 * 2 * (self.pressure * self.pascal_to_psf) ** 0.25
                     * (1 / self.kg_to_pounds))

        return np.round(np.array([mass_wing, mass_h_tail, mass_v_tail, mass_fuselage, mass_gear,
                                  mass_control, mass_electric, mass_misc, self.mass_battery, self.mass_motor]))

    def combine(self):
        data_cessna = self.cessna()
        data_raymer = self.raymer()
        data_torenbeek = self.torenbeek()
        data_usaf = self.usaf()

        data = np.vstack((data_cessna, data_raymer, data_torenbeek, data_usaf))
        average = []

        for i in range(len(data[0])):
            average_column = np.sum(data[:, i]) / len(np.where(data[:, i] > 0)[0])
            average.append(round(average_column))

        average = np.array(average)
        average_relative = np.round((average / np.sum(average[:-1]) * 100))

        data_average = np.vstack((average, average_relative))

        row_labels = np.array(['Cessna', 'Raymer', 'Torenbeek', 'USAF'])
        column_labels = np.array(['Wing', 'H-tail', 'V-tail', 'Fuselage', 'Gear', 'Control',
                                  'Electric', 'Misc', 'Batteries', 'Motor'])

        dataframe = pd.DataFrame(data, columns=column_labels, index=row_labels)
        dataframe_average = pd.DataFrame(data_average, columns=column_labels, index=['Average [kg]', 'Average [%]'])

        print(f'{dataframe} \n\n {dataframe_average} \n\n Expected take-off mass: {np.sum(data_average)} [kg]')


if __name__ == '__main__':
    class_II = MassMethods()
    class_II.combine()

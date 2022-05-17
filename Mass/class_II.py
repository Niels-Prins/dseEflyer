import numpy as np

from class_I import class_I


class MassMethods:

    def __init__(self):
        # Load factor attributes.
        self.safety_margin = 1.5
        self.load_factor = 1
        self.load_factor_ultimate = self.load_factor * self.load_factor

        # Class I attributes.
        self.mass_empty, self.mass_takeoff = class_I()

        # Wing attributes.
        self.wing_aspect_ratio = 1
        self.wing_area = 1
        self.wing_span = 1
        self.wing_sweep_half = 1
        self.wing_sweep_quarter = 1
        self.wing_taper = 1
        self.wing_t_to_c = 1
        self.wing_t_max = 1

        # Horizontal tail attributes.
        self.h_tail_aspect_ratio = 1
        self.h_tail_area = 1
        self.h_tail_span = 1
        self.h_tail_sweep_half = 1
        self.h_tail_sweep_quarter = 1
        self.h_tail_taper = 1
        self.h_tail_t_to_c = 1
        self.h_tail_t_max = 1

        # Vertical tail attributes.
        self.v_tail_aspect_ratio = 1
        self.v_tail_area = 1
        self.v_tail_span = 1
        self.v_tail_sweep_half = 1
        self.v_tail_sweep_quarter = 1
        self.v_tail_taper = 1
        self.v_tail_t_to_c = 1
        self.v_tail_t_max = 1

        # Fuselage attributes.
        self.fuselage_area = 1
        self.fuselage_length = 1
        self.fuselage_height = 1
        self.fuselage_width = 1
        self.fuselage_radius = max(self.fuselage_height, self.fuselage_width) / 2

        # Gear attributes.
        self.gear_length = 0.5
        self.gear_load_factor = 4

        # Empirical attributes.
        self.h_tail_arm = 0.40 * self.wing_span

        # Conversion factors.
        self.kg_to_pounds = 2.2046
        self.meters_to_feet = 3.2808
        self.pascal_to_empirical = self.kg_to_pounds / self.meters_to_feet

        # Flight condition attributes.
        self.velocity = 100
        self.density = 200
        self.pressure = (self.density / 2) * self.velocity ** 2

    def cessna(self):
        mass_wing = ((self.load_factor_ultimate * self.mass_takeoff * self.kg_to_pounds) ** 0.397
                     * (self.wing_area * self.meters_to_feet ** 2) ** 0.360
                     * self.wing_aspect_ratio ** 1.712
                     * (0.04674 / self.kg_to_pounds))

        mass_h_tail = ((self.mass_takeoff * self.kg_to_pounds) ** 0.887
                       * (self.h_tail_area * self.meters_to_feet ** 2) ** 0.101
                       * ((self.h_tail_aspect_ratio ** 0.138) / (self.h_tail_t_max ** 0.223))
                       * (0.0183 / self.kg_to_pounds))

        mass_v_tail = ((self.mass_takeoff * self.kg_to_pounds) ** 0.567
                       * (self.v_tail_area * self.meters_to_feet ** 2) ** 0.125
                       * ((self.v_tail_aspect_ratio ** 0.482) / (self.v_tail_t_max ** 0.747))
                       * (1 / (np.cos(self.v_tail_sweep_quarter) ** 0.882))
                       * (0.0026 / self.kg_to_pounds))

        mass_fuselage = ((self.mass_takeoff * self.kg_to_pounds) ** 0.692
                         * (self.fuselage_radius * self.meters_to_feet) ** 0.374
                         * (self.fuselage_length * self.meters_to_feet) ** 0.590
                         * (0.0468 / self.kg_to_pounds))

        mass_gear = ((6.2 + 0.0143 * self.mass_takeoff * self.kg_to_pounds
                      + 0.3620 * (self.mass_takeoff * self.kg_to_pounds) ** 0.417
                      * self.gear_load_factor ** 0.950 * (self.gear_length * self.kg_to_pounds) ** 0.183
                      + 0.0072 * (self.mass_takeoff * self.kg_to_pounds) ** 0.749
                      * self.load_factor_ultimate * (self.gear_length * self.kg_to_pounds) ** 0.788)
                     * (1 / self.kg_to_pounds))

        mass_control = 0.0168 * self.mass_takeoff

        mass_electric = 0.0268 * self.mass_takeoff

        mass_misc = (0.0911 * (self.mass_takeoff * self.kg_to_pounds) ** 0.489) / self.kg_to_pounds

        mass = np.array([mass_wing, mass_h_tail, mass_v_tail, mass_fuselage, mass_gear, mass_control, mass_electric, mass_misc])

    def raymer(self):
        mass_wing = (0.036 * (self.meters_to_feet ** 2 * self.wing_area) ** 0.758
                     * (self.wing_aspect_ratio / ((np.cos(self.wing_sweep_quarter)) ** 2)) ** 0.6
                     * (0.5 * self.density * self.velocity ** 2 * self.pascal_to_empirical) ** 0.006
                     * self.wing_taper ** 0.04
                     * ((100 * self.wing_t_to_c) / (np.cos(self.wing_sweep_quarter))) ** -0.3
                     * (self.load_factor_ultimate * self.mass_takeoff * self.kg_to_pounds) ** 0.49
                     * 1 / self.kg_to_pounds)

        mass_h_tail = (0.016 * (self.load_factor_ultimate * self.mass_takeoff * self.kg_to_pounds) ** 0.414
                       * (0.5 * self.density * self.velocity ** 2 * self.pascal_to_empirical) ** 0.168
                       * self.h_tail_area ** 0.896
                       * ((100 * self.wing_t_to_c) / (np.cos(self.wing_sweep_quarter))) ** -0.12
                       * (self.h_tail_aspect_ratio / ((np.cos(self.h_tail_sweep_quarter)) ** 2)) ** 0.043
                       * self.h_tail_taper ** -0.12
                       * 1 / self.kg_to_pounds)

        mass_v_tail = (0.073 * (self.load_factor_ultimate * self.mass_takeoff * self.kg_to_pounds) ** 0.376
                       * (0.5 * self.density * self.velocity ** 2 * self.pascal_to_empirical) ** 0.112
                       * self.v_tail_area ** 0.873
                       * ((100 * self.wing_t_to_c) / (np.cos(self.v_tail_sweep_quarter))) ** -0.49
                       * (self.wing_aspect_ratio / (np.cos(self.v_tail_sweep_quarter)) ** 2) ** 0.357
                       * self.v_tail_taper ** 0.039
                       * 1 / self.kg_to_pounds)

        mass_fuselage = (0.052 * self.fuselage_area ** 1.086
                         * (self.load_factor_ultimate * self.mass_takeoff * self.kg_to_pounds) ** 0.177
                         * (self.h_tail_arm * self.meters_to_feet) ** -0.051
                         * ((self.fuselage_length / self.fuselage_width) * self.meters_to_feet) ** -0.072
                         * (0.5 * self.density * self.velocity ** 2 * self.pascal_to_empirical) ** 0.241
                         * 1 / self.kg_to_pounds)
        mass_gear_main = (0.095 * (self.load_factor_ultimate * self.mass_takeoff * self.kg_to_pounds) ** 0.749
                          * (self.gear_length * self.meters_to_feet) ** 0.409
                          * 1 / self.kg_to_pounds)

        mass_gear_nose = (0.125 * (self.gear_load_factor * self.mass_takeoff * self.mass_takeoff) ** 0.566
                          * (self.gear_length * self.meters_to_feet) ** 0.846
                          * 1 / self.kg_to_pounds)

        mass_control = (0.053 * (self.fuselage_length * self.meters_to_feet) ** 1.536
                        * (self.wing_span * self.meters_to_feet) ** 0.371
                        * (self.load_factor_ultimate * self.mass_takeoff * self.kg_to_pounds * 10 ** -4) ** 0.8
                        * 1 / self.kg_to_pounds)

        mass_electric = None

        mass_misc = (((0.0582 * self.mass_takeoff * self.kg_to_pounds) - 65)
                     * 1 / self.kg_to_pounds)

        print(mass_wing, mass_h_tail, mass_v_tail, mass_fuselage, mass_gear_main, mass_gear_nose, mass_control, mass_electric, mass_misc )

        return mass_wing, mass_h_tail, mass_v_tail, mass_fuselage, mass_gear_main, mass_gear_nose, mass_control, mass_electric, mass_misc

    def torenbeek(self):
        mass_wing = None
        mass_h_tail = None
        mass_v_tail = None
        mass_fuselage = None
        mass_gear_main = None
        mass_gear_nose = None
        mass_electric = None
        mass_misc = None

    def usaf(self):
        mass_wing = None
        mass_h_tail = None
        mass_v_tail = None
        mass_fuselage = None
        mass_gear_main = None
        mass_gear_nose = None
        mass_electric = None
        mass_misc = None


if __name__ == '__main__':
    class_II = MassMethods()
    class_II.cessna()
    class_II.raymer()

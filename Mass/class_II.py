import numpy as np


class MassMethods:

    def __init__(self):
        # Load factor attributes.
        self.safety_margin = 1.5
        self.load_factor = 1
        self.load_factor_ultimate = self.load_factor * self.load_factor

        # Class I attributes.
        self.mass_empty = 1
        self.mass_takeoff = 1

        # Wing attributes.
        self.wing_aspect_ratio = 1
        self.wing_area = 1
        self.wing_span = 1
        self.wing_sweep_half = 1
        self.wing_sweep_quarter = 1
        self.wing_taper = 1
        self.wing_t_to_c = 1
        self.wing_t_max = 1

        # Horizontal tail attributes
        self.horizontal_tail_aspect_ratio = 1
        self.horizontal_tail_area = 1
        self.horizontal_tail_span = 1
        self.horizontal_tail_sweep_half = 1
        self.horizontal_tail_sweep_quarter = 1
        self.horizontal_tail_taper = 1
        self.horizontal_tail_t_to_c = 1
        self.horizontal_tail_t_max = 1

        # Conversion factors.
        self.kg_to_pounds = 2.2046
        self.meters_to_feet = 3.2808
        self.pascal_to_empirical = (self.kg_to_pounds / self.meters_to_feet)

        # Flight condition attributes.
        self.velocity = 100
        self.density = 200
        self.pressure = (self.density / 2) * self.velocity ** 2

    def cessna(self):
        mass_wing = ((self.load_factor_ultimate * self.mass_takeoff * self.kg_to_pounds) ** 0.397
                     * (self.wing_area * self.meters_to_feet ** 2) ** 0.360
                     * self.wing_aspect_ratio ** 1.712
                     * (0.04674 / self.kg_to_pounds))

    def raymer(self):
        mass_wing = 0.036 * (self.meters_to_feet**2 * self.wing_area)**0.758 \
                    * (self.wing_aspect_ratio / ((np.cos(self.wing_sweep_quarter))**2))**0.6 \
                    * (0.5 * self.density * self.velocity**2 * self.pascal_to_empirical)**0.006 \
                    * self.wing_taper**0.04 \
                    * ((100 * self.wing_t_to_c) / (np.cos(self.wing_sweep_quarter)))**-0.3 \
                    * (self.load_factor_ultimate * self.mass_takeoff)**0.49

        mass_horizontal_tail = 0.016 * (self.load_factor_ultimate * self.mass_takeoff)
        return mass_wing, mass_horizontal_tail

    def torenbeek(self):
        pass

    def usaf(self):
        pass

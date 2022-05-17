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

        # Flight condition attributes.
        self.velocity = 100
        self.density = 200
        self.pressure = (self.density / 2) * self.velocity ** 2

        # Conversion factors.
        self.kg_to_pounds = 2.2046
        self.meters_to_feet = 3.2808

    def cessna(self):
        asdkandf = 5
        sdkfnsdf = 5

    def raymer(self):
        mass = 122
        return mass

    def torenbeek(self):
        pass

    def usaf(self):
        pass

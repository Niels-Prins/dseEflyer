import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)
pd.set_option('display.width', None)


class Control:

    def __init__(self):
        # Class I attributes.
        self.mass_empty_1, self.mass_takeoff_1, self.mass_battery, self.mass_motor, self.mass_occupants = self.class_I()

        # Load factor attributes.
        self.safety_margin = 1.5
        self.load_factor = 8
        self.load_factor_ultimate = self.load_factor * self.load_factor

        # Wing attributes.
        self.wing_aspect_ratio = 5.80
        self.wing_area = 14.50
        self.wing_span = 9.2
        self.wing_sweep_half = -7.45 * (np.pi / 180)
        self.wing_sweep_quarter = -3.74 * (np.pi / 180)
        self.wing_taper = 0.45
        self.wing_t_to_c = 0.15
        self.wing_t_max = 0.31
        self.wing_chord_root = 2.18
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
        self.fuselage_length = 7.68
        self.fuselage_height = 1.58
        self.fuselage_width = 0.80
        self.fuselage_radius = max(self.fuselage_height, self.fuselage_width) / 2

        # Fuselage attributes extended.
        self.fuselage_circular = False

        self.fuselage_nose_length = 2
        self.fuselage_main_length = 3.70
        self.fuselage_tail_length = 2
        self.fuselage_length = self.fuselage_nose_length + self.fuselage_main_length + self.fuselage_tail_length

        self.fuselage_nose_height_start = 0.50
        self.fuselage_nose_height_end = 1.60
        self.fuselage_nose_height_slope = ((self.fuselage_nose_height_end / self.fuselage_nose_height_start) /
                                           self.fuselage_nose_length)

        self.fuselage_nose_width_start = 0.50
        self.fuselage_nose_width_end = 0.80
        self.fuselage_nose_width_slope = ((self.fuselage_nose_width_end / self.fuselage_nose_width_start) /
                                          self.fuselage_nose_length)

        self.fuselage_main_height = self.fuselage_nose_height_end
        self.fuselage_main_width = self.fuselage_nose_width_end

        self.fuselage_tail_height_begin = self.fuselage_main_height
        self.fuselage_tail_height_end = 0.50
        self.fuselage_tail_width_start = self.fuselage_main_width
        self.fuselage_tail_width_end = 0.50

        self.fuselage_max_height = self.fuselage_main_height
        self.fuselage_max_width = self.fuselage_main_width

        # Gear attributes.
        self.gear_length = 0.5
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
        self.mass_bat_1 = (self.mass_battery * 4) / 10
        self.mass_bat_2 = (self.mass_battery * 6) / 10

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

        # Added attributes ...
        self.mach = 0.23
        self.fuselage_base = self.fuselage_width * self.wing_chord_root
        self.fuselage_frontal = 1.42
        self.fuselage_top_area = 12
        self.alpha = 5 * (np.pi / 180)

    @staticmethod
    def class_I():
        structure = 345
        batteries = 406
        motor = 66
        margin = 0.20
        OEM = (1 + margin) * (structure + batteries + motor)

        occupants = 180
        TOM = OEM + occupants

        # print(f'Class I OEM: {OEM} \nClass I TOM: {TOM}')

        return OEM, TOM, batteries, motor, occupants

    def fuselage(self):
        def fuselage_drag_raymer():
            C_D_0_fus_base = 1              # ????? first term of C_D_0_fus
            C_D_b_fus = (0.029 * ((np.sqrt(4/np.pi * self.fuselage_frontal)) / (self.fuselage_radius * 2)) ** 3
                         / np.sqrt(C_D_0_fus_base * self.wing_area / self.fuselage_area)
                         * (self.fuselage_area / self.wing_area))
            C_D_0_fus = (1.075 * 0.0030759 * (1 + 60 / (self.fuselage_length / self.fuselage_radius * 2) ** 3 + 0.00025
                                              * (self.fuselage_length / self.fuselage_radius * 2)) * self.fuselage_area
                         / self.wing_area + C_D_b_fus)

            n = (0.00001146 * (self.fuselage_length / (self.fuselage_radius * 2)) ** 3 - 0.0008522
                 * (self.fuselage_length / (self.fuselage_radius * 2)) ** 2 + 0.02499
                 * (self.fuselage_length / (self.fuselage_radius * 2)) + 0.507)

            C_d_c = (-2.77 * (self.mach * np.sin(self.alpha)) ** 4 + 3.88 * (self.mach * np.sin(self.alpha)) ** 3 - 0.527
                     * (self.mach * np.sin(self.alpha)) ** 2 + 0.0166 * (self.mach * np.sin(self.alpha)) + 1.2)

            C_D_L_fus = (2 * self.alpha ** 2 * self.fuselage_base / self.wing_area + n * C_d_c * abs(self.alpha ** 3)
                         * self.fuselage_top_area / self.wing_area)

            C_D_fus = C_D_L_fus + C_D_0_fus
            print('CD0_fus', C_D_0_fus)
            print('CDL_fus', C_D_L_fus)
            print('CD_fus', C_D_fus)
            return


        def fuselage_sideforce_sidewash():
            Zw = y_el = self.fuselage_height / 2 - self.wing_t_max / 2
            x_el = np.sqrt((1 - ((y_el ** 2) / (self.fuselage_height / 2) ** 2)) * (self.fuselage_width / 2) ** 2)
            dF2 = np.sqrt(x_el ** 2 + y_el ** 2)
            Ki = 1.6 * Zw / dF2

            x1 = self.fuselage_nose_length + self.fuselage_main_length
            x0 = self.fuselage_length * 0.378 + 0.527 * x1
            a = 0.5 * self.fuselage_tail_height_begin - (((0.5 * self.fuselage_tail_height_begin - 0.5
                                                           * self.fuselage_tail_height_end) / self.fuselage_tail_length)
                                                         * (x0 - x1))
            b = 0.5 * self.fuselage_tail_width_start - (((0.5 * self.fuselage_tail_width_start - 0.5
                                                          * self.fuselage_tail_width_end) / self.fuselage_tail_length)
                                                        * (x0 - x1))
            S0 = np.pi * a * b
            C_Y_beta_f = -2 * Ki * S0 / self.wing_area

            print(C_Y_beta_f)
            return

        fuselage_drag_raymer()
        fuselage_sideforce_sidewash()


if __name__ == '__main__':
    fuselage = Control()
    fuselage.fuselage()

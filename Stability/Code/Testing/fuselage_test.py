import unittest

from Stability.Code.aircraft import Fuselage


class TestFuselage(unittest.TestCase):

    def test_washes(self, margin=0.01):
        downwash_true = None
        downwash_calc = None

        sidewash_true = None
        sidewash_calc = None

        self.assertAlmostEqual(downwash_true, downwash_calc, delta=margin * downwash_true)
        self.assertAlmostEqual(sidewash_true, sidewash_calc, delta=margin * sidewash_true)

    def test_wing_corrections(self, margin=0.01):
        C_L_correction_true = None
        C_L_correction_calc = None

        C_D_0_correction_true = None
        C_D_0_correction_calc = None

        C_M_ac_correction_true = None
        C_M_ac_correction_calc = None

        X_ac_correction_true = None
        X_ac_correction_calc = None

        self.assertAlmostEqual(C_L_correction_true, C_L_correction_calc, delta=margin * C_L_correction_true)
        self.assertAlmostEqual(C_D_0_correction_true, C_D_0_correction_calc, delta=margin * C_D_0_correction_true)
        self.assertAlmostEqual(C_M_ac_correction_true, C_M_ac_correction_calc, delta=margin * C_M_ac_correction_true)
        self.assertAlmostEqual(X_ac_correction_true, X_ac_correction_calc, delta=margin * X_ac_correction_true)

    def test_velocity_corrections(self, margin=0.01):
        # Kelvin tests.
        kelvin_strength_true = None
        kelvin_strength_calc = None

        kelvin_location_true = None
        kelvin_location_calc = None

        kelvin_Y_dot_true = None
        kelvin_Y_dot_calc = None

        kelvin_Z_dot_true = None
        kelvin_Z_dot_calc = None

        self.assertAlmostEqual(kelvin_strength_true, kelvin_strength_calc, delta=margin * kelvin_strength_true)
        self.assertAlmostEqual(kelvin_location_true, kelvin_location_calc, delta=margin * kelvin_location_true)
        self.assertAlmostEqual(kelvin_Y_dot_true, kelvin_Y_dot_calc, delta=margin * kelvin_Y_dot_true)
        self.assertAlmostEqual(kelvin_Z_dot_true, kelvin_Z_dot_calc, delta=margin * kelvin_Z_dot_true)

        # Rankine tests.
        rankine_strength_true = None
        rankine_strength_calc = None

        rankine_location_true = None
        rankine_location_calc = None

        rankine_Y_dot_true = None
        rankine_Y_dot_calc = None

        rankine_Z_dot_true = None
        rankine_Z_dot_calc = None

        self.assertAlmostEqual(rankine_strength_true, rankine_strength_calc, delta=margin * rankine_strength_true)
        self.assertAlmostEqual(rankine_location_true, rankine_location_calc, delta=margin * rankine_location_true)
        self.assertAlmostEqual(rankine_Y_dot_true, rankine_Y_dot_calc, delta=margin * rankine_Y_dot_true)
        self.assertAlmostEqual(rankine_Z_dot_true, rankine_Z_dot_calc, delta=margin * rankine_Z_dot_true)

        # Doublet tests.
        doublet_strength_true = None
        doublet_strength_calc = None

        doublet_location_true = None
        doublet_location_calc = None

        doublet_Y_dot_true = None
        doublet_Y_dot_calc = None

        doublet_Z_dot_true = None
        doublet_Z_dot_calc = None

        self.assertAlmostEqual(doublet_strength_true, doublet_strength_calc, delta=margin * doublet_strength_true)
        self.assertAlmostEqual(doublet_location_true, doublet_location_calc, delta=margin * doublet_location_true)
        self.assertAlmostEqual(doublet_Y_dot_true, doublet_Y_dot_calc, delta=margin * doublet_Y_dot_true)
        self.assertAlmostEqual(doublet_Z_dot_true, doublet_Z_dot_calc, delta=margin * doublet_Z_dot_true)

    def test_drag_corrections(self, margin=0.01):
        drag_correction_true = None
        drag_correction_calc = None

        self.assertAlmostEqual(drag_correction_true, drag_correction_calc, delta=margin * drag_correction_true)


if __name__ == '__main__':
    unittest.main()

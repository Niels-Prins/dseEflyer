"""
Average values for wing and empennage planform according to Roskam
(Single engine propeller aircraft)
"""

# Wing attributes #
wing_area = 13.1                        # m^2
wing_span = 8.72                        # m
wing_chord = 1.5                        # m

# Aileron
Sa_S = 0.079545455                      # -         # Increase?
aileron_area = wing_area * Sa_S         # m^2
aileron_position_in = 0.570909091       # %span
aileron_position_out = 0.941818182      # %span
aileron_position_root = 0.239090909     # %chord
aileron_position_tip = 0.260909091      # %chord


# Horizontal tail attributes #
h_tail_area = 3.18319598                # m^2       # Increase?
h_tail_arm = 4.644043636                # m         # determine after location
h_tail_volume = 0.667272727             # -         # increase (result of larger volume)
h_tail_dihedral = 0                     # deg
h_tail_incidence = 0                    # deg
h_tail_aspect_ratio = 5.15              # -         # ?
h_tail_sweep_quarter = 5                # deg       # ?
h_tail_taper_ratio = 0.7025             # -

# Elevator
Se_Sh = 0.627272727                     # -
elevator_area = h_tail_area * Se_Sh     # m^2
elevator_position_root = 0.42           # %chord
elevator_position_tip = 0.431428571     # %chord


# Vertical tail attributes #
v_tail_area = 1.490671505               # m^2       # Increase?
v_tail_arm = 4.397432727                # m         # determine after location
v_tail_volume = 0.043636364             # -         # increase (result of larger volume)
v_tail_dihedral = 90                    # deg
v_tail_incidence = 0                    # deg
v_tail_aspect_ratio = 1.65              # -
v_tail_sweep_quarter = 27               # deg
v_tail_taper_ratio = 0.45               # -

# Rudder
Sr_Sv = 0.360909091                     # -
rudder_area = v_tail_area * Sr_Sv       # m^2
rudder_position_root = 0.330909091      # %chord
rudder_position_tip = 0.484545455       # %chord

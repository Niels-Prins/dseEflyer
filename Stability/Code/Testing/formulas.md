### Formulas

Derivations (when required) have been done on paper and proofread.

#### 1. area

area = (span / 2) * chord_root * (1 + taper)

#### 2. mean aerodynamic chord

mac = (2 * chord_root / 3) x (( 1 + taper + taper ** 2 ) / (1 + taper))
****
#### 3.1. aerodynamic center X-coordinate

initial_position = X_le - (chord_root / 2)

shift_due_ac = (1 - (2 * %_ac)) * (mac / 2)
shift_due_sweep = - Y_ac * sweep

X_ac = initial_position + shift_due_ac + shift_due_sweep

#### 3.2. aerodynamic center Y-coordinate

symmetric:
Y_ac = Y_le - (span * (1 + (2 * taper))) / (6 * (1 + taper))

asymmetric:
Y_ac = Y_le - (span * (1 + (2 * taper))) / (3 * (1 + taper))

#### 3.3. aerodynamic center Z-coordinate

initial_position = Z_le

shift_due_dihedral = - Y_ac * dihedral

Z_ac = initial_position + shift_due_dihedral

#### 4.1 transformation matrix due to surface rotation (either 1 or -1)

transformation = [[1, 0, 0], [0, 0, rotation], [0, - rotation, 0]]

#### 4.2 transformation matrix due to roll angle

transformation = [[1, 0, 0], [0, cos(roll), sin(roll)], [0, - sin(roll), cos(roll)]]

#### 4.3 transformation matrix due to pitch angle

transformation = [[cos(pitch), 0, - sin(pitch)], [0, 1, 0], [sin(pitch), 0, cos(pitch)]]

#### 4.4 transformation matrix due to yaw angle

transformation = [[cos(yaw), sin(yaw), 0], [- sin(yaw), cos(yaw), 0], [0, 0, 1]]

#### 5.1 downwash gradient

#### 5.2 downwash angle

#### 6.1 roll effect

#### 6.2 pitch effect

#### 6.3 yaw effect
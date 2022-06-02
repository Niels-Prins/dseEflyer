from Stability.Code.solver import main

# Folder in which the aircraft is located.
# path = 'Aircraft/EFlyer'
path = 'Code/Testing/Aircraft/Cessna'

# Analysis with fuselage effects in case True.
fuselage = True

# Run the entire program.
main(path, fuselage=fuselage)

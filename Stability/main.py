from Code.surface import AerodynamicSurface
from Code.solver import main

# TODO: generalise this

# Folder in which the aircraft is located.
path = 'Aircraft/Example'

# All the surfaces belonging to the aircraft (symmetric vertical fin used in example).
wing = AerodynamicSurface(path=f'{path}/Wing')
vertical_left = AerodynamicSurface(path=f'{path}/VerticalLeft', symmetric=True, vertical=True)
vertical_right = AerodynamicSurface(path=f'{path}/VerticalRight', symmetric=True, vertical=True)

# List of all the surfaces, wing ALWAYS first item!
aircraft = [wing, vertical_left, vertical_right]

# Run the entire program.
main(path, aircraft)

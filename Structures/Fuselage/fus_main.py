#%%

from fus_idealization import *
from fus_loads import *
from fus_stress_calc import *

#%%

fuselage = Fuselage_structure().create_fuselage()
load_results =  Fuselage_apply_loads.apply_stress(fuselage.x, fuselage.data)
#%%
import numpy as np
from fus_discr import *

test = Fuselage_structure()
# %%
class Fuselage_apply_loads:
    def __init__(self, moment_distr, force_dist):
        # self.fuselage = Fuselage_structure().create_fuselage()
        pass

    # def shear_stress(self, shear_force):
    #     shearStress = []
    #     for idx in range(self.stringerNum):
    #         shearStress.append(
    #             (-shear_force / self.totalInertia)
    #             * self.areas[idx]
    #             * (self.z[idx] - self.neutralY)
    #         )

    #     residualStress = 0 - sum(shearStress)

    def bending_stress(self):
        pass


#%%

areas = [250, 400, 100, 100, 400, 250, 200, 200]
y = [100, 100, 50, -50, -100, -100, -30, 30]
x = [360, 120, 0, 0, 120, 360, 600, 600]

totalInertia = sum([areas[idx] * y[idx] ** 2 for idx in range(len(y))])

coordinates = (120, 0)


def shear_stress(shear_force, totalInertia, areas, zLoc, yLoc, wetted_area, shearCoor):
    shearStress = [0]
    moments = 0
    for idx in range(1, 8):
        localStress = (-shear_force / totalInertia) * areas[idx] * (
            zLoc[idx]
        ) + shearStress[idx - 1]
        
        shearStress.append(localStress)
        
        moments += -shearStress[idx] * (
            yLoc[idx + 1 if idx != len(yLoc) - 1 else 0]
            - shearCoor[0]
            - (yLoc[idx] - shearCoor[0])
        ) * (zLoc[idx] - shearCoor[1]) + shearStress[idx] * (
            zLoc[idx + 1 if idx != len(yLoc) - 1 else 0]
            - shearCoor[1]
            - (zLoc[idx] - shearCoor[1])
        ) * (
            yLoc[idx] - shearCoor[0]
        )
    residual = -(moments) / (2 * wetted_area)
    shearStress = [i + residual for i in shearStress]
    return shearStress


test = shear_stress(10000, totalInertia, areas, y, x, 97000, coordinates)
# %%

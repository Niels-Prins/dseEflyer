#%%
#Importing for stress calculations
import numpy as np
from fus_idealization import *
from fus_loads import *


# %%
#Initialize class for stress calculations
class Fuselage_apply_loads:
    @staticmethod
    def shear_stress_closed(shear_force, zLoc, totalInertia, areas, thickness, neutralY):
        shearflow = [0]
        # moments = 0
        for idx in range(1, len(zLoc)):
            localflow = (-shear_force / totalInertia) * areas[idx] * (
                zLoc[idx] - neutralY
            ) + shearflow[idx - 1] 

            shearflow.append(localflow)

        # Resultant shear stuff, not applicable rn
        #     moments += -shearStress[idx] * (
        #         yLoc[idx + 1 if idx != len(yLoc) - 1 else 0]
        #         - shearCoor[0]
        #         - (yLoc[idx] - shearCoor[0])
        #     ) * (zLoc[idx] - shearCoor[1]) + shearStress[idx] * (
        #         zLoc[idx + 1 if idx != len(yLoc) - 1 else 0]
        #         - shearCoor[1]
        #         - (zLoc[idx] - shearCoor[1])
        #     ) * (
        #         yLoc[idx] - shearCoor[0]
        #     )
        # # residual = -(moments) / (2 * wetted_area)
        # shearStress = [i for i in shearStress]  #+ residual
        shearStress = [x / thickness for x in shearflow]

        return shearStress

    # @staticmethod
    # def shear_stress_open(shear_force, zLoc, totalInertia, areas):
    #     shearStress = [0]
    #     for idx in range(1, len(zLoc)):
    #         localStress = (-shear_force / totalInertia) * areas[idx] * (
    #             zLoc[idx]
    #         ) + shearStress[idx - 1]

    #         shearStress.append(localStress)

    #     return shearStress

    @staticmethod
    def bending_stress_closed(moment, zLoc, neutralY, totalInertia):
        bendingStress = []
        for idx in range(len(zLoc)):
            stress = (moment / totalInertia) * (zLoc[idx] - neutralY)
            bendingStress.append(stress)

        return bendingStress
    @staticmethod
    def get_forces(x: float = 700, n: int = 12):
        v, m = VM(x, n)
        return v, m
    
    @staticmethod
    def apply_forces(obj, n):
        results = pd.DataFrame()
        # Get load
        for i in obj.x:
            v, m = Fuselage_apply_loads.get_forces(i / 1000, n)
            m = m * 1000

            temp = obj.data[obj.data["xcoor"] == i]
            temp.loc[:,"bendingstress"] = Fuselage_apply_loads.bending_stress_closed(
                m,
                temp["zcoor"],
                temp.loc[0, "neutralY"],
                temp.loc[0, "totalinertia"],
            )

            temp.loc[:,"shearstress"] = Fuselage_apply_loads.shear_stress_closed(
                v, temp["zcoor"], temp.loc[0, "totalinertia"], temp["area"], obj.skinThickness, obj.neutralY
            )
            results = results.append(temp, ignore_index=True)

        return results
        

if __name__ == "__main__":
    test = Fuselage_idealized()
    test.create_fuselage()
    result = Fuselage_apply_loads.apply_forces(test, 12)

# %%

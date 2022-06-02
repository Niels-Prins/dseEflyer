#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from math import *


class Fuselage_structure:
    def __init__(
        self,
        point_spacing_cross: int = 5,
        spacing_3d: int = 80,
        stringer_spacing: float = 140,
        stringer_thickness: float= 2.5, 
        stringer_width: float=30,
        stringer_height: float=30,
        skin_thickness: float = 2
        
    ):
        self.spacingCross = point_spacing_cross
        self.spacing_3d = spacing_3d
        self.stringerSpacing = stringer_spacing
        self.stringerThickness = stringer_thickness
        self.stringerWidth = stringer_width
        self.stringerHeight = stringer_height
        self.skinThickness = skin_thickness
        
    @staticmethod
    def calc_distance(x1, x2, y1, y2):
        distance = np.sqrt(
                    (x1 - x2) ** 2 + (y1 - y2) ** 2
                )
        return distance

    def create_2d(self):
        x_bot = [-400, 0, 400]
        y_bot = [0, -620, 0]
        bottomCurve = interpolate.interp1d(x_bot, y_bot, kind="quadratic")

        x_side = [400, 400 - 25.8141, 400 - 67.348, 400 - 109.619, 400 - 141.198, 250]
        y_side = [0, 88.06, 188.067, 288.067, 402.401, 500]
        side_curve = interpolate.interp1d(x_side, y_side, kind="cubic")

        x_top = np.array([250, 0, -250])
        y_top = [750, 1000, 750]
        topCurve = interpolate.interp1d(x_top, y_top, kind="quadratic")

        xCoorSide = []
        yCoorSide = []

        xCoorTop = []
        yCoorTop = []

        xCoorBot = []
        yCoorBot = []
        # Side curve spacing
        oldySide = 500
        oldyTop = 1000
        oldyBot = -620
        distanceSide = 0
        distanceTop = 0
        distanceBot = 0
        
        xChecking = np.arange(0, 400, 0.01)
        
        for i in range(1, len(xChecking)):
            if 0 <=xChecking[i] <= 250:
                yTop = topCurve(xChecking[i])
                distanceTop += Fuselage_structure.calc_distance(xChecking[i],xChecking[i - 1], yTop, oldyTop)
                if isclose(distanceTop, self.stringerSpacing, rel_tol=0.1):
                    xCoorTop.append(xChecking[i])
                    yCoorTop.append(yTop)
                    distanceTop = 0
                oldyTop = yTop
                
            elif 250 <= xChecking[i] <= 400:
                ySide = side_curve(xChecking[i])
                distanceSide += Fuselage_structure.calc_distance(xChecking[i],xChecking[i - 1], ySide, oldySide)
                if isclose(distanceSide, self.stringerSpacing, rel_tol=0.1):
                    xCoorSide.append(xChecking[i])
                    yCoorSide.append(ySide)
                    distanceSide = 0
                oldySide = ySide
                
            if 0 <= xChecking[i] <= 400:
                yBot = bottomCurve(xChecking[i])
                distanceBot += Fuselage_structure.calc_distance(xChecking[i],xChecking[i - 1], yBot, oldyBot)
                
                if isclose(distanceBot, self.stringerSpacing, rel_tol=0.1):
                    xCoorBot.append(xChecking[i])
                    yCoorBot.append(yBot)
                    distanceBot = 0
                oldyBot = yBot
                
        yCoorSide2 = np.arange(500, 750, 114)
        xCoorSide2 = np.full(len(yCoorSide2), 250)

        self.y = np.concatenate([xCoorBot, xCoorSide, xCoorSide2, xCoorTop])
        self.z = np.concatenate([yCoorBot, yCoorSide, yCoorSide2, yCoorTop])

    def plot_cross_sec(self):
        plt.scatter(self.x, self.y, s=0.5)
        plt.axis("equal")
        plt.show()

    def create_3d(self):
        self.x = np.arange(0, 765 + self.spacing_3d, self.spacing_3d)

    def plot_fuselage(self):
        self.create_3d()
        fig = plt.figure()
        ax = plt.axes(projection="3d")
        for i in range(len(self.x)):
            ax.scatter(self.x[i], self.y, self.z, s=1.0)

        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        ax.set_zlabel("z [mm]")
        ax.set_ylim(0, 1620)
        ax.set_xlim(0, 1620)
        plt.show()
        
    def stringer_area(self):
        self.stringerArea = self.stringerThickness *self.stringerWidth *2 + self.stringerThickness*(self.stringerHeight - 2*self.stringerThickness)
        
    def calc_b(self):
        pass
        
        

#%%
test = Fuselage_structure(10)
test.create_2d()
test.plot_fuselage()



# %%

# angle = atan(500/150)
# x = cos(angle) * 114
# %%

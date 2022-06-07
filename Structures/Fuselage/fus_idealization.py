#%%
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from math import *
import pandas as pd
from tqdm import tqdm


class Fuselage_structure:
    def __init__(
        self,
        point_spacing_cross: int = 10,
        spacing_3d: int = 100,
        stringer_spacing: float = 140,
        stringer_thickness: float = 2.5,
        stringer_width: float = 30,
        stringer_height: float = 30,
        skin_thickness: float = 2,
        density: float=1750 # density, [kg/m^3]
    ):
        self.spacingCross = point_spacing_cross
        self.spacing_3d = spacing_3d
        self.stringerSpacing = stringer_spacing
        self.stringerThickness = stringer_thickness
        self.stringerWidth = stringer_width
        self.stringerHeight = stringer_height
        self.skinThickness = skin_thickness
        self.x = np.arange(3696, 4940 + self.spacing_3d, self.spacing_3d)
        self.data = pd.DataFrame()
        self.density = density

    def create_fuselage(self):  # Function to create the fuselage in one go
        self.create_2d()
        self.calc_b()
        self.calc_inertiaZZ()
        self.create_3d()

    @staticmethod
    def calc_distance(
        x1, x2, y1, y2
    ):  # Basic function to calculate distance between two points
        distance = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
        return distance

    def create_2d(self):  # Create 2d cross-section
        # Data from Catia, interpolate curve through these points
        y_bot = [-400, 0, 400]
        z_bot = [0, -620, 0]
        bottomCurve = interpolate.interp1d(y_bot, z_bot, kind="quadratic")

        y_side = [400, 400 - 25.8141, 400 - 67.348, 400 - 109.619, 400 - 141.198, 250]
        z_side = [0, 88.06, 188.067, 288.067, 402.401, 500]
        side_curve = interpolate.interp1d(y_side, z_side, kind="cubic")

        y_top = np.array([250, 0, -250])
        z_top = [750, 1000, 750]
        topCurve = interpolate.interp1d(y_top, z_top, kind="quadratic")

        # Create empty coordinate lists for stringer/idealization
        yCoorSide, zCoorSide = [], []
        yCoorTop, zCoorTop = [0], [1000]
        yCoorBot, zCoorBot = [0], [-620]

        # Initilalize spacing
        oldZSide = 500
        oldZTop = 1000
        oldZBot = -620
        distanceSide = 0
        distanceTop = 0
        distanceBot = 0

        yChecking = np.arange(0, 400, 0.01)
        self.totalDistance = 0

        # Function to place the stringers equally spaced around the  50% fuselage
        for i in tqdm(range(1, len(yChecking))):
            if 0 <= yChecking[i] <= 250:
                zTop = topCurve(yChecking[i])
                distanceTop += Fuselage_structure.calc_distance(
                    yChecking[i], yChecking[i - 1], zTop, oldZTop
                )
                if isclose(distanceTop, self.stringerSpacing, rel_tol=0.1):
                    yCoorTop.append(yChecking[i])
                    zCoorTop.append(zTop)
                    self.totalDistance += distanceTop
                    distanceTop = 0
                oldZTop = zTop

            elif 250 <= yChecking[i] <= 400:
                zSide = side_curve(yChecking[i])
                distanceSide += Fuselage_structure.calc_distance(
                    yChecking[i], yChecking[i - 1], zSide, oldZSide
                )
                if isclose(distanceSide, self.stringerSpacing, rel_tol=0.1):
                    yCoorSide.append(yChecking[i])
                    zCoorSide.append(zSide)
                    self.totalDistance += distanceSide
                    distanceSide = 0
                oldZSide = zSide

            if 0 <= yChecking[i] <= 400:
                zBot = bottomCurve(yChecking[i])
                distanceBot += Fuselage_structure.calc_distance(
                    yChecking[i], yChecking[i - 1], zBot, oldZBot
                )

                if isclose(distanceBot, self.stringerSpacing, rel_tol=0.1):
                    yCoorBot.append(yChecking[i])
                    zCoorBot.append(zBot)
                    self.totalDistance += distanceTop
                    distanceBot = 0
                oldZBot = zBot

        # Sort coordinates to loop counter clockwise around the fuselage
        zCoorBot = [x for y, x in sorted(zip(yCoorBot, zCoorBot))]
        zCoorSide = [x for y, x in sorted(zip(yCoorSide, zCoorSide))][::-1]
        yCoorSide = yCoorSide[::-1]

        zCoorTop = [x for y, x in sorted(zip(yCoorTop, zCoorTop))][::-1]
        yCoorTop = yCoorTop[::-1]

        zCoorSide2 = np.arange(500, 750, self.stringerSpacing)
        yCoorSide2 = np.full(len(zCoorSide2), 250)

        yHalf = np.concatenate([yCoorBot, yCoorSide, yCoorSide2, yCoorTop])
        zHalf = np.concatenate([zCoorBot, zCoorSide, zCoorSide2, zCoorTop])

        # otherHalfy = yHalf[1:-1]

        # Create full fuselage from 1 half (If wanted)
        # self.y = np.concatenate([yHalf, -yHalf[1:-1][::-1]])
        # self.z = np.concatenate([zHalf, zHalf[1:-1][::-1]])
        self.y = yHalf  
        self.z = zHalf
        self.stringerNum = len(self.z)

    def plot_cross_sec(self):  # Plot 2d cross section
        plt.scatter(self.x, self.y, s=0.5)
        plt.axis("equal")
        plt.show()

    def create_3d(self):  # create 3d spacing of cross section
        for i in self.x:
            temp = pd.DataFrame()
            temp["xcoor"] = np.full((self.stringerNum), i)
            temp["ycoor"] = self.y
            temp["zcoor"] = self.z
            temp["area"] = self.areas
            temp["inertia"] = self.inertias
            temp["totalinertia"] = np.full((self.stringerNum), self.totalInertia)
            temp["neutralY"] = np.full((self.stringerNum), self.neutralY)
            self.data = self.data.append(temp)

    def plot_fuselage(self):  # Plot 3d fuselage
        fig = plt.figure()
        ax = plt.axes(projection="3d")
        for i in range(len(self.x)):
            temp = self.data[self.data["xcoor"] == self.x[i]]
            ax.scatter(self.x[i], temp["ycoor"], temp["zcoor"], s=1.0)

        ax.set_xlabel("x [mm]")
        ax.set_ylabel("y [mm]")
        ax.set_zlabel("z [mm]")
        ax.set_ylim(-810, 810)
        ax.set_xlim(3000, 5000)
        plt.show()

    def calc_stringer_area(self):  # Calculate area of a z-stringer
        self.stringerArea = (
            self.stringerThickness * self.stringerWidth * 2
            + self.stringerThickness
            * (self.stringerHeight - 2 * self.stringerThickness)
        )

    def calc_b(self):  # Calculate areas for idealization
        self.calc_stringer_area()
        areaslst = [self.stringerArea]
        for idx in range(1, int(self.stringerNum)-1):
            b_further = Fuselage_structure.calc_distance(
                self.y[idx], self.y[idx + 1], self.z[idx], self.z[idx + 1]
            )
            b_previous = Fuselage_structure.calc_distance(
                self.y[idx], self.y[idx - 1], self.z[idx], self.z[idx - 1]
            )

            if self.y[idx] == 0:
                area = self.stringerArea

            else:
                area = (
                    self.stringerArea
                    + (self.skinThickness * b_further / 6)
                    * (2 + self.y[idx + 1] / self.y[idx])
                    + (self.skinThickness * b_previous / 6)
                    * (2 + self.y[idx - 1] / self.y[idx])
                )
            areaslst.append(area)
            
        areaslst.append(self.stringerArea)

        # reversed = areaslst[::-1]
        self.areas = areaslst

    def calc_neutraly(
        self,
    ):  # Calculate the location of the neutral axis in z direction
        self.neutralY = sum(self.z * self.areas) / sum(self.areas)

    def calc_inertiaZZ(
        self,
    ):  # Calculate the moment of inertia per boom and total of structure
        self.calc_neutraly()
        self.inertias = self.areas * (self.z - self.neutralY) ** 2
        self.totalInertia = sum(self.inertias)
        
    def calc_weight(self):
        totalstringerweight = self.stringerArea * self.density * (self.x[-1] - self.x[0])*(self.stringerNum*2 - 2) /10**9
        totalskinweight = self.skinThickness * self.totalDistance *2 * (self.x[-1] - self.x[0]) * self.density / 10**9
        self.weight = totalstringerweight + totalskinweight
        

#%%
if __name__ == "__main__":
    test = Fuselage_structure(stringer_spacing=150)
    test.create_fuselage()
    test.plot_fuselage()
    test.calc_weight()
    

# %%

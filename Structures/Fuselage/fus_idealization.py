#%%
""" Program to create the idealization of the fuselage section. Calculates areas, moments of inertias and weight_batteries
Niels Prins"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from math import *
import pandas as pd
from tqdm import tqdm

#%%
# Create Class for structure idealization
class Fuselage_idealized:
    """Class for creating the idealization of the fuselage"""

    def __init__(
        self,
        spacing_3d: int = 10,
        stringer_spacing: float = 140,
        flange_thickness: float = 1.5,
        web_thickness: float = 1,
        stringer_width: float = 30,
        stringer_height: float = 30,
        skin_thickness: float = 1,
        density: float = 1480,  # density, [kg/m^3]
        start_x: float = 3696.0,
        end_x: float = 4940,
        E: float = 42.31906348135664e9,
        v: float = 0.28079904417448404,
        material_comp_stress: float = 533.8656197958434e6,
        material_tens_stress: float = 890.2474375138635e6,
        manual_place: bool = False,
        bot = [0.0, 250, 400],  # From LEFT to Right in y
        side = [250],  # From RIGHT to LEFT in y
        side2 = [750],  ## This one are the z coordinates, Bottom to TOP!
        top = [250, 0.0]  # From right to left in y

    ):
        self.spacing_3d = spacing_3d
        self.stringerSpacing = stringer_spacing
        self.flangeThickness = flange_thickness
        self.webThickness = web_thickness
        self.stringerWidth = stringer_width
        self.stringerHeight = stringer_height
        self.skinThickness = skin_thickness
        self.x = np.arange(start_x, end_x, self.spacing_3d)
        self.data = pd.DataFrame()
        self.density = density
        self.E = E
        self.v = v
        self.materialComp = material_comp_stress
        self.materialTens = material_tens_stress
        self.shear = 102
        self.start = start_x
        self.end = end_x
        self.manualPlace = manual_place
        self.topMan = top
        self.botMan = bot
        self.sideMan = side
        self.side2Man = side2

    def create_fuselage(self):
        """Function to create a fuselage section including properties"""
        self.create_2d()
        self.calc_b()
        self.calc_inertiaZZ()
        self.create_3d()
        self.calc_weight()
        self.crippling_comp()
        self.crippling_tens()

    @staticmethod
    def calc_distance(x1: float, x2: float, y1: float, y2: float):
        """Basic function to calculate distance between two points in 2D

        Args:
            x1 (float): x1 coordinate
            x2 (float): x2 coordinate
            y1 (float): y1 coordinate
            y2 (float): y3 coordinate

        Returns:
            float: distance
        """

        distance = np.sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
        return distance

    def create_2d(self):
        """Function to create one 2d cross section of the fuselage based on the Catia coordiantes"""
        # coordinates Bottom of fuselage
        y_bot = [-400, 0, 400]
        z_bot = [0, -620, 0]
        bottomCurve = interpolate.interp1d(y_bot, z_bot, kind="quadratic")

        # Side part of fuselage
        y_side = [400, 400 - 25.8141, 400 - 67.348, 400 - 109.619, 400 - 141.198, 250]
        z_side = [0, 88.06, 188.067, 288.067, 402.401, 500]
        sideCurve = interpolate.interp1d(y_side, z_side, kind="cubic")

        # Coordinates of top part of the fuselage
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
        self.totalDistance = 250

        # Function to place the stringers equally spaced around the 50% fuselage
        for i in tqdm(range(1, len(yChecking))):

            # Top part of fuselage
            if 0 <= yChecking[i] <= 250:
                zTop = topCurve(yChecking[i])
                distance_between = Fuselage_idealized.calc_distance(
                    yChecking[i], yChecking[i - 1], zTop, oldZTop
                )
                self.totalDistance += distance_between
                distanceTop += distance_between
                if isclose(
                    distanceTop, self.stringerSpacing, rel_tol=0.1
                ):  # If distance is sufficient enough place stringer
                    yCoorTop.append(yChecking[i])
                    zCoorTop.append(zTop)
                    distanceTop = 0
                oldZTop = zTop

            # Side part of fuselage
            elif 250 <= yChecking[i] <= 400:
                zSide = sideCurve(yChecking[i])
                distance_between = Fuselage_idealized.calc_distance(
                    yChecking[i], yChecking[i - 1], zSide, oldZSide
                )
                distanceSide += distance_between
                self.totalDistance += distance_between
                if isclose(
                    distanceSide, self.stringerSpacing, rel_tol=0.1
                ):  # If distance is sufficient enough place stringer
                    yCoorSide.append(yChecking[i])
                    zCoorSide.append(zSide)
                    distanceSide = 0
                oldZSide = zSide

            # Bottom part of fuselage
            if 0 <= yChecking[i] <= 400:
                zBot = bottomCurve(yChecking[i])
                distance_between = Fuselage_idealized.calc_distance(
                    yChecking[i], yChecking[i - 1], zBot, oldZBot
                )
                distanceBot += distance_between
                self.totalDistance += distance_between
                if isclose(
                    distanceBot, self.stringerSpacing, rel_tol=0.1
                ):  # If distance is sufficient enough place stringer
                    yCoorBot.append(yChecking[i])
                    zCoorBot.append(zBot)
                    distanceBot = 0
                oldZBot = zBot

        if self.manualPlace:
            zCoorBot = bottomCurve(self.botMan)
            yCoorBot = self.botMan

            zCoorTop = topCurve(self.topMan)
            yCoorTop = self.topMan

            zCoorSide = sideCurve(self.sideMan)
            yCoorSide = self.sideMan

            zCoorSide2 = self.side2Man
            yCoorSide2 = np.full(len(self.side2Man), 250)

        else:
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


    def create_3d(self):
        """Create 3d version of fuselage by placing multiple cross-sections behind eachother"""
        for i in self.x:
            temp = pd.DataFrame()
            temp["xcoor"] = np.full((self.stringerNum), i)
            temp["ycoor"] = self.y
            temp["zcoor"] = self.z
            temp["area"] = self.areas
            temp["inertia"] = self.inertias
            temp["totalinertia"] = np.full((self.stringerNum), self.totalInertia)
            temp["neutralY"] = np.full((self.stringerNum), self.neutralY)

            self.data = self.data.append(temp)  # Store all data in pandas dataframe

    def plot_fuselage(self):
        """Plot 3d model of fuselage"""
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

    def calc_stringer_area(self):
        """Calculate area of one stringer"""
        self.stringerArea = (
            self.flangeThickness * self.stringerWidth * 2
            + self.webThickness * self.stringerHeight
            + 2 * self.webThickness * self.flangeThickness
        )

    def calc_b(self):
        """Calculate the area of each of the booms in idealization. booms are placed at stringers"""
        self.calc_stringer_area()  # Get stringer area
        areaslst = [self.stringerArea]

        # Loop through all booms, add half of area of skin between neighbouring booms to current boom
        for idx in range(1, int(self.stringerNum) - 1):
            b_further = Fuselage_idealized.calc_distance(
                self.y[idx], self.y[idx + 1], self.z[idx], self.z[idx + 1]
            )
            b_previous = Fuselage_idealized.calc_distance(
                self.y[idx], self.y[idx - 1], self.z[idx], self.z[idx - 1]
            )

            if (
                self.y[idx] == 0
            ):  # Prevent dividing by 0 for booms placed in symmetry axis
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

        self.areas = areaslst

    def calc_neutraly(self):
        """Calculate z coordinate of neutral axis"""
        self.neutralY = sum(self.z * self.areas) / sum(self.areas)

    def calc_inertiaZZ(self):
        """Calculate moment of inertia zz"""
        self.calc_neutraly()
        self.inertias = self.areas * (self.z - self.neutralY) ** 2
        self.totalInertia = sum(self.inertias)

    def calc_weight(self):
        """Calculate weight of the fuselage in Kg"""
        totalstringerweight = (
            self.stringerArea
            * self.density
            * (self.end - self.start)
            * (self.stringerNum * 2 - 2)
            / 10 ** 9
        )
        totalskinweight = (
            self.skinThickness
            * self.totalDistance
            * 2
            * (self.end - self.start)
            * self.density
            / 10 ** 9
        )

        self.weight = totalstringerweight + totalskinweight

    def crippling_tens(self):
        bottopratio = 0.8 * (
            0.425
            * np.pi ** 2
            * self.E
            / (self.materialTens * 12 * (1 - self.v ** 2))
            * (self.flangeThickness / self.stringerWidth) ** 2
        ) ** (1 - 0.6)
        if bottopratio >= 1:
            bottopstress = self.materialTens
        elif bottopratio < 1:
            bottopstress = self.materialTens * bottopratio

        middleratio = 0.8 * (
            4
            * np.pi ** 2
            * self.E
            / (self.materialTens * 12 * (1 - self.v ** 2))
            * (self.webThickness / self.stringerHeight) ** 2
        ) ** (1 - 0.6)
        if middleratio >= 1:
            middlestress = self.materialTens
        elif middleratio < 1:
            middlestress = middleratio * self.materialTens

        self.stringerCripplingTens = (
            (
                2 * bottopstress * self.flangeThickness * self.stringerWidth
                + self.stringerHeight * self.webThickness * middlestress
            )
            / (
                2 * self.flangeThickness * self.stringerWidth
                + self.webThickness * self.stringerHeight
            )
            / 10 ** 6
        )

    def crippling_comp(self):
        bottopratio = 0.8 * (
            0.425
            * np.pi ** 2
            * self.E
            / (self.materialComp * 12 * (1 - self.v ** 2))
            * (self.flangeThickness / self.stringerWidth) ** 2
        ) ** (1 - 0.6)
        if bottopratio >= 1:
            bottopstress = self.materialComp
        else:
            bottopstress = self.materialComp * bottopratio

        middleratio = 0.8 * (
            4
            * np.pi ** 2
            * self.E
            / (self.materialComp * 12 * (1 - self.v ** 2))
            * (self.webThickness / self.stringerHeight) ** 2
        ) ** (1 - 0.6)
        if middleratio >= 1:
            middlestress = self.materialComp
        elif middleratio < 1:
            middlestress = middleratio * self.materialComp

        self.stringerCripplingComp = (
            (
                2 * bottopstress * self.flangeThickness * self.stringerWidth
                + self.stringerHeight * self.webThickness * middlestress
            )
            / (
                2 * self.flangeThickness * self.stringerWidth
                + self.webThickness * self.stringerHeight
            )
            / 10 ** 6
        )
        
    def plot_cross_sec(self):
        """Plot one cross section"""
        
        q = plt.scatter(self.y, self.z, c=self.areas, label="Booms",cmap="turbo", s=75)
        plt.scatter(-self.y, self.z, c=self.areas, cmap="turbo", s=75)
        plt.plot(self.y, self.z, color="deepskyblue", label="Fuselage\ncross-section")
        plt.plot(-self.y, self.z, color="deepskyblue")
        plt.grid(visible=True, which="major", color="#666666", linestyle="-")
        plt.minorticks_on()
        plt.grid(visible=True, which="minor", color="#999999", linestyle="-", alpha=0.2)
        plt.axis("equal")
        plt.xlabel("y [mm]")
        plt.ylabel("z [mm]")
        plt.legend()
        plt.colorbar(q, label="Boom Area [$mm^2$]")
        annotations = np.arange(1, int(len(self.y)/2)+1)
        for i, label in enumerate(annotations):
            plt.annotate(label, # this is the text
                 (self.y[i], self.z[i]), # these are the coordinates to position the label
                 textcoords="offset points", # how to position the text
                 xytext=(10,-7), # distance from text to points (x,y)
                 ha='center')
            
        
        annotations2 = np.arange(int(len(self.y)/2)+1, len(self.y)+1 )
        for i , label in enumerate(annotations2):
            plt.annotate(label, # this is the text
                 (self.y[i+ int(len(self.y)/2)], self.z[i+ int(len(self.y)/2)]), # these are the coordinates to position the label
                 textcoords="offset points", # how to position the text
                 xytext=(10,0), # distance from text to points (x,y)
                 ha='center')
            
        plt.savefig("crossIdea")
        plt.show()

#%%
if __name__ == "__main__":
    test = Fuselage_idealized(stringer_spacing=150)
    test.create_fuselage()
    test.plot_cross_sec()
    # test.plot_fuselage()
    # test.calc_weight()


# %%

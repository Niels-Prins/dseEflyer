# Stability program

### Purpose

### Guide

First, an example aircraft, which is a flying wing with two vertical surfaces at the wing tips, can be found 
in the 'Aircraft' folder under the name 'Example'. Each parameter is self-explanatory,
except the following items under 'control':

- motion, meaning the axis affected by the control surface, 'X' for ailerons, 'Y' for elevator and 'Z' for rudder;
- deflection, True if deflection is symmetric or False if asymmetric (e.g. ailerons);
- alpha, change in alpha due to controller deflection;
- moment, change in aerodynamic center moment due to controller deflection.

A step-by-step guide:

1. start by creating the main aircraft folder and fill in the aircraft.txt file as in the example folder.
2. create sub folders for each surface you would like to add e.g. 'Wing'.
3. for each of these folders fill in the values file as in the example folder and add the polars.txt file from any
aerodynamic program. A similar format has to be used as in the example folder in order for the program to work.
4. fill in the code below. Each surface in the program is created from the AerodynamicSurface class. The defaults are 
as follows, symmetric is True (e.g. wing), define False if needed (e.g. fin), vertical is False (e.g. fin) and 
main_body is False. In case downwash is applicable (e.g. stabilizer) main_body should be set to the object from
which the downwash is generated e.g. the wing in case of the horizontal stabilizer.

### Assumptions

- provided data contains no errors (e.g. negative span) and is in the same format as the example files;
- asymmetric surface coordinates run from - span to 0 instead of - span / 2 to span / 2;
- if a vertical surface is used the inputs are pre-transformed by 90 degrees along the roll axis;
- induced drag constant along span (elliptical distribution), C_L constant along wing span;
- small angle approximations are used except for sweep;
- only the turn rate aerodynamic effects are taken into account;
- incompressible flow, Mach effects neglected;
- X_ac shift along span for pitch rate effect neglected;
- only C_M_alpha and C_N_beta effect of fuselage are taken into account, others are neglected;
- main body wake effects neglected;
- Applied sweep and airfoil are intertwined;

### Limitations

- Scipy could not be used for integration as most methods return tuples, which Scipy cannot handle;

# def class_I():
#     structure = 360
#     batteries = 210
#     motor = 60
#     margin = 0.2
#     OEM = (1 + margin) * (structure + batteries + motor)
#
#     occupants = 180
#     TOM = OEM + occupants
#
#     return OEM, TOM, batteries, motor, occupants

def class_I():
    structure = 345
    batteries = 365
    motor = 60
    margin = 0.3
    OEM = (1 + margin) * (structure + batteries + motor)

    occupants = 180
    TOM = OEM + occupants

    return OEM, TOM, batteries, motor, occupants

def class_I():
    structure = 345
    batteries = 220
    motor = 60
    margin = 0.2
    OEM = (1 + margin) * (structure + batteries + motor)

    occupants = 180
    TOM = OEM + occupants

    print(TOM)

    return OEM, TOM, batteries, motor, occupants

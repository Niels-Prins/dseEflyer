def class_I():
    structure = 350
    batteries = 424
    motor = 66
    margin = 0.20
    OEM = (1 + margin) * (structure + batteries + motor)

    occupants = 180
    TOM = OEM + occupants

    print(f'Class I OEM: {OEM} \nClass I TOM: {TOM}')

    return OEM, TOM, batteries, motor, occupants

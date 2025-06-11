import single_stage_simulation as sim

engine_count = {
    "Rapier": 2,
    "Nerv": 2,
    "Dawn": 30
}

min_engine_count = {
    "Nerv": 1,
    "Rapier": 2,
    "Dawn": 1,
}

# Vall capture/escape = 270 m/s
# Kerbin escape = 850 m/s
# Tylo capture/escape = 800 m/s
flight_profile = [
    "900",

    "100",
    "Minmus 0",
    "100",

    "300",
    "Gilly 5000",
    "400",

    "300",
    "Pol 0",
    "300",

    "200",
    "Bop 20000",
    "200",

    "150",
    "Ike 2500",
    "150",

    "250",
    "Mun 2500",
    "250",

    "500",
    "Dres 2500",
    "500",

    "400",
    "Eeloo 0",
    "400",

    #"LiquidFuel 0.14",
    #"450",

    "LiquidFuel 0.1"
    "600",

    "300",
    "Vall 0",
    "300",

    #"600",
    #"Moho 2500",
    #"600",
]

sim.optimize_engine_layout(engine_count, min_engine_count, flight_profile,
                           allow_fractional = True, allow_new_engines = True, mass_in_orbit = 105000)

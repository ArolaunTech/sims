import single_stage_simulation as sim

engine_count = {
    "Rapier": 1,
    "Nerv": 3,
    "Dawn": 3,
}

min_engine_count = {
    "Rapier": 1,
    "Dawn": 0.1,
}

# Vall capture/escape = 270 m/s
# Kerbin escape = 850 m/s
# Tylo capture/escape = 800 m/s
flight_profile = [
    "850",
    "270",
    "Vall 6000",
    "270",
    "800",
    "Tylo 9100",
    "800",
]

sim.optimize_engine_layout(engine_count, min_engine_count, flight_profile,
                           allow_fractional = True, allow_new_engines = True, mass_in_orbit = 100000)

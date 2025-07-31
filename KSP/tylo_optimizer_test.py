import single_stage
import single_stage_simulation as sim

engine_count = {
	"Rapier": 2,
	"Nerv": 2,
	"Dawn": 30,
	"Wolfhound": 0
}

#engine_count = {'Rapier': 2, 'Nerv': 2.044980127165045, 'Dawn': 30.09975106871209, 'Cheetah': 0.05085751165946217, 'Poodle': 0.003505421848380852, 'Spider': 0.02753955597302982, 'Puff': 0.09994449478536359, 'Terrier': 0.0034910098003195335, 'Dart': 0.014640225194918711, 'Mastodon': 0.07934429121237305, 'Twitch': 0.03149914109509784}

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

	"270",
	"Vall 6000",
	"270",

	"800",
	"Tylo 10000",
	"800",
]

sim.optimize_engine_layout(
	engine_count, 
	min_engine_count, 
	flight_profile,
	allow_fractional = True, 
	allow_new_engines = True, 
	mass_in_orbit = 140000
)
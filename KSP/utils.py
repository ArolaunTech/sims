# Utility file for some functions that I think could be reusable.
import math

# Useful constants and data structures

GRAV_CONSTANT = 6.67 / 10.0 ** 11

# Map of common aliases
aliases = dict()

# Map of celestial body name to some useful quantities.
# Current is a dict of mass in kg and radius in m
bodies = dict()

# Map of fuel type name to mass per unit
fuel_masses = {
    "LiquidFuel": 5,
    "LFOx": 5,
    "Monopropellant": 4,
    "XenonGas": 0.1,
    "SolidFuel": 7.5,
}

# Map of fuel type to tuple of fuel count and mass ratio of best fuel tank
# Tank dry mass = fuel wet mass / (mass ratio)
fuel_tank_info = {
    "LiquidFuel": (50, 10.0),
    "LFOx": (40, 8.0),
    "Monopropellant": (120, 7.5), #Btw the FL-R750 tank has a better monoprop fraction
    "XenonGas": (720, 3),
    "SolidFuel": (0, 2), # idc about solid fuel. Doesn't make sense to use this line anyway
}

# Map of engine name to tuple of (vac isp, fueltype, thrust, fuel/sec, mass)
# Twinboar and SRBS are skipped because the built-in fuel tank requirement is annoying.
engine_info = {
    "Dawn":     (4200, "XenonGas",       2000,    0.486,   250),
    "Nerv":     (800,  "LiquidFuel",     60000,   1.53,    3000),
    "Wolfhound":(380,  "LFOx",           375000,  20.126,  3300),
    "Cheetah":  (355,  "LFOx",           125000,  7.181,   1000),
    "Poodle":   (350,  "LFOx",           250000,  14.568,  1750),
    "Terrier":  (345,  "LFOx",           60000,   3.547,   500),
    "Dart":     (340,  "LFOx",           180000,  10.797,  1000),
    "Rhino":    (340,  "LFOx",           2000000, 119.968, 9000),
    "Skiff":    (330,  "LFOx",           300000,  18.54,   1600),
    "Spark":    (320,  "LFOx",           20000,   1.275,   130),
    "Swivel":   (320,  "LFOx",           215000,  13.703,  1500),
    "Skipper":  (320,  "LFOx",           650000,  41.426,  3000),
    "Ant":      (315,  "LFOx",           2000,    0.129,   20),
    "Vector":   (315,  "LFOx",           1000000, 64.745,  4000),
    "Mammoth":  (315,  "LFOx",           4000000, 258.978, 15000),
    "Reliant":  (310,  "LFOx",           240000,  15.789,  1250),
    "Mainsail": (310,  "LFOx",           1500000, 98.683,  6000),
    "Cub":      (310,  "LFOx",           32000,   2.105,   180),
    "Bobcat":   (310,  "LFOx",           400000,  26.316,  2000),
    "Mastodon": (305,  "LFOx",           1350000, 90.271,  5000),
    "Thud":     (305,  "LFOx",           120000,  8.024,   900),
    "Rapier":   (305,  "LFOx",           180000,  12.036,  2000),
    "Kodiak":   (300,  "LFOx",           260000,  17.675,  1250),
    "Spider":   (290,  "LFOx",           2000,    0.141,   20),
    "Twitch":   (290,  "LFOx",           16000,   1.125,   80),
    "Puff":     (250,  "Monopropellant", 20000,   2.039,   90),
}

aliases["kerbol"] = "the sun"
aliases["mun"] = "the mun"

# Initialize dictionaries. Maybe later we can make this read from a text file instead of it being hardcoded.
bodies["the sun"] = {"mass": 1.7565459 * 10 ** 28, "radius": 261600000, "rotationalperiod": 432000}
bodies["moho"]    = {"mass": 2.5263314 * 10 ** 21, "radius": 250000,    "rotationalperiod": 1210000}
bodies["eve"]     = {"mass": 1.2243980 * 10 ** 23, "radius": 700000,    "rotationalperiod": 80500}
bodies["gilly"]   = {"mass": 1.2420363 * 10 ** 17, "radius": 13000,     "rotationalperiod": 28255}
bodies["kerbin"]  = {"mass": 5.2915158 * 10 ** 22, "radius": 600000,    "rotationalperiod": 21549.4251830898}
bodies["the mun"] = {"mass": 9.7599066 * 10 ** 20, "radius": 200000,    "rotationalperiod": 138984.376574476}
bodies["minmus"]  = {"mass": 2.6457580 * 10 ** 19, "radius": 60000,     "rotationalperiod": 40400}
bodies["duna"]    = {"mass": 4.5154270 * 10 ** 21, "radius": 320000,    "rotationalperiod": 65517.859375}
bodies["ike"]     = {"mass": 2.7821615 * 10 ** 20, "radius": 130000,    "rotationalperiod": 65517.8621348081}
bodies["dres"]    = {"mass": 3.2190937 * 10 ** 20, "radius": 138000,    "rotationalperiod": 34800}
bodies["jool"]    = {"mass": 4.2332127 * 10 ** 24, "radius": 6000000,   "rotationalperiod": 36000}
bodies["laythe"]  = {"mass": 2.9397311 * 10 ** 22, "radius": 500000,    "rotationalperiod": 52980.8790593796}
bodies["vall"]    = {"mass": 3.1087655 * 10 ** 21, "radius": 300000,    "rotationalperiod": 105962.088893924}
bodies["tylo"]    = {"mass": 4.2332127 * 10 ** 22, "radius": 600000,    "rotationalperiod": 211926.35802123}
bodies["bop"]     = {"mass": 3.7261090 * 10 ** 19, "radius": 65000,     "rotationalperiod": 544507.428516654}
bodies["pol"]     = {"mass": 1.0813507 * 10 ** 19, "radius": 44000,     "rotationalperiod": 901902.623531173}
bodies["eeloo"]   = {"mass": 1.1149224 * 10 ** 21, "radius": 210000,    "rotationalperiod": 19460}

for alias in aliases:
    bodies[alias] = bodies[aliases[alias]]

# Returns the fraction of your thrust that has to be upward in order to fight gravity.
def sine_constant_altitude_thrust_angle(thrust, craft_mass, vx, altitude, preset = "", mass = 0, radius = 0):
    if preset in bodies:
        mass = bodies[preset]["mass"]
        radius = bodies[preset]["radius"]

    r = radius + altitude
    vertical_accel = (GRAV_CONSTANT * mass / r - vx * vx) / r
    
    sine_angle = vertical_accel * craft_mass / thrust
    return min(sine_angle, 1)
    
# Returns the angle upwards you have to burn, given parameters, to maintain vertical velocity.
# Used for constant altitude landings, hence the name.
# Thrust should be in Newtons, not kN.
# Mass should be in kg, not tons.
# Returns pitch angle in radians.
def constant_altitude_thrust_angle(thrust, craft_mass, vx, altitude, preset = "", mass = 0, radius = 0):
    a = sine_constant_altitude_thrust_angle(thrust, craft_mass, vx, altitude, preset, mass, radius)
    return math.asin(a)

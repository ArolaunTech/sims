# Utility file for some functions that I think could be reusable.
import math

# Map of common aliases
aliases = dict()

# Map of celestial body name to some useful quantities.
# Current is a dict of mass in kg and radius in m
bodies = dict()

GRAV_CONSTANT = 6.67 / 10.0 ** 11

def initialize():
    aliases["kerbol"] = "the sun"
    aliases["mun"] = "the mun"

    bodies["the sun"] = {"mass": 1.7565459 * 10 ** 28, "radius": 261600000}
    bodies["moho"]    = {"mass": 2.5263314 * 10 ** 21, "radius": 250000}
    bodies["eve"]     = {"mass": 1.2243980 * 10 ** 23, "radius": 700000}
    bodies["gilly"]   = {"mass": 1.2420363 * 10 ** 17, "radius": 13000}
    bodies["kerbin"]  = {"mass": 5.2915158 * 10 ** 22, "radius": 600000}
    bodies["the mun"] = {"mass": 9.7599066 * 10 ** 20, "radius": 200000}
    bodies["minmus"]  = {"mass": 2.6457580 * 10 ** 19, "radius": 60000}
    bodies["duna"]    = {"mass": 4.5154270 * 10 ** 21, "radius": 320000}
    bodies["ike"]     = {"mass": 2.7821615 * 10 ** 20, "radius": 130000}
    bodies["dres"]    = {"mass": 3.2190937 * 10 ** 20, "radius": 138000}
    bodies["jool"]    = {"mass": 4.2332127 * 10 ** 24, "radius": 6000000}
    bodies["laythe"]  = {"mass": 2.9397311 * 10 ** 22, "radius": 500000}
    bodies["vall"]    = {"mass": 3.1087655 * 10 ** 21, "radius": 300000}
    bodies["tylo"]    = {"mass": 4.2332127 * 10 ** 22, "radius": 600000}
    bodies["bop"]     = {"mass": 3.7261090 * 10 ** 19, "radius": 65000}
    bodies["pol"]     = {"mass": 1.0813507 * 10 ** 19, "radius": 44000}
    bodies["eeloo"]   = {"mass": 1.1149224 * 10 ** 21, "radius": 210000}

# Returns gravitational acceleration minus centripetal acceleration. Positive = downward
# Takes in parameters of the parent body, altimeter height, and horizontal velocity.
def effective_gravity(vx, altitude, preset = "", mass = 0, radius = 0):
    if preset in aliases:
        preset = aliases[preset]
    if preset in bodies:
        mass = bodies[preset]["mass"]
        radius = bodies[preset]["radius"]
    else:
        if mass == 0 or radius == 0:
            # Custom planetary body is messed up
            print("No planet parameters found!")
            return 0

    r = radius + altitude
    
    # Centripetal acceleration = v^2 / r
    ac = float(vx * vx) / r
    # Gravitational acceleration = GM / r ^ 2
    ag = GRAV_CONSTANT * mass / (r * r)

    return ag - ac

# Returns the angle upwards you have to burn, given parameters, to maintain vertical velocity.
# Used for constant altitude landings, hence the name.
# Thrust should be in Newtons, not kN.
# Mass should be in kg, not tons.
# Returns pitch angle in radians.
def constant_altitude_thrust_angle(thrust, craft_mass, vx, altitude, preset = "", mass = 0, radius = 0):
    accel = float(thrust) / craft_mass
    vertical_accel = effective_gravity(vx, altitude, preset, mass, radius)

    # vertical_accel should be accel * sin(angle)
    a = vertical_accel / accel
    if (a > 1):
        # Cannot maintain velocity!
        return math.pi / 2

    return math.asin(a)

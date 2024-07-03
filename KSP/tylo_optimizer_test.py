import utils
import single_stage

import math

planet_name = "tylo"
altitude = 9100

utils.initialize()

# Calculate initial orbital velocity
r = utils.bodies[planet_name]["radius"] + altitude
ag = utils.GRAV_CONSTANT * utils.bodies[planet_name]["mass"] / (r * r)

vx = math.sqrt(ag * r)
orbital_velocity = vx

mass_in_orbit = 466650
engines = [
    (4200, "XenonGas",   2000,  0.486),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (800,  "LiquidFuel", 60000, 1.53),
    (380,  "LFOx",       375000, 20.126),
    (380,  "LFOx",       375000, 20.126),
    (380,  "LFOx",       375000, 20.126),
    (380,  "LFOx",       375000, 20.126),
    (380,  "LFOx",       375000, 20.126),
    (380,  "LFOx",       375000, 20.126),
    (305,  "LFOx",       180000, 12.036),
    (305,  "LFOx",       180000, 12.036),
    (305,  "LFOx",       180000, 12.036),
    (305,  "LFOx",       180000, 12.036),
    (305,  "LFOx",       180000, 12.036),
    (305,  "LFOx",       180000, 12.036),
    (305,  "LFOx",       180000, 12.036),
    (305,  "LFOx",       180000, 12.036),
    (305,  "LFOx",       180000, 12.036),
    (305,  "LFOx",       180000, 12.036),
    (305,  "LFOx",       180000, 12.036),
    (305,  "LFOx",       180000, 12.036),
    (305,  "LFOx",       180000, 12.036),
    (305,  "LFOx",       180000, 12.036),
    (305,  "LFOx",       180000, 12.036),
    (305,  "LFOx",       180000, 12.036),
    (305,  "LFOx",       180000, 12.036),
]

fuel_levels = dict()

craft = single_stage.SingleStageCraft()
craft.manual_instantiate(engines, fuel_levels, mass_in_orbit)

craft.toggle_engines_by_isp(4200)
craft.toggle_engines_by_isp(800)
list_of_isp = [4200, 800, 380, 305]
next_engines_index = 2

# Landing.
elapsed_time = 0
iteration_count = 0
tick_time = 0.1
ticks_per_log = 100

while (vx > 0):
    pitch_angle = utils.constant_altitude_thrust_angle(craft.total_thrust,
                                                 craft.wet_mass, vx,
                                                 altitude, preset=planet_name)
    current_effective_isp = craft.get_average_isp() * math.cos(pitch_angle)
    # Check if it's more efficient to turn on some engines
    if next_engines_index < len(list_of_isp):
        craft.toggle_engines_by_isp(list_of_isp[next_engines_index])
        test_pitch_angle = utils.constant_altitude_thrust_angle(craft.total_thrust,
                                                          craft.wet_mass, vx,
                                                          altitude, preset=planet_name)
        test_effective_isp = craft.get_average_isp() * math.cos(test_pitch_angle)
        if test_effective_isp > current_effective_isp:
            pitch_angle = test_pitch_angle
            print(f"Turning on engines with Isp {list_of_isp[next_engines_index]}")
            next_engines_index += 1
            
        else:
            craft.toggle_engines_by_isp(list_of_isp[next_engines_index])
    
    # Update velocity
    vx -= tick_time * math.cos(pitch_angle) * craft.total_thrust / craft.wet_mass 

    # Simulate burning
    craft.simulate_burn(tick_time)
    
    # Print stuff
    pitch_in_degrees = int(pitch_angle * 180 / math.pi)
    if (iteration_count % ticks_per_log == 0):
        print(f"T={iteration_count * tick_time}: \tPitch: {pitch_in_degrees}\tVelocity: {int(vx)}")
    iteration_count += 1

print("Fuel used per type:")
total_fuel_mass = 0
for fuel_type, fuel_level in craft.fuel_levels.items():
    mass = -fuel_level * craft.fuel_masses[fuel_type]
    print(f"{fuel_type}\t{-fuel_level}\tMass: {mass}")
    total_fuel_mass += mass
print(f"Total fuel mass: {total_fuel_mass}")
print("Landed!")
iteration_count = 0
print("Liftoff!")
# Ascent
next_engines_index -= 1
while (vx < orbital_velocity):
    pitch_angle = utils.constant_altitude_thrust_angle(craft.total_thrust,
                                                 craft.wet_mass, vx,
                                                 altitude, preset=planet_name)
    current_effective_isp = craft.get_average_isp() * math.cos(pitch_angle)
    # Check if it's more efficient to turn off some engines
    if next_engines_index >= 1:
        craft.toggle_engines_by_isp(list_of_isp[next_engines_index])
        test_pitch_angle = utils.constant_altitude_thrust_angle(craft.total_thrust,
                                                          craft.wet_mass, vx,
                                                          altitude, preset=planet_name)
        test_effective_isp = craft.get_average_isp() * math.cos(test_pitch_angle)
        if test_effective_isp > current_effective_isp:
            pitch_angle = test_pitch_angle
            print(f"Turning off engines with Isp {list_of_isp[next_engines_index]}")
            next_engines_index -= 1
            
        else:
            craft.toggle_engines_by_isp(list_of_isp[next_engines_index])
    
    # Update velocity
    vx += tick_time * math.cos(pitch_angle) * craft.total_thrust / craft.wet_mass 

    # Simulate burning
    craft.simulate_burn(tick_time)
    
    # Print stuff
    pitch_in_degrees = int(pitch_angle * 180 / math.pi)
    if (iteration_count % ticks_per_log == 0):
        print(f"T={iteration_count * tick_time}: \tPitch: {pitch_in_degrees}\tVelocity: {int(vx)}")
    iteration_count += 1

print("Fuel used per type:")
total_fuel_mass = 0
for fuel_type, fuel_level in craft.fuel_levels.items():
    mass = -fuel_level * craft.fuel_masses[fuel_type]
    print(f"{fuel_type}\t{-fuel_level}\tMass: {mass}")
    total_fuel_mass += mass
print(f"Total fuel mass: {total_fuel_mass}")

import utils
import single_stage

import math
import random
import time

def simulate_landing_with_craft(craft, planet_name, altitude, logging = False):
    r = utils.bodies[planet_name]["radius"] + altitude
    ag = utils.GRAV_CONSTANT * utils.bodies[planet_name]["mass"] / (r * r)

    vx = math.sqrt(ag * r)
    orbital_velocity = vx

    # Build list of specific impulse
    list_of_isp = []
    for engine, count in craft.engines.items():
        if count == 0:
            continue
        isp = utils.engine_info[engine][0]
        if isp not in list_of_isp:
            list_of_isp.append(isp)
    list_of_isp.sort(reverse=True)
    next_engines_index = 1

    # Landing
    lithobrake_speed = 60

    elapsed_time = 0
    iteration_count = 0
    tick_time = 0.2
    ticks_per_log = 50
    craft.toggle_engines_by_isp(list_of_isp[0])

    n_groups = len(list_of_isp)
    
    while (vx > lithobrake_speed):
        sine_pitch_angle = utils.sine_constant_altitude_thrust_angle(craft.total_thrust,
                                                           craft.wet_mass, vx,
                                                           altitude, preset=planet_name)
        thrust_effectiveness = math.sqrt(1 - sine_pitch_angle * sine_pitch_angle)
        current_effective_isp = craft.get_average_isp() * thrust_effectiveness
        # Check if it's more efficient to turn on some engines
        if next_engines_index < n_groups:
            craft.toggle_engines_by_isp(list_of_isp[next_engines_index])
            sine_test_angle = utils.sine_constant_altitude_thrust_angle(craft.total_thrust,
                                                                    craft.wet_mass, vx,
                                                                    altitude, preset=planet_name)
            test_effective_isp = craft.get_average_isp() * math.sqrt(1 - sine_test_angle * sine_test_angle)
            if test_effective_isp > current_effective_isp:
                sine_pitch_angle = sine_test_angle
                if (logging):
                    print(f"Turning on engines with Isp {list_of_isp[next_engines_index]}")
                next_engines_index += 1
            
            else:
                craft.toggle_engines_by_isp(list_of_isp[next_engines_index])
    
        # Because this is tick-based with a pretty coarse simulation, if velocity - target is negative
        # at the end of a tick, the optimizer thinks that more fuel was used and penalizes not reaching
        # the target speed on exactly tick borders.
        # This is a simple linear approximation for the last fraction of a tick.
        
        x_accel = thrust_effectiveness * craft.total_thrust / craft.wet_mass
        burn_time = tick_time
        delta_speed = tick_time * x_accel
        target_delta_speed = vx - lithobrake_speed
        if delta_speed > target_delta_speed:
            burn_time = target_delta_speed / x_accel
        
        # Update velocity    
        vx -= burn_time * x_accel
        
        # Simulate burning
        craft.simulate_burn(burn_time)
    
        # Print stuff
        if (logging and iteration_count % ticks_per_log == 0):
            pitch_in_degrees = int(math.asin(sine_pitch_angle) * 180 / math.pi)
            print(f"T={iteration_count * tick_time}: \tPitch: {pitch_in_degrees}\tVelocity: {int(vx)}")
        iteration_count += 1
    if (logging):
        print("Fuel used per type:")
    total_fuel_mass = 0
    for fuel_type, fuel_level in craft.fuel_used.items():
        mass = fuel_level * utils.fuel_masses[fuel_type]
        if (logging):
            print(f"{fuel_type}\t{fuel_level}\tMass: {mass}")
        total_fuel_mass += mass
    if (logging):
        print(f"Total fuel mass: {total_fuel_mass}")
        print("Landed!")
    iteration_count = 0
    if (logging):
        print("Liftoff!")
    # Ascent
    vx = 0
    next_engines_index -= 1
    while (vx < orbital_velocity):
        sine_pitch_angle = utils.sine_constant_altitude_thrust_angle(craft.total_thrust,
                                                           craft.wet_mass, vx,
                                                           altitude, preset=planet_name)
        thrust_effectiveness = math.sqrt(1 - sine_pitch_angle * sine_pitch_angle)
        current_effective_isp = craft.get_average_isp() * thrust_effectiveness
        # Check if it's more efficient to turn off some engines
        if next_engines_index > 0:
            craft.toggle_engines_by_isp(list_of_isp[next_engines_index])
            sine_test_angle = utils.sine_constant_altitude_thrust_angle(craft.total_thrust,
                                                                    craft.wet_mass, vx,
                                                                    altitude, preset=planet_name)
            test_effective_isp = craft.get_average_isp() * math.sqrt(1 - sine_test_angle * sine_test_angle)
            if test_effective_isp > current_effective_isp:
                sine_pitch_angle = sine_test_angle
                if (logging):
                    print(f"Turning off engines with Isp {list_of_isp[next_engines_index]}")
                next_engines_index -= 1
            
            else:
                craft.toggle_engines_by_isp(list_of_isp[next_engines_index])
    
        x_accel = thrust_effectiveness * craft.total_thrust / craft.wet_mass
        burn_time = tick_time
        delta_speed = tick_time * x_accel
        target_delta_speed = orbital_velocity - vx
        if delta_speed > target_delta_speed:
            burn_time = target_delta_speed / x_accel

        # Update velocity
        vx += burn_time * x_accel

        # Simulate burning
        craft.simulate_burn(burn_time)

        # Print stuff
        if (logging and iteration_count % ticks_per_log == 0):
            pitch_in_degrees = int(math.asin(sine_pitch_angle) * 180 / math.pi)
            print(f"T={iteration_count * tick_time}: \tPitch: {pitch_in_degrees}\tVelocity: {int(vx)}")
        iteration_count += 1
    craft.turn_engines_off()
    if (logging):
        print("Fuel used per type:")
    total_fuel_mass = 0
    for fuel_type, fuel_level in craft.fuel_used.items():
        mass = fuel_level * utils.fuel_masses[fuel_type]
        if (logging):
            print(f"{fuel_type}\t{fuel_level}\tMass: {mass}")
        total_fuel_mass += mass
    if (logging):
        print(f"Total fuel mass: {total_fuel_mass}")


def payload_capacity(craft):
    # Mass - total engine mass - total (used) fuel tank mass
    mass = craft.wet_mass
    for engine, count in craft.engines.items():
        mass -= utils.engine_info[engine][4] * count
    for fuel_type, amount in craft.fuel_used.items():
        fuel_mass = utils.fuel_masses[fuel_type] * amount
        fuel_tank_mass = float(fuel_mass) / utils.fuel_tank_info[fuel_type][1]
        mass -= (fuel_tank_mass)
    return mass


planet_name = "tylo"
altitude = 9100

mass_in_orbit = 146000
engines = []

engine_count = {
    #"Nerv": 8.0,
    #"Wolfhound": 1.0,
    #"Rapier": 4.0,
    "Nerv": 7,
    "Wolfhound": 0.5,
    "Rapier": 4.0,
    "Rhino": 0,
}
min_engine_count = {
    "Nerv": 4.0,
    "Rapier": 4.0,
}


improved_last_iteration = True

best_engine_count = engine_count.copy()
output_ascent_profiles = False
allow_fractional = True
allow_new_engines = False
# Base case
start_time = time.time()
simulation_time = 0
craft = single_stage.SingleStageCraft()
fuel_levels = dict()
craft.manual_instantiate(engine_count, fuel_levels, mass_in_orbit)
simulate_landing_with_craft(craft, "vall", 6000, logging = output_ascent_profiles)
simulate_landing_with_craft(craft, planet_name, altitude, logging = output_ascent_profiles)
max_payload_remaining = payload_capacity(craft)
print(f"Reference payload capacity: {max_payload_remaining}")
print(f"Reference payload fraction: {max_payload_remaining / mass_in_orbit}")

# Simulated annealing.
# Temperature controls the average magnitude of proposed changes to the engine count.
# This will cool down over time until it's equal to one, at which point we do a
# deterministic greedy search.
# Hopefully this can hit a bunch of different local maxima.
temperature = 10.0
cooling_factor = 0.99
n_iterations = 0
anneal_threshold = 1
if (allow_fractional):
    temperature = 0.01
    anneal_threshold = 0.001
while improved_last_iteration or temperature > anneal_threshold:
    n_iterations += 1
    if (n_iterations % 20 == 0):
        print(f"Temperature: {temperature}")
    if (temperature > anneal_threshold and not improved_last_iteration):
        temperature *= cooling_factor
        if (temperature < anneal_threshold):
            temperature = anneal_threshold
    improved_last_iteration = False
    # Each iteration:
    # Pick a random number of engines to try adding or removing, and see what happens
    
    for engine, engine_info in utils.engine_info.items():
        if not allow_new_engines and engine not in engine_count:
            continue
        if allow_fractional:
            add_value = random.random() * temperature
        else:
            add_value = random.randrange(1, int(temperature) + 1)
        # Try adding one or more engines
        test_engine_count = best_engine_count.copy()
        if engine in test_engine_count:
            test_engine_count[engine] += add_value
        else:
            test_engine_count[engine] = add_value
        fuel_levels = dict()
        test_craft = single_stage.SingleStageCraft()
        test_craft.manual_instantiate(test_engine_count, fuel_levels, mass_in_orbit)
        sim_start = time.time()
        simulate_landing_with_craft(test_craft, "vall", 6000, logging = False)
        simulate_landing_with_craft(test_craft, planet_name, altitude, logging = False)
        sim_end = time.time()
        simulation_time += sim_end - sim_start
        test_payload = payload_capacity(test_craft)
        if test_payload > max_payload_remaining:
            print(f"Added {add_value}x {engine}")
            improved_last_iteration = True
            max_payload_remaining = test_payload
            print(f"Payload capacity: {test_payload}")
            print(f"Payload fraction: {test_payload / mass_in_orbit}")
            best_engine_count = test_engine_count

        # Removing an engine
        if allow_fractional:
            subtract_value = random.random() * temperature
        else:
            subtract_value = random.randrange(1, int(temperature) + 1)
        test_engine_count = best_engine_count.copy()
        if engine in test_engine_count:
            min_engines = 0
            if engine in min_engine_count:
                min_engines = min_engine_count[engine]
            subtract_value = min(subtract_value, test_engine_count[engine] - min_engines)
            test_engine_count[engine] -= subtract_value
            if test_engine_count[engine] == 0:
                del test_engine_count[engine]
            fuel_levels = dict()
        else:
            continue
        craft = single_stage.SingleStageCraft()
        craft.manual_instantiate(test_engine_count, fuel_levels, mass_in_orbit)
        sim_start = time.time()
        simulate_landing_with_craft(craft, "vall", 6000, logging = False)
        simulate_landing_with_craft(craft, planet_name, altitude)
        sim_end = time.time()
        simulation_time += sim_end - sim_start
        test_payload = payload_capacity(craft)
        if test_payload > max_payload_remaining:
            print(f"Removed {subtract_value}x {engine}")
            improved_last_iteration = True
            max_payload_remaining = test_payload
            print(f"Payload fraction: {test_payload / mass_in_orbit}")
            print(f"Payload capacity: {test_payload}")
            best_engine_count = test_engine_count

end_time = time.time()
print("Local max reached")
print(f"Iterations: {n_iterations}")
for a, b in best_engine_count.items():
    if b > 0:
        print(f"{a}: {b}")
# Just for fun, print out ascent profile info
craft = single_stage.SingleStageCraft()
craft.manual_instantiate(best_engine_count, dict(), mass_in_orbit)
simulate_landing_with_craft(craft, "vall", 6000, logging=output_ascent_profiles)
simulate_landing_with_craft(craft, planet_name, altitude, logging=output_ascent_profiles)

print(f"Total time: {end_time - start_time}")
print(f"Simulation time: {simulation_time}")

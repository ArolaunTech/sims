import utils
import single_stage

import math

planet_name = "tylo"
altitude = 9100

# Calculate initial orbital velocity

mass_in_orbit = 466650
engines = []

engine_count = {
    "Nerv": 15,
    "Rapier": 15,
    "Dawn": 15,
}
min_engine_count = {
    "Nerv": 15,
    "Rapier": 15,
}

def simulate_landing_with_craft(craft, planet_name, altitude, logging = False):
    r = utils.bodies[planet_name]["radius"] + altitude
    ag = utils.GRAV_CONSTANT * utils.bodies[planet_name]["mass"] / (r * r)

    vx = math.sqrt(ag * r)
    orbital_velocity = vx

    # Build list of specific impulse
    list_of_isp = []
    for engine in craft.engines:
        isp = utils.engine_info[engine][0]
        if isp not in list_of_isp:
            list_of_isp.append(isp)
    list_of_isp.sort(reverse=True)
    next_engines_index = 1

    # Landing

    elapsed_time = 0
    iteration_count = 0
    tick_time = 1
    ticks_per_log = 10
    craft.toggle_engines_by_isp(list_of_isp[0])
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
                if (logging):
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
        if (logging and iteration_count % ticks_per_log == 0):
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
                if (logging):
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
            if (logging):
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


def payload_capacity(craft):
    # Mass - total engine mass - total (used) fuel tank mass
    mass = craft.wet_mass
    for engine in craft.engines:
        mass -= utils.engine_info[engine][4]
    for fuel_type, amount in craft.fuel_used.items():
        fuel_mass = utils.fuel_masses[fuel_type] * amount
        fuel_tank_mass = float(fuel_mass) / utils.fuel_tank_info[fuel_type][1]
        mass -= (fuel_tank_mass)
    return mass

def build_engine_array_from_count(engine_count):
    engine_list = []
    for a, b in engine_count.items():
        for i in range(b):
            engine_list.append(a)
    return engine_list

improved_last_iteration = True

best_engine_count = engine_count.copy()

# Base case

craft = single_stage.SingleStageCraft()
fuel_levels = dict()
engines = build_engine_array_from_count(engine_count)
craft.manual_instantiate(engines, fuel_levels, mass_in_orbit)
simulate_landing_with_craft(craft, "tylo", altitude)
max_payload_remaining = payload_capacity(craft)
print(f"Reference payload capacity (Rapier/Nerv only): {max_payload_remaining}")

while improved_last_iteration:
    improved_last_iteration = False
    # Each iteration:
    # Try adding one engine and removing one engine, and do a simulation each time
    for engine, engine_info in utils.engine_info.items():
        # Try adding an engine
        test_engine_count = best_engine_count.copy()
        if engine in test_engine_count:
            test_engine_count[engine] += 1
        else:
            test_engine_count[engine] = 1
        fuel_levels = dict()
        test_engines = build_engine_array_from_count(test_engine_count)
        test_craft = single_stage.SingleStageCraft()
        test_craft.manual_instantiate(test_engines, fuel_levels, mass_in_orbit)
        simulate_landing_with_craft(test_craft, "tylo", altitude, logging = False)
        test_payload = payload_capacity(test_craft)
        if test_payload > max_payload_remaining:
            print(f"Added {engine}")
            improved_last_iteration = True
            max_payload_remaining = test_payload
            print(f"Payload capacity: {test_payload}")
            best_engine_count = test_engine_count
        # Removing an engine
        test_engine_count = best_engine_count.copy()
        if engine in test_engine_count:
            if engine in min_engine_count and test_engine_count[engine] == min_engine_count[engine]:
                # Can't remove this engine as it's needed for something else
                continue
            else:
                if test_engine_count[engine] > 0:
                    test_engine_count[engine] -= 1
                fuel_levels = dict()
        else:
            continue
        test_engines = build_engine_array_from_count(test_engine_count)
        craft = single_stage.SingleStageCraft()
        craft.manual_instantiate(test_engines, fuel_levels, mass_in_orbit)
        simulate_landing_with_craft(craft, "tylo", altitude)
        test_payload = payload_capacity(craft)
        if test_payload > max_payload_remaining:
            print(f"Removed {engine}")
            improved_last_iteration = True
            max_payload_remaining = test_payload
            print(f"Payload capacity: {test_payload}")
            best_engine_count = test_engine_count

print("Local max reached")
for a, b in best_engine_count.items():
    if b > 0:
        print(f"{a}: {b}")

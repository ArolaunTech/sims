from utils import *
import utils
import single_stage

import math
import random
import time
from copy import deepcopy

def simulate_spring_landing_with_craft(craft, planet_name, altitude, k1, k2, mid1, mid2, starting_altitude, logging=False, refuel=False):
	planetRadius = utils.bodies[planet_name]["radius"]
	planetGravParameter = utils.GRAV_CONSTANT * utils.bodies[planet_name]["mass"]
	planetRotationalPeriod = utils.bodies[planet_name]["rotationalperiod"]

	rotvel = 2 * math.pi * (planetRadius + altitude) / planetRotationalPeriod

	lithobrake_speed = 60

	numEngines = len(craft.engines)

	# State
	x = planetRadius + starting_altitude
	y = 0
	vx = 0
	vy = math.sqrt(planetGravParameter / x)
	time = 0
	timesteps = 0

	vh, vv = vy, vx
	sma, eccentricity, periapsis = 0, 0, 0
	dist = 0

	landingdist, landingvv = 0, 0

	# Simulation
	timestep = 0.1

	# Logging
	timestepsPerLog = round(10 / timestep)

	while vh > rotvel + lithobrake_speed:
		# Calculations
		dist = math.sqrt(x * x + y * y)

		vh = (-y * vx + x * vy) / dist
		vv = (x * vx + y * vy) / dist

		# Gravity
		grav = planetGravParameter / dist / dist
		gx = -x * grav / dist
		gy = -y * grav / dist

		# Thrust
		vaccel = grav - vh * vh / dist
		vaccel -= k1 * (dist - mid1 - planetRadius)

		craft.toggle_engines_by_vaccel(vaccel)

		accel = craft.total_thrust / craft.wet_mass
		bestvaccelclamped = max(min(vaccel, accel), -accel)
		besthaccel = math.sqrt(accel * accel - bestvaccelclamped * bestvaccelclamped)

		ax = (y * besthaccel + x * bestvaccelclamped) / dist
		ay = (-x * besthaccel + y * bestvaccelclamped) / dist

		# Integration
		if timesteps % timestepsPerLog == 0 and logging:
			print(f"T = {round(time, 1)}s:", end="\t")
			print(f"Alt: {round(dist - planetRadius, 2)}m", end="\t")
			print(f"Hvel: {round(vh)} m/s", end="\t")
			print(f"Vvel: {round(vv, 1)} m/s")

		vx += (ax + gx) * timestep
		vy += (ay + gy) * timestep
		x += vx * timestep
		y += vy * timestep
		time += timestep
		timesteps += 1

		craft.simulate_burn(timestep)

		# End conditions
		landingdist = dist
		landingvv = vv

		if craft.wet_mass < 0:
			return 0, -vh, dist - planetRadius, vv, landingdist - planetRadius, landingvv, time, time
		if dist < planetRadius:
			return 1, -vh, dist - planetRadius, vv, landingdist - planetRadius, landingvv, time, time

	if logging:
		print(f"T = {round(time, 1)}s:", end="\t")
		print(f"Alt: {round(dist - planetRadius, 2)}m", end="\t")
		print(f"Hvel: {round(vh)} m/s", end="\t")
		print(f"Vvel: {round(vv, 1)} m/s")

	error = abs(dist - altitude - planetRadius)
	if error > 1:
		if dist - planetRadius < altitude:
			return 2, craft.wet_mass, dist - planetRadius, vv, landingdist - planetRadius, landingvv, time, time
		return 3, craft.wet_mass, dist - planetRadius, vv, landingdist - planetRadius, landingvv, time, time

	if logging:
		print("\nFuel used per type:")
	
		total_fuel_mass = 0
		for fuel_type, fuel_level in craft.fuel_used.items():
			mass = fuel_level * utils.fuel_masses[fuel_type]
			
			print(f"{fuel_type}\t{fuel_level}\tMass: {mass}")
			total_fuel_mass += mass
		
		print(f"Total fuel mass: {total_fuel_mass}")
		print(f"\nSuccessfully landed, wet mass = {craft.wet_mass} kg")
		print("Liftoff!\n")

	if refuel:
		craft.refuel()

	# Successfully landed
	x = planetRadius + altitude
	y = 0
	vx = 0
	vy = rotvel

	landingtime = time

	time = 0
	timesteps = 0

	while periapsis < planetRadius:
		# Calculations
		dist = math.sqrt(x * x + y * y)

		vh = (-y * vx + x * vy) / dist
		vv = (x * vx + y * vy) / dist

		# Gravity
		grav = planetGravParameter / dist / dist
		gx = -x * grav / dist
		gy = -y * grav / dist

		# Thrust
		vaccel = grav - vh * vh / dist
		vaccel -= k2 * (dist - mid2 - planetRadius)

		#print(dist, vh, vv, vaccel)

		craft.toggle_engines_by_vaccel(vaccel)

		accel = craft.total_thrust / craft.wet_mass
		bestvaccelclamped = max(min(vaccel, accel), -accel)
		besthaccel = math.sqrt(accel * accel - bestvaccelclamped * bestvaccelclamped)

		ax = (-y * besthaccel + x * bestvaccelclamped) / dist
		ay = (x * besthaccel + y * bestvaccelclamped) / dist

		# Integration
		if timesteps % timestepsPerLog == 0 and logging:
			print(f"T = {round(time, 1)}s:", end="\t")
			print(f"Alt: {round(dist - planetRadius, 2)}m", end="\t")
			print(f"Hvel: {round(vh)} m/s", end="\t")
			print(f"Vvel: {round(vv, 1)} m/s")
			#print(sum(craft.throttles.values()))

		vx += (ax + gx) * timestep
		vy += (ay + gy) * timestep
		x += vx * timestep
		y += vy * timestep
		time += timestep
		timesteps += 1

		craft.simulate_burn(timestep)

		# End conditions
		dist = math.sqrt(x * x + y * y)
		vh = (-y * vx + x * vy) / dist
		sma = 1 / (2 / dist - (vx * vx + vy * vy) / planetGravParameter)

		if (dist * vh) ** 2 > (sma * planetGravParameter):
			eccentricity = 0
		else:
			eccentricity = math.sqrt(1 - (dist * vh) ** 2 / (sma * planetGravParameter))
		
		periapsis = sma * (1 - eccentricity)

		if craft.wet_mass < 0:
			return 4, vh, dist - planetRadius, vv, landingdist - planetRadius, landingvv, time, landingtime
		if dist < planetRadius:
			return 5, vh, dist - planetRadius, vv, landingdist - planetRadius, landingvv, time, landingtime

	if logging:
		print(f"T = {round(time, 1)}s:", end="\t")
		print(f"Alt: {round(dist - planetRadius, 2)}m", end="\t")
		print(f"Hvel: {round(vh)} m/s", end="\t")
		print(f"Vvel: {round(vv, 1)} m/s")

	if logging:
		print("\nFuel used per type:")
	
		total_fuel_mass = 0
		for fuel_type, fuel_level in craft.fuel_used.items():
			mass = fuel_level * utils.fuel_masses[fuel_type]
			
			print(f"{fuel_type}\t{fuel_level}\tMass: {mass}")
			total_fuel_mass += mass
		
		print(f"Total fuel mass: {total_fuel_mass}")
		print(f"\nSuccessfully back in orbit, wet mass = {craft.wet_mass} kg\n")

	return 6, craft.wet_mass, dist - planetRadius, vv, landingdist - planetRadius, landingvv, time, landingtime

def simulate_landing_with_craft(craft, planet_name, altitude, logging = False, refuel = False):
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
		list_of_isp.append((isp, engine))
	list_of_isp.sort(reverse=True)
	next_engines_index = 1

	# Landing
	lithobrake_speed = 60

	elapsed_time = 0
	iteration_count = 0
	tick_time = 0.2
	ticks_per_log = 50
	throttle_resolution = 0.1 # might cause some sluggishness but I doubt it
	# Turn on most efficient engines (always best, I think)
	craft.set_throttle_by_name(list_of_isp[0][1], 1)

	n_groups = len(list_of_isp)

	landing_time = 0
	ascent_time = 0
	
	while (vx > lithobrake_speed):
		pitch_angle = utils.constant_altitude_thrust_angle(craft.total_thrust,
														   craft.wet_mass, vx,
														   altitude, preset=planet_name)
		thrust_effectiveness = math.cos(pitch_angle)
		current_effective_isp = craft.get_average_isp() * thrust_effectiveness
		# Check if it's more efficient to increase throttle of the next engine
		for next_engines_index in range(len(list_of_isp)):
			engine_name = list_of_isp[next_engines_index][1]
			current_throttle = craft.throttles[engine_name]
			new_throttle = max(0, current_throttle - throttle_resolution)
			craft.set_throttle_by_name(engine_name, new_throttle)
			test_angle = utils.constant_altitude_thrust_angle(craft.total_thrust,
																	craft.wet_mass, vx,
																	altitude, preset=planet_name)
			test_effective_isp = craft.get_average_isp() * math.cos(test_angle)
			if test_effective_isp > current_effective_isp:
				pitch_angle = test_angle
				thrust_effectiveness = math.cos(pitch_angle)
				if (logging):
					print(f"Setting {engine_name} throttle to {new_throttle}")
				if (new_throttle <= 0):
					next_engines_index -= 1
			
			else:
				craft.set_throttle_by_name(engine_name, current_throttle)

		for next_engines_index in range(len(list_of_isp)):
			engine_name = list_of_isp[next_engines_index][1]
			current_throttle = craft.throttles[engine_name]
			new_throttle = min(1, current_throttle + throttle_resolution)
			craft.set_throttle_by_name(engine_name, new_throttle)
			test_angle = utils.constant_altitude_thrust_angle(craft.total_thrust,
																	craft.wet_mass, vx,
																	altitude, preset=planet_name)
			test_effective_isp = craft.get_average_isp() * math.cos(test_angle)
			if test_effective_isp > current_effective_isp:
				pitch_angle = test_angle
				thrust_effectiveness = math.cos(pitch_angle)
				if (logging):
					print(f"Setting {engine_name} throttle to {new_throttle}")
				if (new_throttle <= 0):
					next_engines_index -= 1
			
			else:
				craft.set_throttle_by_name(engine_name, current_throttle)
	
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
			pitch_in_degrees = int(pitch_angle * 180 / math.pi)
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

	if refuel:
		craft.refuel()

	landing_time = iteration_count * tick_time
	iteration_count = 0
	if (logging):
		print("Liftoff!")
	# Ascent
	vx = 0

	# If the currently pointed engine is completely off (or doesn't exist), go back to considering the
	# lowest-isp active engine
	if next_engines_index == len(list_of_isp) or craft.throttles[list_of_isp[next_engines_index][1]] == 0:
		next_engines_index -= 1

	while (vx < orbital_velocity):
		pitch_angle = utils.constant_altitude_thrust_angle(craft.total_thrust,
																	 craft.wet_mass, vx,
																	 altitude, preset=planet_name)
		thrust_effectiveness = math.cos(pitch_angle)
		current_effective_isp = craft.get_average_isp() * thrust_effectiveness
		# Check if it's more efficient to decrease throttle of the next engine
		for next_engines_index in range(len(list_of_isp)):
			engine_name = list_of_isp[next_engines_index][1]
			current_throttle = craft.throttles[engine_name]
			new_throttle = max(0, current_throttle - throttle_resolution)
			craft.set_throttle_by_name(engine_name, new_throttle)
			test_angle = utils.constant_altitude_thrust_angle(craft.total_thrust,
																	craft.wet_mass, vx,
																	altitude, preset=planet_name)
			test_effective_isp = craft.get_average_isp() * math.cos(test_angle)
			if test_effective_isp > current_effective_isp:
				pitch_angle = test_angle
				thrust_effectiveness = math.cos(pitch_angle)
				if (logging):
					print(f"Setting {engine_name} throttle to {new_throttle}")
				if (new_throttle <= 0):
					next_engines_index -= 1
			
			else:
				craft.set_throttle_by_name(engine_name, current_throttle)

		for next_engines_index in range(len(list_of_isp)):
			engine_name = list_of_isp[next_engines_index][1]
			current_throttle = craft.throttles[engine_name]
			new_throttle = min(1, current_throttle + throttle_resolution)
			craft.set_throttle_by_name(engine_name, new_throttle)
			test_angle = utils.constant_altitude_thrust_angle(craft.total_thrust,
																	craft.wet_mass, vx,
																	altitude, preset=planet_name)
			test_effective_isp = craft.get_average_isp() * math.cos(test_angle)
			if test_effective_isp > current_effective_isp:
				pitch_angle = test_angle
				thrust_effectiveness = math.cos(pitch_angle)
				if (logging):
					print(f"Setting {engine_name} throttle to {new_throttle}")
				if (new_throttle <= 0):
					next_engines_index -= 1
			
			else:
				craft.set_throttle_by_name(engine_name, current_throttle)
			
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
			pitch_in_degrees = int(pitch_angle * 180 / math.pi)
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
		print(f"Final craft mass: {craft.wet_mass}")

	ascent_time = iteration_count * tick_time
	return landing_time, ascent_time

def optimize_spring_landing(craft, planet_name, altitude, logging = False, refuel = False):
	landing_time, ascent_time = simulate_landing_with_craft(deepcopy(craft), planet_name, altitude)

	# Initial guess for k parameters based on assumption that
	# vertical motion is harmonic
	k1 = (math.pi / landing_time) ** 2
	k2 = (math.pi / ascent_time) ** 2

	mid1 = 0.5 * altitude
	mid2 = 0.5 * altitude

	starting_altitude = 10
	final_altitude = 50

	# Iteratively improve guesses for k1, mid1, k2, and mid2
	# The convergence is very rapid so 10 iterations is more than enough
	# Feel free to reduce this number to save computation time
	iterations = 0
	finished = False
	while iterations < 100 and not finished:
		code, score, alt, vv, landalt, landvv, time2, time1 = simulate_spring_landing_with_craft(
			deepcopy(craft), 
			planet_name, 
			altitude, 
			k1, 
			k2, 
			mid1, 
			mid2, 
			starting_altitude,
			False,
			refuel
		)

		accelmaxalt1 = (altitude - mid1) * k1
		timediff = landvv / accelmaxalt1
		maxalt = math.sqrt((landalt - mid1) ** 2 + landvv ** 2 / k1) + mid1

		#print(maxalt)

		mid1 += max(-100, min(100, 0.5 * (altitude - maxalt)))
		k1 = min((math.pi / (math.pi / math.sqrt(k1) - timediff)) ** 2, 2 * k1)

		k1finalest = (math.pi / time1) ** 2

		if math.sqrt(k1 / k1finalest) > 1.5:
			k1 = k1finalest

		if code > 3:
			accelminalt2 = (final_altitude - mid2) * k2
			timediff = vv / accelminalt2
			minalt = -math.sqrt((alt - mid2) ** 2 + vv ** 2 / k2) + mid2

			#print(minalt)

			mid2 += max(-100, min(100, 0.5 * (final_altitude - minalt)))
			k2 = min((math.pi / (math.pi / math.sqrt(k2) - timediff)) ** 2, 2 * k2)

			k2finalest = (math.pi / time2) ** 2

			if math.sqrt(k2 / k2finalest) > 1.5:
				k2 = k2finalest

		#print(code, score, alt, vv, landalt, landvv)
		#print(k1, k2, mid1, mid2)

		finished = (code == 6) and (abs(vv) < 1) and (abs(landvv) < 1) and (abs(alt - final_altitude) < 1)
		iterations += 1

	simulate_spring_landing_with_craft(
		craft, 
		planet_name, 
		altitude, 
		k1, 
		k2, 
		mid1, 
		mid2, 
		starting_altitude,
		logging,
		refuel
	)

def payload_capacity(craft):
	# Mass - total engine mass - total (used) fuel tank mass
	mass = craft.wet_mass
	for engine, count in craft.engines.items():
		mass -= utils.engine_info[engine][4] * count
	for fuel_type, amount in craft.fuel_used.items():
		fuel_mass = utils.fuel_masses[fuel_type] * amount
		fuel_tank_mass = float(fuel_mass) / utils.fuel_tank_info[fuel_type][1]
		mass -= fuel_tank_mass
	return mass

"""
flight_profile is a list of landings and vacuum burns.
Landings are of the format 'name altitude' and vacuum burns are of the format 'dv'

"""
def simulate_flight_profile(craft, flight_profile, log=False, refuel = False):
	for item in flight_profile:
		if " " in item:
			# This is a landing
			name = item.split(" ")[0]
			altitude = float(item.split(" ")[1])
			if name.lower() in bodies:
				#simulate_landing_with_craft(craft, name.lower(), altitude, logging=log, refuel = refuel)
				optimize_spring_landing(craft, name.lower(), altitude, logging=log, refuel = refuel)
			elif name == "LiquidFuel" or name == "LFOx" or name == "Monopropellant" or name == "XenonGas":
				craft.use_fuel(name, altitude * craft.wet_mass)
			else:
				print(name)
		else:
			# This is a low twr vacuum burn
			dv = int(item)
			craft.use_dv(dv)

#planet_name = "tylo"
#altitude = 9100

mass_in_orbit = 40650
engines = []

engine_count = {'Rapier': 1.0, 'Nerv': 2.082927573786803, 'Dawn': 1.0002454085822676, 'Rhino': 0.05251314180328612, 'Dart': 1.194070695619644e-05}

min_engine_count = {
	"Nerv": 1,
	"Rapier": 1,
	"Dawn": 1,
}

def optimize_engine_layout(engine_count = dict(), min_engine_count = dict(), flight_profile = [], output_ascent_profiles = True, allow_fractional = False, allow_new_engines = True, mass_in_orbit = 0, refuel = False):
	improved_last_iteration = True
	best_engine_count = engine_count.copy()
	
	# Base case
	start_time = time.time()
	simulation_time = 0
	craft = single_stage.SingleStageCraft()
	fuel_levels = dict()
	craft.manual_instantiate(engine_count, fuel_levels, mass_in_orbit)
	simulate_flight_profile(craft, flight_profile, output_ascent_profiles, refuel)
	max_payload_remaining = (payload_capacity(craft) - 6000)/best_engine_count["Dawn"]
	print(f"Reference payload capacity: {max_payload_remaining}")
	print(f"Reference payload fraction: {max_payload_remaining / mass_in_orbit}")
	
	"""
	Simulated annealing.
	Temperature controls the average magnitude of proposed changes to the engine count.
	This will cool down over time until it's equal to one, at which point we do a
	deterministic greedy search.
	Hopefully this can hit a bunch of different local maxima.
	"""
	temperature = 10.0
	cooling_factor = 0.95
	n_iterations = 0
	anneal_threshold = 1
	if (allow_fractional):
		temperature = 0.1
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

			simulate_flight_profile(test_craft, flight_profile, False, refuel)

			sim_end = time.time()
			simulation_time += sim_end - sim_start
			test_payload = payload_capacity(test_craft)

			test_payload -= 6000
			test_payload /= test_engine_count["Dawn"]

			if test_payload > max_payload_remaining:
				print(f"Added {add_value}x {engine}")
				improved_last_iteration = True
				max_payload_remaining = test_payload
				print(f"Payload capacity: {test_payload}")
				print(f"Payload fraction: {test_payload / mass_in_orbit}")
				print(test_engine_count)
				#print(test_engine_count)
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
			test_craft2 = single_stage.SingleStageCraft()
			test_craft2.manual_instantiate(test_engine_count, fuel_levels, mass_in_orbit)
			sim_start = time.time()

			simulate_flight_profile(test_craft2, flight_profile, False, refuel)

			sim_end = time.time()
			simulation_time += sim_end - sim_start
			test_payload = payload_capacity(test_craft2)

			test_payload -= 6000
			test_payload /= test_engine_count["Dawn"]

			if test_payload > max_payload_remaining:
				print(f"Removed {subtract_value}x {engine}")
				improved_last_iteration = True
				max_payload_remaining = test_payload
				print(f"Payload capacity: {test_payload}")
				print(f"Payload fraction: {test_payload / mass_in_orbit}")
				print(test_engine_count)
				#print(test_engine_count)
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
	
	simulate_flight_profile(craft, flight_profile, output_ascent_profiles, refuel)

	print(f"Total time: {end_time - start_time}")
	print(f"Simulation time: {simulation_time}")

	return best_engine_count

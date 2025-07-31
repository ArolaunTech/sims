from utils import *
import utils #Eh
import math

class SingleStageCraft:
	def __init__(self):
		# List of engines. Map name to count
		self.engines = {}
		# List of throttle values (0 to 1) for each engine
		self.throttles = {}
		
		# Track fuel consumption per sec per type.
		# We're going to assume vacuum landings for now.
		self.engine_activations = dict()
		self.fuel_consumption = {
			"LiquidFuel": 0,
			"LFOx": 0,
			"XenonGas": 0,
			"Monopropellant": 0,
			"SolidFuel": 0,
		}
		
		# Map of fueltype name to fuel (in "Units")
		self.fuel_levels = {
			"LiquidFuel": 0,
			"LFOx": 0,
			"XenonGas": 0,
			"Monopropellant": 0,
			"SolidFuel": 0,
		}
		self.fuel_used = self.fuel_levels.copy()
		
		#self.dry_mass = 0
		
		self.wet_mass = 0
	
		self.total_thrust = 0
		self.mass_flow = 0
		self.mass_tank_flow = 0
	# Arguments:
	# List of tuples of (isp, fueltype, thrust)
	# Dictionary of fuel type and fuel level
	# Number for dry mass
	def manual_instantiate(self, _engines, _fuel_levels, _wet_mass, respect_tank_mass=False):
		self.engines = _engines.copy()
		self.wet_mass = _wet_mass
		for engine in self.engines:
			self.throttles[engine] = 0
		self.total_thrust = 0

	def instantiate_from_craft_file(self, craft_file_path):
		# Not implemented yet!
		print("WIP")

	def set_throttle_by_name(self, name, value):
		# Subtract current fuel consumption, thrust, and mass flow
		engine_info = utils.engine_info[name]
		throttle = self.throttles[name]
		fuel_type = engine_info[1]
		count = self.engines[name]
		fuel_drain = engine_info[3]
		thrust = engine_info[2]

		self.total_thrust -= thrust * count * throttle
		self.fuel_consumption[fuel_type] -= fuel_drain * count * throttle
		self.mass_flow -= fuel_drain * count * throttle * utils.fuel_masses[fuel_type]
		self.mass_tank_flow -= fuel_drain * count * throttle * utils.fuel_masses[fuel_type] * (1 + 1 / utils.fuel_tank_info[fuel_type][1])
		# Set new throttle
		throttle = value
		self.throttles[name] = throttle
		self.total_thrust += thrust * count * throttle
		self.fuel_consumption[fuel_type] += fuel_drain * count * throttle
		self.mass_flow += fuel_drain * count * throttle * utils.fuel_masses[fuel_type]
		self.mass_tank_flow += fuel_drain * count * throttle * utils.fuel_masses[fuel_type] * (1 + 1 / utils.fuel_tank_info[fuel_type][1])

	def get_throttle_by_name(self, name):
		return self.throttles[name]

	def toggle_engines_by_isp(self, target_isp):
		for name in self.engines:
			engine_info = utils.engine_info[name]
			if engine_info[0] == target_isp:
				engineOn = self.throttles[name] > 0
				if not engineOn:
					self.set_throttle_by_name(name, 1)
				else:
					self.set_throttle_by_name(name, 0)

	def toggle_engines_by_total_thrust(self, total_thrust):
		# Build Isp list
		ispList = []
		for engine, count in self.engines.items():
			if count == 0:
				continue

			fueltype = utils.engine_info[engine][1]
			fuelratio = utils.fuel_tank_info[fueltype][1]

			# Include fuel tanks in isp to maximize payload fraction,
			# not just fuel usage
			isp = utils.engine_info[engine][0] / (1 + 1 / fuelratio)
			ispList.append((isp, engine))
		ispList.sort(reverse=True)

		# Clamp total_thrust
		total_thrust_adjusted = total_thrust
		if total_thrust_adjusted < 0:
			total_thrust_adjusted = 0
		if total_thrust_adjusted > len(ispList):
			total_thrust_adjusted = len(ispList)

		# Set throttles
		self.turn_engines_off()
		for i in range(int(total_thrust_adjusted)):
			self.set_throttle_by_name(ispList[i][1], 1)

		if int(total_thrust_adjusted) < len(ispList):
			self.set_throttle_by_name(
				ispList[int(total_thrust_adjusted)][1],
				total_thrust_adjusted % 1
			)

	def toggle_engines_by_vaccel(self, vaccel):
		# Build Isp list
		ispList = []
		for engine, count in self.engines.items():
			if count == 0:
				continue

			fueltype = utils.engine_info[engine][1]
			fuelratio = utils.fuel_tank_info[fueltype][1]

			# Include fuel tanks in isp to maximize payload fraction,
			# not just fuel usage
			isp = utils.engine_info[engine][0] / (1 + 1 / fuelratio)
			ispList.append((isp, engine))
		ispList.sort(reverse=True)

		# Calculate best throttle
		for name in self.engines:
			self.set_throttle_by_name(name, 1)

		# All engines on
		bestthrottle = len(ispList)

		accel = self.total_thrust / self.wet_mass
		#print("accel", accel, vaccel)
		if accel < abs(vaccel):
			return

		haccel = math.sqrt(accel * accel - vaccel * vaccel)
		#print("haccel", haccel)
		bestscore = haccel / accel * self.get_average_isp_tank()

		#print(bestscore)

		# Some engines off
		numgroups = len(ispList)
		for i in range(numgroups - 1, -1, -1):
			engine_name = ispList[i][1]
			self.set_throttle_by_name(engine_name, 0)

			accel = self.total_thrust / self.wet_mass

			#print(accel, vaccel)

			if accel >= abs(vaccel):
				haccel = math.sqrt(accel * accel - vaccel * vaccel)

				score = haccel / accel * self.get_average_isp_tank()

				if score > bestscore:
					bestscore = score
					bestthrottle = i

			t1 = self.total_thrust
			f1 = self.mass_tank_flow

			t0 = self.engines[engine_name] * utils.engine_info[engine_name][2]
			f0 = t0 / ispList[i][0] / 9.80665

			t0 /= self.wet_mass
			t1 /= self.wet_mass

			optimalthrottle = (f0 * (t1 * t1 - vaccel * vaccel) - t0 * t1 * f1) / (t0 * (t0 * f1 - t1 * f0))
			if optimalthrottle < 0:
				continue
			if optimalthrottle > 1:
				continue

			self.set_throttle_by_name(engine_name, optimalthrottle)

			accel = self.total_thrust / self.wet_mass
			if accel >= abs(vaccel):
				haccel = math.sqrt(accel * accel - vaccel * vaccel)

				score = haccel / accel * self.get_average_isp_tank()

				if score > bestscore:
					bestscore = score
					bestthrottle = i + optimalthrottle

			self.set_throttle_by_name(engine_name, 0)
			accel = self.total_thrust / self.wet_mass

			if accel < abs(vaccel):
				break

		self.toggle_engines_by_total_thrust(bestthrottle)

	def turn_engines_off(self):
		for name in self.engines:
			self.set_throttle_by_name(name, 0)
	
	def simulate_burn(self, delta_time):
		# just reduce fuel levels and wet mass
		# Assumes infinite of each fuel type (for optimizing relative amounts).
		for fuel_type, consumption_rate in self.fuel_consumption.items():
			fuel_used_tick = consumption_rate * delta_time
			self.fuel_used[fuel_type] += fuel_used_tick
		self.wet_mass -= delta_time * self.mass_flow

	def get_average_isp(self):
		if self.mass_flow <= 0:
			return 0
		# Thrust in newtons divided by mass flow rate in kg/s
		# Technically this is "effective exhaust velocity" but whatever
		return self.total_thrust / self.mass_flow

	def get_average_isp_tank(self):
		if self.mass_tank_flow <= 0:
			return 0
		return self.total_thrust / self.mass_tank_flow

	def use_dv(self, dv):
		# Decreases fuel levels by minimum possible to use given dV
		# Useful for simulating orbital transfers
		best_isp = 0
		best_fuel = ""
		for engine in self.engines:
			engine_info = utils.engine_info[engine]
			engine_isp = engine_info[0]
			engine_fueltype = engine_info[1]

			if engine_isp > best_isp:
				best_isp = engine_isp
				best_fuel = engine_fueltype

		fuel_ratio = math.exp(-dv/(best_isp * 9.81))
		new_mass = self.wet_mass * fuel_ratio

		#print(self.wet_mass, new_mass, (self.wet_mass - new_mass)/fuel_masses[best_fuel], self.fuel_consumption)

		self.fuel_used[best_fuel] += (self.wet_mass - new_mass)/fuel_masses[best_fuel]
		self.wet_mass = new_mass

	def use_fuel(self, fueltype, fuel):
		self.fuel_used[fueltype] += fuel
		self.wet_mass -= fuel

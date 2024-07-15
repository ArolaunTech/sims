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
    # Arguments:
    # List of tuples of (isp, fueltype, thrust)
    # Dictionary of fuel type and fuel level
    # Number for dry mass
    def manual_instantiate(self, _engines, _fuel_levels, _wet_mass, respect_tank_mass=False):
        self.engines = _engines.copy()
        self.wet_mass = _wet_mass
        for engine in self.engines:
            self.throttles[engine] = 0

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
        # Set new throttle
        throttle = value
        self.throttles[name] = throttle
        self.total_thrust += thrust * count * throttle
        self.fuel_consumption[fuel_type] += fuel_drain * count * throttle
        self.mass_flow += fuel_drain * count * throttle * utils.fuel_masses[fuel_type]
        
    def toggle_engines_by_isp(self, target_isp):
        for name in self.engines:
            engine_info = utils.engine_info[name]
            if engine_info[0] == target_isp:
                engineOn = self.throttles[name] > 0
                if not engineOn:
                    self.set_throttle_by_name(name, 1)
                else:
                    self.set_throttle_by_name(name, 0)

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
        if(self.mass_flow == 0):
            return 0
        # Thrust in newtons divided by mass flow rate in kg/s
        return self.total_thrust / self.mass_flow

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

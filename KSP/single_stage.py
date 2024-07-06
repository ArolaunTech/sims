import utils

class SingleStageCraft:
    def __init__(self):
        # List of engines. String name
        self.engines = []
        
        # Track whether each engine is on or off + fuel consumption per sec per type.
        # We're going to assume vacuum landings for now.
        self.engine_activations = []
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
        
        self.dry_mass = 0
        
        self.wet_mass = 0
    
        self.total_thrust = 0

    # Arguments:
    # List of tuples of (isp, fueltype, thrust)
    # Dictionary of fuel type and fuel level
    # Number for dry mass
    def manual_instantiate(self, _engines, _fuel_levels, _wet_mass, respect_tank_mass=False):
        self.engines = _engines.copy()
        #self.fuel_levels = _fuel_levels
        self.wet_mass = _wet_mass
        for engine in self.engines:
            self.engine_activations.append(False)

    def instantiate_from_craft_file(self, craft_file_path):
        # Not implemented yet!
        print("WIP")

    def toggle_engines_by_isp(self, target_isp):
        for i in range(len(self.engines)):
            engine_name = self.engines[i]
            engine_info = utils.engine_info[engine_name]
            if engine_info[0] == target_isp:
                engineOn = self.engine_activations[i]
                if not engineOn:
                    # Turn engine on
                    self.engine_activations[i] = True
                    
                    self.total_thrust += engine_info[2]
                    self.fuel_consumption[engine_info[1]] += engine_info[3]
                else:
                    # Turn engine off
                    self.engine_activations[i] = False
                    self.total_thrust -= engine_info[2]
                    self.fuel_consumption[engine_info[1]] -= engine_info[3]

    def simulate_burn(self, delta_time):
        # just reduce fuel levels and wet mass
        # Assumes infinite of each fuel type (for optimizing relative amounts).
        for fuel_type, consumption_rate in self.fuel_consumption.items():
            fuel_used_tick = consumption_rate * delta_time
            self.fuel_used[fuel_type] += fuel_used_tick
            self.wet_mass -= fuel_used_tick * utils.fuel_masses[fuel_type]

    def get_average_isp(self):
        # Thrust in newtons divided by mass flow rate in kg/s
        mass_flow = 0
        for fuel_type, consumption_rate in self.fuel_consumption.items():
            mass_flow += consumption_rate * utils.fuel_masses[fuel_type]
        return self.total_thrust / mass_flow

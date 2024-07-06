import utils

class SingleStageCraft:
    def __init__(self):
        # List of engines. Map name to count
        self.engines = {}
        
        # Track whether each engine group is on or off + fuel consumption per sec per type.
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
            self.engine_activations[engine] = False

    def instantiate_from_craft_file(self, craft_file_path):
        # Not implemented yet!
        print("WIP")

    def toggle_engines_by_isp(self, target_isp):
        for name, count in self.engines.items():
            engine_info = utils.engine_info[name]
            if engine_info[0] == target_isp:
                engineOn = self.engine_activations[name]
                fuel_type = engine_info[1]
                fuel_added = engine_info[3] * count
                if not engineOn:
                    # Turn engine on
                    self.engine_activations[name] = True
                    self.total_thrust += engine_info[2] * count
                    self.fuel_consumption[fuel_type] += fuel_added
                    self.mass_flow += fuel_added * utils.fuel_masses[fuel_type]
                else:
                    # Turn engine off
                    self.engine_activations[name] = False
                    self.total_thrust -= engine_info[2] * count
                    self.mass_flow -= fuel_added * utils.fuel_masses[fuel_type]
                    self.fuel_consumption[fuel_type] -= fuel_added

    def simulate_burn(self, delta_time):
        # just reduce fuel levels and wet mass
        # Assumes infinite of each fuel type (for optimizing relative amounts).
        for fuel_type, consumption_rate in self.fuel_consumption.items():
            fuel_used_tick = consumption_rate * delta_time
            self.fuel_used[fuel_type] += fuel_used_tick
        self.wet_mass -= delta_time * self.mass_flow

    def get_average_isp(self):
        # Thrust in newtons divided by mass flow rate in kg/s
        return self.total_thrust / self.mass_flow

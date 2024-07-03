class SingleStageCraft:
    # List of engines. Tuple of (isp, fueltype, thrust, fuel/sec)
    engines = []

    # Track whether each engine is on or off + fuel consumption per sec per type.
    # We're going to assume vacuum landings for now.
    engine_activations = []
    fuel_consumption = {
        "LiquidFuel": 0,
        "LFOx": 0,
        "XenonGas": 0,
        "Monopropellant": 0,
        "SolidFuel": 0,
    }

    # Map of fueltype name to fuel (in "Units")
    fuel_levels = {
        "LiquidFuel": 0,
        "LFOx": 0,
        "XenonGas": 0,
        "Monopropellant": 0,
        "SolidFuel": 0,
    }
    original_fuel_levels = fuel_levels.copy()

    # Map of fueltype to mass in units
    fuel_masses = {
        "LiquidFuel": 5,
        "LFOx": 5,
        "Monopropellant": 4,
        "XenonGas": 0.1,
        "SolidFuel": 7.5,
    }

    # Sum of engine mass and fuel tank mass.
    dry_mass = 0

    wet_mass = 0

    total_thrust = 0
    # Arguments:
    # List of tuples of (isp, fueltype, thrust)
    # Dictionary of fuel type and fuel level
    # Number for dry mass
    def manual_instantiate(self, _engines, _fuel_levels, _wet_mass, respect_tank_mass=False):
        self.engines = _engines
        #self.original_fuel_levels = _fuel_levels
        #self.fuel_levels = _fuel_levels
        self.wet_mass = _wet_mass
        #for fuel_type, fuel_count in _fuel_levels.items():
        #    self.wet_mass += fuel_count * self.fuel_masses(fuel_type)
        for engine in self.engines:
            self.engine_activations.append(False)

    def instantiate_from_craft_file(self, craft_file_path):
        # Not implemented yet!
        print("WIP")

    def toggle_engines_by_isp(self, target_isp):
        for i in range(len(self.engines)):
            if self.engines[i][0] == target_isp:
                engineOn = self.engine_activations[i]
                if not engineOn:
                    # Turn engine on
                    self.engine_activations[i] = True
                    self.total_thrust += self.engines[i][2]
                    self.fuel_consumption[self.engines[i][1]] += self.engines[i][3]
                else:
                    # Turn engine off
                    self.engine_activations[i] = False
                    self.total_thrust -= self.engines[i][2]
                    self.fuel_consumption[self.engines[i][1]] -= self.engines[i][3]

    def simulate_burn(self, delta_time):
        # just reduce fuel levels and wet mass/
        # Assumes infinite of each fuel type (for optimizing relative amounts).
        for fuel_type, consumption_rate in self.fuel_consumption.items():
            fuel_used = consumption_rate * delta_time
            self.fuel_levels[fuel_type] -= fuel_used
            self.wet_mass -= fuel_used * self.fuel_masses[fuel_type]

    def get_average_isp(self):
        # Thrust in newtons divided by mass flow rate in kg/s
        mass_flow = 0
        for fuel_type, consumption_rate in self.fuel_consumption.items():
            mass_flow += consumption_rate * self.fuel_masses[fuel_type]
        return self.total_thrust / mass_flow

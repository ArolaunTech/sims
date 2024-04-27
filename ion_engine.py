#Imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colormaps as cm

pilot_mass = 0.07
engine_mass = 0.25
seat_mass = 0.05
reaction_mass = 0.05

dry_mass = pilot_mass + engine_mass + seat_mass + reaction_mass

G = 6.67408e-11

moon_dv = 1676 + 900 + 200
moon_return_dv = 200

moon_g = 0.166
moon_radius = 200000
moon_mass = 9.81*moon_g*moon_radius*moon_radius/G
print(moon_mass)

orbital_altitude = 10000
moon_orbital_velocity = np.sqrt(G * moon_mass / (moon_radius + orbital_altitude))
print(moon_orbital_velocity)

landing_altitude = 5000
second_burn = moon_orbital_velocity - np.sqrt(G * moon_mass * (2/(moon_radius + orbital_altitude) - 2/(2 * moon_radius + landing_altitude + orbital_altitude)))
moon_dv += second_burn
moon_return_dv += second_burn

first_burn = np.sqrt(G * moon_mass * (2/(moon_radius + landing_altitude) - 2/(2 * moon_radius + landing_altitude + orbital_altitude)))
moon_rotational_period = 138984.38
rotational_velocity = 2 * np.pi * (moon_radius + landing_altitude)/moon_rotational_period
first_burn -= rotational_velocity
surface_orbital_velocity = np.sqrt(G * moon_mass/(moon_radius + landing_altitude))
print(first_burn, surface_orbital_velocity)

#Functions
def calcXenonCapacity(targetXenon):
	if targetXenon <= 0:
		return 0
	out = 1000000000
	for container in [405, 720, 5700]:
		out = min(out, container + calcXenonCapacity(targetXenon - container))
	return out

def calcXenonDryMass(targetXenon):
	return calcXenonCapacity(targetXenon)/30000

def calcBatteryMass(batteryCapacity):
	actualBatteryCapacity = np.ceil(batteryCapacity/100)
	return 0.005*actualBatteryCapacity

def calcSolarMass(solarGen):
	if solarGen <= 0:
		return 0
	availablePanels = [ #Power/mass pairs
		[0.35, 0.005], #OX-STAT
		[1.64, 0.0175], #OX-4L 1x6 Photovoltaic Panels
		[2.8, 0.04], #OX-STAT-XL
		[8.25, 0.09] #OX-10L 1x5 Photovoltaic Panels
	]
	out = 1000000000
	for panel in availablePanels:
		out = min(out, panel[1] + calcSolarMass(solarGen - panel[0]))
	return out

def calcMass(batteryCapacity, solarGen, xenon):
	return dry_mass + calcXenonDryMass(xenon) + calcBatteryMass(batteryCapacity) + calcSolarMass(solarGen) + xenon/10000

def calcPlausibility(batteryCapacity, solarGen, xenon):
	wetMass = calcMass(batteryCapacity, solarGen, xenon)
	dryMass = wetMass - xenon/10000
	#print(wetMass - pilot_mass)

	dV = 4200 * 9.81 * np.log(wetMass/dryMass)

	if dV < moon_dv:
		return False
	#print(dV)
	moonOrbitdV = dV - moon_dv
	moonOrbitMass = dryMass * np.exp(moonOrbitdV/(4200*9.81))

	#Moon landing and ascent
	currVel = rotational_velocity + first_burn

	currentEC = batteryCapacity
	timestep = 1

	while currVel > rotational_velocity:
		orb_frac = currVel/surface_orbital_velocity
		verticalAcceleration = 9.81 * moon_g * (orb_frac * orb_frac - 1)
		totalAcceleration = 2/moonOrbitMass

		#print(verticalAcceleration, totalAcceleration, orb_frac)
		if currentEC < 0:
			totalAcceleration *= 0.5 + 0.5*solarGen/8.74
		if totalAcceleration < verticalAcceleration:
			return False
		if moonOrbitMass < dryMass:
			return False
		horizontalAcceleration = np.sqrt(totalAcceleration * totalAcceleration - verticalAcceleration*verticalAcceleration)

		currVel -= horizontalAcceleration * timestep
		currentEC -= 8.74 * timestep
		currentEC += solarGen * timestep
		moonOrbitMass -= 0.486*0.0001
		#print(currVel, currentEC, moonOrbitMass)

	currentEC = batteryCapacity
	currVel = rotational_velocity
	while currVel < rotational_velocity + first_burn:
		orb_frac = currVel/surface_orbital_velocity
		verticalAcceleration = 9.81 * moon_g * (orb_frac * orb_frac - 1)
		totalAcceleration = 2/moonOrbitMass

		#print(verticalAcceleration, totalAcceleration, orb_frac)
		if currentEC < 0:
			totalAcceleration *= 0.5 + 0.5*solarGen/8.74
		if totalAcceleration < verticalAcceleration:
			return False
		if moonOrbitMass < dryMass:
			return False
		horizontalAcceleration = np.sqrt(totalAcceleration * totalAcceleration - verticalAcceleration*verticalAcceleration)

		currVel += horizontalAcceleration * timestep
		currentEC -= 8.74 * timestep
		currentEC += solarGen * timestep
		moonOrbitMass -= 0.486*0.0001

		#print(currVel, currentEC, moonOrbitMass)
	#print(currVel, currentEC, moonOrbitMass - pilot_mass)
	dV = 4200 * 9.81 * np.log(moonOrbitMass/dryMass)
	#print(dV)
	return dV > moon_return_dv


def calc_battery_capacity(xenon, solarGen):
	for optbattery in range(100, 2000, 100):
		if calcPlausibility(optbattery, solarGen, xenon):
			return optbattery
	return -1

#Optimization
#print(calcMass(1600, 0.34, 490))
xenonTests = np.linspace(400, 700, num=100)
solarTests = np.linspace(3.28, 5, num=100)
xenonTests, solarTests = np.meshgrid(xenonTests, solarTests)
print(xenonTests[0], solarTests)
massTests = np.zeros(xenonTests.shape)

#massTests[0][0] = np.nan
minMass = 1000
for r in range(len(massTests)):
	for c in range(len(massTests[0])):
		minimumBatteries = calc_battery_capacity(xenonTests[r][c], solarTests[r][c])
		if minimumBatteries == -1:
			massTests[r][c] = np.nan
		else:
			massTests[r][c] = -pilot_mass + calcMass(minimumBatteries, solarTests[r][c], xenonTests[r][c])
			if massTests[r][c] < minMass:
				minMass = massTests[r][c]
				print(minMass, minimumBatteries, solarTests[r][c], xenonTests[r][c])
	print(r)

plottingcolormap = plt.get_cmap('inferno').copy()
plottingcolormap.set_bad('white')

plt.imshow(massTests, cmap=plottingcolormap, extent=[400, 700, 5, 3.28], aspect=320/7, interpolation=None)
plt.xlabel('Xenon in ion stage (units)')
plt.ylabel('Solar generation (EC/s)')
plt.colorbar()
plt.show()
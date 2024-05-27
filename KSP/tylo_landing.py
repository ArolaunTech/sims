import numpy as np

g = 9.81
G = 6.67408e-11

moon_radius = 600000
moon_g = 0.8
moon_mass = g*moon_g*moon_radius*moon_radius/G

print(moon_mass, g*moon_g)

orbital_height = 15000
landing_height = 0

rotational_period = 211926.36

orbital_vel = np.sqrt(G*moon_mass/(moon_radius + orbital_height))
surface_orbital_vel = np.sqrt(G*moon_mass/(moon_radius + landing_height))
pre_orbital_vel = np.sqrt(G*moon_mass*(2/(moon_radius + orbital_height) - 2/(moon_radius*2 + orbital_height + landing_height)))
pre_landing_vel = np.sqrt(G*moon_mass*(2/(moon_radius + landing_height) - 2/(moon_radius*2 + orbital_height + landing_height)))
rotational_vel = 2*np.pi*moon_radius/(rotational_period)

print(orbital_vel, pre_orbital_vel, pre_landing_vel, rotational_vel)

availableTanks = np.array([ #Fuel, dry
	[0.11, 0.01375], #Dumpling
	[0.2, 0.025], #Oscar-B
	[0.3, 0.0375], #Doughnut
	[0.27, 0.03375] #Baguette
])

def calcMinMass(fuel_mass):
	if fuel_mass <= 0:
		return 0
	out = 1000000
	for tank in availableTanks:
		out = min(out, calcMinMass(fuel_mass - tank[0]) + tank[1])
	return out

print(calcMinMass(0.1))

availableEngines = np.array([ #Mass, Thrust, Isp, Gimbal
	[0.02, 2, 290, True], #Spider
	[0.02, 2, 315, False], #Ant
	[0.08, 16, 290, True], #Twitch
	[0.13, 20, 320, True], #Spark
])

#Set engine[2] to fuel consumption
for i in range(len(availableEngines)):
	availableEngines[i][2] = availableEngines[i][1]/(availableEngines[i][2] * g)

pilot_mass = 0.09
seat_mass = 0.05
core_mass = 0.001
couple_mass = 0.01

working_crafts = []
for opt in range(1000000):
	stage_number = np.random.randint(1, 4)

	if opt%1000 == 0:
		print(opt, len(working_crafts))
	stages = []
	for i in range(stage_number):
		fuel_mass = np.random.random() * 1
		stages.append([[fuel_mass, calcMinMass(fuel_mass)],availableEngines[np.random.choice(len(availableEngines), size=np.random.randint(3))], 0, 0, 0, 0]) #Tanks, engines, wet mass, dry mass, thrust, fuel consumption
	thrust = 0
	for i in range(stage_number):
		for engine in stages[i][1]:
			thrust += engine[1]
		stages[i][4] = thrust

	mass = pilot_mass + seat_mass
	for i in range(stage_number):
		mass += stages[i][0][0] + stages[i][0][1]
		for engine in stages[i][1]:
			mass += engine[0]
		mass += core_mass

		#print(mass)
		stages[i][2] = mass

		mass += couple_mass

	fuelConsump = 0
	for i in range(stage_number):
		for engine in stages[i][1]:
			fuelConsump += engine[2]
		stages[i][5] = fuelConsump

	for i in range(stage_number):
		stages[i][3] = stages[i][2]
		stages[i][3] -= stages[i][0][0]

	if stages[-1][2] > 1.5: #Max mass
		continue

	gimbal = False #Need control
	for engine in stages[0][1]:
		if engine[3]:
			gimbal = True
			break
	if not gimbal:
		continue

	#Simulate stage
	currVel = pre_landing_vel
	currStage = stage_number - 1
	currMass = stages[currStage][2]
	currThrust = stages[currStage][4]
	currFuelConsumption = stages[currStage][5]

	viable = True
	timestep = 1

	while currVel > rotational_vel:
		acceleration = currThrust/currMass
		orbfrac = currVel/surface_orbital_vel
		verticalAcceleration = g*moon_g*(1 - orbfrac*orbfrac)
		#print(verticalAcceleration, acceleration)
		if np.abs(verticalAcceleration) > acceleration:
			viable = False
			break
		horizontalAcceleration = np.sqrt(acceleration*acceleration - verticalAcceleration*verticalAcceleration)

		currVel -= horizontalAcceleration * timestep
		currMass -= currFuelConsumption * timestep
		if currMass < stages[currStage][3]:
			#stage
			currStage -= 1
			if currStage == -1:
				viable = False
				break
			currMass = stages[currStage][2]
			currThrust = stages[currStage][4]
			currFuelConsumption = stages[currStage][5]
		#print(currVel)

	#currMass -= 0.04

	if viable:
		currVel = rotational_vel
		targetVelocity = 1600
		while currVel < targetVelocity:
			acceleration = currThrust/currMass
			orbfrac = currVel/surface_orbital_vel
			verticalAcceleration = g*moon_g*(1 - orbfrac*orbfrac)
			#print(verticalAcceleration, acceleration)
			if np.abs(verticalAcceleration) > acceleration:
				viable = False
				break
			horizontalAcceleration = np.sqrt(acceleration*acceleration - verticalAcceleration*verticalAcceleration)

			currVel += horizontalAcceleration * timestep
			currMass -= currFuelConsumption * timestep
			if currMass < stages[currStage][3]:
				#stage
				currStage -= 1
				if currStage == -1:
					viable = False
					break
				currMass = stages[currStage][2]
				currThrust = stages[currStage][4]
				currFuelConsumption = stages[currStage][5]
			#print(currVel)

	if not viable:
		continue

	working_crafts.append(stages)

working_crafts.sort(reverse=True, key=lambda x: x[-1][2])
for craft in working_crafts:
	for stage in craft:
		print(stage)
	print('\n')
print(len(working_crafts))
#New Moho SSTO optimizer, goal is less than 6,080 kg
import numpy as np
from copy import deepcopy as dc

#Constants
g = 9.81
ssto_lf_fraction = 0.07963 #Empirical from Brad SSTO mission

#Define Moho
moho_radius = 250000
moho_g = 0.275 * g
moho_stdgp = moho_g * moho_radius**2
moho_rotp = 1.21e+6

landing_alt = 5000
orb_alt = 10000

moho_rv = 2*np.pi*(moho_radius + landing_alt)/moho_rotp
moho_orbvel = np.sqrt(moho_stdgp/(moho_radius + landing_alt))
first_burn = np.sqrt(moho_stdgp * (2/(moho_radius + landing_alt) - 2/(2 * moho_radius + landing_alt + orb_alt)))
second_burn = np.sqrt(moho_stdgp * (2/(moho_radius + orb_alt) - 2/(2*moho_radius + landing_alt + orb_alt)))
second_burn = np.sqrt(moho_stdgp/(moho_radius + orb_alt)) - second_burn

print('Rotational Velocity of Moho: ' + str(moho_rv) + 'm/s')
print('Orbital Velocity of Moho: ' + str(moho_orbvel) + 'm/s')
print('First burn: ' + str(first_burn) + 'm/s')
print('Second burn: ' + str(second_burn) + 'm/s')
print('\n')

#Define Kerbin
kerbin_radius = 600000
kerbin_g = g
kerbin_stdgp = kerbin_g * kerbin_radius**2
kerbin_rotp = 21549.425
kerbin_atmo_alt = 70000
max_kerbin_circularization_altitude = 200000
krkaaws = (kerbin_radius + kerbin_atmo_alt)**2
krmkcaws = (kerbin_radius + max_kerbin_circularization_altitude)**2
ska = 72000

space_kerbin_orbvel = np.sqrt(kerbin_stdgp/(kerbin_radius + ska))
min_circularization = space_kerbin_orbvel - np.sqrt(kerbin_stdgp*(2/(kerbin_radius + ska) - 1/(kerbin_radius + 0.5*(ska+18500))))
print(min_circularization)

#Functions
def simulate_ion_circularization(apoapsis, periapsis, code, mass, battery, solar, ionThrust, chemThrust):
	#Simulates a given ion burning code.
	#Initialization
	starting_altitude = code[0]
	sma = 0.5*(apoapsis+periapsis) + kerbin_radius
	starting_velocity = np.sqrt(kerbin_stdgp*(2/(kerbin_radius + starting_altitude) - 1/sma))
	#print(starting_velocity)

	eccentricity = (apoapsis+kerbin_radius)/sma - 1
	#print(eccentricity)
	horizontal_velocity = np.sqrt((1 - eccentricity**2)*sma*kerbin_stdgp)/(starting_altitude + kerbin_radius)
	vertical_velocity = np.sqrt(starting_velocity**2 - horizontal_velocity**2)

	x = 0
	y = starting_altitude + kerbin_radius
	vx = horizontal_velocity
	vy = vertical_velocity
	currMass = mass
	currEC = battery
	time = 0
	nextCodeIndex = 1
	currCodeAngle = 0
	currCodeIonThrust = 0
	currCodeChemThrust = 0
	totalXenonUsed = 0
	totalLFOXUsed = 0

	timestep = 0.5 #seconds (everything is in SI units to save a lot of headaches)

	#Simulation
	for i in range(int(3000/timestep)):
		#Calculate gravity forces
		distance = np.sqrt(x**2 + y**2)
		gravForce = kerbin_stdgp/(distance * distance)
		gx = -x*gravForce/distance
		gy = -y*gravForce/distance
		vx += gx * timestep
		vy += gy * timestep

		#Read code
		if nextCodeIndex < len(code):
			if time >= code[nextCodeIndex][0]:
				currCodeAngle = code[nextCodeIndex][1]
				currCodeIonThrust = code[nextCodeIndex][2]
				currCodeChemThrust = code[nextCodeIndex][3]
				nextCodeIndex += 1
		#Calculate thrust
		ionPush = ionThrust * currCodeIonThrust
		ionXenonConsumption = ionPush/(4200*9.81)
		ionECConsumption = ionXenonConsumption * 180000
		#print(ionPush, ionXenonConsumption, ionECConsumption)
		netECGen = solar - ionECConsumption
		if (currEC <= 0) and (netECGen <= 0):
			ionPush *= 0.5 + 0.5*(solar/ionECConsumption) #Ion engine electricity starvation

		chemPush = chemThrust * currCodeChemThrust
		chemLFOXConsumption = chemPush/(305*9.81)
		totalPush = chemPush + ionPush

		#Calculate craft direction
		zeroAngleDirection = [y, -x] #Horizon
		zeroAngleDirectionLength = np.sqrt(zeroAngleDirection[0]**2 + zeroAngleDirection[1]**2)
		zeroAngleDirection[0] /= zeroAngleDirectionLength
		zeroAngleDirection[1] /= zeroAngleDirectionLength

		craftDirection = [
			zeroAngleDirection[0]*np.cos(currCodeAngle) - zeroAngleDirection[1]*np.sin(currCodeAngle),
			zeroAngleDirection[1]*np.cos(currCodeAngle) + zeroAngleDirection[0]*np.sin(currCodeAngle)
		]

		#Update velocities according to thrust from engines
		vx += craftDirection[0]*totalPush*timestep/currMass
		vy += craftDirection[1]*totalPush*timestep/currMass

		x += vx * timestep
		y += vy * timestep
		totalXenonUsed += ionXenonConsumption * timestep
		totalLFOXUsed += chemLFOXConsumption * timestep
		currMass -= (ionXenonConsumption + chemLFOXConsumption) * timestep
		currEC += netECGen * timestep
		currEC = max(currEC, 0)
		currEC = min(currEC, battery)
		time += timestep

		#Calculate Orbital Parameters
		currSMA = 1/(2/distance - (vx**2 + vy**2)/kerbin_stdgp)
		currHVel = (vx * y - vy * x)/distance
		currEccentricity = np.sqrt(1 - (currHVel*distance)**2/(currSMA*kerbin_stdgp))
		currPeri = currSMA * (1 - currEccentricity)
		currApo = currSMA * (1 + currEccentricity)

		if currApo > 9e+6: #The craft was lost in deep space. The kerbal died. The end.
			break
		if currApo < 0: #The craft was DEFINITELY lost in deep space. The kerbal died. The end.
			break
		if currMass <= 0: #The craft dipped into the atmosphere. The kerbal died. The end.
			break

		if currPeri >= kerbin_radius + ska: #The craft succeeded! Farewell and good luck on your further missions!
			periapsisVel = np.sqrt(kerbin_stdgp * (2/currPeri - 1/currSMA))
			periapsisVel -= np.sqrt(kerbin_stdgp/(kerbin_radius + ska))
			massRewarded = currMass * np.exp(periapsisVel/(4200 * 9.81)) #Higher apoapsis means that we'll need less fuel for the rest of the mission
			return currMass, massRewarded, totalXenonUsed, totalLFOXUsed, currApo - kerbin_radius

		#print(x, y, vx, vy, distance, mass, currEC, ionPush, currPeri, currApo)
		if x**2 + y**2 > krmkcaws:
			break #The craft was lost in deep space. The kerbal died. The end.
		if x**2 + y**2 < krkaaws:
			break
	return None #The craft dipped into the atmosphere. The kerbal died. The end.

print(
	simulate_ion_circularization(
		80000,
		-300000,
		[70000, [10, 0.5, 1, 0]],
		5.5,
		1000,
		3,
		12,
		180
	)
)

def calc_optimal_ion_circularization(apoapsis, periapsis, mass, battery, solar, ionThrust, chemThrust):
	#Calculates optimal ion ascent such that a final orbit can be reached.
	print(apoapsis, periapsis, mass, battery, solar, ionThrust, chemThrust)
	bestRewardMass = -1
	bestCode = []
	bestLFOX = -1
	bestXenon = -1
	countObjs = 0

	for i in range(20):
		#Create random codes until one works
		bestRandomRewardMass = -1
		bestRandomCode = []
		bestRandomXenon = -1
		bestRandomLFOX = -1

		for k in range(100):
			randomCode = [np.random.randint(kerbin_atmo_alt, apoapsis)]
			for j in range(np.random.randint(2,4)):
				randomCode.append([0, np.random.random()*np.pi/2, np.random.random(), np.random.random()])
			times = list(np.random.rand(len(randomCode) - 1) * 200)
			times.sort()
			for j in range(1, len(randomCode)):
				randomCode[j][0] = times[j-1]

			#print(randomCode)
			res = simulate_ion_circularization(apoapsis, periapsis, randomCode, mass, battery, solar, ionThrust, chemThrust)
			#print(res)
			if res != None:
				#print(randomCode, res)
				#print(res[1])
				#print('\n')

				if res[1] > bestRandomRewardMass:
					bestRandomRewardMass = res[1]
					bestRandomCode = randomCode
					bestRandomXenon = res[2]
					bestRandomLFOX = res[3]
		print(bestRandomRewardMass, bestRandomCode, bestRandomXenon, bestRandomLFOX)

		if bestRandomRewardMass == -1:
			continue

		#Slowly change the code to make it even better
		oldCode = bestRandomCode
		oldRewardMass = bestRandomRewardMass
		oldXenon = bestRandomXenon
		oldLFOX = bestRandomLFOX
		prevWorse = True

		objCode = bestRandomCode
		objRewardMass = bestRandomRewardMass
		objXenon = bestRandomXenon
		objLFOX = bestRandomLFOX

		r1 = [-1 for t in oldCode]
		r2 = [-1 for t in oldCode]
		r3 = [-1 for t in oldCode]
		r4 = [-1 for t in oldCode]

		for j in range(5000):
			if j % 200 == 0:
				oldCode = dc(objCode)
				oldRewardMass = objRewardMass
				oldXenon = objXenon
				oldLFOX = objLFOX
				print(j)

			temperature = 0.5*(1 - j/8600)
			a = 25
			b = 0.1
			#Calculate random offset
			if prevWorse:
				r5 = np.random.randn() * a

			testCode = dc(oldCode) 
			testCode[0] += r5
			testCode[0] = max(testCode[0], kerbin_atmo_alt)
			testCode[0] = min(testCode[0], apoapsis)

			for k in range(1, len(testCode)):
				if prevWorse:
					r1[k] = np.random.randn() * a
					r2[k] = np.random.randn() * b
					r3[k] = np.random.randn() * b
					r4[k] = np.random.randn() * b
				testCode[k][0] += r1[k]
				testCode[k][0] = max(testCode[k][0], 0)
				if k + 1 < len(testCode):
					testCode[k][0] = min(testCode[k][0], testCode[k+1][0])

				testCode[k][1] += r2[k]
				testCode[k][1] = max(testCode[k][1], 0)
				testCode[k][1] = min(testCode[k][1], np.pi/2)

				testCode[k][2] += r3[k]
				testCode[k][2] = max(testCode[k][2], 0)
				testCode[k][2] = min(testCode[k][2], 1)

				testCode[k][3] += r4[k]
				testCode[k][3] = max(testCode[k][3], 0)
				testCode[k][3] = min(testCode[k][3], 1)

			#Evaluate
			prevWorse = True
			res = simulate_ion_circularization(apoapsis, periapsis, testCode, mass, battery, solar, ionThrust, chemThrust)
			if res == None:
				continue
			if (res[1] > oldRewardMass) or ((np.random.random() < temperature) and (res[1] > oldRewardMass * 0.998)):
				oldCode = testCode
				oldRewardMass = res[1]
				oldXenon = res[2]
				oldLFOX = res[3]
				#prevWorse = False
				#print(oldRewardMass, j)
			if (res[1] > objRewardMass):
				objCode = testCode
				objRewardMass = res[1]
				objXenon = res[2]
				objLFOX = res[3]
				countObjs += 1
				print(objRewardMass, objXenon, objLFOX, apoapsis, periapsis, mass, j)
		print(oldCode, oldRewardMass)
		print(objCode, objRewardMass)

		if objRewardMass > bestRewardMass:
			bestRewardMass = objRewardMass
			bestCode = objCode
			bestLFOX = objLFOX
			bestXenon = objXenon
	print(countObjs)
	return bestRewardMass, bestXenon, bestLFOX, 0, bestCode

"""
print(calc_optimal_ion_circularization(
	80000,
	-300000,
	5.5,
	1000,
	3,
	4,
	180
))
"""

#Simulation
def simulate(mass, battery, solar, ionThrust, chemThrust):
	print('Starting Mass: ' + str(mass) + 't')
	xe_used = 0
	lfox_used = 0
	
	lf_used = ssto_lf_fraction * mass
	curr_mass = mass - lf_used
	print('Mass at RAPIER closed-cycle: ' + str(curr_mass) + 't')

	"""for ion_circ in range(0, 500, 10):
		space_vel = space_kerbin_orbvel - ion_circ
		space_a = 1/(2/(kerbin_radius + ska) - space_vel*space_vel/kerbin_stdgp)
		high_alt_vel = np.sqrt(kerbin_stdgp*(2/(kerbin_radius + 18500) - 1/(space_a)))
		high_alt_vel -= 2*np.pi*(kerbin_radius + 18500)/(kerbin_rotp)
		high_alt_vel -= 1650

		rapier_ratio = np.exp(high_alt_vel/(305*9.81))

		calc_mass = curr_mass
		calc_xe = xe_used
		calc_lf = lf_used
		calc_lfox = lfox_used

		calc_lfox += calc_mass - calc_mass/(rapier_ratio)
		calc_mass /= rapier_ratio
		print(calc_mass, calc_xe, calc_lf, calc_lfox)"""
	#apostep = 10000

	lowPeriapsis = 0
	highPeriapsis = kerbin_radius + 18500
	lowApoapsis = kerbin_radius + kerbin_atmo_alt
	highApoapsis = kerbin_radius + 140000

	periapsis = kerbin_radius - 300000
	apoapsis = kerbin_radius + 80000

	sma = 0.5*(apoapsis + periapsis)
	highvel = np.sqrt(kerbin_stdgp*(2/(kerbin_radius+18500) - 1/sma))
	highvel -= 2*np.pi*(kerbin_radius+18500)/kerbin_rotp
	highvel -= 1650
	if highvel < 0:
		highvel = 0
			
	calcMass = curr_mass
	calcLF = lf_used
	calcLFOX = lfox_used
	calcXenon = xe_used

	massRatio = np.exp(highvel/(305*9.81))
	calcLFOX += calcMass - calcMass/massRatio
	calcMass /= massRatio

	res = calc_optimal_ion_circularization(apoapsis-kerbin_radius, periapsis-kerbin_radius, calcMass, battery, solar, ionThrust, chemThrust)
	if res[0] != -1:
		print(res)
		calcMass = res[0]
		calcLFOX += res[2]
		calcXenon += res[1]
	else:
		print(periapsis, apoapsis, sma, highvel, calcMass)

	"""
	peristep = highPeriapsis - lowPeriapsis
	peristep /= 5

	for apoapsis in range(kerbin_radius + kerbin_atmo_alt+1, kerbin_radius + 130000, apostep):
		for i in range(5):
			periapsisMasses = []
			for periapsis in range(int(lowPeriapsis), int(highPeriapsis), int(peristep)):
				sma = 0.5*(apoapsis + periapsis)
				highvel = np.sqrt(kerbin_stdgp*(2/(kerbin_radius+18500) - 1/sma))
				highvel -= 2*np.pi*(kerbin_radius+18500)/kerbin_rotp
				highvel -= 1650
				if highvel < 0:
					highvel = 0
			
				calcMass = curr_mass
				calcLF = lf_used
				calcLFOX = lfox_used
				calcXenon = xe_used

				massRatio = np.exp(highvel/(305*9.81))
				calcLFOX += calcMass - calcMass/massRatio
				calcMass /= massRatio

				res = calc_optimal_ion_circularization(apoapsis-kerbin_radius, periapsis-kerbin_radius, calcMass, battery, solar, ionThrust, chemThrust)
				if res[0] != -1:
					print(res)
					calcMass = res[0]
					calcLFOX += res[2]
					calcXenon += res[1]
					periapsisMasses.append(res[0])
				else:
					periapsisMasses.append(10**10)
				print(periapsis, apoapsis, sma, highvel, calcMass, lowPeriapsis, highPeriapsis, peristep)

			if min(periapsisMasses) == 10**10:
				break
			lowPeriapsis = max(0, (periapsisMasses.index(max(periapsisMasses)) - 2)*peristep)
			highPeriapsis = min(kerbin_radius+18500, (periapsisMasses.index(max(periapsisMasses)) + 2)*peristep)
			peristep /= 5
	"""
	for i in range(20):
		oldApoapsis, oldPeriapsis = apoapsis, periapsis

		apoapsis += 10000*np.random.randn()
		periapsis += 100000*np.random.randn()

		apoapsis = min(apoapsis, highApoapsis)
		apoapsis = max(apoapsis, lowApoapsis)

		periapsis = min(periapsis, highPeriapsis)
		periapsis = max(periapsis, lowPeriapsis)

		sma = 0.5*(apoapsis + periapsis)
		highvel = np.sqrt(kerbin_stdgp*(2/(kerbin_radius+18500) - 1/sma))
		highvel -= 2*np.pi*(kerbin_radius+18500)/kerbin_rotp
		highvel -= 1650
		if highvel < 0:
			highvel = 0
			
		calcMass = curr_mass
		calcLF = lf_used
		calcLFOX = lfox_used
		calcXenon = xe_used

		massRatio = np.exp(highvel/(305*9.81))
		calcLFOX += calcMass - calcMass/massRatio
		calcMass /= massRatio

		res = calc_optimal_ion_circularization(apoapsis-kerbin_radius, periapsis-kerbin_radius, calcMass, battery, solar, ionThrust, chemThrust)
		if res[0] != -1:
			print(res)
			calcMass = res[0]
			calcLFOX += res[2]
			calcXenon += res[1]
		else:
			apoapsis, periapsis = oldApoapsis, oldPeriapsis
			print(periapsis, apoapsis, sma, highvel, calcMass)
	return xe_used, lf_used, lfox_used


print(simulate(6, 1000, 3, 4, 180))

#Optimization
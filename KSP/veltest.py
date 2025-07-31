import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

#Constants
generator = np.random.default_rng()

timestep = 3
g0 = 9.80665

planetRadius = 6e5
planetGravity = 0.8 * g0
planetRotationalPeriod = 211926.35802123

planetGravParameter = planetRadius * planetRadius * planetGravity

#Boundary conditions
initDist = 605000
initHV = 2 * np.pi * initDist / planetRotationalPeriod
engineAccels = [8]
engineIsps = [320]

finalDist = 602500
finalHV = np.sqrt(planetGravParameter / finalDist)
finalTime = 500

#Calculations
STEPS = int(finalTime / timestep)

def getConstAltitudeInitialization():
	mass = 1

	gravity = planetGravParameter / initDist / initDist
	rotvel = 2 * np.pi * initDist / planetRotationalPeriod

	vel = np.sqrt(planetGravParameter / initDist)
	sumthrottles = 1

	vels = [vel]

	while vel > rotvel:
		#Effective gravity
		vaccel = gravity - vel * vel / initDist

		#Throttles
		bestfuelconsumption = 0
		bestratio = 0
		besthaccel = 0
		
		newthrottle = sumthrottles - 0.001
		improvedPrev = True
		while improvedPrev:
			newthrottle += 0.001
			if newthrottle > len(engineAccels):
				break

			thrust = 0
			fuelconsumption = 0
			for i in range(int(newthrottle)):
				thrust += engineAccels[i]
				fuelconsumption += engineAccels[i] / engineIsps[i] / g0
		
			if int(newthrottle) < len(engineAccels):
				thrust += engineAccels[int(newthrottle)] * (newthrottle % 1)
				fuelconsumption += \
					engineAccels[int(newthrottle)] * \
					(newthrottle % 1) / \
					engineIsps[int(newthrottle)] / \
					g0

			accel = thrust / mass
			if accel > vaccel:
				haccel = np.sqrt(accel * accel - vaccel * vaccel)
			else:
				haccel = 0

			if haccel / fuelconsumption > bestratio:
				bestratio = haccel / fuelconsumption
				besthaccel = haccel
				bestfuelconsumption = fuelconsumption

				sumthrottles = newthrottle

			improvedPrev = haccel / fuelconsumption >= bestratio

		#Integrate
		vel -= haccel * timestep
		mass -= fuelconsumption * timestep

		vels.append(vel)

	vels[-1] = rotvel
	vels = [np.sqrt(planetGravParameter / initDist) for i in range(STEPS - len(vels))] + vels
	
	vel = rotvel

	while vel < np.sqrt(planetGravParameter / initDist):
		#Effective gravity
		vaccel = gravity - vel * vel / initDist

		#Throttles
		bestfuelconsumption = 0
		bestratio = 0
		besthaccel = 0
		
		newthrottle = sumthrottles + 0.001
		improvedPrev = True
		while improvedPrev:
			newthrottle -= 0.001
			if newthrottle <= 0.5:
				break

			thrust = 0
			fuelconsumption = 0
			for i in range(int(newthrottle)):
				thrust += engineAccels[i]
				fuelconsumption += engineAccels[i] / engineIsps[i] / g0
		
			if int(newthrottle) < len(engineAccels):
				thrust += engineAccels[int(newthrottle)] * (newthrottle % 1)
				fuelconsumption += \
					engineAccels[int(newthrottle)] * \
					(newthrottle % 1) / \
					engineIsps[int(newthrottle)] / \
					g0

			accel = thrust / mass
			if accel > vaccel:
				haccel = np.sqrt(accel * accel - vaccel * vaccel)
			else:
				haccel = 0

			if haccel / fuelconsumption > bestratio:
				bestratio = haccel / fuelconsumption
				besthaccel = haccel
				bestfuelconsumption = fuelconsumption

				sumthrottles = newthrottle

			improvedPrev = haccel / fuelconsumption >= bestratio

		#Integrate
		vel += haccel * timestep
		mass -= fuelconsumption * timestep

		vels.append(vel)

	vels[-1] = np.sqrt(planetGravParameter / initDist)
	vels = vels + [np.sqrt(planetGravParameter / initDist) for i in range(2 * STEPS - len(vels) - 1)]

	radiis = np.full((2 * STEPS - 1,), initDist)

	return np.array(vels), radiis

def simulate(vels, radiis, retAccels = False):
	thetaprime = vels / radiis
	rprime = (radiis - np.roll(radiis, 1)) / timestep
	rprime[0] = rprime[1]

	thetaprimeprime = (thetaprime - np.roll(thetaprime, 1)) / timestep
	thetaprimeprime[0] = thetaprimeprime[1]
	rprimeprime = (rprime - np.roll(rprime, 1)) / timestep
	rprimeprime[0] = rprimeprime[1]

	av = rprimeprime - radiis * thetaprime * thetaprime + planetGravParameter / radiis / radiis
	ah = 2 * rprime * thetaprime + radiis * thetaprimeprime

	accels = np.sqrt(av * av + ah * ah)

	costs = 0
	mass = 1
	maxAccels = []
	for i in range(len(vels)):
		thrust = accels[i] * mass

		engineIdx = 0
		while thrust > 0 and engineIdx < len(engineAccels):
			groupThrust = min(thrust, engineAccels[engineIdx])
			
			costs += groupThrust * timestep / engineIsps[engineIdx] / g0
			mass -= groupThrust * timestep / engineIsps[engineIdx] / g0

			thrust -= groupThrust
			engineIdx += 1

		if thrust > 0:
			costs += thrust * timestep

		maxAccels.append(sum(engineAccels) / mass)

	if retAccels:
		return accels, ah, av, maxAccels

	return costs, 1 - mass

radiis = np.linspace(0, 1, num=2 * STEPS - 1)
radiis = np.abs(radiis - 0.5) * 2

vels = initHV + radiis * (finalHV - initHV)
radiis = radiis * radiis * (3 - 2 * radiis)
radiis = initDist + (finalDist - initDist) * radiis

"""vels, radiis = getConstAltitudeInitialization()

radiis = np.linspace(0, 1, num=2 * STEPS - 1)
radiis = np.abs(radiis - 0.5) * 2
radiis = radiis * radiis * (3 - 2 * radiis)
radiis = initDist + (finalDist - initDist) * radiis"""

#Optimize
bestvels, bestradiis = deepcopy((vels, radiis))

score, mass = simulate(vels, radiis)

bestscore = score

scores = []
anyOrbit = True
improvedPrev = False
for i in range(200000):
	newvels, newradiis = deepcopy((vels, radiis))

	if not improvedPrev:
		offsetloc = int(generator.uniform(2, 2 * STEPS - 2))
		offsetheight = abs(generator.normal(loc=0, scale=50))
		offsetwidth = np.exp(generator.normal(loc=4, scale=4))
		offsetangle = generator.uniform(0, 2 * np.pi)
		offsettype = generator.integers(2)
	minRange = 2
	maxRange = 2 * STEPS - 3
	if anyOrbit:
		minRange = 0
		maxRange += 2
	for j in range(minRange, maxRange):
		if j == STEPS - 2 or j == STEPS - 1 or j == STEPS:
			continue

		if offsettype == 0:
			magnitude = offsetheight * np.exp(-((j - offsetloc) / offsetwidth) ** 2)
		else:
			magnitude = 0
			if (offsetloc < STEPS - 1) == (j < STEPS - 1):
				cycles = int((STEPS - 1) / offsetwidth)
				magnitude = offsetheight * (np.cos(2 * np.pi * cycles * j / (STEPS - 1)) - 1)
			
		newvels[j] += np.cos(offsetangle) * magnitude
		newradiis[j] += np.sin(offsetangle) * magnitude

		if newradiis[j] < planetRadius:
			newradiis[j] = planetRadius

	newvels[0] = np.sqrt(planetGravParameter / newradiis[0])
	newvels[-1] = np.sqrt(planetGravParameter / newradiis[-1])

	newscore, newmass = simulate(newvels, newradiis)

	if newscore < score or generator.random() < 0.0:
		score = newscore

		vels, radiis = deepcopy((newvels, newradiis))

	if newscore < bestscore:
		bestscore = newscore

		bestvels, bestradiis = deepcopy((newvels, newradiis))

	improvedPrev = newscore < score

	if i % 100 == 0:
		print(i, round(bestscore, 5))

	if i % 1000 == 0:
		score = bestscore
		vels, radiis = deepcopy((bestvels, bestradiis))

	scores.append(score)

#plt.plot(scores)
fig, axs = plt.subplots(5)

accels, ah, av, maxAccels = simulate(bestvels, bestradiis, True)
score, mass = simulate(bestvels, bestradiis)

grav = planetGravParameter / bestradiis / bestradiis - bestvels * bestvels / bestradiis

print(accels / (finalHV - bestvels))
print(mass)

axs[0].plot(bestvels)
axs[1].plot(bestradiis)
axs[2].plot(accels)
axs[2].plot(maxAccels)
axs[3].plot(ah)
axs[3].plot(av)
axs[3].plot(grav)
axs[3].plot(av - grav)
axs[4].plot(accels / (finalHV - bestvels))
plt.show()
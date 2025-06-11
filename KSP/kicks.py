#Optimize the timing of periapsis kicks
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy

#Functions
def plotOrbit(a, e, argp, col='white'):
	v = np.linspace(0, 2 * np.pi, num=100)
	dist = a * (1 - e * e)/(1 + e * np.cos(v))
	x = dist * np.cos(v + argp)
	y = dist * np.sin(v + argp)
	plt.plot(x, y, c=col)

def plotOrbitPoint(a, e, argp, v, col='bo'):
	dist = a * (1 - e * e)/(1 + e * np.cos(v))
	x = dist * np.cos(v + argp)
	y = dist * np.sin(v + argp)
	plt.plot(x, y, col)

def getTransferTime(waits, burns):
	vel = initVel
	sma = 1/(2/(planetR + initPeriapsis) - (vel * vel)/planetStdGP)
	period = 2 * np.pi * sma * np.sqrt(sma/planetStdGP)
	time = timeUntilPeriapsis + period * (waits[0] + 1)
	for i in range(len(burns)):
		vel += burns[i]
		sma = 1/(2/(planetR + initPeriapsis) - (vel * vel)/planetStdGP)
		period = 2 * np.pi * sma * np.sqrt(sma/planetStdGP)
		time += period * (waits[i + 1] + 1)
	vel = initVel + burnRequired
	sma = 1/(2/(planetR + initPeriapsis) - (vel * vel)/planetStdGP)
	period = 2 * np.pi * sma * np.sqrt(sma/planetStdGP)
	time += period * waits[-1]
	time += finalTransferTime
	return time

def eval(waits, burns):
	time = getTransferTime(waits, burns)
	err = (time - requiredTime)/targetPeriod%1
	err = 1.01 - 2 * abs(err - 0.5)
	return err * max(time, targetPeriod)

#Planet parameters
planetR = 600000
planetG = 9.81
planetStdGP = planetG * planetR * planetR

#Initial orbit
initPeriapsis = 72500
initApoapsis = 191200
timeUntilPeriapsis = 1708

#Target body
targetSMA = 12000000
targetEccentricity = 0
targetArgP = 0
targetV = 2

#Calculations
initSMA = planetR + 0.5 * (initPeriapsis + initApoapsis)
initEccentricity = (planetR + initApoapsis)/initSMA - 1
initVel = np.sqrt(planetStdGP * (2/(planetR + initPeriapsis) - 1/initSMA))

e = targetEccentricity
p = initPeriapsis + planetR
d = targetSMA * (1 - targetEccentricity * targetEccentricity)
c = 0
f = targetArgP
q1 = p * p * (e * e - 1)
q2 = 2 * p * d
qa = q1 + d * d - q2 * e * np.cos(c-f)
qb = 2 * q1 - q2 * e * np.cos(c-f) + q2
qc = q1 - d * d + q2
determinant = qb*qb - 4*qa*qc
finalEccentricity = (-qb + np.sqrt(determinant))/(2 * qa)
finalSMA = (initPeriapsis + planetR)/(1 - finalEccentricity)
finalVel = np.sqrt(planetStdGP * (2/(planetR + initPeriapsis) - 1/finalSMA))
burnRequired = finalVel - initVel
a = (1 + finalEccentricity) * (planetR + initPeriapsis)

c1 = finalEccentricity * d * np.cos(c) - a * e * np.cos(f)
c2 = finalEccentricity * d * np.sin(c) - a * e * np.sin(f)
cd = np.sqrt(c1 * c1 + c2 * c2)
s = np.arctan2(c1/cd, c2/cd)
intersect = np.pi/2 - s if a > d else 3*np.pi/2 - s
targetIntersectTrueAnomaly = intersect - f
finalIntersectTrueAnomaly = intersect - c

targetPeriod = 2 * np.pi * targetSMA * np.sqrt(targetSMA/planetStdGP)
finalPeriod = 2 * np.pi * finalSMA * np.sqrt(finalSMA/planetStdGP)
targetEccentricAnomaly = np.arctan(np.tan(targetV/2)/np.sqrt((1 + targetEccentricity)/(1 - targetEccentricity))) * 2
targetMeanAnomaly = targetEccentricAnomaly - targetEccentricity * np.sin(targetEccentricAnomaly)
targetMeanAnomaly *= targetPeriod/(2 * np.pi)

targetIntersectEccentricAnomaly = np.arctan(np.tan(targetIntersectTrueAnomaly/2)/np.sqrt((1 + targetEccentricity)/(1 - targetEccentricity))) * 2
targetIntersectMeanAnomaly = targetIntersectEccentricAnomaly - targetEccentricity * np.sin(targetIntersectEccentricAnomaly)
targetIntersectMeanAnomaly *= targetPeriod/(2 * np.pi)
requiredTime = targetIntersectMeanAnomaly - targetMeanAnomaly

finalIntersectEccentricAnomaly = np.arctan(np.tan(finalIntersectTrueAnomaly/2)/np.sqrt((1 + targetEccentricity)/(1 - targetEccentricity))) * 2
finalTransferTime = finalIntersectEccentricAnomaly - finalEccentricity * np.sin(finalIntersectEccentricAnomaly)
finalTransferTime *= finalPeriod/(2 * np.pi)

#Optimization settings
numBurns = 0
maxBurn = 50

if maxBurn * numBurns < burnRequired:
	print("Not enough burns provided, bumping it up to " + str(int(np.ceil(burnRequired/maxBurn))) + " burns")
	numBurns = int(np.ceil(burnRequired/maxBurn))

#Optimization
oldWaits = [0 for i in range(numBurns + 1)]
oldBurns = [burnRequired/numBurns for i in range(numBurns - 1)]
oldEval = eval(oldWaits, oldBurns)
for i in range(100000):
	newWaits = deepcopy(oldWaits)
	newBurns = deepcopy(oldBurns)
	for j in range(numBurns + 1):
		if np.random.random() < 0.4:
			newWaits[j] -= 1
		if np.random.random() < 0.2:
			newWaits[j] += 1
		if newWaits[j] < 0:
			newWaits[j] = 0
	for j in range(numBurns - 1):
		delta = 30 * np.random.normal()
		newBurns[j] += delta
		for k in range(numBurns - 1):
			newBurns[k] -= delta/numBurns
		fails = False
		for k in range(numBurns - 1):
			if newBurns[k] < 0 or newBurns[k] > maxBurn:
				fails = True
				break
		if fails:
			newBurns[j] -= delta
			for k in range(numBurns - 1):
				newBurns[k] += delta/numBurns
	s = sum(newBurns)
	if s > burnRequired:
		for j in range(numBurns - 1):
			newBurns[j] *= burnRequired/s
	if s < burnRequired - maxBurn:
		for j in range(numBurns - 1):
			newBurns[j] *= (burnRequired - maxBurn)/s

	newEval = eval(newWaits, newBurns)
	if newEval < oldEval:
		oldWaits = newWaits
		oldBurns = newBurns
		oldEval = newEval

print("Best evaluation: " + str(oldEval))
print("Time to first encounter: " + str(getTransferTime(oldWaits, oldBurns)) + "s")
print("List of burns (m/s): " + str(oldBurns + [burnRequired - sum(oldBurns)]))
print("Number of orbits to skip before/after each burn: " + str(oldWaits))

#Plot
planet = plt.Circle((0,0), planetR, color='blue')

ax = plt.gca()
ax.set_aspect('equal')
ax.add_patch(planet)
ax.set_facecolor('black')
vel = initVel
sma = 1/(2/(planetR + initPeriapsis) - (vel * vel)/planetStdGP)
eccen = 1 - (planetR + initPeriapsis)/sma
plotOrbit(sma, eccen, 0)
for i in range(numBurns - 1):
	vel += oldBurns[i]
	sma = 1/(2/(planetR + initPeriapsis) - (vel * vel)/planetStdGP)
	eccen = 1 - (planetR + initPeriapsis)/sma
	plotOrbit(sma, eccen, 0)
plotOrbit(finalSMA, finalEccentricity, 0)
plotOrbit(targetSMA, targetEccentricity, targetArgP, 'orange')
plotOrbitPoint(targetSMA, targetEccentricity, targetArgP, targetV, 'ro')
plt.show()
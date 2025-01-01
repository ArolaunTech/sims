import numpy as np
import matplotlib.pyplot as plt

#Planet parameters
planetSMA = 5263138304
planetEccentricity = 0.2
planetSTDGP = 250000*250000*2.7
planetSOI = 9646663

#Sun parameters
solarRadius = 261600000
solarSTDGP = 1.7 * 9.81 * solarRadius * solarRadius #G*M_sun
#planetSOI = planetSMA * (planetSTDGP/solarSTDGP)**0.4

#Ship parameters
initRV = 30
initAngle = np.random.random() * 2 * np.pi
planetPeriapsis = 270000
maxAssists = 200
endingRV = 2110
print(initAngle)

def getMaxDeflection(relvel):
	if relvel < np.sqrt(planetSTDGP/planetSOI):
		return np.pi
	sma = 1/(2/planetSOI - relvel * relvel/planetSTDGP)
	eccentricity = np.sqrt(1 - (planetPeriapsis**2) * (2/planetPeriapsis - 1/sma)/sma)
	q = (sma * (1 - eccentricity * eccentricity)/planetSOI - 1)/eccentricity
	hvel = planetPeriapsis * np.sqrt(planetSTDGP * (2/planetPeriapsis - 1/sma))/planetSOI
	vvel = np.sqrt(relvel * relvel - hvel * hvel)
	return 2 * np.arcsin((hvel * np.sqrt(1 - q * q) - vvel * q)/relvel)

def getDist(a, e, p, v):
	return a * (1 - e * e) / (1 + e * np.cos(v - p))

def getVelocityHV(a, e, p, v):
	#Find horizontal and vertical velocity
	dist = getDist(a, e, p, v)
	vel = np.sqrt(solarSTDGP * (2/dist - 1/a))
	hvel = np.sqrt(a * solarSTDGP * (1 - e * e))/dist
	vvel = np.sqrt(vel * vel - hvel * hvel)

	#Find if before or after apoapsis
	if (v - p)%(2 * np.pi) > np.pi:
		vvel *= -1

	return hvel, vvel

def getVelocityXY(a, e, p, v):
	hvel, vvel = getVelocityHV(a, e, p, v)
	hx = -np.sin(v)
	hy = np.cos(v)
	vx = hy
	vy = -hx

	return hx * hvel + vx * vvel, hy * hvel + vy * vvel

def plotOrbit2D(a, e, p, col):
	orbitAngles = np.linspace(0, 2*np.pi, num=100)
	dist = getDist(a, e, p, orbitAngles)
	x = dist * np.cos(orbitAngles)
	y = dist * np.sin(orbitAngles)
	plt.plot(x, y, color=col)

def getNewIntersection(a, e, p, initAngle, initRV, ejectionAngle, plotOrbits=False):
	initDist = getDist(a, e, p, initAngle)
	xvel, yvel = getVelocityXY(a, e, p, initAngle)
	rx = initRV * np.cos(ejectionAngle)
	ry = initRV * np.sin(ejectionAngle)

	nx = xvel + rx
	ny = yvel + ry

	hx = -np.sin(initAngle)
	hy = np.cos(initAngle)
	vx = hy
	vy = -hx

	nhvel = nx * hx + ny * hy
	nvvel = nx * vx + ny * vy

	#print(nhvel, nvvel)
	na = 1/(2/initDist - (nx * nx + ny * ny)/solarSTDGP)
	ne = np.sqrt(1 - ((nhvel * initDist)**2)/(na * solarSTDGP))
	nAngleApoapsis = np.arccos((na * (1 - ne * ne)/initDist - 1)/ne)

	if nvvel < 0:
		nperi = initAngle + nAngleApoapsis
	else:
		nperi = initAngle - nAngleApoapsis

	c1o = a * (1 - e * e)
	c1n = na * (1 - ne * ne)

	c1 = e*c1n*np.cos(p) - ne*c1o*np.cos(nperi)
	c2 = e*c1n*np.sin(p) - ne*c1o*np.sin(nperi)
	c3 = np.sqrt(c1 * c1 + c2 * c2)
	c4 = np.arctan2(c1, c2)
	s1 = np.arcsin((c1o-c1n)/c3)-c4
	s2 = np.pi-np.arcsin((c1o-c1n)/c3)-c4

	solution = s1
	if (np.abs(initAngle - s1))%(2 * np.pi) < (np.abs(initAngle - s2))%(2 * np.pi):
		solution = s2

	solutionDist = getDist(a, e, p, solution)

	solutionHOld, solutionVOld = getVelocityXY(a, e, p, solution)
	solutionHNew, solutionVNew = getVelocityXY(na, ne, nperi, solution)
	intersectionAngle = np.arctan2(solutionVNew - solutionVOld, solutionHNew - solutionHOld)

	nrelvel = np.sqrt((solutionHOld - solutionHNew)**2 + (solutionVOld - solutionVNew)**2)

	solutionX = solutionDist * np.cos(solution)
	solutionY = solutionDist * np.sin(solution)

	#print(solution, solutionDist)

	if plotOrbits:
		plotOrbit2D(na, ne, nperi, 'white')
		plt.plot(solutionX, solutionY, 'ro')

	return nrelvel, solution, intersectionAngle

def findBestSequence(a, e, p, initAngle, initRV, numAssists, targetRelVel, plotting):
	assistData = []

	currRV = initRV
	currAngle = initAngle
	incomingAngle = None
	assistsDone = 0
	for i in range(numAssists):
		bestAngle = 0
		bestRelvel = 0
		bestNewAngle = 0
		bestIncomingAngle = 0

		for j in range(1000):
			if j > 500:
				angle = 0.01*(np.random.random() * 2 - 1) + bestAngle
			else:
				angle = j * np.pi * 0.004
			relvel, newAngle, newincomingAngle = getNewIntersection(planetSMA, planetEccentricity, 0, currAngle, currRV, angle, False)
			if relvel > bestRelvel:
				bestRelvel = relvel
				bestAngle = angle
				bestNewAngle = newAngle
				bestIncomingAngle = newincomingAngle
		if incomingAngle is None:
			requiredDeflection = 0
		else:
			incomingAngle = incomingAngle % (2 * np.pi)
			bestAngle = bestAngle % (2 * np.pi)
			requiredDeflection = abs(incomingAngle - bestAngle)
			if requiredDeflection > np.pi:
				requiredDeflection = 2 * np.pi - requiredDeflection
		maxDeflection = getMaxDeflection(currRV)
		neededAssists = int(np.ceil(requiredDeflection/maxDeflection))
		if neededAssists == 0:
			neededAssists = 1
		if plotting:
			print("Iteration:", i, "# of assists before iteration:", assistsDone, "Relative velocity to planet during iteration:", currRV, "m/s")
		assistData.append([i, assistsDone, currRV])
		assistsDone += neededAssists
		if bestRelvel > targetRelVel:
			if plotting:
				print("Done after " + str(assistsDone) + " assists, Relative velocity: " + str(bestRelvel) + " m/s")
			assistData.append([i+1, assistsDone, bestRelvel])
			break

		if plotting:
			getNewIntersection(planetSMA, planetEccentricity, 0, currAngle, currRV, bestAngle, True)

		currAngle = bestNewAngle
		currRV = bestRelvel
		incomingAngle = bestIncomingAngle
	return assistData

fig, ax = plt.subplots()
ax.set_facecolor('black')
ax.set_aspect('equal')

sun = plt.Circle((0, 0), radius=solarRadius, color=(1.0, 1.0, 0.8, 1.0))
ax.add_patch(sun)

plotOrbit2D(planetSMA, planetEccentricity, 0, 'purple')

initDist = getDist(planetSMA, planetEccentricity, 0, initAngle)
initX = initDist * np.cos(initAngle)
initY = initDist * np.sin(initAngle)

initHV, initVV = getVelocityHV(planetSMA, planetEccentricity, 0, initAngle)
initXV, initYV = getVelocityXY(planetSMA, planetEccentricity, 0, initAngle)
print(initX, initY)
print("Hvel: " + str(initHV) + "m/s, Vvel: " + str(initVV) + "m/s")
print("Xvel: " + str(initXV) + "m/s, Yvel: " + str(initYV) + "m/s")

savingsPerAssist = []
assists = []

for i in range(10):
	assistData = findBestSequence(planetSMA, planetEccentricity, 0, initAngle, initRV, maxAssists, endingRV, False)
	numAssists = assistData[-1][1]
	finalRV = assistData[-1][2]
	finalSMA = 1/(2/planetSOI - finalRV * finalRV/planetSTDGP)
	finalPeriapsisVelocity = np.sqrt(planetSTDGP * (2/planetPeriapsis - 1/finalSMA))
	for point in assistData:
		relvel = point[2]
		sma = 1/(2/planetSOI - relvel * relvel/planetSTDGP)
		periapsisVelocity = np.sqrt(planetSTDGP * (2/planetPeriapsis - 1/sma))
		savingsPerAssist.append(finalPeriapsisVelocity - periapsisVelocity - 0.1 * (numAssists - point[1]))
		assists.append(numAssists - point[1])

findBestSequence(planetSMA, planetEccentricity, 0, initAngle, initRV, maxAssists, endingRV, True)

plt.plot(initX, initY, 'bx')
plt.show()

plt.plot(assists, savingsPerAssist, "rx")
plt.axhline(y=1040, color="black", linestyle="dashed", label="Messenger Assists (est.)")
plt.legend()
plt.title("Δv savings of repeated Moho assists over direct Moho-Eve transfer")
plt.xlabel("Number of Moho gravity assists done")
plt.ylabel("Δv savings over direct transfer (m/s)")
plt.show()
import math
import matplotlib.pyplot as plt
import random

g0 = 9.80665

planetG = 1.628
planetR = 200000
rotPeriod = 138984.376574476

isp = 320
veff = isp * g0

landingAlt = 5000

timestep = 0.1

planetSTDGP = planetG * planetR * planetR
landingOrbVel = math.sqrt(planetSTDGP/(planetR + landingAlt))
landingGrav = planetSTDGP/(planetR + landingAlt)**2
rotVel = 2 * math.pi * (planetR + landingAlt)/rotPeriod

def testConstAltitude(accel, log):
	mass = 1
	thrust = accel
	vel = landingOrbVel
	time = 0

	flightLog = []

	while vel > rotVel:
		grav = landingGrav - vel * vel/(planetR + landingAlt)
		accel = thrust/mass

		haccel = 0
		if accel > grav:
			haccel = math.sqrt(accel * accel - grav * grav)

		vel -= haccel * timestep
		mass -= thrust/veff * timestep
		time += timestep

		if log:
			flightLog.append((time, vel, mass))

	vel = rotVel

	while vel < landingOrbVel:
		grav = landingGrav - vel * vel/(planetR + landingAlt)
		accel = thrust/mass

		haccel = 0
		if accel > grav:
			haccel = math.sqrt(accel * accel - grav * grav)

		vel += haccel * timestep
		mass -= thrust/veff * timestep
		time += timestep

		if log:
			flightLog.append((time, vel, mass))

	return flightLog

def testConstHaccelDescentScore(accel, start, haccel):
	mass = 1
	thrust = accel

	x = 0
	y = planetR + start
	vx = math.sqrt(planetSTDGP/(planetR + start))
	vy = 0
	vh = vx
	vv = 0
	while vh > rotVel:
		dist = math.sqrt(x * x + y * y)
		if dist < planetR:
			return 1e9, dist, vv, mass
		grav = planetSTDGP/(dist * dist)
		gx = -x * grav/dist
		gy = -y * grav/dist

		accel = thrust/mass
		if abs(haccel) < abs(accel):
			vaccel = math.sqrt(accel * accel - haccel * haccel)
			actualhaccel = haccel
		else:
			vaccel = 0
			actualhaccel = accel

		ax = (vaccel * x - actualhaccel * y)/dist
		ay = (vaccel * y + actualhaccel * x)/dist

		vx += (gx + ax) * timestep
		vy += (gy + ay) * timestep
		x += vx * timestep
		y += vy * timestep
		mass -= thrust/veff * timestep
		if mass < 0:
			return 1e9, dist, vv, mass

		vh = (y * vx - x * vy)/dist
		vv = (x * vx + y * vy)/dist

	if (accel < landingGrav) or (dist < planetR + landingAlt):
		return 1e9, dist, vv, mass
	height = dist - planetR - landingAlt
	if 2 * (accel - landingGrav) * height < vv * vv:
		return 1e9, dist, vv, mass
	landingBurn = thrust/veff * math.sqrt((2 * height * landingGrav + vv * vv)/(accel * (accel - landingGrav)))
	if landingBurn > mass:
		return 1e9, dist, vv, mass
	mass -= landingBurn

	return -mass, dist, vv, mass

def testConstHaccel1(accel, log, start, haccel):
	mass = 1
	thrust = accel

	flightLog = []

	x = 0
	y = planetR + start
	vx = math.sqrt(planetSTDGP/(planetR + start))
	vy = 0
	vh = vx
	vv = 0

	time = 0
	while vh > rotVel:
		dist = math.sqrt(x * x + y * y)
		grav = planetSTDGP/(dist * dist)
		gx = -x * grav/dist
		gy = -y * grav/dist

		accel = thrust/mass
		if haccel < accel:
			vaccel = math.sqrt(accel * accel - haccel * haccel)
			actualhaccel = haccel
		else:
			vaccel = 0
			actualhaccel = accel

		ax = (vaccel * x - actualhaccel * y)/dist
		ay = (vaccel * y + actualhaccel * x)/dist

		vx += (gx + ax) * timestep
		vy += (gy + ay) * timestep
		x += vx * timestep
		y += vy * timestep
		time += timestep
		mass -= thrust/veff * timestep
		if mass < 0:
			break

		vh = (y * vx - x * vy)/dist
		vv = (x * vx + y * vy)/dist

		if log:
			flightLog.append((time, vh, mass))

	if (accel < landingGrav) or (dist < planetR + landingAlt):
		return []
	height = dist - planetR - landingAlt
	if 2 * (accel - landingGrav) * height < vv * vv:
		return []
	landingBurn = thrust/veff * math.sqrt((2 * height * landingGrav + vv * vv)/(accel * (accel - landingGrav)))
	if landingBurn > mass:
		return []
	mass -= landingBurn

	vel = rotVel
	while vel < landingOrbVel:
		grav = landingGrav - vel * vel/(planetR + landingAlt)
		accel = thrust/mass

		haccel = 0
		if accel > grav:
			haccel = math.sqrt(accel * accel - grav * grav)

		vel += haccel * timestep
		mass -= thrust/veff * timestep
		time += timestep

		if log:
			flightLog.append((time, vel, mass))

	return flightLog

accel = 2.53

starts = []
haccels = []

vstarts = []
vhaccels = []

bestScores = []
bestStarts = []
bestHaccels = []

bestStart = 0
bestHaccel = 0
bestScore = 1e9
particles = 100
for i in range(particles):
	start = random.random() * 2 * landingAlt
	haccel = random.random() * 2 * accel

	score, dist, vv, mass = testConstHaccelDescentScore(accel, start, haccel)
	if score < bestScore:
		bestStart = start
		bestHaccel = haccel
		bestScore = score

	starts.append(start)
	bestStarts.append(start)
	haccels.append(haccel)
	bestHaccels.append(haccel)
	bestScores.append(score)
	vstarts.append(random.gauss(0, 10))
	vhaccels.append(random.gauss(0, 0.1))
	print(bestScore)

i = 0
for i in range(0):
	for j in range(particles):
		vstarts[j] =\
			0.9 * vstarts[j] +\
			0.5 * random.random() * (bestStarts[j] - starts[j]) +\
			0.5 * random.random() * (bestStart - starts[j])
		vhaccels[j] =\
			0.9 * vhaccels[j] +\
			0.5 * random.random() * (bestHaccels[j] - haccels[j]) +\
			0.5 * random.random() * (bestHaccel - haccels[j])
		
		starts[j] += vstarts[j]
		if starts[j] < 0:
			starts[j] = 0
		haccels[j] += vhaccels[j]
		score, dist, vv, mass = testConstHaccelDescentScore(accel, starts[j], haccels[j])

		if score < bestScores[j]:
			bestScores[j] = score
			bestStarts[j] = starts[j]
			bestHaccels[j] = haccels[j]
		if score < bestScore:
			bestScore = score
			bestStart = starts[j]
			bestHaccel = haccels[j]
			print(i, bestScore, bestStart, bestHaccel, dist, vv, mass)
	i += 1
log = testConstHaccel1(accel, True, bestStart, bestHaccel)
for elem in log:
	print(elem)

n = 100
tests, results = [], []
for i in range(n):
	test = 5 * i/n
	if test < 0.1:
		continue
	mass = testConstAltitude(test, True)[-1][2]
	print(str(test)+", "+str(mass))

	tests.append(test)
	results.append(1-mass)

minUsed = 1 - math.exp(-2 * landingOrbVel/veff)
print(minUsed)

plt.plot(tests, results, label="Constant Altitude Landing")
plt.xlabel("Initial acceleration (m/s/s)")
plt.ylabel("Fraction of craft used as fuel")
plt.title("Fuel used by Mun landers with 320s Isp")

plt.plot(
	[2, 1.75, 1.5, 1.25], 
	[0.329368336792, 0.343041615525, 0.363318424963, 0.408181694272], 
	"rx", 
	label="Const. Acceleration Landing"
)

"""
plt.plot(
	[2, 1.5, 1],
	[0.355, 0.374, 0.572],
	"ro",
	label="Spring landing"
)
"""

plt.xlim(0, 5)
plt.ylim(0, 1)

plt.hlines(minUsed, 0, 5, colors="black", linestyles="dashed", label="Theoretical minimum")
plt.legend()

plt.show()
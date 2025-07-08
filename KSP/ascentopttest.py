import numpy as np
import matplotlib.pyplot as plt
import bisect
from copy import deepcopy

g0 = 9.80665
pi = 3.141592653589

generator = np.random.default_rng()

planetRadius = 600000
planetGees = 0.8
planetRotPeriod = 138984

gravParameter = planetRadius ** 2 * planetGees * g0

class ValueCurve:
	def __init__(self, xs=[], ys=[]):
		self.xs = xs
		self.ys = ys

	def evaluate(self, x):
		idx = bisect.bisect(self.xs, x) - 1

		if idx < 0:
			return self.ys[0]
		if idx >= len(self.ys):
			return self.ys[-1]

		return self.ys[idx]

	def sort(self):
		self.ys = [y for _, y in sorted(zip(self.xs, self.ys))]
		self.xs.sort()

	def mutate(self, xCutoff = -1, variancex = 0.1, variance = 0.1):
		for i in range(len(self.xs)):
			if self.xs[i] < xCutoff:
				continue
			self.xs[i] += variancex * generator.normal()
			self.ys[i] += variance * generator.normal()

		self.xs[0] = 0

		self.sort()

def randomValueCurve(n, xscale=10, yscale=1):
	xs = []
	ys = []

	currx = 0
	curry = yscale * generator.normal()
	for i in range(n):
		xs.append(currx)
		ys.append(curry)

		currx += xscale * generator.random()
		curry += 0.5 * yscale * generator.normal()

	return ValueCurve(xs, ys)

def simulateAscentConstantAltitude(thrust, isp, alt):
	vel = 2 * pi * (planetRadius + alt) / planetRotPeriod
	mass = 1

	fuelConsumption = thrust / (isp * g0)
	landingGravity = gravParameter / (planetRadius + alt) ** 2

	timestep = 0.1
	orbVel = np.sqrt(gravParameter / (planetRadius + alt))
	while vel < orbVel:
		accel = thrust / mass
		vaccel = landingGravity - vel * vel / (planetRadius + alt)

		if accel < vaccel:
			accel = vaccel
			mass = thrust / vaccel

		haccel = np.sqrt(accel * accel - vaccel * vaccel)

		mass -= fuelConsumption * timestep
		vel += haccel * timestep

	return mass

def simulateAscentCode(thrust, isp, alt, ax, ay):
	x = planetRadius + alt
	y = 0
	vx = 0
	vy = 2 * pi * (planetRadius + alt) / planetRotPeriod
	mass = 1

	timestep = 0.1
	time = 0

	fuelConsumption = thrust / (isp * g0)

	sma = 0
	eccentricity = 0
	periapsis = 1e-6
	while periapsis < planetRadius:
		dist = np.sqrt(x * x + y * y)

		vv = (x * vx + y * vy) / dist

		if dist < planetRadius or eccentricity > 1 or time > 1e4:
			#print(105, round(periapsis, 2), sma, round(eccentricity, 2), mass)
			return 1e-8 * (periapsis - planetRadius) + 1e-6 * mass, sma, periapsis, mass

		grav = gravParameter / dist / dist
		gx = -x * grav / dist
		gy = -y * grav / dist

		accel = max(0, ax.evaluate(sma))
		angle = ay.evaluate(sma)

		haccel = abs(accel * np.cos(angle))
		vaccel = accel * np.sin(angle)

		xaccel = (-y * haccel + x * vaccel) / dist
		yaccel = ( x * haccel + y * haccel) / dist

		if accel > thrust / mass:
			xaccel *= thrust / (mass * accel)
			yaccel *= thrust / (mass * accel)
			accel = thrust / mass

		vx += (xaccel + gx) * timestep
		vy += (yaccel + gy) * timestep
		x += vx * timestep
		y += vy * timestep
		mass -= accel * mass / (isp * g0) * timestep
		time += timestep

		vh = (-y * vx + x * vy) / dist

		sma = 1 / (2 / dist - (vx * vx + vy * vy) / gravParameter)
		eccentricity = np.sqrt(1 - (planetRadius * vh) ** 2 / (sma * gravParameter))
		periapsis = sma * (1 - eccentricity)

		#print(x, y, vx, vy, vh, vv, periapsis, angle % (2 * pi))

	apoapsis = sma * (1 + eccentricity)

	velPeriapsis = np.sqrt(gravParameter * (2 / periapsis - 1 / sma))
	velPeriapsisAfterBurn1 = np.sqrt(gravParameter * (2 / periapsis - 2 / (periapsis + planetRadius + alt)))

	dvUsed = velPeriapsis - velPeriapsisAfterBurn1

	massFraction = np.exp(-dvUsed / (isp * g0))

	mass *= massFraction

	print(144, mass)

	return mass, sma, periapsis, mass

def optimizeAscentCode(thrust, isp, alt):
	bestx = randomValueCurve(30, planetRadius / 15, 10)
	besty = randomValueCurve(30, planetRadius / 15, 1)
	bestScore, bestSMA, bestPeri, bestMass = simulateAscentCode(thrust, isp, alt, bestx, besty)

	oldx, oldy, oldScore, oldMass = deepcopy(bestx), deepcopy(besty), deepcopy(bestScore), deepcopy(bestMass)

	for i in range(1000000):
		currx = deepcopy(oldx)
		curry = deepcopy(oldy)

		if generator.random() < 0.9:
			currx.mutate(bestSMA - 150000, variancex=5000, variance=2 * generator.random())
			curry.mutate(bestSMA - 150000, variancex=5000, variance=0.2 * generator.random())
		else:
			currx.mutate(variancex=10000, variance=2 * generator.random())
			curry.mutate(variancex=10000, variance=0.2 * generator.random())

		score, sma, peri, mass = simulateAscentCode(thrust, isp, alt, currx, curry)

		if score >= bestScore:
			bestx = deepcopy(currx)
			besty = deepcopy(curry)
			bestScore = score
			bestSMA = sma
			bestPeri = peri
			bestMass = mass

			print([(round(x), round(y, 2)) for x, y in zip(bestx.xs, bestx.ys)])
			print([(round(x), round(y % (2 * pi), 2)) for x, y in zip(besty.xs, besty.ys)])
			print(i, bestScore, bestSMA, bestPeri, bestMass, isp * g0 * np.log(1 / bestMass))

		if score >= oldScore or generator.random() < 0.1:
			oldx = deepcopy(currx)
			oldy = deepcopy(curry)
			oldScore = score
			oldMass = mass

		if i % 100 == 0:
			oldx = deepcopy(bestx)
			oldy = deepcopy(besty)
			oldScore = bestScore
			oldMass = bestMass

print(optimizeAscentCode(10, 320, 5000))


twrs = np.linspace(0.5, 2, num=200)
dvs = []
minthrusts = []

alt = 5000
isp = 320

orbVel = np.sqrt(gravParameter / (planetRadius + alt))

for twr in twrs:
	thrust = twr * gravParameter / (planetRadius + alt) ** 2

	dvs.append(-isp * g0 * np.log(simulateAscentConstantAltitude(thrust, isp, alt)))
	minthrusts.append(-isp * g0 * np.log(twr))

plt.plot(twrs, dvs, 'b-', label="Constant Altitude Ascent")
plt.plot(twrs, minthrusts, 'b--', label="TWR = 1")
plt.hlines(orbVel - 2 * pi * (planetRadius + alt) / planetRotPeriod, min(twrs), max(twrs), colors="black", linestyles="dashed", label="Theoretical Limit")

plt.xlim(min(twrs), max(twrs))
plt.ylim(0, max(dvs))
plt.xlabel("Initial TWR")
plt.ylabel("dV used (m/s)")

plt.legend()

plt.show()
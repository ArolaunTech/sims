import math
import random
import matplotlib.pyplot as plt
from copy import deepcopy as dc

#All units in SI
#Consts
g0 = 9.80665

#Planet Stats
planetG = 1.63
planetR = 200000
rotationalPeriod = 138984.38

#Flight Stats
landingAltitude = 5000

#Craft Stats
wetMass = 500
ionThrust = 2000
ionIsp = 4200
maxEC = 100
genEC = 3.28

#Sim stats
timestep = 0.5

#Calculations
planetSTDGP = planetR * planetR * planetG

orbitVel = math.sqrt(planetSTDGP/(planetR + landingAltitude))
rotVel = 2 * math.pi * (planetR + landingAltitude)/rotationalPeriod
landingGrav = planetSTDGP/(planetR + landingAltitude)**2

ionFuel = ionThrust/(ionIsp * g0)
ionEC = ionFuel * 180

def calcIonThrust(ec, solar, throttle):
	#thrust, ec, xe
	demandEC = ionEC * throttle
	if (ec > 0) or (solar > demandEC): #Full EC
		return ionThrust * throttle, demandEC, ionFuel * throttle
	ecFraction = solar/demandEC
	if ecFraction > 0.1: #>10% EC starvation
		return ionThrust * throttle * (0.5 + 0.5 * ecFraction), solar, ionFuel * throttle
	return ionThrust * throttle * 5.5 * ecFraction, solar, ionFuel * throttle * 10 * ecFraction

def experimentalIonThrust(v, mass):
	grav = landingGrav - v * v/(planetR + landingAltitude)
	return min(1, max(genEC/ionEC, 4*grav*grav*ionEC*mass*mass/(genEC*ionThrust*ionThrust)-genEC/ionEC))

def simulate(code, log, experiment):
	landLog = []
	mass = wetMass
	ec = maxEC

	v = orbitVel

	l = len(code)
	time = 0
	throttle = code[0][1]
	nextSwitch = code[0][0]
	i = 0
	while v > rotVel:
		if (time > nextSwitch) and (i < l - 1):
			i += 1
			nextSwitch += code[i][0]
			throttle = code[i][1]
		if experiment:
			throttle = experimentalIonThrust(v, mass)

		grav = landingGrav - v * v/(planetR + landingAltitude)
		thrust, ecUsage, xeUsage = calcIonThrust(ec, genEC, throttle)
		accel = thrust/mass
		if accel < grav:
			return -1, ()
		haccel = math.sqrt(accel * accel - grav * grav)

		v -= haccel * timestep
		mass -= xeUsage * timestep
		if mass < 0:
			return -1, ()
		ec += (genEC - ecUsage) * timestep
		if ec < 0:
			ec = 0
		if ec > maxEC:
			ec = maxEC
		time += timestep
		if log:
			landLog.append([time, v, mass, ec, throttle, 0])
			print(time, v, mass, ec, throttle)

	v = rotVel
	ec = maxEC
	while v < orbitVel:
		if (time > nextSwitch) and (i < l - 1):
			i += 1
			nextSwitch += code[i][0]
			throttle = code[i][1]
		if experiment:
			throttle = experimentalIonThrust(v, mass)

		grav = landingGrav - v * v/(planetR + landingAltitude)
		thrust, ecUsage, xeUsage = calcIonThrust(ec, genEC, throttle)
		accel = thrust/mass
		if accel < grav:
			return -1, ()
		haccel = math.sqrt(accel * accel - grav * grav)

		v += haccel * timestep
		mass -= xeUsage * timestep
		if mass < 0:
			return -1, ()
		ec += (genEC - ecUsage) * timestep
		if ec < 0:
			ec = 0
		if ec > maxEC:
			ec = maxEC
		time += timestep
		if log:
			landLog.append([time, v, mass, ec, throttle, 1])
			print(time, v, mass, ec, throttle)
	return mass, tuple(landLog)

xMass, xLog = simulate([[10000,1]], True, True)
print(xMass)

bestCode = []
bestMass = -5000
bestLog = ()
for i in range(5000000):
	if (i == 0) or (random.random() < 0.1):
		code = []
		for j in range(random.randint(10, 20)):
			code.append([100*random.random(), random.randint(0,4)/4])
	elif (random.random() < 0.1):
		code = []
		time = 0
		for j in range(random.randint(20,30)):
			delta = 100 * random.random()
			if time > xLog[-1][0]:
				break
			code.append([delta, round(4*xLog[int(time/timestep)-1][4])/4])
			time += delta
	else:
		code = dc(bestCode)
		for j in range(len(code)):
			code[j][0] += random.gauss(mu=0, sigma=0.01)
			if code[j][0] < 0:
				code[j][0] = 0
			code[j][1] += random.gauss(mu=0, sigma=0.2)
			if code[j][1] < 0:
				code[j][1] = 0
			if code[j][1] > 1:
				code[j][1] = 1
			code[j][1] = round(4*code[j][1])/4
	m, log = simulate(code, False, False)
	if i%1000 == 0:
		print(i)
	if m < bestMass:
		continue
	bestCode = dc(code)
	bestMass = m
	print(bestCode)
	print(code == bestCode)
	print(m)
	print(i)

bestMass, bestLog = simulate(bestCode, True, False)
defMass, defLog = simulate([[10000, 1]], True, False)

saved = bestMass - defMass

print(saved, defMass, bestMass, xMass)
print(bestCode)

plt.plot([i[1] for i in bestLog], [i[4] for i in bestLog], label="optimizer (uses " + str(round(wetMass-bestMass,2)) + "kg Xe)")
plt.plot([i[1] for i in defLog],  [i[4] for i in defLog], label="full throttle (uses " + str(round(wetMass-defMass,2)) + "kg Xe)")
plt.plot([i[1] for i in defLog],  [experimentalIonThrust(i[1], i[2]) for i in defLog], label="experiment (uses " + str(round(wetMass-xMass, 2)) + "kg Xe)")

plt.legend()

plt.xlabel("Velocity (m/s)")
plt.ylabel("Throttle")

plt.title("Velocity-throttle plot")

plt.show()
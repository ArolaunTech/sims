import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy

g0 = 9.80665
timestep = 0.1
isp = 320

generator = np.random.default_rng()

veff = isp * g0

planetRadius = 6e5
planetGravity = 0.8 * g0
planetRotationalPeriod = 211926.35802123

planetGravParameter = planetRadius * planetRadius * planetGravity

def simulate(initDist, thrust, law, coeffs=[], logging=True):
	x = initDist
	y = 0
	vx = 0
	vy = 2 * np.pi * initDist / planetRotationalPeriod
	time = 0

	mass = 1

	log = []

	periapsis = 0

	valid = True
	while periapsis < planetRadius:
		dist = np.sqrt(x * x + y * y)
		vh = (-y * vx + x * vy) / dist
		vv = (x * vx + y * vy) / dist

		grav = planetGravParameter / dist / dist
		gx = -x * grav / dist
		gy = -y * grav / dist

		accel = thrust / mass

		if law == 0:
			#Constant Altitude control law
			effgrav = grav - vh * vh / dist
			haccel = np.sqrt(accel * accel - effgrav * effgrav)
			ax = haccel * gy / grav - gx * effgrav / grav
			ay = -haccel * gx / grav - gy * effgrav / grav
		elif law == 1:
			#Polynomial acceleration
			effgrav = grav - vh * vh / dist -\
				0.58025796 +\
				0.00501925044 * time +\
				2.556e-8 * time * time +\
				6.416e-9 * time * time * time
			if effgrav > accel:
				effgrav = accel
			if effgrav < -accel:
				effgrav = -accel

			haccel = np.sqrt(accel * accel - effgrav * effgrav)
			vaccel = effgrav

			ax = haccel * gy / grav - gx * vaccel / grav
			ay = -haccel * gx / grav - gy * vaccel / grav
		elif law == 2:
			#Const acceleration
			haccel = 5.66
			vaccel = np.sqrt(accel * accel - haccel * haccel)
			ax = haccel * gy / grav - gx * vaccel / grav
			ay = -haccel * gx / grav - gy * vaccel / grav

			if vv > 0:
				effgrav = grav - vh * vh / dist
				haccel = np.sqrt(accel * accel - effgrav * effgrav)
				ax = haccel * gy / grav - gx * effgrav / grav
				ay = -haccel * gx / grav - gy * effgrav / grav
		elif law == 3:
			#Const vertical acceleration
			effgrav = 7.101
			if vv > 0:
				effgrav = grav - vh * vh / dist

			haccel = np.sqrt(accel * accel - effgrav * effgrav)
			ax = haccel * gy / grav - gx * effgrav / grav
			ay = -haccel * gx / grav - gy * effgrav / grav
		elif law == 4:
			#Sproing
			effgrav = grav - vh * vh / dist - (dist - 602500) * 19.6e-5
			if effgrav > accel:
				effgrav = accel
			if effgrav < -accel:
				effgrav = -accel

			haccel = np.sqrt(accel * accel - effgrav * effgrav)
			ax = haccel * gy / grav - gx * effgrav / grav
			ay = -haccel * gx / grav - gy * effgrav / grav

		vx += (ax + gx) * timestep
		vy += (ay + gy) * timestep
		x += vx * timestep
		y += vy * timestep
		mass -= thrust * timestep / veff
		time += timestep

		log.append([time, x, y, vx, vy, vh, vv, ax, ay])
		#print(x, y, vx, vy)

		vh = (-y * vx + x * vy) / dist
		sma = 1 / (2 / dist - (vx * vx + vy * vy) / planetGravParameter)
		eccentricity = np.sqrt(1 - (dist * vh) ** 2 / (sma * planetGravParameter))
		periapsis = sma * (1 - eccentricity)

		if mass <= 0:
			if logging:
				print(law)
			valid = False
			break
		if eccentricity > 1:
			if logging:
				print(law)
			valid = False
			break
		if dist < planetRadius:
			if logging:
				print(law)
			valid = False
			break

	if not valid:
		return log, 0

	#mass += (periapsis - planetRadius) / 5000000

	return log, mass

#Optimize polynomial
cs = [0.416169, -1.067822e-4, -1.550779e-8, -5.152017e-13]

logs = []
masses = []
for i in range(5):
	log, mass = simulate(605000, 8, i, cs)
	logs.append(log)
	masses.append(mass)

	print(f"Law #{i}: End mass = {round(100 * mass, 3)}%")

print("")
print(cs)
#for logEntry in log:
	#print(logEntry)

#plt.plot([i[1] for i in log], [i[2] for i in log])
for idx, log in enumerate(logs):
	plt.plot([i[0] for i in log], [np.sqrt(i[1] * i[1] + i[2] * i[2]) for i in log], label=f"Law #{idx}: {round(100 * masses[idx], 3)}%")
#plt.gca().set_aspect('equal')
#plt.plot([i[0] for i in log], [i[2] for i in log])

plt.legend()
plt.show()
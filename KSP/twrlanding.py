import numpy as np
import matplotlib.pyplot as plt

#Constants
g = 9.81

#Planet
planetRadius = 600000
planetGrav = 0.8 * g
planetRotationalPeriod = 211926.36

planetSTDGP = planetRadius * planetRadius * planetGrav

#Craft
landingAltitude = 5000
Isp = 305

#Simulation
timestep = 0.01

#Calculations
orbVel = np.sqrt(planetSTDGP/(planetRadius + landingAltitude))
rotVel = 2 * np.pi * (planetRadius + landingAltitude)/planetRotationalPeriod
landingGrav = planetSTDGP/(planetRadius + landingAltitude)**2

print(orbVel, rotVel, landingGrav)

def calc_dv(acceleration):
	mass = 1
	fuelConsumption = acceleration/(Isp * g)
	currvel = orbVel
	while currvel > rotVel:
		accel = acceleration/mass
		vaccel = (1 - (currvel/orbVel)**2) * landingGrav
		if vaccel > accel:
			return -1

		haccel = np.sqrt(accel * accel - vaccel * vaccel)

		mass -= fuelConsumption * timestep
		currvel -= haccel * timestep

	currvel = rotVel
	while currvel < orbVel:
		accel = acceleration/mass
		vaccel = (1 - (currvel/orbVel)**2) * landingGrav
		if vaccel > accel:
			return -1

		haccel = np.sqrt(accel * accel - vaccel * vaccel)

		mass -= fuelConsumption * timestep
		currvel += haccel * timestep

	return -Isp * g * np.log(mass)

xs = np.linspace(0.1, 100, num=1000)
ys = []
i = 0
for x in xs:
	ys.append(calc_dv(x) - 2 * orbVel + 2 * rotVel)
	i += 1
	print(str(x) + ", " + str(ys[-1]))

plt.plot(xs, ys)
plt.show()
import numpy as np
from itertools import permutations as permute

maxAngles = np.array([90, 32, 37, 11, 10])
maxAngles = np.sin(maxAngles*np.pi/180)

def calcVacuumAscent(accel, Isp, rad, stdgp, lalt, oalt, rp, chemIsp, chemThrust, nervisp, nervthrust, moonindex):
	currvel = 2*np.pi*(rad+lalt)/rp
	orbvel = np.sqrt(stdgp/(rad+lalt))
	necessary = np.sqrt(stdgp*(2/(rad+lalt)-2/(2*rad+lalt+oalt)))
	grav = stdgp/(rad+lalt)**2

	out = 0
	mass = 1
	timestep = 0.5
	xenonUsed = 0
	lfUsed = 0
	lfoxUsed = 0
	while currvel < necessary:
		orbfrac = currvel/orbvel
		vthrust = grav * (1 - orbfrac*orbfrac)
		baccel = accel
		if accel/mass * maxAngles[moonindex] < vthrust:
			baccel = vthrust*mass/maxAngles[moonindex]
		caccel = baccel - accel

		daccel = 0
		if caccel > nervthrust:
			daccel = caccel - nervthrust
			caccel = nervthrust
		if daccel > chemThrust:
			return None
		hthrust = np.sqrt((baccel/mass)**2 - (vthrust)**2)
		currvel += hthrust * timestep

		xenonUsed += accel * timestep/(Isp * 9.81)
		lfUsed += caccel * timestep/(nervisp * 9.81)
		lfoxUsed += daccel * timestep/(chemIsp * 9.81)

		mass -= accel * timestep/(Isp * 9.81)
		mass -= caccel * timestep/(nervisp * 9.81)
		mass -= daccel * timestep/(chemIsp * 9.81)
		out += baccel * timestep/(mass)
		#print(mass, currvel)
	return out, mass, xenonUsed, lfUsed, lfoxUsed

def calcVacuumLanding(accel, Isp, rad, stdgp, lalt, oalt, rp, chemIsp, chemThrust, nervisp, nervthrust, moonindex):
	rotvel = 2*np.pi*(rad+lalt)/rp
	orbvel = np.sqrt(stdgp/(rad+lalt))
	currvel = np.sqrt(stdgp*(2/(rad+lalt)-2/(2*rad+lalt+oalt)))
	grav = stdgp/(rad+lalt)**2

	out = 0
	mass = 1
	timestep = 0.5
	xenonUsed = 0
	lfUsed = 0
	lfoxUsed = 0
	while currvel > rotvel:
		orbfrac = currvel/orbvel
		vthrust = grav * (1 - orbfrac*orbfrac)
		baccel = accel
		if accel/mass * maxAngles[moonindex] < vthrust:
			baccel = vthrust*mass/maxAngles[moonindex]
		caccel = baccel - accel

		daccel = 0
		if caccel > nervthrust:
			daccel = caccel - nervthrust
			caccel = nervthrust
		if daccel > chemThrust:
			return None
		hthrust = np.sqrt((baccel/mass)**2 - (vthrust)**2)
		currvel -= hthrust * timestep

		xenonUsed += accel * timestep/(Isp * 9.81)
		lfUsed += caccel * timestep/(nervisp * 9.81)
		lfoxUsed += daccel * timestep/(chemIsp * 9.81)

		mass -= accel * timestep/(Isp * 9.81)
		mass -= caccel * timestep/(nervisp * 9.81)
		mass -= daccel * timestep/(chemIsp * 9.81)
		out += baccel * timestep/(mass)
		#print(mass, currvel)
	return out, mass, xenonUsed, lfUsed, lfoxUsed

def lowToHighOrbit(moonIndex):
	return np.sqrt(stdgps[moonIndex]*(2/(radii[moonIndex]+lowOrbAlt[moonIndex]) - 2/(radii[moonIndex] + lowOrbAlt[moonIndex] + sois[moonIndex]))) - np.sqrt(stdgps[moonIndex]/(radii[moonIndex] + lowOrbAlt[moonIndex]))

moons = ['Laythe', 'Vall', 'Tylo', 'Bop', 'Pol']
radii = [500000, 300000, 600000, 65000, 44000]
grav =  [0.8, 0.235, 0.8, 0.06, 0.038]
rotPeriod = [52980.879, 105962.09, 211926.36, 544507.43, 901902.62]
orbAlt = [27184000, 43152000, 68500000, 128500000, 179890000]
lowOrbAlt = [55000, 10000, 15000, 25000, 10000]
landingAlt = [0, 7000, 8000, 22000, 5000]
lth = []

primarystdgp = 0.8*9.81*3.6e+13
sois = []
for i in range(len(grav)):
	grav[i] *= 9.81
stdgps = []
for i in range(len(grav)):
	stdgps.append(grav[i]*radii[i]*radii[i])

for i in range(len(grav)):
	sois.append(orbAlt[i]*(((stdgps[i])/(primarystdgp))**0.4))

for i in range(len(grav)):
	lth.append(lowToHighOrbit(i))

print(moons, '\n', radii, '\n', grav, '\n', stdgps, '\n', sois, '\n', lth)

locations = []
for moon in moons:
	locations.append(moon)
	locations.append('Low ' + moon + ' Orbit')
	locations.append('High ' + moon + ' Orbit')

print(locations)

#print(calcVacuumLanding(1, 4200, 65000, stdgps[3], 22000, 25000, rotPeriod[3], 412))

startingMass = 1
print('startingMass: ' + str(startingMass))
joolDv = 900
rocketIsp = 412
rocketThrust = 7
ionThrust = 0.001
NERVIsp = 800
NERVThrust = 2

startingMass /= np.exp(joolDv/(4200*9.81))
startingXenon = 1 - startingMass
scores = [-1 for i in range(120)]
xUse = [-1 for i in range(120)]
lUse = [-1 for i in range(120)]
oUse = [-1 for i in range(120)]
print('mass at jool: ' + str(startingMass))
for j, order in enumerate(permute(moons)):
	currMass = startingMass

	currLF = 0.5 #Proxima suggestion
	currLFOX = 0

	if order[0] == 'Bop':
		currMass /= np.exp(50/(4200*9.81))
	elif order[0] == 'Pol':
		currMass /= np.exp(140/(4200*9.81))

	currXenon = startingXenon + startingMass - currMass
	#print(currXenon, currLF, currLFOX)
	#print('mass at ' + order[0] + ' (High Orbit): ' + str(currMass))
	failed = False
	for i, destination in enumerate(order):
		if destination != 'Laythe':
			currXenon += currMass - currMass/np.exp(lth[moons.index(destination)]/(4200*9.81))
			currMass /= np.exp(lth[moons.index(destination)]/(4200*9.81))
		#print('mass at ' + destination + ' (Low Orbit): ' + str(currMass))
		moonIndex = moons.index(destination)
		if destination == 'Bop':
			#print(calcVacuumLanding(ionThrust/currMass, 4200, radii[3], stdgps[3], 22000, 25000, rotPeriod[3]))
			landingRes = calcVacuumLanding(ionThrust/currMass, 4200, radii[3], stdgps[3], 22000, 25000, rotPeriod[3], rocketIsp, rocketThrust/currMass, NERVIsp, NERVThrust/currMass, moonIndex)
			if landingRes == None:
				failed = True
				break

			currXenon += landingRes[2]*currMass
			currLF += landingRes[3]*currMass
			currLFOX += landingRes[4]*currMass

			currMass *= landingRes[1]
		elif destination == 'Pol':
			landingRes = calcVacuumLanding(ionThrust/currMass, 4200, radii[4], stdgps[4], 5000, 10000, rotPeriod[4], rocketIsp, rocketThrust/currMass, NERVIsp, NERVThrust/currMass, moonIndex)
			if landingRes == None:
				failed = True
				break

			currXenon += landingRes[2]*currMass
			currLF += landingRes[3]*currMass
			currLFOX += landingRes[4]*currMass

			currMass *= landingRes[1]
		elif destination != 'Laythe':
			moonIndex = moons.index(destination)
			landingRes = calcVacuumLanding(ionThrust/currMass, 4200, radii[moonIndex], stdgps[moonIndex], landingAlt[moonIndex], lowOrbAlt[moonIndex], rotPeriod[moonIndex], rocketIsp, rocketThrust/currMass, NERVIsp, NERVThrust/currMass, moonIndex) 
			if landingRes == None:
				failed = True
				break

			currXenon += landingRes[2]*currMass
			currLF += landingRes[3]*currMass
			currLFOX += landingRes[4]*currMass

			currMass *= landingRes[1]
		#print('mass on surface: ' + str(currMass), destination)
		if destination == 'Bop':
			#print(calcVacuumLanding(ionThrust/currMass, 4200, radii[3], stdgps[3], 22000, 25000, rotPeriod[3]))
			landingRes = calcVacuumAscent(ionThrust/currMass, 4200, radii[3], stdgps[3], 22000, 25000, rotPeriod[3], rocketIsp, rocketThrust/currMass, NERVIsp, NERVThrust/currMass, moonIndex)
			if landingRes == None:
				failed = True
				break
			currXenon += landingRes[2]*currMass
			currLF += landingRes[3]*currMass
			currLFOX += landingRes[4]*currMass

			currMass *= landingRes[1]
		elif destination == 'Pol':
			landingRes = calcVacuumAscent(ionThrust/currMass, 4200, radii[4], stdgps[4], 5000, 10000, rotPeriod[4], rocketIsp, rocketThrust/currMass, NERVIsp, NERVThrust/currMass, moonIndex)
			if landingRes == None:
				failed = True
				break
			currXenon += landingRes[2]*currMass
			currLF += landingRes[3]*currMass
			currLFOX += landingRes[4]*currMass

			currMass *= landingRes[1]
		elif destination == 'Laythe':
			calcMass = currMass
			currMass *= (1-0.0694) #LF portion of ascent with RAPIERS
			currMass /= np.exp(192/(NERVIsp*9.81))

			currLF += calcMass - currMass
		else:
			moonIndex = moons.index(destination)
			landingRes = calcVacuumAscent(ionThrust/currMass, 4200, radii[moonIndex], stdgps[moonIndex], landingAlt[moonIndex], lowOrbAlt[moonIndex], rotPeriod[moonIndex], rocketIsp, rocketThrust/currMass, NERVIsp, NERVThrust/currMass, moonIndex)
			if landingRes == None:
				failed = True
				break
			currXenon += landingRes[2]*currMass
			currLF += landingRes[3]*currMass
			currLFOX += landingRes[4]*currMass
			currMass *= landingRes[1]
		#print('mass in Low Orbit: ' + str(currMass))

		currXenon += currMass - currMass/np.exp(lth[moons.index(destination)]/(4200*9.81))
		currMass /= np.exp(lth[moons.index(destination)]/(4200*9.81))
		#print('mass in High Orbit: ' + str(currMass))

		#Calculate transfer costs
		moonIndex = moons.index(destination)
		returncosts = [0, 0, 0, 50, 140]
		currXenon += currMass - currMass/np.exp(returncosts[moonIndex]/(4200*9.81))
		currMass /= np.exp(returncosts[moonIndex]/(4200*9.81))

		if i < 4:
			nextDestination = order[i+1]
			nextMoonIndex = moons.index(nextDestination)
			currXenon += currMass - currMass/np.exp(returncosts[nextMoonIndex]/(4200*9.81))
			currMass /= np.exp(returncosts[nextMoonIndex]/(4200*9.81))
	if not failed:
		scores[j] = currMass
		xUse[j] = currXenon
		lUse[j] = currLF
		oUse[j] = currLFOX
	print('Final Mass: ' + str(currMass))
	#print(currXenon, currLF, currLFOX)
	#print('\n')

print(scores)
m = max(scores)
p = permute(moons)
for o, s, x, l, ox in zip(p, scores, xUse, lUse, oUse):
	if s >= 0.995*m:
		print(o, s, x/3, l/10, ox/8, s-x/3-l/10-ox/8)
print(m)
m = max(scores)
for o, s, x, l, ox in zip(p, scores, xUse, lUse, oUse):
	if s >= 0.999999*m:
		print(o, s, x, l, ox)
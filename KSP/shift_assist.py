import numpy as np
import matplotlib.pyplot as plt

#Functions
def eccentricityNewton(x, pEjection):
	return 2 * np.sqrt(1 + x*np.cos(pEjection)) + (1 - x*x)/(1 + x*np.cos(pEjection)) - (3 - (vRelative/planetVelocity)**2)

def eccentricityNewtonDerivative(x, pEjection):
	return np.cos(pEjection)/np.sqrt(1 + x*np.cos(pEjection)) + (x*x-1)*(np.cos(pEjection)/(1+x*np.cos(pEjection))**2) - 2*x/(1+x*np.cos(pEjection))

def calcOrbitParameters(startAngles):
	pEjection = np.pi + 0.5*startAngles
	eEjection = 0.5

	for i in range(10):
		eEjection -= eccentricityNewton(eEjection, pEjection)/eccentricityNewtonDerivative(eEjection, pEjection) #Newton's method

	aEjection = planetAltitude*(1 + eEjection*np.cos(pEjection))/(1 - eEjection*eEjection)
	return aEjection, eEjection, pEjection

#Constants
G = 6.67408e-11

#Solar parameters
solarRadius = 261600000
solarGravASL = 1.7
solarMass = 9.81*solarGravASL*solarRadius**2/G

#Planet parameters
planetAltitude = 9832684544
planetVelocity = (G*solarMass/planetAltitude)**0.5
planetOrbitalPeriod = 2*np.pi*planetAltitude/planetVelocity

#Simulation parameters
vRelative = 3000
maxPlanetOrbits = 10
targetShiftSum = 20 * np.pi/180
maxShifts = 5

#Generate allowed phase shift angles
possibleDelta = np.array([])
for numPlanetOrbits in range(maxPlanetOrbits):

	startAngles = np.linspace(0, 4*np.pi, num=100)

	for i in range(20):
		planetStartDelay = (startAngles)%(2*np.pi) * planetOrbitalPeriod/(2*np.pi)
		aEjection, eEjection, pEjection = calcOrbitParameters(startAngles)
		PEjection = 2*np.pi*aEjection**(3/2)/np.sqrt(G*solarMass)

		EEjection = np.arccos((1-planetAltitude/aEjection)/eEjection)
		delayEjection = (np.where(np.sin(pEjection)>0, 1, -1) * 2 * (EEjection - eEjection*np.sin(EEjection)))%(2*np.pi)
		delayEjection /= 2*np.pi
		delayEjection *= PEjection

		delayPlanet = planetStartDelay + numPlanetOrbits*planetOrbitalPeriod
		idealEjectionOrbits = (delayPlanet - delayEjection)/PEjection
		idealEjectionOrbits = np.where(idealEjectionOrbits < 0, 0, np.round(idealEjectionOrbits))
		delayEjection = idealEjectionOrbits * PEjection + delayEjection
		delayEjection -= delayPlanet
		delayEjection1 = delayEjection

		startAngles += 0.001
		planetStartDelay = (startAngles)%(2*np.pi) * planetOrbitalPeriod/(2*np.pi)
		aEjection, eEjection, pEjection = calcOrbitParameters(startAngles)
		PEjection = 2*np.pi*aEjection**(3/2)/np.sqrt(G*solarMass)

		EEjection = np.arccos((1-planetAltitude/aEjection)/eEjection)
		delayEjection = (np.where(np.sin(pEjection)>0, 1, -1) * 2 * (EEjection - eEjection*np.sin(EEjection)))%(2*np.pi)
		delayEjection /= 2*np.pi
		delayEjection *= PEjection

		delayPlanet = planetStartDelay + numPlanetOrbits*planetOrbitalPeriod
		idealEjectionOrbits = (delayPlanet - delayEjection)/PEjection
		idealEjectionOrbits = np.where(idealEjectionOrbits < 0, 0, np.round(idealEjectionOrbits))
		delayEjection = idealEjectionOrbits * PEjection + delayEjection
		delayEjection -= delayPlanet
		delayEjection -= delayEjection1

		delayEjection /= 0.001

		startAngles -= 0.001
		startAngles -= delayEjection1/delayEjection
	startAngles %= 4*np.pi
	startAngles = np.sort(startAngles)
	startAngles = np.unique(startAngles.round(decimals=6))

	possibleDelta = np.concatenate((possibleDelta, startAngles))

possibleDelta = np.sort(possibleDelta)
possibleDelta = np.unique(possibleDelta.round(decimals=6))
print(len(possibleDelta))
print(possibleDelta)
print(targetShiftSum)

#print(aEjection, eEjection, pEjection)
#print(PEjection, PEjection/planetOrbitalPeriod)
#print(delayEjection, delayEjection1)

#print(startAngles)

#Optimize the targetShiftSum
minError = 10
minErrorChoice = np.array([])
for i in range(10000):
	choice = np.random.choice(possibleDelta, size=(1+np.random.randint(maxShifts)))
	error = np.fabs(np.sum(choice)%(2*np.pi) - targetShiftSum)
	if error < minError:
		minError = error
		minErrorChoice = np.copy(choice)
print(minError)
print(minErrorChoice)
aEjection, eEjection, pEjection = calcOrbitParameters(minErrorChoice)
prefix = 0
for i in range(len(aEjection)):
	pEjection[i] += prefix
	prefix += minErrorChoice[i]

#Plot
fig, ax = plt.subplots()

sun = plt.Circle((0,0), solarRadius, color='yellow')
ax.set_aspect(1)
orbitAngles = np.linspace(0, 2*np.pi, num=100)
plt.plot(planetAltitude * np.cos(orbitAngles), planetAltitude * np.sin(orbitAngles))

for a, e, p in zip(aEjection, eEjection, pEjection):
	orbitAngles = np.linspace(0, 2*np.pi, num=100)
	ejectionRadius = a*(1-e*e)/(1+e*np.cos(orbitAngles-p))
	plt.plot(ejectionRadius * np.cos(orbitAngles), ejectionRadius * np.sin(orbitAngles))
ax.add_patch(sun)
#plt.plot(startAngles, planetStartDelay)
plt.show()
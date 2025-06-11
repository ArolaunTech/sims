import matplotlib.pyplot as plt
import math
import numpy as np
import bisect
from copy import deepcopy

from utils import *

KERBIN_STDGP = 3531600000000
KERBIN_ROT = 21549.425
KERBIN_ADIABATIC = 1.4
KERBIN_MOLARMASS = 0.0289644002914429
KERBIN_RADIUS = 6e5

VALL_STDGP = 207481499473.751
VALL_ROT = 105962.088893924
VALL_RADIUS = 3e5

TYLO_STDGP = 2825280042099.95
TYLO_ROT = 211926.35802123
TYLO_RADIUS = 6e5

R = 8.31446261815324
phi = (1 + np.sqrt(5))/2
generator = np.random.default_rng()

g0 = 9.80665

timestep = 0.5

kerbinPressureKeys = [ #Kerbin pressure keys from Kittopia Dumps, x, y, y', y'
	[    0.0,   101.325,        0.0,          -0.01501631  ],
	[ 1241.025,  84.02916,     -0.01289846,   -0.01289826  ],
	[ 2439.593,  69.68138,     -0.01107876,   -0.01107859  ],
	[ 3597.11,   57.78001,     -0.009515483,  -0.009515338 ],
	[ 4714.942,  47.90862,     -0.00817254,   -0.008172415 ],
	[ 5794.409,  39.72148,     -0.00701892,   -0.007018813 ],
	[ 6836.791,  32.93169,     -0.006027969,  -0.006027877 ],
	[ 7843.328,  27.30109,     -0.005176778,  -0.0051767   ],
	[ 8815.22,   22.63206,     -0.004445662,  -0.004445578 ],
	[10786.42,   15.3684,      -0.003016528,  -0.00301646  ],
	[12101.4,    11.87313,     -0.002329273,  -0.00232922  ],
	[13417.05,    9.172798,    -0.001798594,  -0.001798554 ],
	[16678.47,    4.842261,    -0.0009448537, -0.0009448319],
	[21143.1,     2.050097,    -0.0003894095, -0.0003894005],
	[26977.92,    0.6905929,   -0.0001252565, -0.0001252534],
	[33593.82,    0.2201734,   -3.626878e-5,  -3.626788e-5 ],
	[42081.87,    0.05768469,  -9.063159e-6,  -9.062975e-6 ],
	[49312.13,    0.01753794,  -3.029397e-6,  -3.029335e-6 ],
	[56669.95,    0.004591824, -8.827175e-7,  -8.826996e-7 ],
	[62300.84,    0.001497072, -3.077091e-7,  -3.077031e-7 ],
	[70000.0,     0.0,          0.0,           0.0         ]
]

kerbinTempBaseKeys = [
	[    0.0,  288.15,   0.0,         -0.008125   ],
	[ 8815.22, 216.65,  -0.008096968,  0.0        ],
	[16050.39, 216.65,   0.0,          0.001242164],
	[25729.23, 228.65,   0.001237475,  0.003464929],
	[37879.44, 270.65,   0.00344855,   0.0        ],
	[41129.24, 270.65,   0.0,         -0.003444189],
	[57440.13, 214.65,  -0.003422425, -0.002444589],
	[68797.88, 186.946, -0.002433851,  0.0        ],
	[70000.0,  186.946,  0.0,          0.0        ]
]

kerbinSunMultKeys = [
	[    0.00, 1.0,  0.0,          0.0         ],
	[ 8815.22, 0.3, -5.91316e-5,  -5.91316e-5  ],
	[16050.39, 0.0,  0.0,          0.0         ],
	[25729.23, 0.0,  0.0,          0.0         ],
	[37879.44, 0.2,  0.0,          0.0         ],
	[57440.13, 0.2,  0.0,          0.0         ],
	[63902.72, 1.0,  0.0001012837, 0.0001012837],
	[70000.00, 1.2,  0.0,          0.0         ]
]

kerbinLatitudeBiasKeys = [
	[ 0.0,  17.0,      0.0,       -0.3316494],
	[10.0,  12.0,     -0.65,      -0.65     ],
	[18.0,   6.36371, -0.4502313, -0.4502313],
	[30.0,   0.0,     -1.3,       -1.3      ],
	[35.0, -10.0,     -1.65,      -1.65     ],
	[45.0, -23.0,     -1.05,      -1.05     ],
	[55.0, -31.0,     -0.6,       -0.6      ],
	[70.0, -37.0,     -0.6689383, -0.6689383],
	[90.0, -50.0,     -0.02418368, 0.0      ]
]

kerbinLatitudeSunMultKeys = [
	[ 0.0,  9.0,       0.0,          0.1554984  ],
	[40.0, 14.2,       0.08154097,   0.08154097 ],
	[55.0, 14.9,      -0.006055089, -0.006055089],
	[68.0, 12.16518,  -0.2710912,   -0.2710912  ],
	[76.0,  8.582909, -0.6021729,   -0.6021729  ],
	[90.0,  5.0,       0.0,          0.0        ]
]

rapierPressureKeys = [
	[0.0,   0.0,  0.0,       0.0],
	[0.018, 0.09, 7.914787,  7.914787],
	[0.08,  0.3,  1.051923,  1.051923],
	[0.35,  0.5,  0.3927226, 0.3927226],
	[1.0,   1.0,  1.055097,  0.0]
]

rapierMachKeys = [
	[0.0,  1.0,  0.0,        0.08333334],
	[0.2,  0.98, 0.42074,    0.42074   ],
	[0.7,  1.8,  2.290406,   2.290406  ],
	[1.4,  4.0,  3.887193,   3.887193  ],
	[3.75, 8.5,  0.0,        0.0       ],
	[4.5,  7.3, -2.831749,  -2.831749  ],
	[5.5,  3.0, -5.260566,  -5.260566  ],
	[6.0,  0.0, -0.02420209, 0.0       ]
]

nervIspKeys = [
	[0.0, 800.0,  0.0, -615.0],
	[1.0, 185.0, -400, -400.0],
	[2.0, 0.0,   -185,    0.0]
]

liftAoAKeys = [ #Thanks to Lt_Duckweed
	[0.0,       0.0,        0.0,        1.965926 ],
	[0.258819,  0.5114774,  1.990092,   1.905806 ],
	[0.5,       0.9026583,  0.7074468, -0.7074468],
	[0.7071068, 0.5926583, -2.087948,  -1.990095 ],
	[1.0,       0.0,       -2.014386,  -2.014386 ]
]

dragAoAKeys = [
	[0.0,       0.01, 0.0,       0.0      ],
	[0.3420201, 0.06, 0.1750731, 0.1750731],
	[0.5,       0.24, 2.60928,   2.60928  ],
	[0.7071068, 1.7,  3.349777,  3.349777 ],
	[1.0,       2.4,  1.387938,  0.0      ]
]

liftMachKeys = [
	[ 0.0, 1.0,    0.0,           0.0       ],
	[ 0.3, 0.5,   -1.671345,     -0.8273422 ],
	[ 1.0, 0.125, -0.0005291355, -0.02625772],
	[ 5.0, 0.0625, 0.0,           0.0       ],
	[25.0, 0.05,   0.0,           0.0       ]
]

dragMachKeys = [
	[ 0.0,  0.35,  0.0,         -0.8463008 ],
	[ 0.15, 0.125, 0.0,          0.0       ],
	[ 0.9,  0.275, 0.541598,     0.541598  ],
	[ 1.1,  0.75,  0.0,          0.0       ],
	[ 1.4,  0.4,  -0.3626955,   -0.3626955 ],
	[ 1.6,  0.35, -0.1545923,   -0.1545923 ],
	[ 2.0,  0.3,  -0.09013031,  -0.09013031],
	[ 5.0,  0.22,  0.0,          0.0       ],
	[25.0,  0.3,   0.0006807274, 0.0       ]
]

class ValueCurve:
	def __init__(self, xs, ys):
		self.xs = xs
		self.ys = ys

	def evaluate(self, x):
		return self.ys[bisect.bisect(self.xs, x) - 1]

	def sort(self):
		sortedKeys = sorted(zip(self.xs, self.ys))
		self.xs = [x for x,_ in sortedKeys]
		self.ys = [y for _,y in sortedKeys]

	def __repr__(self):
		return "ValCurve {" + str(self.xs) + ", " + str(self.ys) + "}"

def readHermiteSpline(x, spline):
	minX = spline[0][0]
	maxX = spline[-1][0]
	if x <= minX:
		return spline[0][1] + (x - minX) * spline[0][2]
	if x >= maxX:
		return spline[-1][1] + (x - maxX) * spline[-1][3]

	for i in range(len(spline)):
		if x <= spline[i + 1][0]:
			break

	lx = spline[i][0]
	ux = spline[i + 1][0]

	ly = spline[i][1]
	uy = spline[i + 1][1]

	ld = spline[i][3]
	ud = spline[i + 1][2]

	interval = ux - lx
	t = (x - lx)/interval

	t2 = t * t
	t3 = t2 * t

	h0 =  2.0 * t3 - 3.0 * t2 + 1.0
	h1 =        t3 - 2.0 * t2 + t
	h2 =        t3 -       t2
	h3 = -2.0 * t3 + 3.0 * t2

	h1 *= interval
	h2 *= interval

	return h0 * ly + h1 * ld + h2 * ud + h3 * uy

def pressureAlt(alt):
	return readHermiteSpline(alt, kerbinPressureKeys)

def temperatureAlt(alt):
	return readHermiteSpline(alt, kerbinTempBaseKeys) + readHermiteSpline(alt, kerbinSunMultKeys) * readHermiteSpline(0, kerbinLatitudeBiasKeys)

densityASL = 1000 * pressureAlt(0)/temperatureAlt(0) * KERBIN_MOLARMASS/R

def flowcap(x):
	if x <= 3:
		return x
	return 6 * x/(x + 3)

def RAPIERThrust(mach, alt):
	if alt > 27000:
		return 0
	density = 1000 * pressureAlt(alt)/temperatureAlt(alt) * KERBIN_MOLARMASS/R
	pressureMult = readHermiteSpline(density/densityASL, rapierPressureKeys)
	machmult = readHermiteSpline(mach, rapierMachKeys)

	return 105000.0 * flowcap(pressureMult * machmult)

def nervIsp(alt):
	return readHermiteSpline(pressureAlt(alt)/101.325, nervIspKeys)

def calculateLDRatio(mach, aoa):
	liftAoAMult = readHermiteSpline(np.sin(aoa), liftAoAKeys)
	dragAoAMult = readHermiteSpline(np.sin(aoa), dragAoAKeys)
	liftMachMult = readHermiteSpline(mach, liftMachKeys)
	dragMachMult = readHermiteSpline(mach, dragMachKeys)

	lift = liftAoAMult * liftMachMult
	drag = dragAoAMult * dragMachMult/2.4

	truelift = lift * np.cos(aoa)
	lid = lift * np.sin(aoa)

	return truelift/(drag + lid)

def findBestLD(mach):
	xmin = np.pi/180
	xmax = np.pi/36
	while xmax - xmin > 1e-3:
		xi1 = xmin/phi + xmax*(2-phi)
		xi2 = xmax/phi + xmin*(2-phi)

		vals = [calculateLDRatio(mach, xmin), calculateLDRatio(mach, xi1), calculateLDRatio(mach, xi2), calculateLDRatio(mach, xmax)]
		maxld = max(vals)

		if vals[0] == maxld:
			xmax = xi1
		elif vals[1] == maxld:
			xmax = xi2
		elif vals[2] == maxld:
			xmin = xi1
		else:
			xmin = xi2

	return max(vals)

machs = []
ldratios = []
for i in np.linspace(0, 25, num=2500):
	machs.append(i)
	ldratios.append(findBestLD(i))

ldratio = ValueCurve(machs, ldratios)

alt = 16500
vallAlt = 7500
tyloAlt = 11000
density = 1000 * pressureAlt(alt)/temperatureAlt(alt) * KERBIN_MOLARMASS/R
soundspeed = np.sqrt(KERBIN_ADIABATIC * R * temperatureAlt(alt)/KERBIN_MOLARMASS)

gravity = KERBIN_STDGP/(KERBIN_RADIUS + alt)**2
vallGravity = VALL_STDGP/(VALL_RADIUS + vallAlt)**2
tyloGravity = TYLO_STDGP/(TYLO_RADIUS + tyloAlt)**2

vallrotVel = 2 * np.pi * (VALL_RADIUS + vallAlt)/VALL_ROT
tylorotVel = 2 * np.pi * (TYLO_RADIUS + tyloAlt)/TYLO_ROT

def calcThrustFuel(engines, throttles, rapierclosedcycle, drag, mach, alt):
	thrust = 0
	fuelconsumption = 1e-6 #To prevent divide-by-zero errors
	lf = 0
	lfox = 0
	mp = 0
	for engine in engines:
		if engine == "Rapier":
			if rapierclosedcycle:
				thrust += 180000 * engines[engine] * throttles[engine]
				fuelconsumption += 180000/305/g0 * 9/8 * engines[engine] * throttles[engine]
				lfox += 180000/305/g0 * engines[engine] * throttles[engine]
			else:
				rapierthrust = RAPIERThrust(mach, alt)
				thrust += rapierthrust * engines[engine] * throttles[engine]
				fuelconsumption += rapierthrust/3200/g0 * 11/10 * engines[engine] * throttles[engine]
				lf += rapierthrust/3200/g0 * engines[engine] * throttles[engine]
			continue
		if engine == "Dawn":
			continue
		count = engines[engine]
		thrust += engine_info[engine][2] * engines[engine] * throttles[engine]
		engineconsumption = engine_info[engine][2] * engines[engine] * throttles[engine]/engine_info[engine][0]/g0
		if engine_info[engine][1] == "LiquidFuel":
			lf += engineconsumption
			engineconsumption *= 11/10
		if engine_info[engine][1] == "LFOx":
			lfox += engineconsumption
			engineconsumption *= 9/8
		if engine_info[engine][1] == "Monopropellant":
			mp += engineconsumption
			engineconsumption *= 8.5/7.5
		fuelconsumption += engineconsumption
	return thrust, fuelconsumption, (thrust - drag)/fuelconsumption, lf, lfox, mp


def simulate(initmass, engines):
	vel = 1000
	mass = initmass * np.exp(-2 * vel/3200/g0) #The 2 is not a constant, it can vary based on ascent profile
	thrust = 1000
	drag = 1
	time = 0
	totalConsumption = [initmass-mass, 0, 0, 0]

	timestep = 1

	throttles = {engine:0 for engine in engines}
	throttles["Rapier"] = 1

	closedcycle = False

	highestIsp = 0
	highestFuelIdx = 0
	highestName = ""
	for engine in engines:
		if engine_info[engine][0] > highestIsp:
			highestIsp = engine_info[engine][0]
			highestFuelIdx = ["LiquidFuel", "LFOx", "Monopropellant", "XenonGas"].index(engine_info[engine][1])
			highestName = engine

	while vel < 2243:
		machnumber = vel/soundspeed
		effectivegravity = mass * abs(gravity - (vel + 2*np.pi*(KERBIN_RADIUS+alt)/21549.425)**2/(KERBIN_RADIUS + alt))
		drag = effectivegravity/ldratio.evaluate(machnumber)

		bestThrust, bestFuel, bestVeff, bestLF, bestLFOX, bestMP = calcThrustFuel(engines, throttles, closedcycle, drag, machnumber, alt)
		for i in range(10):
			for engine in engines:
				throttles[engine] = 1-throttles[engine]
				thrust, fuel, veff, lf, lfox, mp = calcThrustFuel(engines, throttles, closedcycle, drag, machnumber, alt)
				if veff > bestVeff:
					bestThrust = thrust
					bestFuel = fuel
					bestVeff = veff
					bestLF = lf
					bestLFOX = lfox
					bestMP = mp
				else:
					throttles[engine] = 1-throttles[engine]
			closedcycle = not closedcycle
			thrust, fuel, veff, lf, lfox, mp = calcThrustFuel(engines, throttles, closedcycle, drag, machnumber, alt)
			if veff > bestVeff:
				bestThrust = thrust
				bestFuel = fuel
				bestVeff = veff
				bestLF = lf
				bestLFOX = lfox
				bestMP = mp
			else:
				closedcycle = not closedcycle
		
		accel = (bestThrust - drag)/mass

		#print(time, round(vel), round(mass), throttles, closedcycle)

		totalConsumption[0] += bestLF * timestep
		totalConsumption[1] += bestLFOX * timestep
		totalConsumption[2] += bestMP * timestep
		mass -= bestLF * timestep
		mass -= bestLFOX * timestep
		mass -= bestMP * timestep
		vel += accel * timestep
		time += timestep

	mult = np.exp(-50/800/g0)
	totalConsumption[0] += mass * (1 - mult)
	mass *= mult

	mult = np.exp(-1300/highestIsp/g0)
	totalConsumption[highestFuelIdx] += mass * (1 - mult)
	mass *= mult

	#Landing on Vall
	vel = np.sqrt(VALL_STDGP/(VALL_RADIUS + vallAlt))
	throttles = {engine: 0 for engine in engines}
	throttles[highestName] = 1
	time = 0
	burnDir = -1
	while vel < np.sqrt(VALL_STDGP/(VALL_RADIUS + vallAlt)) + 1:
		effectivegravity = vallGravity - vel * vel/(VALL_RADIUS + vallAlt)

		thrust = 0
		fuelconsumption = 0
		for engine in engines:
			thrust += throttles[engine] * engine_info[engine][2] * engines[engine]
			engineconsumption = engine_info[engine][2] * engines[engine] * throttles[engine]/engine_info[engine][0]/g0

			if engine_info[engine][1] == "LiquidFuel":
				engineconsumption *= 11/10

			if engine_info[engine][1] == "LFOx":
				engineconsumption *= 9/8

			if engine_info[engine][1] == "Monopropellant":
				engineconsumption *= 8.5/7.5

			if engine_info[engine][1] == "XenonGas":
				engineconsumption *= 4/3

			fuelconsumption += engineconsumption
		for i in range(10):
			oldveff = 0
			if thrust > mass * effectivegravity:
				oldveff = np.sqrt(thrust * thrust - (mass * effectivegravity)**2)/fuelconsumption
			oldthrottles = deepcopy(throttles)

			oldthrust = thrust
			oldconsumption = fuelconsumption

			for engine in engines:
				if engine == highestName:
					continue

				throttles[engine] += generator.normal(loc=0.0, scale=0.1)
				if throttles[engine] < 0:
					throttles[engine] = 0
				if throttles[engine] > 1:
					throttles[engine] = 1
				
			thrust = 0
			fuelconsumption = 0
			for engine in engines:
				thrust += throttles[engine] * engine_info[engine][2] * engines[engine]
				engineconsumption = engine_info[engine][2] * engines[engine] * throttles[engine]/engine_info[engine][0]/g0

				if engine_info[engine][1] == "LiquidFuel":
					engineconsumption *= 11/10

				if engine_info[engine][1] == "LFOx":
					engineconsumption *= 9/8

				if engine_info[engine][1] == "Monopropellant":
					engineconsumption *= 8.5/7.5

				if engine_info[engine][1] == "XenonGas":
					engineconsumption *= 4/3

				fuelconsumption += engineconsumption

			newveff = 0
			if thrust > mass * effectivegravity:
				newveff = np.sqrt(thrust * thrust - (mass * effectivegravity)**2)/fuelconsumption
			if oldveff > newveff:
				throttles = deepcopy(oldthrottles)

				thrust = oldthrust
				fuelconsumption = oldconsumption


		fuelconsumptions = [0, 0, 0, 0]
		for engine in engines:
			engineconsumption = engine_info[engine][2] * engines[engine] * throttles[engine]/engine_info[engine][0]/g0

			if engine_info[engine][1] == "LiquidFuel":
				fuelconsumptions[0] += engineconsumption

			if engine_info[engine][1] == "LFOx":
				fuelconsumptions[1] += engineconsumption

			if engine_info[engine][1] == "Monopropellant":
				fuelconsumptions[2] += engineconsumption

			if engine_info[engine][1] == "XenonGas":
				fuelconsumptions[3] += engineconsumption

		accel = thrust/mass
		if accel < effectivegravity:
			return -1
		haccel = np.sqrt(accel * accel - effectivegravity * effectivegravity)

		#print(time, vel, throttles, burnDir, mass)

		mass -= sum(fuelconsumptions) * timestep
		for i in range(4):
			totalConsumption[i] += fuelconsumptions[i] * timestep
		vel += burnDir * haccel * timestep
		time += timestep

		if vel < vallrotVel:
			vel = vallrotVel
			burnDir = 1
	#print(vel, vallGravity, vallrotVel)

	mult = np.exp(-1200/highestIsp/g0)
	totalConsumption[highestFuelIdx] += mass * (1 - mult)
	mass *= mult

	#Landing on Tylo (scary)
	vel = np.sqrt(TYLO_STDGP/(TYLO_RADIUS + tyloAlt))
	throttles = {engine: 0 for engine in engines}
	throttles[highestName] = 1
	time = 0
	burnDir = -1
	while vel < np.sqrt(TYLO_STDGP/(TYLO_RADIUS + tyloAlt)) + 1:
		effectivegravity = tyloGravity - vel * vel/(TYLO_RADIUS + tyloAlt)

		thrust = 0
		fuelconsumption = 0
		for engine in engines:
			thrust += throttles[engine] * engine_info[engine][2] * engines[engine]
			engineconsumption = engine_info[engine][2] * engines[engine] * throttles[engine]/engine_info[engine][0]/g0

			if engine_info[engine][1] == "LiquidFuel":
				engineconsumption *= 11/10

			if engine_info[engine][1] == "LFOx":
				engineconsumption *= 9/8

			if engine_info[engine][1] == "Monopropellant":
				engineconsumption *= 8.5/7.5

			if engine_info[engine][1] == "XenonGas":
				engineconsumption *= 4/3

			fuelconsumption += engineconsumption
		for i in range(10):
			oldveff = 0
			if thrust > mass * effectivegravity:
				oldveff = np.sqrt(thrust * thrust - (mass * effectivegravity)**2)/fuelconsumption
			oldthrottles = deepcopy(throttles)

			oldthrust = thrust
			oldconsumption = fuelconsumption

			for engine in engines:
				if engine == highestName:
					continue

				throttles[engine] += generator.normal(loc=0.0, scale=0.1)
				if throttles[engine] < 0:
					throttles[engine] = 0
				if throttles[engine] > 1:
					throttles[engine] = 1
				
			thrust = 0
			fuelconsumption = 0
			for engine in engines:
				thrust += throttles[engine] * engine_info[engine][2] * engines[engine]
				engineconsumption = engine_info[engine][2] * engines[engine] * throttles[engine]/engine_info[engine][0]/g0

				if engine_info[engine][1] == "LiquidFuel":
					engineconsumption *= 11/10

				if engine_info[engine][1] == "LFOx":
					engineconsumption *= 9/8

				if engine_info[engine][1] == "Monopropellant":
					engineconsumption *= 8.5/7.5

				if engine_info[engine][1] == "XenonGas":
					engineconsumption *= 4/3

				fuelconsumption += engineconsumption

			newveff = 0
			if thrust > mass * effectivegravity:
				newveff = np.sqrt(thrust * thrust - (mass * effectivegravity)**2)/fuelconsumption
			if oldveff > newveff:
				throttles = deepcopy(oldthrottles)

				thrust = oldthrust
				fuelconsumption = oldconsumption


		fuelconsumptions = [0, 0, 0, 0]
		for engine in engines:
			engineconsumption = engine_info[engine][2] * engines[engine] * throttles[engine]/engine_info[engine][0]/g0

			if engine_info[engine][1] == "LiquidFuel":
				fuelconsumptions[0] += engineconsumption

			if engine_info[engine][1] == "LFOx":
				fuelconsumptions[1] += engineconsumption

			if engine_info[engine][1] == "Monopropellant":
				fuelconsumptions[2] += engineconsumption

			if engine_info[engine][1] == "XenonGas":
				fuelconsumptions[3] += engineconsumption

		accel = thrust/mass
		if accel < effectivegravity:
			return -1
		haccel = np.sqrt(accel * accel - effectivegravity * effectivegravity)

		#print(time, vel, throttles, burnDir, mass)

		mass -= sum(fuelconsumptions) * timestep
		for i in range(4):
			totalConsumption[i] += fuelconsumptions[i] * timestep
		vel += burnDir * haccel * timestep
		time += timestep

		if vel < tylorotVel:
			vel = tylorotVel
			burnDir = 1

	mult = np.exp(-900/highestIsp/g0)
	totalConsumption[highestFuelIdx] += mass * (1 - mult)
	mass *= mult

	return mass, totalConsumption

bestLayout = {'Nerv': 2.963931979119301, 'Rapier': 0.7, 'Dawn': 0.11043690460225665, 'Rhino': 0.015865908795365435}
bestMass, bestConsumption = simulate(70000, bestLayout)
numEngines = len(engine_info)
for engine in bestLayout:
	bestMass -= engine_info[engine][4] * bestLayout[engine]
bestMass -= bestConsumption[0]/10
bestMass -= bestConsumption[1]/8
bestMass -= bestConsumption[2]/7.5
bestMass -= bestConsumption[3]/3
layout = deepcopy(bestLayout)
for i in range(100000):
	if generator.uniform() < 0.1:
		layout = deepcopy(bestLayout)
	for engine in layout:
		layout[engine] *= generator.normal(loc=1, scale=0.01)

	if generator.uniform() < 0.25:
		layout[list(engine_info.keys())[generator.integers(numEngines)]] = np.exp(generator.uniform(-10, 0))
	if generator.uniform() < 0.25:
		if len(list(layout.keys())) == 0:
			layout = deepcopy(bestLayout)
		del layout[list(layout.keys())[generator.integers(len(list(layout.keys())))]]

	if "Rapier" not in layout:
		continue
	if "Nerv" not in layout:
		continue

	if layout["Rapier"] < 0.7:
		layout["Rapier"] = 0.7

	mass, consumption = simulate(70000, layout)
	for engine in layout:
		mass -= engine_info[engine][4] * layout[engine]
	mass -= consumption[0]/10
	mass -= consumption[1]/8
	mass -= consumption[2]/7.5
	mass -= consumption[3]/3

	if mass > bestMass:
		bestMass = mass
		bestLayout = deepcopy(layout)
		print(i, bestMass, bestLayout)
		continue
	if i % 100 == 0:
		print(i, mass)
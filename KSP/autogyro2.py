import matplotlib.pyplot as plt
import math
import numpy as np
import bisect
from copy import deepcopy

KERBIN_STDGP = 3531600000000
KERBIN_ROT = 21549.425
KERBIN_ADIABATIC = 1.4
KERBIN_MOLARMASS = 0.0289644002914429
R = 8.31446261815324
phi = (1 + np.sqrt(5))/2
generator = np.random.default_rng()

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

def temperatureAlt(alt):
	return readHermiteSpline(alt, kerbinTempBaseKeys) + readHermiteSpline(alt, kerbinSunMultKeys) * readHermiteSpline(0, kerbinLatitudeBiasKeys)

soundspeed = np.sqrt(KERBIN_ADIABATIC * R * temperatureAlt(0)/KERBIN_MOLARMASS)

def vectorLength(a):
	sqdist = 0
	for elem in a:
		sqdist += elem * elem
	return np.sqrt(sqdist)

def vectorNormalize(a):
	length = vectorLength(a)
	b = []
	for elem in a:
		b.append(elem/length)
	return b

def cross(a, b):
	return [a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2], a[0]*b[1] - a[1]*b[0]]

def dot(a, b):
	out = 0
	for e1, e2 in zip(a, b):
		out += e1 * e2
	return out

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

#Calculate autogyro parameters for multi-blade autogyro rings

radius = 10
radialvelocity = 25
timestep = 0.02

tipvelocity = radius * radialvelocity
oeoffset = radialvelocity * timestep

torqueOffset = 50

def calculate(bladeTop, topVector, planevelocity, torquethreshold, allowspeedup):
	#bladeTop = vectorNormalize([-0.01, 1, 0.0]) #Positive X = forward, Positive Y = up, Positive Z = outwards

	#topVector = vectorNormalize([-0.07, 1, 0.0]) #Positive X = prograde/east, Positive Y = radial out/up, Positive Z = antinormal/south
	initVector = vectorNormalize([-topVector[1], topVector[0], 0])
	orthogonalVector = cross(topVector, initVector)

	totalForce = [0, 0, 0]
	totalTorque = [0,0,0]

	for rotation in np.linspace(0, 2 * np.pi, num=1000):
		position = [
			(initVector[0] * np.cos(rotation) + orthogonalVector[0] * np.sin(rotation)) * radius,
			(initVector[1] * np.cos(rotation) + orthogonalVector[1] * np.sin(rotation)) * radius,
			(initVector[2] * np.cos(rotation) + orthogonalVector[2] * np.sin(rotation)) * radius,
		]

		position2 = [
			(initVector[0] * np.cos(rotation - oeoffset) + orthogonalVector[0] * np.sin(rotation - oeoffset)) * radius,
			(initVector[1] * np.cos(rotation - oeoffset) + orthogonalVector[1] * np.sin(rotation - oeoffset)) * radius,
			(initVector[2] * np.cos(rotation - oeoffset) + orthogonalVector[2] * np.sin(rotation - oeoffset)) * radius,
		]

		velocity = [
			(position[0] - position2[0])/timestep + planevelocity,
			(position[1] - position2[1])/timestep,
			(position[2] - position2[2])/timestep,
		]

		currBladeTop = [
			(-initVector[0] * np.sin(rotation) + orthogonalVector[0] * np.cos(rotation)) * bladeTop[0] + topVector[0] * bladeTop[1] + position[0]/radius * bladeTop[2],
			(-initVector[1] * np.sin(rotation) + orthogonalVector[1] * np.cos(rotation)) * bladeTop[0] + topVector[1] * bladeTop[1] + position[1]/radius * bladeTop[2],
			(-initVector[2] * np.sin(rotation) + orthogonalVector[2] * np.cos(rotation)) * bladeTop[0] + topVector[2] * bladeTop[1] + position[2]/radius * bladeTop[2],
		]

		aoadotproduct = dot(vectorNormalize(velocity), vectorNormalize(currBladeTop))
		aoa = np.arccos(aoadotproduct) - 0.5 * np.pi

		mach = vectorLength(velocity)/soundspeed

		Cl = readHermiteSpline(np.sin(abs(aoa)), liftAoAKeys) * readHermiteSpline(mach, liftMachKeys) * mach
		Cd = readHermiteSpline(np.sin(abs(aoa)), dragAoAKeys) * readHermiteSpline(mach, dragMachKeys)/2.4 * mach

		if aoa < 0:
			Cl *= -1

		normalizedVelocity = vectorNormalize(velocity)
		normalizedCurrBladeTop = vectorNormalize(currBladeTop)

		aero = [
			normalizedCurrBladeTop[0] * Cl - normalizedVelocity[0] * Cd,
			normalizedCurrBladeTop[1] * Cl - normalizedVelocity[1] * Cd,
			normalizedCurrBladeTop[2] * Cl - normalizedVelocity[2] * Cd,
		]

		torque = cross(position, aero)
		
		totalTorque[0] += torque[0]
		totalTorque[1] += torque[1]
		totalTorque[2] += torque[2]

		totalForce[0] += aero[0]
		totalForce[1] += aero[1]
		totalForce[2] += aero[2]

	abslift = np.sqrt(totalForce[1]*totalForce[1] + totalForce[2]*totalForce[2])
	finalTorque = vectorLength(totalTorque)
	if (dot(totalTorque, topVector) > 0):
		finalTorque = 0.0001
	if (dot(totalTorque, topVector) > 0) and (not allowspeedup):
		return -1e9, totalForce, finalTorque
	if abslift < torquethreshold * (finalTorque + torqueOffset) * (1 - (planevelocity/2500)**2):
		return -1e9, totalForce, finalTorque

	return totalForce[0]/abslift, totalForce, finalTorque

bestBladeTop = vectorNormalize([-1, 0.1, 0.0])
bestTopVector = vectorNormalize([0.1, 0.01, 0.0])
currThreshold = 0.0
bestScore, bestForce, bestTorque = calculate(bestBladeTop, bestTopVector, 400, currThreshold, True)

bladeTop = deepcopy(bestBladeTop)
topVector = deepcopy(bestTopVector)

for i in range(10000):
	if generator.uniform() < 0.02:
		bladeTop = deepcopy(bestBladeTop)
		topVector = deepcopy(bestTopVector)

	bladeTop[0] += generator.normal(loc=0, scale=0.1)
	bladeTop[2] += generator.normal(loc=0, scale=0.1)
	topVector[0] += generator.normal(loc=0, scale=0.1)
	topVector[0] = abs(topVector[0])
	#topVector[2] += generator.normal(loc=0, scale=0.1)
	bladeTop = vectorNormalize(bladeTop)
	topVector = vectorNormalize(topVector)

	score, force, torque = calculate(bladeTop, topVector, 400, currThreshold, True)
	#print(i, -1/score, -1/bestScore, bestForce, bestTorque, np.sqrt(force[1] * force[1] + force[2] * force[2])/(torque + torqueOffset)/currThreshold, currThreshold)
	if score < -1e8:
		continue
	if (score < bestScore) and ((generator.uniform() < 0.99) or (np.sqrt(force[1] * force[1] + force[2] * force[2])/(torque + torqueOffset)/currThreshold < 1)):
		continue

	bestScore = score
	bestBladeTop = deepcopy(bladeTop)
	bestTopVector = deepcopy(topVector)
	bestForce = deepcopy(force)
	bestTorque = torque

	currThreshold = max(currThreshold, min(i/10000, np.sqrt(force[1] * force[1] + force[2] * force[2])/(torque+torqueOffset)))

	print(i, score, bladeTop, topVector)
	print(i, -1/score, -1/bestScore, bestForce, bestTorque, np.sqrt(force[1] * force[1] + force[2] * force[2])/(torque + torqueOffset)/currThreshold, currThreshold)

print(-1/bestScore, bestForce, bestBladeTop, bestTopVector)

topVectors = [
	vectorNormalize([0.9990556970014287, 0.0434478341116056, 0.0]),
	vectorNormalize([0.9957765295551383, 0.09181014751717251, 0.0]),
	vectorNormalize([0.9914819529041644, 0.13024414407313734, 0.0])
]
bladeTops = [
	vectorNormalize([-0.8155251469963354, 0.1760634267263142, -0.5512897644487819]), 
	vectorNormalize([-0.7744845139695163, 0.08772025507162573, -0.6264813600352137]),
	vectorNormalize([-0.7851573422672639, 0.4356136232623523, -0.4401916845103414])
]

vels = []
old = []
news = [[] for top in topVectors]
for i in np.linspace(0, 2000, num=100):
	vels.append(i)
	old.append(ldratio.evaluate(i/soundspeed))
	for j in range(len(news)):
		score, force, torque = calculate(bladeTops[j], topVectors[j], i, 0, True)
		#print(score, force, torque)
		news[j].append(-1/score)

		#print(-1/score, force[1])
	print(i, ldratio.evaluate(i/soundspeed))

plt.plot(vels, old, label="Conventional wing")
for i in range(len(news)):
	plt.plot(vels, news[i], label="Config " + str(i))
plt.legend()
plt.title("Upwards force/Backwards force ratio for various wing configurations")
plt.xlabel("Velocity (m/s)")
plt.ylim(0, 30)
plt.xlim(0, 2000)
plt.show()
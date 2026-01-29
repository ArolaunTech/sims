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

	return 105.0 * flowcap(pressureMult * machmult)

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

alt = 22500
density = 1000 * pressureAlt(alt)/temperatureAlt(alt) * KERBIN_MOLARMASS/R
soundspeed = np.sqrt(KERBIN_ADIABATIC * R * temperatureAlt(alt)/KERBIN_MOLARMASS)

gravity = KERBIN_STDGP/(600000 + alt)**2
nervs = 1

def simulate(initmass):
	vel = 1400
	mass = initmass * np.exp(-2 * vel/3200/9.80665) #The 2 is not a constant, it can vary based on ascent profile
	thrust = 1000
	drag = 1
	time = 0

	timestep = 0.5

	rocketIsp = nervIsp(alt)

	while drag/thrust < 1-rocketIsp/3200:
		machnumber = vel/soundspeed
		thrust = 1000 * RAPIERThrust(machnumber, alt)
		drag = mass * (gravity - ((vel + 2 * np.pi * (600000 + alt)/KERBIN_ROT)**2)/(600000 + alt))/ldratio.evaluate(machnumber)

		fuelConsumption = thrust/3200/9.80665

		vel += (thrust-drag)/mass * timestep
		mass -= fuelConsumption * timestep
		time += timestep
	while vel < 2243:
		machnumber = vel/soundspeed
		thrust = 1000 * RAPIERThrust(machnumber, alt) + rocketIsp/800 * nervs * 6e4
		drag = mass * (gravity - vel * vel/(600000 + alt))/ldratio.evaluate(machnumber)

		fuelConsumption = thrust/3200/9.80665 + nervs * 6e4/800/9.80665

		vel += (thrust-drag)/mass * timestep
		mass -= fuelConsumption * timestep
		time += timestep

	return mass, time

if __name__ == '__main__':
	xs = []
	ys = []
	ts = []
	ds = []
	for initmass in np.linspace(10000, 120000, num=200):
		finalmass, finaltime = simulate(initmass)
		finalmass *= np.exp(-45/800/9.80665)
		print(initmass, finalmass, finaltime, 1000*finalmass/(0.1 * (initmass - finalmass) + 2000 + (initmass)/70), (0.1 * (initmass - finalmass) + 2000 + (initmass)/70))

		xs.append(initmass)
		ys.append(finalmass)
		ts.append(finaltime)
		ds.append(1000*finalmass/(0.1 * (initmass - finalmass) + 2000 + (initmass)/70))

	plt.plot(xs, ys)
	#plt.plot(xs, ts)
	plt.plot(xs, ds)
	plt.gca().set_xlim(10000, 150000)
	plt.gca().set_ylim(0, 70000)
	plt.title("Maximum mass of a RAPIER-NERV SSTO after reaching 1800 m/s")
	plt.xlabel("Mass per RAPIER on takeoff")
	plt.ylabel("Mass per RAPIER after reaching 1800 m/s")
	plt.show()
import matplotlib.pyplot as plt
import numpy as np

#Constants
molarMass = 0.0289644002914429
gamma = 1.4
R = 8.31446261815324

#Lift and drag curves
liftAoAKeys = [
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

#Kerbin atmospheric curves
kerbinPressureKeys = [
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

#Hermite
def hermite(x, spline):
	"""
	Reads a curve defined by cubic Hermite splines.

	Returns:
	  - y - y-value of curve

	Parameters:
	  - x - x-value of curve
	  - spline[][4] - curve, defined as a set of cubic Hermite splines.
	    - Each point contains:
	      - x - x position
	      - y - y position
	      - y'_1 - left derivative
	      - y'_2 - right derivative
	"""
	nKeys = len(spline)
	minX = spline[0][0] #Edge cases
	maxX = spline[-1][0]
	if (x <= minX):
		return spline[0][1] + (x - minX) * spline[0][2] #Linear extension of spline
	if (x >= maxX):
		return spline[-1][1] + (x - maxX) * spline[-1][3] #Linear extension of spline

	#Determining which region of function we are in
	for i in range(nKeys):
		if (spline[i][0] >= x):
			break

	lx = spline[i - 1][0] #Finding values
	ux = spline[i][0]

	ly = spline[i - 1][1]
	uy = spline[i][1]

	ld = spline[i - 1][3]
	ud = spline[i][2]

	interval = ux - lx
	t = (x - lx)/interval

	t2 = t * t
	t3 = t2 * t

	h0 =  2.0 * t3 - 3.0 * t2 + 1.0 #Hermite splines
	h1 =        t3 - 2.0 * t2 + t
	h2 =        t3 -       t2
	h3 = -2.0 * t3 + 3.0 * t2

	h1 *= interval
	h2 *= interval

	y = h0 * ly + h1 * ld + h2 * ud + h3 * uy
	return y

def temp(alt, sundotnormalized):
	baseTemp = hermite(alt, kerbinTempBaseKeys)
	sunMult = hermite(alt, kerbinSunMultKeys)
	latitudeBias = 17
	latSunMult = 9
	return baseTemp + sunMult * (latitudeBias + latSunMult * sundotnormalized)

def soundSpeed(alt, sundotnormalized):
	return np.sqrt(gamma * R * temp(alt, sundotnormalized)/molarMass)

def density(alt, sundotnormalized):
	pressure = 1000*hermite(alt, kerbinPressureKeys)
	return pressure * molarMass/(R * temp(alt, sundotnormalized))

def dragCoefficient(aoa, mach):
	sinAoA = np.sin(aoa)
	return hermite(sinAoA, dragAoAKeys) * hermite(mach, dragMachKeys)

def liftCoefficient(aoa, mach):
	sinAoA = np.sin(aoa)
	return hermite(sinAoA, liftAoAKeys) * hermite(mach, liftMachKeys)

def LDratio(aoa, mach):
	cD = dragCoefficient(aoa, mach)
	cL = liftCoefficient(aoa, mach)
	return cL * np.cos(aoa)/(cL * np.sin(aoa) + cD)

def test(alt, sundot, vy, r, va, aoa):
	vx = r * va
	v = np.sqrt(vx*vx + vy*vy)

	propAngle = np.arctan2(vy, vx)
	actualAoA = aoa - propAngle
	speedOfSound = soundSpeed(alt, sundot)
	atmDensity = density(alt, sundot)
	mach = v/speedOfSound

	cD = dragCoefficient(actualAoA, mach)
	cL = liftCoefficient(actualAoA, mach)
	dx = -vx/v
	dy = -vy/v
	lx = -np.sin(aoa)
	ly = np.cos(aoa)
	fx = 0.5 * atmDensity * v * v * bladeArea * (36 * cL * lx + 15 * cD * dy)
	fy = 0.5 * atmDensity * v * v * bladeArea * (36 * cL * ly + 15 * cD * dy)

	resistance = -fx * r
	resistance *= np.cos(va / 50) # DUMP correction
	return resistance, fy

def findva(alt, sundot, vy, r, aoa):
	resistance, fy = test(alt, sundot, vy, r, 460 * np.pi / 30, aoa)
	if resistance <= torque:
		return 460 * np.pi / 30

	lo = 0
	hi = 460 * np.pi / 30
	while (hi - lo) > 1e-3:
		middle = (lo + hi)/2
		resistance, fy = test(alt, sundot, vy, r, middle, aoa)
		if resistance <= torque:
			lo = middle
		else:
			hi = middle
	return middle

def findvavy(alt, sundot, r, aoa):
	lo = 0
	hi = 1000
	while (hi - lo) > 1e-3:
		middle = (lo + hi)/2
		va = findva(alt, sundot, middle, r, aoa)
		resistance, fy = test(alt, sundot, middle, r, va, aoa)
		if fy <= lift:
			hi = middle
		else:
			lo = middle
	return va, middle

def canAscend(alt, sundot, r, aoa):
	va = findva(alt, sundot, 0, r, aoa)
	resistance, fy = test(alt, sundot, 0, r, va, aoa)
	return fy > lift

def findmaxalt(sundot, r, aoa):
	print(sundot, r, aoa)
	if not canAscend(0, sundot, r, aoa):
		return 0
	lo = 0
	hi = 70000
	while (hi - lo) > 1e-3:
		middle = (lo + hi)/2
		if canAscend(middle, sundot, r, aoa):
			lo = middle
		else:
			hi = middle
	return middle

"""
def findmaxalt(sundot, r, aoa):
	return 3000*(aoa + r)
"""

#Craft
lift = 500 * 9.81 #Weight (N)
torque = 15000 #Ns
bladeArea = 0.24

rs = np.linspace(0, 20, num=101)
aoas = np.linspace(90, 0, num=91)
rs, aoas = np.meshgrid(rs, aoas)
alt = np.vectorize(findmaxalt)(0.707, rs, aoas * np.pi/180)
maxalt = 0
maxi = 0
maxj = 0
for i in range(len(alt)):
	for j in range(len(alt[0])):
		if alt[i][j] > maxalt:
			maxalt = alt[i][j]
			maxi = i
			maxj = j
plt.title("Max. alt. for 500 kg prop craft with 15 kNm torque")
plt.xlabel("Prop radius (m)")
plt.ylabel("Angle of attack (degrees)")
plt.imshow(alt, extent=(0, 20, 0, 90), aspect=2/9, cmap="inferno")
plt.plot(rs[maxi][maxj], aoas[maxi][maxj], 'kx')
plt.annotate("Best alt: " + str(int(maxalt)) + "m", (rs[maxi][maxj], aoas[maxi][maxj]))
plt.colorbar()
cs = plt.contour(rs, aoas, alt, levels=np.arange(0, 30000, 2500), colors="white")
plt.clabel(cs)
plt.show()
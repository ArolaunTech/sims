#include <iostream>
#include <cmath>
#include <vector>
#include <random>

#define pi 3.14159
#define KERBIN_STDGP (KERBIN_RADIUS * KERBIN_RADIUS * 9.81)
#define KERBIN_RADIUS 6.0e5
#define KERBIN_ROT 21549.425
#define KERBIN_ATMO_HEIGHT 7.0e4
#define RAPIER_MAXALT 2.7e4
#define KERBIN_MOLARMASS 0.0289644002914429
#define KERBIN_ADIABATIC 1.4
#define R 8.31446261815324

#define kerbinPKeys 21
#define kerbinTKeys 9
#define kerbinSMKeys 8
#define kerbinLBKeys 9
#define kerbinLSMKeys 6
const double kerbinPressureKeys[kerbinPKeys][4] = { //Kerbin pressure keys from Kittopia Dumps, x, y, y', y'
	{    0.0,   101.325,        0.0,          -0.01501631  },
	{ 1241.025,  84.02916,     -0.01289846,   -0.01289826  },
	{ 2439.593,  69.68138,     -0.01107876,   -0.01107859  },
	{ 3597.11,   57.78001,     -0.009515483,  -0.009515338 },
	{ 4714.942,  47.90862,     -0.00817254,   -0.008172415 },
	{ 5794.409,  39.72148,     -0.00701892,   -0.007018813 },
	{ 6836.791,  32.93169,     -0.006027969,  -0.006027877 },
	{ 7843.328,  27.30109,     -0.005176778,  -0.0051767   },
	{ 8815.22,   22.63206,     -0.004445662,  -0.004445578 },
	{10786.42,   15.3684,      -0.003016528,  -0.00301646  },
	{12101.4,    11.87313,     -0.002329273,  -0.00232922  },
	{13417.05,    9.172798,    -0.001798594,  -0.001798554 },
	{16678.47,    4.842261,    -0.0009448537, -0.0009448319},
	{21143.1,     2.050097,    -0.0003894095, -0.0003894005},
	{26977.92,    0.6905929,   -0.0001252565, -0.0001252534},
	{33593.82,    0.2201734,   -3.626878e-5,  -3.626788e-5 },
	{42081.87,    0.05768469,  -9.063159e-6,  -9.062975e-6 },
	{49312.13,    0.01753794,  -3.029397e-6,  -3.029335e-6 },
	{56669.95,    0.004591824, -8.827175e-7,  -8.826996e-7 },
	{62300.84,    0.001497072, -3.077091e-7,  -3.077031e-7 },
	{70000.0,     0.0,          0.0,           0.0         }
};

const double kerbinTempBaseKeys[kerbinTKeys][4] = {
	{    0.0,  288.15,   0.0,         -0.008125   },
	{ 8815.22, 216.65,  -0.008096968,  0.0        },
	{16050.39, 216.65,   0.0,          0.001242164},
	{25729.23, 228.65,   0.001237475,  0.003464929},
	{37879.44, 270.65,   0.00344855,   0.0        },
	{41129.24, 270.65,   0.0,         -0.003444189},
	{57440.13, 214.65,  -0.003422425, -0.002444589},
	{68797.88, 186.946, -0.002433851,  0.0        },
	{70000.0,  186.946,  0.0,          0.0        }
};

const double kerbinSunMultKeys[kerbinSMKeys][4] = {
	{    0.00, 1.0,  0.0,          0.0         },
	{ 8815.22, 0.3, -5.91316e-5,  -5.91316e-5  },
	{16050.39, 0.0,  0.0,          0.0         },
	{25729.23, 0.0,  0.0,          0.0         },
	{37879.44, 0.2,  0.0,          0.0         },
	{57440.13, 0.2,  0.0,          0.0         },
	{63902.72, 1.0,  0.0001012837, 0.0001012837},
	{70000.00, 1.2,  0.0,          0.0         }
};

const double kerbinLatitudeBiasKeys[kerbinLBKeys][4] = {
	{ 0.0,  17.0,      0.0,       -0.3316494},
	{10.0,  12.0,     -0.65,      -0.65     },
	{18.0,   6.36371, -0.4502313, -0.4502313},
	{30.0,   0.0,     -1.3,       -1.3      },
	{35.0, -10.0,     -1.65,      -1.65     },
	{45.0, -23.0,     -1.05,      -1.05     },
	{55.0, -31.0,     -0.6,       -0.6      },
	{70.0, -37.0,     -0.6689383, -0.6689383},
	{90.0, -50.0,     -0.02418368, 0.0      }
};

const double kerbinLatitudeSunMultKeys[kerbinLSMKeys][4] = {
	{ 0.0,  9.0,       0.0,          0.1554984  },
	{40.0, 14.2,       0.08154097,   0.08154097 },
	{55.0, 14.9,      -0.006055089, -0.006055089},
	{68.0, 12.16518,  -0.2710912,   -0.2710912  },
	{76.0,  8.582909, -0.6021729,   -0.6021729  },
	{90.0,  5.0,       0.0,          0.0        }
};

const double rapierIspKeys[3][4] = { //RAPIER curves from KSP Wiki
	{0.0, 305.0, 0.0, -30.0},
	{1.0, 275.0, -305.0/9.0, -305.0/9.0},
	{9.0, 0.001, -275.0/8.0, 0.0}
};

const double rapierPressureKeys[5][4] = {
	{0.0,   0.0,  0.0,       0.0},
	{0.018, 0.09, 7.914787,  7.914787},
	{0.08,  0.3,  1.051923,  1.051923},
	{0.35,  0.5,  0.3927226, 0.3927226},
	{1.0,   1.0,  1.055097,  0.0}
};

const double rapierMachKeys[8][4] = {
	{0.0,  1.0,  0.0,        0.08333334},
	{0.2,  0.98, 0.42074,    0.42074   },
	{0.7,  1.8,  2.290406,   2.290406  },
	{1.4,  4.0,  3.887193,   3.887193  },
	{3.75, 8.5,  0.0,        0.0       },
	{4.5,  7.3, -2.831749,  -2.831749  },
	{5.5,  3.0, -5.260566,  -5.260566  },
	{6.0,  0.0, -0.02420209, 0.0       }
};

const double liftAoAKeys[5][4] = { //Thanks to Lt_Duckweed
	{0.0,       0.0,        0.0,        1.965926 },
	{0.258819,  0.5114774,  1.990092,   1.905806 },
	{0.5,       0.9026583,  0.7074468, -0.7074468},
	{0.7071068, 0.5926583, -2.087948,  -1.990095 },
	{1.0,       0.0,       -2.014386,  -2.014386 }
};

const double dragAoAKeys[5][4] = {
	{0.0,       0.01, 0.0,       0.0      },
	{0.3420201, 0.06, 0.1750731, 0.1750731},
	{0.5,       0.24, 2.60928,   2.60928  },
	{0.7071068, 1.7,  3.349777,  3.349777 },
	{1.0,       2.4,  1.387938,  0.0      }
};

const double liftMachKeys[5][4] = {
	{ 0.0, 1.0,    0.0,           0.0       },
	{ 0.3, 0.5,   -1.671345,     -0.8273422 },
	{ 1.0, 0.125, -0.0005291355, -0.02625772},
	{ 5.0, 0.0625, 0.0,           0.0       },
	{25.0, 0.05,   0.0,           0.0       }
};

const double dragMachKeys[9][4] = {
	{ 0.0,  0.35,  0.0,         -0.8463008 },
	{ 0.15, 0.125, 0.0,          0.0       },
	{ 0.9,  0.275, 0.541598,     0.541598  },
	{ 1.1,  0.75,  0.0,          0.0       },
	{ 1.4,  0.4,  -0.3626955,   -0.3626955 },
	{ 1.6,  0.35, -0.1545923,   -0.1545923 },
	{ 2.0,  0.3,  -0.09013031,  -0.09013031},
	{ 5.0,  0.22,  0.0,          0.0       },
	{25.0,  0.3,   0.0006807274, 0.0       }
};

struct vector2D {
	double x;
	double y;
};

struct simulationResult {
	double finalMass;
	double punishMass;
	double totalLF;
	double totalLFOX;
	double apoapsis;
	double periapsis;
};

//Randomness settings
std::mt19937_64 generator(time(NULL));
std::uniform_real_distribution<double> uniform(0.0, 1.0);
std::normal_distribution<double> normal(0.0, 1.0);

double readHermiteSpline(double x, int nKeys, const double spline[][4]) {
	/*
	Reads a curve defined by cubic Hermite splines.

	Returns:
	  - double y - y-value of curve

	Parameters:
	  - double x - x-value of curve
	  - double spline[][4] - curve, defined as a set of cubic Hermite splines.
	    - Each point contains:
	      - x - x position
	      - y - y position
	      - y'_1 - left derivative
	      - y'_2 - right derivative
	  - int nKeys - number of keys in curve
	*/
	double minX = spline[0][0]; //Edge cases
	double maxX = spline[nKeys - 1][0];
	if (x <= minX) {
		return spline[0][1] + (x - minX) * spline[0][2]; //Linear extension of spline
	}
	if (x >= maxX) {
		return spline[nKeys - 1][1] + (x - maxX) * spline[nKeys - 1][3]; //Linear extension of spline
	}

	int i; //Determining which region of function we are in
	for (i = 0; i < nKeys; i++) {
		if (spline[i][0] >= x) {
			break;
		}
	}

	double lx = spline[i - 1][0]; //Finding values
	double ux = spline[i][0];

	double ly = spline[i - 1][1];
	double uy = spline[i][1];

	double ld = spline[i - 1][3];
	double ud = spline[i][2];

	double interval = ux - lx;
	double t = (x - lx)/interval;

	double t2 = t * t;
	double t3 = t2 * t;

	double h0 =  2.0 * t3 - 3.0 * t2 + 1.0; //Hermite splines
	double h1 =        t3 - 2.0 * t2 + t;
	double h2 =        t3 -       t2;
	double h3 = -2.0 * t3 + 3.0 * t2;

	h1 *= interval;
	h2 *= interval;

	double y = h0 * ly +
	                  h1 * ld +
	                  h2 * ud +
	                  h3 * uy;
	return y;
}

inline double pressureAlt(double alt) {
	/*
	Calculates pressure in Kerbin's atmosphere at a given altitude.

	Returns:
	  - double pressure - pressure of atmosphere in kilopascals.

	Parameters:
	  - double altitude - altitude above sea level in meters.
	*/
	return readHermiteSpline(alt, kerbinPKeys, kerbinPressureKeys);
}

double temperatureAlt(double alt, double x, double y) {
	/*
	Calculates temperature in Kerbin's atmosphere at a given altitude.

	Returns:
	  - double temperature - temperature of atmosphere in kelvins.

	Parameters:
	  - double altitude - altitude above sea level in meters.
	  - double x, double y - position in meters
	*/
	double baseTemp = readHermiteSpline(alt, kerbinTKeys, kerbinTempBaseKeys);
	double sunMult = readHermiteSpline(alt, kerbinSMKeys, kerbinSunMultKeys);
	double latitudeBias = readHermiteSpline(0.0, kerbinLBKeys, kerbinLatitudeBiasKeys);
	double latitudeSunMult = readHermiteSpline(0.0, kerbinLSMKeys, kerbinLatitudeSunMultKeys);
	double dist = std::sqrt(x * x + y * y);
	double sunDotNormalized = std::sqrt(2) * 0.25 * (x - y)/dist + 0.5;

	//std::cout << sunMult << "\n\n";

	return baseTemp + sunMult * (latitudeBias + latitudeSunMult * sunDotNormalized);
}

double flowcap(double x) {
	if (x <= 3.0) return x;
	return 6.0 * x/(x + 3.0);
}

double RAPIERThrust(double vx, double vy, double alt, bool mode, double x, double y, double temperature, double pressure) {
	/*
	Calculates R.A.P.I.E.R. thrust for a given velocity, altitude, and mode.
	
	Returns:
	  - double thrust - thrust from engine in kilonewtons.
	
	Parameters:
	  - double vx, vy - velocity relative to planet in m/s.
	  - double alt - altitude above sea level in meters.
	  - bool mode - determines whether air-breathing or closed-cycle.
	    - Possible values:
	  	  - true: closed-cycle mode
	  	  - false: open-cycle mode
	  - double x, y - position in meters
	*/

	//double pressure = pressureAlt(alt)/101.325; //Pressure in atmospheres.
	//std::cout << pressure << "\n";

	if (mode) {
		//Closed-cycle - rocket engine
		double rapierIsp = readHermiteSpline(pressure/101.325, 3, rapierIspKeys);
		return 180.0 * rapierIsp/305.0;
	} else {
		//Open-cycle - jet engine
		if (alt > RAPIER_MAXALT) {
			return 0.0;
		}
		double pressureMult = readHermiteSpline(pressure/101.325, 5, rapierPressureKeys);

		//double temperature = temperatureAlt(alt, x, y);
		double soundspeed = std::sqrt(KERBIN_ADIABATIC * R * temperature/KERBIN_MOLARMASS);

		double sx = 2.0 * pi * y/KERBIN_ROT;
		double sy = -2.0 * pi * x/KERBIN_ROT;
		double surfacespeed = std::sqrt((vx - sx)*(vx - sx) + (vy - sy)*(vy - sy));
		double machnumber = surfacespeed/soundspeed;
		double machmult = readHermiteSpline(machnumber, 8, rapierMachKeys);
		//std::cout << pressureMult * machmult << "\n" << soundspeed << "\n" << machnumber << "\n\n";
		return 105.0 * flowcap(pressureMult * machmult);
	}
	//return 100.0;
}

vector2D calculateWingForce(double x, double y, double vx, double vy, double cx, double cy, double wingarea, double temperature, double pressure, double clMach, double cdMach) {
	/*
	Calculates force from wing parts.
	
	Returns:
	  - vector2D force - force from wing.
	
	Parameters:
	  - double x, y - position in meters
	  - double vx, vy - velocity in meters
	  - double cx, cy - craft orientation
	  - double wingarea - area of wing
	*/

	double dist = std::sqrt(x*x + y*y);
	double alt = dist - KERBIN_RADIUS;

	double sx = 2.0 * pi * y/KERBIN_ROT;
	double sy = -2.0 * pi * x/KERBIN_ROT;

	double surfacespeed = std::sqrt((vx - sx)*(vx - sx) + (vy - sy)*(vy - sy));

	//double temperature = temperatureAlt(alt, x, y);
	double soundspeed = std::sqrt(KERBIN_ADIABATIC * R * temperature/KERBIN_MOLARMASS);

	double machnumber = surfacespeed/soundspeed;
	//std::cout << "Mach: " << machnumber << "\n";

	double cosAoA = ((vx-sx) * cx + (vy-sy) * cy)/surfacespeed;
	//std::cout << cosAoA << " " << cx << " " << cy << "\n";
	double sinAoA = std::sqrt(1 - std::pow(cosAoA, 2.0));
	double clAoA = readHermiteSpline(sinAoA, 5, liftAoAKeys);
	double cdAoA = readHermiteSpline(sinAoA, 5, dragAoAKeys);

	//double clMach = readHermiteSpline(machnumber, 5, liftMachKeys);
	//double cdMach = readHermiteSpline(machnumber, 9, dragMachKeys);

	//std::cout << clAoA << " " << clMach << "\n";
	//std::cout << cdAoA << " " << cdMach << "\n";

	double cl = clAoA * clMach;
	double cd = cdAoA * cdMach;

	//double pressure = pressureAlt(alt);
	double density = pressure * 1000.0 * KERBIN_MOLARMASS/(R * temperature);

	double drag = cd * wingarea * density * surfacespeed * surfacespeed * 0.5 * 0.001 * 25.0;
	double lift = cl * wingarea * density * surfacespeed * surfacespeed * 0.5 * 0.001 * 25.0;

	double dx = drag * (sx - vx)/surfacespeed;
	double dy = drag * (sy - vy)/surfacespeed;

	double lx = cy;
	double ly = -cx;

	if (((vx - sx) * lx + (vy - sy) * ly) > 0) {
		lx *= -1.0;
		ly *= -1.0;
	}

	lx *= lift;
	ly *= lift;

	//std::cout << "sin(AoA) = " << sinAoA << ", L/D = " << std::fabs(cosAoA) * cl/(std::fabs(sinAoA) * cl + cd) << ", Lift (kN) = " << lift << ", Drag (kN) = " << drag << "\n";

	vector2D res;
	res.x = lx + dx;
	res.y = ly + dy;

	return res;
}

simulationResult simulate(const std::vector<std::vector<double> > code, double incidence, double startingMass, double targetApo, double targetPeri, double wingarea) {
	double x = 0.0;
	double y = KERBIN_RADIUS;
	double vx = 2.0 * pi * y/KERBIN_ROT + 200.0;
	double vy = 0.0;

	double temperature = 0.0;
	double pressure = 0.0;

	double time = 0.0;

	int currCodeIndex = -1;
	double nextCodeTime = code[0][0];

	double currWingMag = 1.0;
	double currAirThrottle = 1.0;
	double currLFOXThrottle = 0.0;
	double currOrientation = 0.4 * pi;
	double currProgradeFrac = 0.0;

	double cosIncidence = std::cos(incidence);
	double sinIncidence = std::sin(incidence);

	double dist, sqdist = 0.0;
	double gravf, gravx, gravy = 0.0;
	double thrust = 0.0;
	double ccthrust = 0.0;
	double octhrust = 0.0;

	double currMass = startingMass;

	double totalLF = 0.0;
	double totalLFOX = 0.0;

	double cx = std::cos(currOrientation);
	double cy = std::sin(currOrientation);

	double cxp = cx;
	double cyp = cy;
	double cpl = 0.0;

	double surfacespeed = 0.0;
	double sx = 2.0 * pi * KERBIN_RADIUS/KERBIN_ROT;
	double sy = 0.0;

	double clMach;
	double cdMach;

	double soundspeed;
	double machnumber;

	vector2D wingForceUp;
	vector2D wingForceDown;

	const double timestep = 0.5;

	for (int i = 0; i < 10000; i++) {
		//std::cout << x << " " << y << " " << vx << " " << vy << " " << dist - KERBIN_RADIUS << " " << currMass << "\n";
		if ((time > nextCodeTime) && ((currCodeIndex + 1) < code.size())) {
			currCodeIndex++;
			if ((currCodeIndex + 1) < code.size()) {
				nextCodeTime += code[currCodeIndex + 1][0];
			}

			currWingMag = code[currCodeIndex][1];
			currAirThrottle = code[currCodeIndex][2];
			currLFOXThrottle = code[currCodeIndex][3];
			currOrientation = code[currCodeIndex][4];
			currProgradeFrac = code[currCodeIndex][5];
		}

		sx = 2.0 * pi * y/KERBIN_ROT;
		sy = -2.0 * pi * x/KERBIN_ROT;

		//std::cout << sx << " " << sy << "\n";

		surfacespeed = std::sqrt((vx - sx)*(vx - sx) + (vy - sy)*(vy - sy));

		//std::cout << "Surface speed: " << surfacespeed << "\n";

		cxp = std::cos(currOrientation);
		cyp = std::sin(currOrientation);

		if (surfacespeed > 1e-10) {
			cxp = cxp * (1.0 - currProgradeFrac) + (vx - sx) * currProgradeFrac/surfacespeed;
			cyp = cyp * (1.0 - currProgradeFrac) + (vy - sy) * currProgradeFrac/surfacespeed;
		}

		cpl = std::sqrt(cxp*cxp + cyp*cyp);

		cx = cxp/cpl;
		cy = cyp/cpl;

		if (cpl < 1e-10) { //Prevent divide-by-zero
			cx = 0.0;
			cy = 0.0;
		}

		sqdist = x * x + y * y;
		dist = std::sqrt(sqdist);

		gravf = KERBIN_STDGP/sqdist;
		gravx = -gravf * x/dist;
		gravy = -gravf * y/dist;

		vx += gravx * timestep;
		vy += gravy * timestep;

		temperature = temperatureAlt(dist - KERBIN_RADIUS, x, y);
		pressure = pressureAlt(dist - KERBIN_RADIUS);

		soundspeed = std::sqrt(KERBIN_ADIABATIC * R * temperature/KERBIN_MOLARMASS);
		machnumber = surfacespeed/soundspeed;

		clMach = readHermiteSpline(machnumber, 5, liftMachKeys);
		cdMach = readHermiteSpline(machnumber, 9, dragMachKeys);

		wingForceUp = calculateWingForce(x, y, vx, vy, cx*cosIncidence - cy*sinIncidence, cy*cosIncidence + cx*sinIncidence, wingarea, temperature, pressure, clMach, cdMach);
		wingForceDown = calculateWingForce(x, y, vx, vy, cx*cosIncidence + cy*sinIncidence, cy*cosIncidence - cx*sinIncidence, wingarea, temperature, pressure, clMach, cdMach);
		//std::cout << "Wing Force: " << wingForce.x << " " << wingForce.y << "\n";
		vx += ((wingForceUp.x * currWingMag) + (wingForceDown.x * (1.0 - currWingMag))) * timestep/currMass;
		vy += ((wingForceUp.y * currWingMag) + (wingForceDown.y * (1.0 - currWingMag))) * timestep/currMass;

		octhrust = RAPIERThrust(vx, vy, dist - KERBIN_RADIUS, false, x, y, temperature, pressure) * currAirThrottle;
		ccthrust = RAPIERThrust(vx, vy, dist - KERBIN_RADIUS, true, x, y, temperature, pressure) * currLFOXThrottle;
		if (ccthrust > octhrust) {
			thrust = ccthrust;

			totalLFOX += timestep * currLFOXThrottle * 180.0/(305.0 * 9.81);
			currMass -= timestep * currLFOXThrottle * 180.0/(305.0 * 9.81);
		} else {
			thrust = octhrust;

			totalLF += timestep * thrust/(3200.0 * 9.81);
			currMass -= timestep * thrust/(3200.0 * 9.81);
		}
		//std::cout << " " << thrust << " " << totalLF << " " << totalLFOX << "\n";
		vx += thrust * timestep * cx/currMass;
		vy += thrust * timestep * cy/currMass;

		//currMass -= timestep * thrust/(3200.0 * 9.81);

		x += vx * timestep;
		y += vy * timestep;

		if (dist < KERBIN_RADIUS) { //Crashed into ground
			simulationResult res;
			res.finalMass = -1.0;
			res.punishMass = -1.0;
			res.totalLF = -1.0;
			res.totalLFOX = -1.0;
			res.apoapsis = -1.0;
			res.periapsis = -1.0;

			return res;
		}
		if (dist > KERBIN_RADIUS + KERBIN_ATMO_HEIGHT) {
			break;
		}

		time += timestep; //Step time
	}
	//std::cout << x << " " << y << " " << vx << " " << vy << "\n";

	dist = std::sqrt(x*x + y*y);

	sx = y/dist;
	sy = -x/dist;

	double sma = 1.0/(2.0/dist - (vx*vx + vy*vy)/KERBIN_STDGP);
	double hvel = sx * vx + sy * vy;
	double eccentricity = std::sqrt(1.0 - (dist * dist * hvel * hvel)/(sma * KERBIN_STDGP));

	double apoapsis = sma * (1.0 + eccentricity) - KERBIN_RADIUS;
	double periapsis = sma * (1.0 - eccentricity) - KERBIN_RADIUS;

	double targetSMA = KERBIN_RADIUS + 0.5 * (targetApo + targetPeri);
	double targetEccentricity = (KERBIN_RADIUS + targetApo)/targetSMA - 1.0;

	double targetVelocity = std::sqrt(KERBIN_STDGP * (2.0/dist - 1.0/targetSMA));
	double targetHorizontal = std::sqrt(KERBIN_STDGP * targetSMA * (1.0 - targetEccentricity * targetEccentricity))/dist;
	double targetVertical = std::sqrt(targetVelocity * targetVelocity - targetHorizontal * targetHorizontal);

	double targetX = targetHorizontal * sx - targetVertical * sy;
	double targetY = targetHorizontal * sy + targetVertical * sx;

	double error = std::sqrt((vx - targetX) * (vx - targetX) + (vy - targetY) * (vy - targetY));

	//std::cout << vx << " " << vy << " " << targetX << " " << targetY << " " << apoapsis << " " << periapsis << "\n";

	simulationResult res;
	res.finalMass = currMass;
	res.punishMass = currMass * std::exp(-error/(305.0 * 9.81));
	res.totalLFOX = totalLFOX + currMass * (1.0 - std::exp(-error/(305.0 * 9.81)));
	res.totalLF = totalLF;
	res.apoapsis = apoapsis;
	res.periapsis = periapsis;

	return res;
}

double findBestAscent(int nKeys, double incidence, double startingMass, double targetApo, double targetPeri, double wingarea) {
	std::vector<std::vector<double> > bestCode;
	std::vector<std::vector<double> > bestRandomCode;
	std::vector<std::vector<double> > randomCode;
	std::vector<std::vector<double> > newCode;
	std::vector<double> key;

	simulationResult randomRes;
	simulationResult bestRandomRes;
	simulationResult newRes;
	simulationResult bestRes;

	bool setRandomBest = false;

	for (int i = 0; i < 10; i++) {
		//Generate random code sequences
		setRandomBest = false;
		for (int j = 0; j < 5000; j++) {
			randomCode.clear();
			for (int k = 0; k < nKeys; k++) {
				key.clear();
				key.push_back(uniform(generator) * 300.0); //Time
				key.push_back(0.5 + 0.5*uniform(generator)); //wingmag
				key.push_back(0.5 + 0.5*uniform(generator)); //LF Throttle
				key.push_back(0.5*uniform(generator)); //LFOX Throttle
				key.push_back((uniform(generator) - 0.5) * pi); //Orientation (not relative to Kerbin's surface but to its center)
				key.push_back(0.5 + 1.0*uniform(generator)); //Prograde fraction

				randomCode.push_back(key);
			}

			randomRes = simulate(randomCode, incidence, startingMass, targetApo, targetPeri, wingarea);
			if (randomRes.punishMass < 0.0) 
				continue;
			if (std::isnan(randomRes.punishMass))
				continue;
			if (!setRandomBest) {
				bestRandomRes = randomRes;
				bestRandomCode = randomCode;
				setRandomBest = true;
			}
			if (randomRes.punishMass < bestRandomRes.punishMass)
				continue;
			bestRandomRes = randomRes;
			bestRandomCode = randomCode;
		}

		if ((i==0) || (bestRandomRes.punishMass > bestRes.punishMass)) {
			bestRes = bestRandomRes;
			bestCode = bestRandomCode;
		}

		//Printout
		/*
		std::cout << "Simulation Results: " << bestRandomRes.punishMass << ", Total LF Used: " << bestRandomRes.totalLF << "t, Total LFOX Used: " << bestRandomRes.totalLFOX << "t\n";
		std::cout << "Random Code: {\n";
		for (int j = 0; j < nKeys; j++) {
			std::cout << "\t{";
			for (int k = 0; k < 6; k++) {
				std::cout << bestRandomCode[j][k];
				if (k < 5) std::cout << ", ";
			}
			std::cout << "}";
			if (j + 1 < nKeys) {
				std::cout << ",";
			}
			std::cout << "\n";
		}
		std::cout << "}\n\n";
		*/

		for (int j = 0; j < 10000; j++) {
			newCode = bestRandomCode;
			for (int row = 0; row < nKeys; row++) {
				newCode[row][0] += 10.0 * normal(generator);
				newCode[row][1] += 0.2 * normal(generator);
				newCode[row][2] += 0.2 * normal(generator);
				newCode[row][3] += 0.2 * normal(generator);
				newCode[row][4] += 0.2 * normal(generator);
				newCode[row][5] += 0.2 * normal(generator);

				//Put stuff in bounds!
				if (newCode[row][0] < 0.0) newCode[row][0] = 0.0;
				if (newCode[row][1] < 0.0) newCode[row][1] = 0.0;
				if (newCode[row][2] < 0.0) newCode[row][2] = 0.0;
				if (newCode[row][3] < 0.0) newCode[row][3] = 0.0;
				if (newCode[row][5] < 0.0) newCode[row][5] = 0.0;
				
				if (newCode[row][1] > 0.0) newCode[row][1] = 1.0;
				if (newCode[row][2] > 1.0) newCode[row][2] = 1.0;
				if (newCode[row][3] > 1.0) newCode[row][3] = 1.0;
			}
			newRes = simulate(newCode, incidence, startingMass, targetApo, targetPeri, wingarea);

			if (newRes.punishMass < 0.0) continue;
			if (std::isnan(newRes.punishMass)) continue;
			if (newRes.punishMass < bestRandomRes.punishMass) continue;

			bestRandomRes = newRes;
			bestRandomCode = newCode;

			if (newRes.punishMass < bestRes.punishMass) continue;

			bestRes = newRes;
			bestCode = newCode;

			std::cout << "Simulation Results: " << newRes.punishMass << ", Total LF Used: " << newRes.totalLF << "t, Total LFOX Used: " << newRes.totalLFOX << "t\n";
			std::cout << "Ascent Code: {\n";
			for (int j = 0; j < nKeys; j++) {
				std::cout << "\t{";
				for (int k = 0; k < 6; k++) {
					std::cout << bestRandomCode[j][k];
					if (k < 5) std::cout << ", ";
				}
				std::cout << "}";
				if (j + 1 < nKeys) {
					std::cout << ",";
				}
				std::cout << "\n";
			}
			std::cout << "}\n\n";
		}
	}
	return 0.0;
}

int main() {
	/*
	std::vector<std::vector<double> > code(1, std::vector<double>(6));
	code[0][0] = 0.0;
	code[0][1] = -1.0;
	code[0][2] = 1.0;
	code[0][3] = 0.0;
	code[0][4] = pi*0.35;
	code[0][5] = 0.2;

	simulationResult res = simulate(code, 0.1, 6.4, 97782.0, -371211.336852, 1.36);
	std::cout << "Simulation Results: " << res.punishMass << ", Total LF Used: " << res.totalLF << "t, Total LFOX Used: " << res.totalLFOX << "t\n";
	*/

	double optimizerRes = findBestAscent(5, 0.0, 50.0, 75000.0, 0.0, 20.0);
	std::cout << "Optimizer Results: " << optimizerRes << "\n";
	return 0;
}
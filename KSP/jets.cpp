#include <iostream>
#include <vector>
#include <cmath>

const double KERBIN_MOLAR_MASS = 0.0289644002914429;
const double KERBIN_ADIABATIC_INDEX = 1.4;
const double KERBIN_RADIUS = 600000;
const double KERBIN_ROTATIONAL_PERIOD = 21549.425;

const double R = 8.31446261815324;
const double pi = 3.1415926535897932384626433832795028841972;
const double g0 = 9.81;

const double inversePhi = (std::sqrt(5) - 1) / 2;

const double KERBIN_ROT_VEL = 2 * pi * KERBIN_RADIUS / KERBIN_ROTATIONAL_PERIOD;
const double KERBIN_GRAVITY = 9.81;
const double KERBIN_STDGP = KERBIN_RADIUS * KERBIN_RADIUS * KERBIN_GRAVITY;

const std::vector<std::vector<double> > kerbinPressureKeys = { //Kerbin pressure keys from Kittopia Dumps, x, y, y', y'
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

const std::vector<std::vector<double> > kerbinTempBaseKeys = {
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

const std::vector<std::vector<double> > kerbinSunMultKeys = {
	{    0.00, 1.0,  0.0,          0.0         },
	{ 8815.22, 0.3, -5.91316e-5,  -5.91316e-5  },
	{16050.39, 0.0,  0.0,          0.0         },
	{25729.23, 0.0,  0.0,          0.0         },
	{37879.44, 0.2,  0.0,          0.0         },
	{57440.13, 0.2,  0.0,          0.0         },
	{63902.72, 1.0,  0.0001012837, 0.0001012837},
	{70000.00, 1.2,  0.0,          0.0         }
};

const std::vector<std::vector<double> > kerbinLatitudeBiasKeys = {
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

const std::vector<std::vector<double> > kerbinLatitudeSunMultKeys = {
	{ 0.0,  9.0,       0.0,          0.1554984  },
	{40.0, 14.2,       0.08154097,   0.08154097 },
	{55.0, 14.9,      -0.006055089, -0.006055089},
	{68.0, 12.16518,  -0.2710912,   -0.2710912  },
	{76.0,  8.582909, -0.6021729,   -0.6021729  },
	{90.0,  5.0,       0.0,          0.0        }
};

const std::vector<std::vector<double> > liftAoAKeys = { //Thanks to Lt_Duckweed
	{0.0,       0.0,        0.0,        1.965926 },
	{0.258819,  0.5114774,  1.990092,   1.905806 },
	{0.5,       0.9026583,  0.7074468, -0.7074468},
	{0.7071068, 0.5926583, -2.087948,  -1.990095 },
	{1.0,       0.0,       -2.014386,  -2.014386 }
};

const std::vector<std::vector<double> > dragAoAKeys = {
	{0.0,       0.01, 0.0,       0.0      },
	{0.3420201, 0.06, 0.1750731, 0.1750731},
	{0.5,       0.24, 2.60928,   2.60928  },
	{0.7071068, 1.7,  3.349777,  3.349777 },
	{1.0,       2.4,  1.387938,  0.0      }
};

const std::vector<std::vector<double> > liftMachKeys = {
	{ 0.0, 1.0,    0.0,           0.0       },
	{ 0.3, 0.5,   -1.671345,     -0.8273422 },
	{ 1.0, 0.125, -0.0005291355, -0.02625772},
	{ 5.0, 0.0625, 0.0,           0.0       },
	{25.0, 0.05,   0.0,           0.0       }
};

const std::vector<std::vector<double> > dragMachKeys = {
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

const std::vector<std::vector<double> > junoAtm = {
	{0,     0,     0,         0.7448742},
	{0.072, 0.13,  2.075459,  2.075459 },
	{0.16,  0.28,  1.464173,  1.464173 },
	{0.42,  0.578, 0.93687,   0.93687  },
	{1,     1,     0.5529748, 0        }
};

const std::vector<std::vector<double> > junoMach = {
	{0,    1,     0,         0        },
	{0.44, 0.897, 0,         0        },
	{1,    1,     0.1988732, 0.1988732},
	{1.3,  1.03,  0,         0        },
	{2,    0.68,  -1.065708, -1.065708},
	{2.4,  0,     0,         0        }
};

const std::vector<std::vector<double> > pantherAtm = {
	{0,          0,         1.666667,  1.666667},
	{0.07066164, 0.1397133, 1.961396,  1.961396},
	{0.34,       0.56,      1.084002,  1.084002},
	{1,          1,         0.5302638, 0.5302638}
};

const std::vector<std::vector<double> > pantherMach = {
	{0,    1,    0,          0        },
	{0.18, 0.97, 0,          0        },
	{0.43, 1,    0.202683,   0.202683 },
	{1,    1.42, 1.280302,   1.280302 },
	{2.5,  3.63, 0,          0        },
	{3,    0.58, -2.708558,  -2.708558},
	{3.35, 0,    -0.6150925, 0        }
};

double lambertw(double x) {
	double guess = 0.43 * std::log(0.15 * x * x + 1.95 * x + 0.78) + 0.16;

	for (int i = 0; i < 15; i++) {
		guess += (x - guess * std::exp(guess)) / (guess + 1) / std::exp(guess);
	}

	return guess;
}

struct JetEngine {
	double mass;
	double thrust;
	double isp;
	double cap;

	std::vector<std::vector<double> > machCurve;
	std::vector<std::vector<double> > atmCurve;
};

struct RocketEngine {
	double mass;
	double thrust;
	double fuelconsumption;

	RocketEngine(double newmass, double newthrust, double isp): 
		mass(newmass), 
		thrust(newthrust), 
		fuelconsumption(newthrust / isp / g0) {}
};

const std::vector<RocketEngine> engines = {
	RocketEngine(500, 60e3, 345), // max = 
	RocketEngine(130, 20e3, 320), // max = 
	RocketEngine(20,  2e3,  315), // max = 6
	RocketEngine(180, 32e3, 310), // max = 
	RocketEngine(80,  16e3, 290)  // max = 
};

struct ValueCurve {
	std::vector<double> xs;
	std::vector<double> ys;

	double operator()(double x) const {
		std::size_t numXs = xs.size();

		if (x <= xs[0]) return ys[0];
		if (x >= xs[numXs - 1]) return ys[numXs - 1];

		std::size_t low = 0;
		std::size_t high = numXs - 1;
		while (high - low > 1) {
			std::size_t mid = (low + high) / 2;

			if (xs[mid] >= x) high = mid;
			else if (xs[mid] <= x) low = mid;
		}

		double lowx = xs[low];
		double lowy = ys[low];
		double highx = xs[high];
		double highy = ys[high];

		return lowy + (highy - lowy) * (x - lowx) / (highx - lowx);
	}
};

double readHermiteSpline(double x, const std::vector<std::vector<double> > spline) {
	/*
	Reads a curve defined by cubic Hermite splines.

	Returns:
	  - double y - y-value of curve

	Parameters:
	  - double x - x-value of curve
	  - std::vector<std::vector<double> > spline - curve, defined as a set of cubic Hermite splines.
	    - Each point contains:
	      - x - x position
	      - y - y position
	      - y'_1 - left derivative
	      - y'_2 - right derivative
	*/

	std::size_t nKeys = spline.size();
	double minX = spline[0][0]; //Edge cases
	double maxX = spline[nKeys - 1][0];
	if (x <= minX) {
		return spline[0][1] + (x - minX) * spline[0][2]; //Linear extension of spline
	}
	if (x >= maxX) {
		return spline[nKeys - 1][1] + (x - maxX) * spline[nKeys - 1][3]; //Linear extension of spline
	}

	std::size_t i; //Determining which region of function we are in
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
	return 1000.0 * readHermiteSpline(alt, kerbinPressureKeys);
}

double maxTemperatureAlt(double alt) {
	/*
	Calculates temperature in Kerbin's atmosphere at a given altitude.

	Returns:
	  - double temperature - temperature of atmosphere in kelvins.

	Parameters:
	  - double altitude - altitude above sea level in meters.
	*/
	double baseTemp = readHermiteSpline(alt, kerbinTempBaseKeys);
	double sunMult = readHermiteSpline(alt, kerbinSunMultKeys);
	double latitudeBias = readHermiteSpline(0.0, kerbinLatitudeBiasKeys);
	double latitudeSunMult = readHermiteSpline(0.0, kerbinLatitudeSunMultKeys);

	return baseTemp + sunMult * (latitudeBias + latitudeSunMult);
}

double densityAlt(double alt) {
	return KERBIN_MOLAR_MASS * pressureAlt(alt) / (R * maxTemperatureAlt(alt));
}

double soundSpeedAlt(double alt) {
	return std::sqrt(KERBIN_ADIABATIC_INDEX * R * maxTemperatureAlt(alt) / KERBIN_MOLAR_MASS);
}

double LDratio(double mach, double aoa) {
	double sinAoa = std::sin(aoa);
	double cosAoa = std::cos(aoa);

	double dragAoaMult = readHermiteSpline(sinAoa, dragAoAKeys) / 2.4;
	double liftAoaMult = readHermiteSpline(sinAoa, liftAoAKeys);
	double dragMachMult = readHermiteSpline(mach, dragMachKeys);
	double liftMachMult = readHermiteSpline(mach, liftMachKeys);

	double dragMult = dragAoaMult * dragMachMult;
	double liftMult = liftAoaMult * liftMachMult;

	return liftMult * cosAoa / (liftMult * sinAoa + dragMult);
}

double maxLDratio(double mach) {
	double a = 0;
	double d = 0.17;

	while (d - a > 1e-3) {
		double b = a + (d - a) * (1 - inversePhi);
		double c = a + (d - a) * inversePhi;

		if (LDratio(mach, b) > LDratio(mach, c)) {
			d = c;
		} else {
			a = b;
		}
	}

	return LDratio(mach, 0.5 * (a + d));
}

double flowcap(double x, double cap) {
	if (x <= cap) return x;
	return 2.0 * x * cap / (x + cap);
}

double jetEngineThrust(double vel, double alt, double thrust, double cap, std::vector<std::vector<double> > machCurve, std::vector<std::vector<double> > atmCurve) {
	double soundspeed = soundSpeedAlt(alt);
	double density = densityAlt(alt);

	double mach = vel / soundspeed;
	double standarddensity = KERBIN_MOLAR_MASS * pressureAlt(0) / (R * readHermiteSpline(0, kerbinTempBaseKeys));

	double machMult = readHermiteSpline(mach, machCurve);
	double densityMult = readHermiteSpline(density / standarddensity, atmCurve);

	double totalMult = machMult * densityMult;
	return thrust * flowcap(totalMult, cap);
}

struct JetResult {
	double mass;
	double vel;
};

JetResult simulateJet(double initmass, JetEngine jet, const ValueCurve& ldcurve) {
	double mass = initmass;
	double vel = 0;
	double effectiveisp;

	const double timestep = 0.1;

	const double soundspeed = soundSpeedAlt(0);

	do {
		double orbitalvel = vel + KERBIN_ROT_VEL;
		double effectivegravity = KERBIN_GRAVITY - orbitalvel * orbitalvel / KERBIN_RADIUS;

		double mach = vel / soundspeed;
		double ldratio = ldcurve(mach);

		double thrust = jetEngineThrust(vel, 0, jet.thrust, jet.cap, jet.machCurve, jet.atmCurve);
		double fuelconsumption = thrust / jet.isp / g0;

		double accel = thrust / mass;

		//First definition is optimal direction, second is horizontal burn
		//double netaccel = (ldratio * ldratio * accel - std::abs(accel - effectivegravity * std::sqrt(ldratio * ldratio + 1))) / ldratio / std::sqrt(ldratio * ldratio + 1);
		double netaccel = accel - effectivegravity / ldratio;
		effectiveisp = netaccel * mass / fuelconsumption / g0;

		vel += netaccel * timestep;
		mass -= fuelconsumption * timestep;
	} while (effectiveisp > 300); //Rough approximation for rocket isp

	double liquidfuel = mass - initmass;
	double tankmass = std::ceil(liquidfuel / 250) * 25;
	mass -= tankmass + jet.mass + 10; //10 - Intake

	JetResult out;
	out.mass = mass;
	out.vel = vel;

	return out;
}

double simulateRocket(double initmass, double initvel, const ValueCurve& ldcurve, double finalapo, double finalperi, double finalalt, double rocketthrust, double rocketisp, double rocketmass) {
	double mass = initmass;

	const double soundspeed = soundSpeedAlt(0);

	double sma = KERBIN_RADIUS + 0.5 * (finalapo + finalperi);
	double eccentricity = (KERBIN_RADIUS + finalapo) / sma - 1;

	double finalvel = std::sqrt(KERBIN_STDGP * (2 / (KERBIN_RADIUS + finalalt) - 1 / sma));
	double finalhvel = std::sqrt((1 - eccentricity * eccentricity) * sma * KERBIN_STDGP) / (KERBIN_RADIUS + finalalt);
	double finalvvel = std::sqrt(finalvel * finalvel - finalhvel * finalhvel);

	finalhvel -= 2 * pi * (KERBIN_RADIUS + finalalt) / KERBIN_ROTATIONAL_PERIOD;
	finalvel = std::sqrt(finalhvel * finalhvel + finalvvel * finalvvel);

	double dv = finalvel - initvel;
	double cosfinalangle = finalhvel / finalvel;
	double sinfinalangle = finalvvel / finalvel;

	double totaldrag = KERBIN_GRAVITY * (sinfinalangle + cosfinalangle / ldcurve(finalvel / soundspeed));

	if (rocketthrust / mass < totaldrag) return 0;

	double finalmass = -rocketthrust * lambertw(-mass * totaldrag / (rocketthrust * std::exp(dv / rocketisp / g0 + mass * totaldrag / rocketthrust))) / totaldrag;

	double fuel = mass - finalmass;
	mass = finalmass - fuel / 8 - rocketmass - 50 - 150; //50 - Wing, 150 - Fairing

	return mass;
}

double optimize(double initmass, JetEngine jet, ValueCurve ldcurve, double finalapo, double finalperi, double finalalt) {
	JetResult rocketphasestate = simulateJet(initmass, jet, ldcurve);

	//Simple brute force
	std::size_t numengines = engines.size();

	std::vector<int> maxes = {3, 3, 7, 3, 3};

	int numiters = 1;
	for (std::size_t i = 0; i < numengines; i++) {
		numiters *= maxes[i];
	}

	double bestmass = 0;
	for (int i = 0; i < numiters; i++) {
		double rocketthrust = 0;
		double rocketmass = 0;
		double rocketfuelconsumption = 0;

		std::size_t div = 1;
		for (std::size_t j = 0; j < numengines; j++) {
			std::size_t enginemult = (i / div) % maxes[j];

			rocketthrust += enginemult * engines[j].thrust;
			rocketmass += enginemult * engines[j].mass;
			rocketfuelconsumption += enginemult * engines[j].fuelconsumption;

			div *= maxes[j];
		}

		double mass = simulateRocket(
			rocketphasestate.mass,
			rocketphasestate.vel,
			ldcurve,
			finalapo, 
			finalperi, 
			finalalt, 
			rocketthrust, 
			rocketthrust / rocketfuelconsumption / g0, 
			rocketmass
		);

		if (mass > bestmass) {
			//std::cout << i << " " << mass << " opt\n";
			bestmass = mass;
		}
	}

	return bestmass;
}

int main() {
	/*=== Initialization ===*/
	ValueCurve ldcurve;

	for (double mach = 0; mach < 25; mach += 0.01) {
		ldcurve.xs.push_back(mach);
		ldcurve.ys.push_back(maxLDratio(mach));
	}

	JetEngine juno;
	juno.mass = 250;
	juno.thrust = 20000;
	juno.isp = 6400;
	juno.cap = 10000;
	juno.machCurve = junoMach;
	juno.atmCurve = junoAtm;

	JetEngine panther;
	panther.mass = 1200;
	panther.thrust = 130000;
	panther.isp = 4000;
	panther.cap = 1.1;
	panther.machCurve = pantherMach;
	panther.atmCurve = pantherAtm;

	/*=== Test ===*/
	for (double mass = 1000; mass < 10000; mass += 100) {
		double finalmass = optimize(mass, juno, ldcurve, 9e4, -3e5, 15e3);
		double finalmasspanther = optimize(mass, panther, ldcurve, 9e4, -3e5, 15e3);

		std::cout << mass << " " << finalmass << " " << finalmasspanther << " fin\n";
	}
}
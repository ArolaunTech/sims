#include <cmath>

#include "curves.h"
#include "configreader.h"

#ifndef AERO_H
#define AERO_H

const HermiteSpline liftAoACurve(std::vector<std::array<double, 4> >{ //Thanks to Lt_Duckweed
	{0.0,       0.0,        0.0,        1.965926 },
	{0.258819,  0.5114774,  1.990092,   1.905806 },
	{0.5,       0.9026583,  0.7074468, -0.7074468},
	{0.7071068, 0.5926583, -2.087948,  -1.990095 },
	{1.0,       0.0,       -2.014386,  -2.014386 }
});

const HermiteSpline dragAoACurve(std::vector<std::array<double, 4> >{
	{0.0,       0.01, 0.0,       0.0      },
	{0.3420201, 0.06, 0.1750731, 0.1750731},
	{0.5,       0.24, 2.60928,   2.60928  },
	{0.7071068, 1.7,  3.349777,  3.349777 },
	{1.0,       2.4,  1.387938,  0.0      }
});

const HermiteSpline liftMachCurve(std::vector<std::array<double, 4> >{
	{ 0.0, 1.0,    0.0,           0.0       },
	{ 0.3, 0.5,   -1.671345,     -0.8273422 },
	{ 1.0, 0.125, -0.0005291355, -0.02625772},
	{ 5.0, 0.0625, 0.0,           0.0       },
	{25.0, 0.05,   0.0,           0.0       }
});

const HermiteSpline dragMachCurve(std::vector<std::array<double, 4> >{
	{ 0.0,  0.35,  0.0,         -0.8463008 },
	{ 0.15, 0.125, 0.0,          0.0       },
	{ 0.9,  0.275, 0.541598,     0.541598  },
	{ 1.1,  0.75,  0.0,          0.0       },
	{ 1.4,  0.4,  -0.3626955,   -0.3626955 },
	{ 1.6,  0.35, -0.1545923,   -0.1545923 },
	{ 2.0,  0.3,  -0.09013031,  -0.09013031},
	{ 5.0,  0.22,  0.0,          0.0       },
	{25.0,  0.3,   0.0006807274, 0.0       }
});

double LDratio(double mach, double aoa) {
	if (aoa < 0) return -LDratio(mach, -aoa);

	double sinAoA = std::sin(aoa);
	double cosAoA = std::cos(aoa);

	double dragAoAMult = dragAoACurve.evaluate(sinAoA) / 2.4;
	double liftAoAMult = liftAoACurve.evaluate(sinAoA);
	double dragMachMult = dragMachCurve.evaluate(mach);
	double liftMachMult = liftMachCurve.evaluate(mach);

	double dragMult = dragAoAMult * dragMachMult;
	double liftMult = liftAoAMult * liftMachMult;

	return liftMult * cosAoA / (liftMult * sinAoA + dragMult);
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

double lift(
	const Body& body,
	double alt,
	double vel,
	double wingarea,
	double aoa
) {
	// Lift assuming 3pm (maximum temperature) and default settings.
	if (!body.atmospheric) return 0;
	if (aoa < 0) return -lift(body, alt, vel, wingarea, -aoa);

	double sinAoA = std::sin(aoa);
	double cosAoA = std::cos(aoa);

	double density = body.atmosphere.mindensity(alt);
	double dynamicpressure = 0.5 * density * vel * vel;
	
	double soundspeed = body.atmosphere.maxsoundspeed(alt);
	double mach = vel / soundspeed;

	double liftMachMult = liftMachCurve.evaluate(mach);
	double liftAoAMult = liftAoACurve.evaluate(sinAoA);

	double liftMult = liftMachMult * liftAoAMult;

	return dynamicpressure * liftMult * wingarea * 36 * cosAoA;
}

double drag(
	const Body& body,
	double alt,
	double vel,
	double wingarea,
	double aoa
) {
	// Drag assuming 3pm (maximum temperature), default settings, and no body drag.
	if (!body.atmospheric) return 0;

	double sinAoA = std::sin(aoa);

	double density = body.atmosphere.mindensity(alt);
	double dynamicpressure = 0.5 * density * vel * vel;

	double soundspeed = body.atmosphere.maxsoundspeed(alt);
	double mach = vel / soundspeed;

	double liftMachMult = liftMachCurve.evaluate(mach);
	double liftAoAMult = liftAoACurve.evaluate(sinAoA);

	double dragMachMult = dragMachCurve.evaluate(mach);
	double dragAoAMult = dragAoACurve.evaluate(sinAoA);

	double liftMult = liftMachMult * liftAoAMult;
	double dragMult = dragMachMult * dragAoAMult;

	return dynamicpressure * liftMult * wingarea * 36 * sinAoA + dynamicpressure * dragMult * wingarea * 15;
}

#endif
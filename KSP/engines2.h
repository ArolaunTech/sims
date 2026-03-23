#include "curves.h"
#include "consts.h"

#ifndef ENGINES2
#define ENGINES2

struct JetEngine {
	int id;

	double thrust;
	double isp;
	double flowmultcap;
	double mass;

	HermiteSpline velCurve;
	HermiteSpline atmCurve;

	double flowcap(double x) const {
		if (x <= flowmultcap) return x;
		return flowmultcap + (x - flowmultcap) / (2 + (x - flowmultcap) / flowmultcap);
	}

	double getthrust(const Body& body, double vel, double alt) const {
		if (!body.atmospheric) return 0;
		if (!body.atmosphere.oxygen) return 0;

		double soundspeed = body.atmosphere.maxsoundspeed(alt);
		double density = body.atmosphere.mindensity(alt);

		double mach = vel / soundspeed;

		double machMult = velCurve.evaluate(mach);
		double densityMult = atmCurve.evaluate(density / standarddensity);

		double totalMult = machMult * densityMult;
		return thrust * flowcap(totalMult);
	}

	double getfuelconsumption(const Body& body, double vel, double alt) const {
		return getthrust(body, vel, alt) / isp / g0;
	}
};

#endif
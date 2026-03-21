#include <cmath>

#include "curves.h"

#ifndef SYSTEM_H
#define SYSTEM_H

struct Atmosphere {
	bool oxygen;

	double adiabaticindex;
	double molarmass;

	bool normalizedpressurecurve;
	bool normalizedtemperaturecurve;

	double maxalt;

	HermiteSpline pressurecurve;
	HermiteSpline temperaturecurve;
	HermiteSpline tempsunmultcurve;
	HermiteSpline templatitudebiascurve;
	HermiteSpline templatitudesunmultcurve;
	HermiteSpline tempaxialsunbiascurve;
	HermiteSpline tempaxialsunmultcurve;
	HermiteSpline tempeccentricitybiascurve;

	double pressure(double alt) const {
		if (normalizedpressurecurve) return 1000 * pressurecurve.evaluate(alt / maxalt);
		return 1000 * pressurecurve.evaluate(alt);
	}

	double maxtemp(double alt) const {
		double x = alt;

		if (normalizedtemperaturecurve) x /= maxalt;

		return temperaturecurve.evaluate(x) + tempsunmultcurve.evaluate(x) * (
			templatitudebiascurve.evaluate(0) +
			templatitudesunmultcurve.evaluate(0) +
			tempaxialsunbiascurve.evaluate(0) * tempaxialsunmultcurve.evaluate(0) +
			tempeccentricitybiascurve.evaluate(0)
		);
	}

	double mindensity(double alt) const {
		return molarmass * pressure(alt) / (R * maxtemp(alt));
	}

	double maxsoundspeed(double alt) const {
		return std::sqrt(adiabaticindex * R * maxtemp(alt) / molarmass);
	}
};

struct Body {
	std::string name;

	double radius;
	double GM;

	double rotperiod;

	bool atmospheric;
	Atmosphere atmosphere;
};

struct System {
	std::unordered_map<std::string, Body> bodies;
};

#endif
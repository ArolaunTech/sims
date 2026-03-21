#include <iostream>
#include <vector>
#include <array>

#include "systemloader.h"
#include "aero.h"
#include "curves.h"
#include "consts.h"
#include "partloader.h"

// Juno
const HermiteSpline jetAtm(std::vector<std::array<double, 4> >{
	{0,     0,     0,         0.7448742},
	{0.072, 0.13,  2.075459,  2.075459 },
	{0.16,  0.28,  1.464173,  1.464173 },
	{0.42,  0.578, 0.93687,   0.93687  },
	{1,     1,     0.5529748, 0        }
});

const HermiteSpline jetMach(std::vector<std::array<double, 4> >{
	{0,    1,     0,         0        },
	{0.44, 0.897, 0,         0        },
	{1,    1,     0.1988732, 0.1988732},
	{1.3,  1.03,  0,         0        },
	{2,    0.68,  -1.065708, -1.065708},
	{2.4,  0,     0,         0        }
});

double flowcap(double x, double cap) {
	if (x <= cap) return x;
	return cap + (x - cap) / (2 + (x - cap) / cap);
}

double jetEngineThrust(
	const Body& body, 
	double vel, 
	double alt, 
	double thrust, 
	double cap, 
	const HermiteSpline& machCurve, 
	const HermiteSpline& atmCurve
) {
	double soundspeed = body.atmosphere.maxsoundspeed(alt);
	double density = body.atmosphere.mindensity(alt);

	double mach = vel / soundspeed;

	double machMult = machCurve.evaluate(mach);
	double densityMult = atmCurve.evaluate(density / standarddensity);

	double totalMult = machMult * densityMult;
	return thrust * flowcap(totalMult, cap);
}

struct Result {
	double score;
};

Result simulate(const Body& body, const std::vector<std::array<double, 2> >& code) {
	Result out;

	out.score = 0;

	return out;
}

int main() {
	PartCategories catalog = loadParts("./parts");

	System stock = loadSystem("./configs");
	Body kerbin = stock.bodies["Kerbin"];

	Result curr = simulate(kerbin, std::vector<std::array<double, 2> >{{0, 0}});
}
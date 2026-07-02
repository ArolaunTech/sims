#include <iostream>
#include <cmath>

#include "aero.h"
#include "consts.h"
#include "partloader.h"
#include "systemloader.h"
#include "engines.h"
#include "vector.h"

Vector2D calcForce(
	const Body& body,
	const JetEngine& jet, 
	double mass, 
	double vel, 
	double alt, 
	double wingarea, 
	double incidence, 
	double aoa
) {
	double hforce = 
		jet.getthrust(body, vel, alt) * std::cos(aoa) - 
		drag(body, alt, vel, wingarea, incidence + aoa);
	double vforce = 
		lift(body, alt, vel, wingarea, incidence + aoa) +
		std::pow(vel + 2 * pi * (body.radius + alt) / body.rotperiod, 2) * mass / (body.radius + alt) +
		jet.getthrust(body, vel, alt) * std::sin(aoa) -
		body.GM * mass / std::pow(body.radius + alt, 2);

	Vector2D out;

	out.x = hforce;
	out.y = vforce;

	return out;
}

struct aoanetforce {
	double aoa;
	double netforce;
};

aoanetforce calcNetForce(
	const Body& body,
	const JetEngine& jet, 
	double mass, 
	double vel, 
	double alt, 
	double wingarea, 
	double incidence
) {
	double lo = -pi / 6;
	double hi = pi / 6;

	for (int it = 0; it < 30; it++) {
		double mid = (lo + hi) / 2;

		Vector2D force = calcForce(body, jet, mass, vel, alt, wingarea, incidence, mid);

		if (std::abs(force.y) < 1) break;
		if (force.y < 0) lo = mid;
		if (force.y > 0) hi = mid;	
	}

	double mid = (lo + hi) / 2;

	Vector2D force = calcForce(body, jet, mass, vel, alt, wingarea, incidence, mid);

	aoanetforce out;

	out.aoa = mid;
	out.netforce = force.x;

	return out;
}

int main() {
	PartCategories catalog = loadParts("./parts");
	System stock = loadSystem("./configs");
	Body kerbin = stock.bodies["Kerbin"];

	JetEngine jet;

	for (const auto& jetpair : catalog.jets) {
		if (jetpair.first == "Juno") jet = jetpair.second;
	}

	const double initmass = 1100;
	double mass = initmass;
	double finalt = 5000;
	double vel = 100;
	double wingarea = .12;
	double incidence = -5 * pi / 180;

	double angle = 0;
	double alt = 0;

	const double dt = 1;

	for (int i = 0; i < 100000; i++) {
		aoanetforce net = calcNetForce(kerbin, jet, mass, vel, alt, wingarea, incidence);
		double accel = net.netforce / mass;

		mass -= jet.getfuelconsumption(kerbin, vel, alt) * dt;
		vel += accel * dt;
		angle += vel / (kerbin.radius + alt) * dt;

		if (vel > 600) alt = finalt;
		if (angle > 2 * pi) break;

		if (i % 25 == 0) {
			std::cout << i * dt << " " << vel << " " << net.aoa * 180 / pi << " " << mass << " " << angle << "\n";
		}
	}

	std::cout << initmass << " " << mass << "\n";
	std::cout << .2 * (initmass - mass) << "\n";
}
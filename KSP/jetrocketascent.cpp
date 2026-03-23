#include <iostream>
#include <vector>
#include <array>
#include <cmath>

#include "systemloader.h"
#include "aero.h"
#include "curves.h"
#include "consts.h"
#include "partloader.h"
#include "engines.h"
#include "random.h"

struct InputParams {
	std::vector<std::array<double, 3> > code;

	double mass;
	double takeoffspeed;
	double wingarea;
	double friction;
	double wingincidence;
};

struct Result {
	double mass;
	double score;
};

Result simulate(const Body& body, const JetEngine& jet, const Engine& rocket, const InputParams& inputs, bool log = false) {
	const double dt = 0.1;

	double mass = inputs.mass;
	double vel = 0;

	// Takeoff
	while(vel < inputs.takeoffspeed) {
		double thrust = jet.getthrust(body, vel, 100);
		double accel = (thrust - inputs.friction) / mass;
		double fc = thrust / jet.isp / g0;

		if (accel < 0) {
			std::cerr << "Not enough thrust!\n";
			break;
		}

		vel += accel * dt;
		mass -= fc * dt;
	}

	// Flight
	double time = 0;
	double x = 0;
	double y = body.radius + 100;
	double vx = 2 * pi * (body.radius + 100) / body.rotperiod + inputs.takeoffspeed;
	double vy = 0;

	int currCodeIndex = 0;
	double nextSwitch = inputs.code[0][1];

	double score;

	double sright, sup, svel;
	double alt, upangle;
	double apoapsis, periapsis;

	for (int i = 0; i < 3000; i++) {
		// Gravity
		double dist = std::sqrt(x * x + y * y);
		double grav = body.GM / dist / dist;

		double gx = -x * grav / dist;
		double gy = -y * grav / dist;

		// Aero
		if (time > nextSwitch && currCodeIndex < inputs.code.size() - 1) {
			currCodeIndex++;
			nextSwitch += inputs.code[currCodeIndex][1];
		}

		upangle = inputs.code[currCodeIndex][0] + inputs.wingincidence;

		alt = dist - body.radius;

		double surfx = 2 * pi * y / body.rotperiod;
		double surfy = -2 * pi * x / body.rotperiod;

		double sx = vx - surfx;
		double sy = vy - surfy;

		svel = std::sqrt(sx * sx + sy * sy);
		sup = (sx * x + sy * y) / dist;
		sright = (sx * y - sy * x) / dist;

		double sangle = std::atan2(sup, sright);

		double liftforce = lift(body, alt, svel, inputs.wingarea, upangle - sangle);
		double dragforce = drag(body, alt, svel, inputs.wingarea, upangle - sangle);

		double aerox = -sx * dragforce / (mass * svel) - sy * liftforce / (mass * svel);
		double aeroy = -sy * dragforce / (mass * svel) + sx * liftforce / (mass * svel);

		// Thrust
		double thrust = jet.getthrust(body, svel, alt) + rocket.thrust * inputs.code[currCodeIndex][2];
		double accel = thrust / mass;
		double fc = jet.getthrust(body, svel, alt) / jet.isp / g0 + rocket.thrust * inputs.code[currCodeIndex][2] / rocket.isp / g0;

		double ax = accel * (std::sin(upangle) * x + std::cos(upangle) * y) / dist;
		double ay = accel * (std::sin(upangle) * y - std::cos(upangle) * x) / dist;

		// Orbit
		double energy = (vx * vx + vy * vy) / 2 - body.GM / dist;
		double sma = -body.GM / (2 * energy);

		double h = vx * y - vy * x;
		double e = std::sqrt(1 - h * h / body.GM / sma);

		apoapsis = sma * (1 + e) - body.radius;
		periapsis = sma * (1 - e) - body.radius;

		if (alt < 0) break;
		if (apoapsis > 90000 && periapsis > -380000) {
			score = mass;
			break;
		}

		// Log
		if (i % 100 == 0 && log) {
			std::cout << "Time " << time << " s:\n";
			std::cout << " - Pos: <" << x << ", " << y << "> m\n";
			std::cout << " - Vel: <" << vx << ", " << vy << "> m/s\n";
			std::cout << " - Surf. vel: <" << sright << ", " << sup << "> m/s (" << svel << " m/s total)\n";
			std::cout << " - Alt: " << alt << " m\n";
			std::cout << " - Mass: " << mass << " kg\n";
			std::cout << " - Upangle: " << upangle * 180 / pi << " deg\n";
			std::cout << " - Apoapsis: " << apoapsis << " m, periapsis: " << periapsis << " m\n";
			std::cout << "\n";
		}

		// Update
		time += dt;
		vx += (gx + aerox + ax) * dt;
		vy += (gy + aeroy + ay) * dt;
		x += vx * dt;
		y += vy * dt;
		mass -= fc * dt;

		score = periapsis;
	}

	if (log) {
		std::cout << "Time " << time << " s:\n";
		std::cout << " - Pos: <" << x << ", " << y << "> m\n";
		std::cout << " - Vel: <" << vx << ", " << vy << "> m/s\n";
		std::cout << " - Surf. vel: <" << sright << ", " << sup << "> m/s (" << svel << " m/s total)\n";
		std::cout << " - Alt: " << alt << " m\n";
		std::cout << " - Mass: " << mass << " kg\n";
		std::cout << " - Upangle: " << upangle * 180 / pi << " deg\n";
		std::cout << " - Apoapsis: " << apoapsis << " m, periapsis: " << periapsis << " m\n";
		std::cout << "\n";
	}

	Result out;

	out.score = score;
	out.mass = mass;

	return out;
}

int main() {
	PartCategories catalog = loadParts("./parts");
	JetEngine juno;

	for (const auto& jetpair : catalog.jets) {
		if (jetpair.first == "Juno") juno = jetpair.second;
	}

	Engine rocket = engines.find("Spark")->second;
	rocket.thrust *= 2;

	System stock = loadSystem("./configs");
	Body kerbin = stock.bodies["Kerbin"];

	InputParams inputs, best;

	inputs.code = std::vector<std::array<double, 3> >{{0.25, 20, 0}, {0.15, 20, 0}, {0.095, 40, 0}, {0.055, 110, 0}, {0.1, 5, 1}, {0.15, 5, 1}, {0.2, 5, 1}, {0.25, 5, 1}, {0.3, 5, 1}, {0.35, 5, 1}, {0.4, 5, 1}, {0.5, 5, 1}};
	inputs.mass = 3000;
	inputs.takeoffspeed = 50;
	inputs.wingarea = 0.5;
	inputs.friction = 15000;
	inputs.wingincidence = 0 * pi / 180;

	best = inputs;

	Result bestout = simulate(kerbin, juno, rocket, inputs);

	std::cout << bestout.score << "\n";

	for (int i = 0; i < 1000000; i++) {
		inputs = best;
		for (std::size_t j = 0; j < inputs.code.size(); j++) {
			inputs.code[j][0] += normal(generator) * 0.001;
			inputs.code[j][1] += normal(generator) * 1;

			if (inputs.code[j][1] < 0) inputs.code[j][1] = 0;

			inputs.code[j][2] += normal(generator) * 0.01;

			if (inputs.code[j][2] < 0) inputs.code[j][2] = 0;
			if (inputs.code[j][2] > 1) inputs.code[j][2] = 1;
		}

		Result currout = simulate(kerbin, juno, rocket, inputs);

		if (currout.score > bestout.score) {
			std::cout << i << " " << currout.score << "\n";

			bestout = currout;
			best = inputs;

			for (std::size_t j = 0; j < best.code.size(); j++) {
				std::cout << best.code[j][0] << " " << best.code[j][1] << " " << best.code[j][2] << "\n";
			}

			simulate(kerbin, juno, rocket, inputs, true);

			std::cout << "\n";
		}
	}
}
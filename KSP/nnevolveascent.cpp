#include <iostream>
#include <string>
#include <cmath>
#include <utility>
#include <algorithm>

#include "consts.h"
#include "craft.h"
#include "engines.h"
#include "random.h"
#include "body.h"
#include "nn.h"

const double timestep = 0.1;

double simulate_ascent(const Body body, const Craft init_craft, double init_alt, const NN agent) {
	/*=== Constants ===*/
	double total_thrust = 0;
	for (auto const & iterator : init_craft.engines) {
		total_thrust += iterator.first.thrust * iterator.second;
	}

	/*=== Initial conditions ===*/
	double x = body.radius + init_alt;
	double y = 0;
	double vx = 0;
	double vy = 2 * pi * x / body.rotational_period;
	Craft craft = init_craft;

	/*=== Iterate ===*/
	double periapsis = 0;
	double dist = 0;
	while (periapsis < body.radius) {
		/*=== NN inputs ===*/
		dist = std::sqrt(x * x + y * y);
		double vh = (-y * vx + x * vy) / dist;
		double vv = ( x * vx + y * vy) / dist;

		std::vector<double> inputs = {
			(dist - body.radius) * 1e-3,
			vh * 1e-3,
			vv * 1e-3,
			std::log(craft.mass),
			craft.mass * 1e-6
		};

		/*=== Action ===*/
		std::vector<double> action = agent.evaluate(inputs);

		double ah = action[0] / std::sqrt(action[0] * action[0] + action[1] * action[1]);
		double av = action[1] / std::sqrt(action[0] * action[0] + action[1] * action[1]);

		double thrust = total_thrust / (1 + std::exp(-0.1 * action[2]));
		double accel = thrust / craft.mass;

		double tick_consumption = 0;
		for (auto const & iterator : craft.engines) {
			double group_thrust = iterator.first.thrust * iterator.second;

			if (group_thrust > thrust) {
				tick_consumption += thrust / iterator.first.isp / g0;
				break;
			}

			tick_consumption += group_thrust / iterator.first.isp / g0;
			thrust -= group_thrust;
		}

		/*=== Calculate accelerations ===*/
		ah *= accel;
		av *= accel;

		double ax = (-y * ah + x * av) / dist;
		double ay = ( x * ah + y * av) / dist;

		double grav = body.grav_parameter / dist / dist;
		double gx = -x * grav / dist;
		double gy = -y * grav / dist;

		/*=== Integrate ===*/
		vx += (ax + gx) * timestep;
		vy += (ay + gy) * timestep;
		x += vx * timestep;
		y += vy * timestep;
		craft.mass -= tick_consumption * timestep;

		/*=== Calculate orbit ===*/
		double sma = 1 / (2 / dist - (vx * vx + vy * vy) / body.grav_parameter);
		double eccentricity = std::sqrt(1 - (vh * vh * dist * dist) / (sma * body.grav_parameter));

		periapsis = sma * (1 - eccentricity);

		//std::cout << x << " " << y << " " << vx << " " << vy << " " << dist << " " << vh << " " << sma << " " << eccentricity << " " << periapsis << "\n";

		/*=== Failure conditions ===*/
		if (dist < body.radius) break;
		if (eccentricity > 1.0) break;
	}

	//std::cout << x << " " << y << " " << vx << " " << vy << " " << dist << " " << periapsis << "\n";

	return (periapsis < body.radius) ? (periapsis / body.radius) : (1 + std::exp(craft.payload_fraction()));
}

int main() {
	/*=== Constants ===*/
	std::vector<int> layer_num_nodes = {5, 32, 32, 32, 32, 32, 2};
	int num_agents = 1000;

	Craft craft;
	craft.mass = 100000;
	craft.engines[engines.at("Skipper")] = 1.5;

	Body mun {65138397520.7807, 2e5, 138984.376574476};

	/*=== Population ===*/
	std::vector<std::pair<double, NN> > population;

	for (int i = 0; i < num_agents; i++) {
		population.push_back(std::pair<double, NN>(0, get_random_nn(layer_num_nodes, Activation::ACTIVATION_LOGISTIC)));
	}

	/*=== Optimize ===*/
	for (int i = 0; i < 10; i++) {
		//Scoring population
		for (int j = 0; j < num_agents; j++) {
			double fitness = simulate_ascent(mun, craft, 5000, population[i].second);

			//std::cout << fitness << "\n";

			population[j].first = fitness;
		}
		//std::cout << std::endl;

		std::sort(population.begin(), population.end());

		for (int j = 0; j < num_agents; j++) {
			std::cout << population[j].first << "\n";
		}
		std::cout << std::endl;

		for (int j = 0; j < num_agents / 10; j++) {
			population[j] = population[num_agents - 1 - j];
			population[j].second.mutate();
		}

		//std::cout << craft.fuel_used[FuelType::FUEL_LFOX] << "\n";
	}
}
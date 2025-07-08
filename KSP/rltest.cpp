#include <iostream>
#include <vector>
#include <cmath>
#include <map>

#include "consts.h"
#include "craft.h"
#include "engines.h"
#include "random.h"
#include "body.h"
#include "nn.h"

//Constants
const double timestep = 0.1;

/*=== Simulation ===*/
struct SARGroup {
	std::vector<std::vector<double> > states;
	std::vector<std::vector<double> > actions;
	std::vector<double> rewards;
};

struct Result {
	bool valid;
	Craft state;
	double score;

	std::vector<SARGroup> critic_data;
};

Result simulate_ascent(Body body, Craft craft, double landing_alt, double orbit_alt, bool fuel_reverse, NN agent) {
	/*=== Constants ===*/
	double total_thrust = 0;
	for (auto const & iterator : craft.engines) {
		total_thrust += iterator.first.thrust * iterator.second;
	}

	/*=== Initial conditions ===*/
	double x = body.radius + landing_alt;
	double y = 0;
	double vx = 0;
	double vy = 2 * pi * x / body.rotational_period;

	/*=== NN inputs ===*/
	std::vector<double> nn_inputs = {
		0,
		0,
		0,
		0,
		body.radius * 1e-5, //Scaled to reduce magnitude for network
		std::log(body.grav_parameter),
		std::sqrt(body.grav_parameter / body.radius) * 1e-3,
		body.grav_parameter / body.radius / body.radius,
		2 * pi * body.radius / body.rotational_period,
		fuel_reverse ? -1.0 : 1.0
	};

	//Add engine inputs
	for (auto const & iterator : engines) {
		if (craft.engines.count(iterator.second) == 0) {
			nn_inputs.push_back(0);
			continue;
		}

		nn_inputs.push_back(craft.engines.at(iterator.second));
	}

	/*=== Recording ===*/
	SARGroup sar;

	/*=== Iterate ===*/
	double periapsis = 0;
	double sma = 0;
	while (periapsis < body.radius) {
		/*=== Fill NN inputs ===*/
		double dist = std::sqrt(x * x + y * y);
		double vh = (-y * vx + x * vy) / dist;
		double vv = ( x * vx + y * vy) / dist;

		nn_inputs[0] = dist;
		nn_inputs[1] = vh;
		nn_inputs[2] = vv;
		nn_inputs[3] = std::log(craft.mass);

		/*=== Find Action ===*/
		std::vector<double> action = agent.evaluate(nn_inputs);

		sar.states.push_back(nn_inputs);
		sar.actions.push_back(action);

		/*=== End conditions ===*/
		if (dist < body.radius) {
			//Crash
			sar.rewards.push_back(-1e6);
			break;
		}

		/*=== Calculate thrust ===*/
		double ah = action[0] / std::sqrt(action[0] * action[0] + action[1] * action[1]);
		double av = action[1] / std::sqrt(action[0] * action[0] + action[1] * action[1]);

		double thrust = total_thrust / (1 + std::exp(-action[2]));
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

		craft.mass -= (fuel_reverse ? -1.0 : 1.0) * tick_consumption * timestep;

		/*=== Calculate periapsis ===*/
		sma = 1 / (2 / dist - (vx * vx + vy * vy) / body.grav_parameter);
		double eccentricity = std::sqrt(1 - (dist * dist * vh * vh) / (sma * body.grav_parameter));

		if (eccentricity > 1.0) {
			//Escape trajectory
			sar.rewards.push_back(-1e6);
			break;
		}

		periapsis = sma * (1 - eccentricity);

		/*=== Calculate rewards ===*/
		sar.rewards.push_back(-tick_consumption);
	}

	/*=== Use delta-v to raise/lower apoapsis ===*/
	double periapsis_vel =            std::sqrt(body.grav_parameter * (2 / periapsis - 1 / sma));
	double periapsis_vel_after_burn = std::sqrt(body.grav_parameter * (2 / periapsis - 2 / (periapsis + body.radius + orbit_alt)));

	craft.use_dv(periapsis_vel_after_burn - periapsis_vel);

	Result out;
	out.valid = (periapsis > body.radius);
	out.state = craft;
	out.score = craft.payload_fraction();

	out.critic_data.push_back(sar);

	return out;
}

Result simulate_landing(Body body, Craft craft, double landing_alt, double orbit_alt, NN agent) {
	/*=== Binary search to simulate landing ===*/
	double low = 0;
	double high = craft.mass;

	Craft state = craft;
	Result result;

	while (high - low > 1e-6) {
		double middle = (low + high) / 2;

		if (middle / craft.mass < 0.1) {
			return result;
		}

		state.mass = middle;

		result = simulate_ascent(body, state, landing_alt, orbit_alt, true, agent);

		if (!result.valid) {
			high = middle;
			continue;
		}

		if (result.state.mass > craft.mass) {
			high = middle;
			continue;
		}

		low = middle;
		continue;
	}

	/*=== Ascent ===*/
	state = result.state;
	Result ascent_result = simulate_ascent(body, state, landing_alt, orbit_alt, false, agent);

	ascent_result.critic_data.push_back(result.critic_data[0]);

	return ascent_result;
}

NN train_network(std::vector<int> layer_node_counts_agent) {
	/*=== Initialize NNs ===*/
	NN agent = get_random_nn(layer_node_counts_agent, Activation::ACTIVATION_RELU);

	layer_node_counts_agent[0] += 2;
	NN critic = get_random_nn(layer_node_counts_agent, Activation::ACTIVATION_RELU);

	/*=== Training ===*/
	for (int i = 0; i < 1000; i++) {
		/*=== Random body / craft ===*/
		Body body;

		body.radius = 850000 * std::pow(uniform(generator), 2.0);
		body.grav_parameter = std::abs(
			body.radius * 
			body.radius * 
			(1 + 0.2 * normal(generator)) * 
			((body.radius < 450000 ? 0.000009 : 0.000044) * (body.radius - 450000) + 4)
		);

		body.rotational_period = 150000 * std::pow(
			std::abs(std::log(uniform(generator))), 
			1.5
		);

		//Filters
		double surface_velocity = 2 * pi * body.radius / body.rotational_period;
		if (surface_velocity * surface_velocity * body.radius / body.grav_parameter > 0.25) {
			continue;
		}

		//Craft
		Craft craft;
		craft.mass = 100000;

		for (auto const & iterator : engines) {
			if (uniform(generator) > 0.1) continue;


		}

		std::cout << 
			body.radius << " " << 
			body.grav_parameter << " " << 
			body.rotational_period << " " << 
			body.grav_parameter / body.radius / body.radius << " " <<
			i << "\n";
	}

	return agent;
}

int main() {
	int num_engines = (int)engines.size();
	std::vector<int> layer_node_counts_agent = {num_engines + 10, 32, 32, 3};

	train_network(layer_node_counts_agent);
}
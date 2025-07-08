#include <iostream>
#include <utility>
#include <cmath>
#include <map>

#include "consts.h"
#include "craft.h"
#include "engines.h"
#include "random.h"
#include "body.h"
#include "nn.h"

const double timestep = 0.1;

class StateValueCollector {
	private:

	public:
		std::vector<std::pair<std::vector<double>, double> > stateValuePairs;
};

struct Result {
	bool valid;
	Craft state;
	double score;
};

Result simulate(Body body, Craft craft, double alt, double randomness, NN estimator, StateValueCollector& collector) {
	/*=== Constants ===*/
	double total_thrust = 0;
	for (auto const & iterator : craft.engines) {
		total_thrust += iterator.first.thrust * iterator.second;
	}

	/*=== State ===*/
	double x = body.radius + alt;
	double y = 0;
	double vx = 0;
	double vy = 2 * pi * x / body.rotational_period;
	Craft state = craft;
	double time = 0;

	/*=== Initialize state vector ===*/
	std::vector<double> stateVector(engines.size() + 6);

	/*=== Collector ===*/
	std::vector<std::pair<std::vector<double>, double> > stateValuePairs;

	/*=== End condition ===*/
	double periapsis = 0;
	bool valid = true;

	/*=== Iterate ===*/
	while (periapsis < body.radius) {
		/*=== Calculate state vector ===*/
		double dist = std::sqrt(x * x + y * y);
		double vh = (-y * vx + x * vy) / dist;
		double vv = ( x * vx + y * vy) / dist;

		stateVector[0] = dist;
		stateVector[1] = vh;
		stateVector[2] = vv;
		stateVector[3] = state.mass / craft.mass;

		int i = 6;
		for (auto const & iterator : engines) {
			if (state.engines.count(iterator.second) > 0) {
				stateVector[i] = state.engines.at(iterator.second) * iterator.second.mass / state.mass;
			} else {
				stateVector[i] = 0;
			}

			i++;
		}

		/*=== Calculate action ===*/
		double value = estimator.evaluate(stateVector)[0];
		for (int i = 0; i < 100; i++) {
			double prevAh = stateVector[4];
			double prevAv = stateVector[5];

			stateVector[4] += 0.1 * normal(generator);
			stateVector[5] += 0.1 * normal(generator);

			double newValue = estimator.evaluate(stateVector)[0];

			if (newValue <= value) {
				stateVector[4] = prevAh;
				stateVector[5] = prevAv;
			}
		}

		/*=== Calculate acceleration ===*/
		double ah = stateVector[4] + randomness * normal(generator); //Add variability for training purposes
		double av = stateVector[5] + randomness * normal(generator);

		double accel = std::sqrt(ah * ah + av * av);
		if (accel > total_thrust / state.mass) {
			ah *= total_thrust / (accel * state.mass);
			av *= total_thrust / (accel * state.mass);
			accel = total_thrust / state.mass;
		}

		double thrust = state.mass * accel;
		double tick_consumption = 0;
		for (auto const & iterator : state.engines) {
			double group_thrust = iterator.first.thrust * iterator.second;

			if (group_thrust > thrust) {
				tick_consumption += thrust / iterator.first.isp / g0;
				break;
			}

			tick_consumption += group_thrust / iterator.first.isp / g0;
			thrust -= group_thrust;
		}

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
		state.mass -= tick_consumption * timestep;
		time += timestep;

		/*=== Calculate orbit ===*/
		double sma = 1 / (2 / dist - (vx * vx + vy * vy) / body.grav_parameter);
		double eccentricity = std::sqrt(1 - (vh * vh * dist * dist) / (sma * body.grav_parameter));

		periapsis = sma * (1 - eccentricity);

		/*=== Collect data ===*/
		stateValuePairs.push_back(std::pair<std::vector<double>, double>(stateVector, -state.mass / craft.mass));

		/*=== Failure conditions ===*/
		if ((dist < body.radius) || (eccentricity > 1.0) || (time > 10000)) {
			valid = false;
			break;
		}
	}

	std::size_t numStateValuePairs = stateValuePairs.size();
	for (std::size_t i = 0; i < numStateValuePairs; i++) {
		collector.stateValuePairs.push_back(
			std::pair<std::vector<double>, double>(
				stateVector, 
				stateValuePairs[i].second + state.mass / craft.mass + (periapsis > body.radius ? 0 : 1e3 * (periapsis - body.radius) / body.radius)
			)
		);
	}

	Result out;
	out.state = state;
	out.valid = valid;
	out.score = state.mass / craft.mass;
	return out;
}

void train_nn(NN& estimator) {
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

			craft.engines[iterator.second] = std::pow(uniform(generator), 2) * craft.mass / iterator.second.mass;
		}

		/*=== Simulate ===*/
		StateValueCollector collector;
		simulate(body, craft, 10000 * uniform(generator), 0.1, estimator, collector);

		/*std::size_t numStateValuePairs = collector.stateValuePairs.size();
		for (std::size_t j = 0; j < numStateValuePairs; j++) {
			std::size_t stateSize = collector.stateValuePairs[j].first.size();
			for (std::size_t k = 0; k < stateSize; k++) {
				std::cout << collector.stateValuePairs[j].first[k] << " ";
			}

			std::cout << collector.stateValuePairs[j].second << "\n";
		}*/

		std::size_t numStateValuePairs = collector.stateValuePairs.size();

		/*=== Get loss ===*/
		double totalLoss = 0;
		for (std::size_t j = 0; j < numStateValuePairs; j++) {
			double actual = estimator.evaluate(collector.stateValuePairs[j].first)[0];

			std::cout << actual << " " << collector.stateValuePairs[j].second << "\n";

			totalLoss += 0.5 * std::pow(actual - collector.stateValuePairs[j].second, 2.0);
		}

		double averageLoss = totalLoss / numStateValuePairs;

		std::cout << averageLoss << "\n";

		/*=== Train network ===*/
		for (std::size_t j = 0; j < numStateValuePairs; j++) {
			estimator.backpropagate(collector.stateValuePairs[j].first, std::vector<double>{collector.stateValuePairs[j].second});
		}
	}
}

int main() {
	/*=== Settings ===*/
	double craftMass = 1000;

	std::map<Engine, double, std::greater<Engine> > craftEngines;
	craftEngines[engines.at("Spark")] = 1;

	/*=== Initialize NN ===*/
	std::vector<int> layerNeuronCounts = {6 + (int)engines.size(), 16, 16, 16, 1};
	NN valueEstimator = get_random_nn(layerNeuronCounts, Activation::ACTIVATION_RELU);

	/*=== Initialize craft ===*/
	Craft craft;

	craft.mass = craftMass;
	craft.engines = craftEngines;

	/*=== Train NN ===*/
	train_nn(valueEstimator);
}
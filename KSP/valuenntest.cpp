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
#include "debug.h"

const double timestep = 0.1;

struct SARGroup {
	std::vector<double> state;
	std::vector<double> actionRandomness;
	double reward;

	SARGroup(std::vector<double> newState, std::vector<double> newActionRandomness, double newReward) :
		state(newState), 
		actionRandomness(newActionRandomness), 
		reward(newReward) {}
};

struct StateValueCollector {
	std::vector<SARGroup> stateValuePairs;
};

struct Result {
	bool valid;
	Craft state;
	double score;
};

Result simulate(Body body, Craft craft, double alt, double randomness, NN agent, StateValueCollector& collector) {
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
	std::vector<double> stateVector(engines.size() + 8);

	/*=== Collector ===*/
	std::vector<SARGroup> stateValuePairs;

	/*=== End condition ===*/
	double periapsis = 0;
	double vh, vv, dist;
	bool valid = true;

	/*=== Iterate ===*/
	while (periapsis < body.radius) {
		/*=== Calculate state vector ===*/
		dist = std::sqrt(x * x + y * y);
		vh = (-y * vx + x * vy) / dist;
		vv = ( x * vx + y * vy) / dist;

		stateVector[0] = dist * 1e-5;
		stateVector[1] = (dist - body.radius) * 1e-3;
		stateVector[2] = body.radius * 1e-5;
		stateVector[3] = std::sqrt(body.grav_parameter / body.radius) * 1e-3;
		stateVector[4] = body.grav_parameter / body.radius / body.radius;
		stateVector[5] = vh * 1e-3;
		stateVector[6] = vv * 1e-3;
		stateVector[7] = state.mass / craft.mass;

		int i = 8;
		for (auto const & iterator : engines) {
			if (state.engines.count(iterator.second) > 0) {
				stateVector[i] = state.engines.at(iterator.second) * iterator.second.mass / state.mass;
			} else {
				stateVector[i] = 0;
			}

			i++;
		}

		/*=== Calculate acceleration ===*/
		std::vector<double> action = agent.evaluate(stateVector);
		double rh = randomness * normal(generator), rv = randomness * normal(generator);

		double ah = action[0] + rh; //Add variability for training purposes
		double av = action[1] + rv;

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
		stateValuePairs.push_back(SARGroup{stateVector, std::vector<double>{rh, rv}, 0});

		/*=== Failure conditions ===*/
		if ((dist < body.radius) || (eccentricity > 1.0) || (time > 10000)) {
			valid = false;
			break;
		}
	}

	std::size_t numStateValuePairs = stateValuePairs.size();
	for (std::size_t i = 0; i < numStateValuePairs; i++) {
		double score = state.mass / craft.mass + (periapsis > body.radius ? 1e6 : vh + time);

		if (std::isnan(score)) {
			continue;
		}

		stateValuePairs[i].reward = score;
		collector.stateValuePairs.push_back(stateValuePairs[i]);
	}

	Result out;
	out.state = state;
	out.valid = valid;
	out.score = state.mass / craft.mass;
	return out;
}

void train_nn(NN& estimator, NN& agent) {
	estimator.log();

	double maxReward = -1000;
	for (int i = 0; i < 100000; i++) {
		/*=== body / craft ===*/
		Body body;

		body.radius = 3e5;
		body.grav_parameter = 3e5 * 3e5 * 0.275 * g0;
		body.rotational_period = 105962.088893924;

		//Filters
		double surface_velocity = 2 * pi * body.radius / body.rotational_period;
		if (surface_velocity * surface_velocity * body.radius / body.grav_parameter > 0.25) {
			continue;
		}

		//Craft
		Craft craft;
		craft.mass = 100000;

		std::map<Engine, double, std::greater<Engine> > craftEngines;
		craftEngines[engines.at("Spark")] = 1;

		craft.engines = craftEngines;

		/*for (auto const & iterator : engines) {
			if (uniform(generator) > 0.1) continue;

			craft.engines[iterator.second] = std::pow(uniform(generator), 2) * craft.mass / iterator.second.mass;
		}*/

		/*=== Simulate ===*/
		StateValueCollector collector;
		simulate(body, craft, 5000, 1, agent, collector);

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
		double reward = 0;
		for (std::size_t j = 0; j < numStateValuePairs; j++) {
			double actual = estimator.evaluate(collector.stateValuePairs[j].state)[0];

			//if (j == 0) std::cout << actual << " " << collector.stateValuePairs[j].reward << "\n";

			reward = collector.stateValuePairs[j].reward;

			if (reward > maxReward) {
				maxReward = reward;
			}
		}

		//std::cout << reward << ", " << maxReward << " 1\n";

		/*=== Train network ===*/
		for (std::size_t j = 0; j < numStateValuePairs; j++) {
			SARGroup sartriad = collector.stateValuePairs[j];

			estimator.backpropagate(sartriad.state, std::vector<double>{sartriad.reward}, 1e-1, 0.001);

			double estimatedReward = estimator.evaluate(sartriad.state)[0];
			double value = sartriad.reward - estimatedReward;

			value *= 20;

			double gx = value * sartriad.actionRandomness[0];
			double gy = value * sartriad.actionRandomness[1];

			std::vector<double> prevAction = agent.evaluate(sartriad.state);

			agent.backpropagate_gradient(sartriad.state, std::vector<double>{gx, gy}, 1e-6, 0.00001);

			std::vector<double> newAction = agent.evaluate(sartriad.state);

			if (i % 100 == 0) {
				std::cout << sartriad.actionRandomness[0] << " " << sartriad.actionRandomness[1] << " " << gx << " " << gy << " " << sartriad.reward << " " << estimatedReward << "\n";
				std::cout << i << " Prev: ";
				logDoubleVector(prevAction);
				std::cout << "\nNew: ";
				logDoubleVector(newAction);
				std::cout << "\n\n";
			}
		}
	}
}

int main() {
	/*=== Settings ===*/
	double craftMass = 1000;

	std::map<Engine, double, std::greater<Engine> > craftEngines;
	craftEngines[engines.at("Spark")] = 1;

	/*=== Initialize NN ===*/
	std::vector<int> agentLayerNeuronCounts = {8 + (int)engines.size(), 16, 16, 2};
	std::vector<int> estimatorLayerNeuronCounts = {8 + (int)engines.size(), 16, 16, 1};
	NN agent = get_random_nn(agentLayerNeuronCounts, Activation::ACTIVATION_LOGISTIC);
	NN valueEstimator = get_random_nn(estimatorLayerNeuronCounts, Activation::ACTIVATION_LOGISTIC);

	/*=== Initialize craft ===*/
	Craft craft;

	craft.mass = craftMass;
	craft.engines = craftEngines;

	/*=== Train NN ===*/
	train_nn(valueEstimator, agent);
}
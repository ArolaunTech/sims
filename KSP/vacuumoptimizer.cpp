#include <iostream>
#include <cmath>
#include <vector>
#include <random>

//Physical constants
const double pi = 3.14159;
const double g = 9.81;

const double fuelTankMult[4] = {1.1, 1.125, 4.0/3.0, 17.0/15.0};

//Planetary constants
const double radius = 6.0e+5;
const double gravity = 0.8;
const double rotationalPeriod = 211926.36;

//Flight constants
const double orbitAltitude = 10000.0;
const double landingAltitude = 5000.0;

//Constant calculations
const double GM = gravity * g * radius * radius;

const double landingRadius = radius + landingAltitude;
const double orbitalRadius = radius + orbitAltitude;

const double rotationalVelocity = 2.0 * pi * landingRadius/rotationalPeriod;
const double landingOrbitalVelocity = std::sqrt(GM/landingRadius);
const double landingTarget = std::sqrt(GM*(2.0/landingRadius - 2.0/(landingRadius + orbitalRadius)));
const double landingGravity = GM/(landingRadius*landingRadius);
const double orbitOrbitalVelocity = std::sqrt(GM/orbitalRadius);
const double circularization = orbitOrbitalVelocity - std::sqrt(GM*(2.0/orbitalRadius - 2.0/(landingRadius + orbitalRadius)));

//Structs
struct engine {
	double isp;
	int fueltype; //0 - Liquid fuel, 1 - Liquid fuel+Oxidiser, 2 - Xenon, 3 - Monopropellant
	double thrust;
	double fuelConsumption;
	double mass;
	engine(double isp, int fueltype, double thrust, double mass) : isp(isp), fueltype(fueltype), thrust(thrust), mass(mass) {
		fuelConsumption = thrust/(isp*g);
	}
};

struct simulationResult {
	double finalMass;
	double fuelUsed[4];
};

//Randomness settings
std::mt19937_64 generator(time(NULL));
std::uniform_real_distribution<double> uniform(0.0, 1.0);
std::normal_distribution<double> normal(0.0, 1.0);

//Simulation
double calculateEffectiveIsp(double mass, std::vector<engine> engineList, std::vector<int> engineNum, std::vector<double> throttles, double velocity) {
	double orbFrac = velocity/landingOrbitalVelocity;
	double effectiveGrav = landingGravity * (1.0 - orbFrac * orbFrac);
	double weight = mass * effectiveGrav;
	double thrust = 0.0;
	double fuelConsumption = 0.0;
	int numEngines = engineNum.size();

	for (int i = 0; i < numEngines; i++) {
		thrust += engineList[i].thrust * (double)engineNum[i] * throttles[i];
		fuelConsumption += engineList[i].fuelConsumption * fuelTankMult[engineList[i].fueltype] * (double)engineNum[i] * throttles[i];
	}
	if (thrust < weight) {
		return 0.0;
	}
	return std::sqrt(thrust * thrust - weight * weight)/(fuelConsumption * g);
}

simulationResult simulateLandingAscent(double initialMass, std::vector<engine> engineList, std::vector<int> engineNum) {
	double fuelUsed[4] = {0.0};

	int bestEngineFuelType;
	double bestEngineIsp = 0.0;

	const int numEngines = engineNum.size();
	for (int i = 0; i < numEngines; i++) {
		if ((engineList[i].isp > bestEngineIsp) && (engineNum[i] > 0)) {
			bestEngineIsp = engineList[i].isp;
			bestEngineFuelType = engineList[i].fueltype;
			//std::cout << bestEngineIsp << " " << bestEngineFuelType << "\n";
		}
	}

	if (bestEngineIsp == 0.0) {
		simulationResult res;
		res.finalMass = -1.0;
		for (int i = 0; i < 4; i++) {
			res.fuelUsed[i] = 0.0;
		}
		return res;
	}

	double currMass = initialMass/std::exp(circularization/(bestEngineIsp * g));
	fuelUsed[bestEngineFuelType] += initialMass - currMass;
	double currVel = landingTarget;

	double orbitalFraction;
	double verticalAcceleration;
	double horizontalAcceleration;

	double time = 0.0;
	const double timestep = 0.5;

	double thrust;
	double acceleration;
	double engineGroupFuel;

	double numEnginesFiring;

	bool valid = true;

	std::vector<double> throttles;
	for (int i = 0; i < numEngines; i++) {
		throttles.push_back(1.0);
	}
	double effectiveIsp;
	double oldThrottle;
	double neweffectiveIsp;

	double thrustWE;
	double fuelWE;
	double ispWithoutEngine;
	double weight;
	double wantedThrust;

	while (currVel > rotationalVelocity) {
		//Optimize throttles
		orbitalFraction = currVel/landingOrbitalVelocity;
		for (int i = 0; i < 2; i++) {
			effectiveIsp = calculateEffectiveIsp(currMass, engineList, engineNum, throttles, currVel);
			for (int j = 0; j < numEngines; j++) {
				oldThrottle = throttles[j];

				thrustWE = 0.0;
				fuelWE = 0.0;
				for (int k = 0; k < numEngines; k++) {
					if (j == k) {
						continue;
					}
					thrustWE += engineList[k].thrust * (double)engineNum[k] * throttles[k];
					fuelWE += engineList[k].thrust * (double)engineNum[k] * throttles[k] * fuelTankMult[engineList[k].fueltype]/(engineList[k].isp * g);
				}

				if (thrustWE == 0.0) {
					continue;
				}

				weight = currMass * landingGravity * (1.0 - orbitalFraction * orbitalFraction);

				ispWithoutEngine = thrustWE/(fuelWE * g);
				double newIsp = engineList[j].isp;

				wantedThrust = (thrustWE * thrustWE)/newIsp - (weight*weight)/newIsp - (thrustWE*thrustWE)/ispWithoutEngine;
				wantedThrust /= thrustWE/ispWithoutEngine - thrustWE/newIsp;
				wantedThrust /= engineList[j].thrust * (double)engineNum[j];

				if (wantedThrust < 0.0) {
					wantedThrust = 0.0;
				}
				if (wantedThrust > 1.0) {
					wantedThrust = 1.0;
				}

				throttles[j] = wantedThrust;
				//if (throttles[j] < 0.0) {
				//	throttles[j] = 0.0;
				//}
				//if (throttles[j] > 1.0) {
				//	throttles[j] = 1.0;
				//}

				neweffectiveIsp = calculateEffectiveIsp(currMass, engineList, engineNum, throttles, currVel);
				if (neweffectiveIsp < effectiveIsp) {
					throttles[j] = oldThrottle;
				} else {
					effectiveIsp = neweffectiveIsp;
				}
			}
		}
		//for (int i = 0; i < numEngines; i++) {
		//	std::cout << throttles[i] << " " << engineNum[i] << "\n";
		//}
		//std::cout << effectiveIsp << "\n";

		//Calculate thrust and fuel consumption
		thrust = 0.0;
		for (int i = 0; i < numEngines; i++) {
			numEnginesFiring = throttles[i] * (double)engineNum[i];

			thrust += engineList[i].thrust * numEnginesFiring; //Thrust of each engine times number of engines
			engineGroupFuel = engineList[i].fuelConsumption * numEnginesFiring;
			
			fuelUsed[engineList[i].fueltype] += engineGroupFuel * timestep;
			currMass -= engineGroupFuel * timestep;
		}

		//std::cout << thrust << "\n";

		//Calculate acceleration
		acceleration = thrust/currMass;
		verticalAcceleration = landingGravity * (orbitalFraction*orbitalFraction - 1.0);

		if (acceleration < std::fabs(verticalAcceleration)) { //Not enough thrust to stay at a constant altitude
			valid = false;
			break;
		}

		horizontalAcceleration = std::sqrt(acceleration * acceleration - verticalAcceleration * verticalAcceleration);

		//Integrate
		time += timestep;
		currVel -= horizontalAcceleration * timestep;

		//Printout
		//std::cout << verticalAcceleration << "\t" << horizontalAcceleration << "\t" << currVel << "\n";
	}

	currVel = rotationalVelocity;
	while (valid && (currVel < landingTarget)) {
		//Optimize throttles
		orbitalFraction = currVel/landingOrbitalVelocity;
		for (int i = 0; i < 2; i++) {
			effectiveIsp = calculateEffectiveIsp(currMass, engineList, engineNum, throttles, currVel);
			for (int j = 0; j < numEngines; j++) {
				oldThrottle = throttles[j];

				thrustWE = 0.0;
				fuelWE = 0.0;
				for (int k = 0; k < numEngines; k++) {
					if (j == k) {
						continue;
					}
					thrustWE += engineList[k].thrust * (double)engineNum[k] * throttles[k];
					fuelWE += engineList[k].thrust * (double)engineNum[k] * throttles[k] * fuelTankMult[engineList[k].fueltype]/(engineList[k].isp * g);
				}

				if (thrustWE == 0.0) {
					continue;
				}

				weight = currMass * landingGravity * (1.0 - orbitalFraction * orbitalFraction);

				ispWithoutEngine = thrustWE/(fuelWE * g);
				double newIsp = engineList[j].isp;

				wantedThrust = (thrustWE * thrustWE)/newIsp - (weight*weight)/newIsp - (thrustWE*thrustWE)/ispWithoutEngine;
				wantedThrust /= thrustWE/ispWithoutEngine - thrustWE/newIsp;
				wantedThrust /= engineList[j].thrust * (double)engineNum[j];

				if (wantedThrust < 0.0) {
					wantedThrust = 0.0;
				}
				if (wantedThrust > 1.0) {
					wantedThrust = 1.0;
				}

				throttles[j] = wantedThrust;
				//if (throttles[j] < 0.0) {
				//	throttles[j] = 0.0;
				//}
				//if (throttles[j] > 1.0) {
				//	throttles[j] = 1.0;
				//}

				neweffectiveIsp = calculateEffectiveIsp(currMass, engineList, engineNum, throttles, currVel);
				if (neweffectiveIsp < effectiveIsp) {
					throttles[j] = oldThrottle;
				} else {
					effectiveIsp = neweffectiveIsp;
				}
			}
		}

		//Calculate thrust and fuel consumption
		thrust = 0.0;
		for (int i = 0; i < numEngines; i++) {
			numEnginesFiring = throttles[i] * (double)engineNum[i];

			thrust += engineList[i].thrust * numEnginesFiring; //Thrust of each engine times number of engines
			engineGroupFuel = engineList[i].fuelConsumption * numEnginesFiring;
			
			fuelUsed[engineList[i].fueltype] += engineGroupFuel * timestep;
			currMass -= engineGroupFuel * timestep;
		}

		//Calculate acceleration
		acceleration = thrust/currMass;
		verticalAcceleration = landingGravity * (orbitalFraction*orbitalFraction - 1.0);

		if (acceleration < std::fabs(verticalAcceleration)) { //Not enough thrust to stay at a constant altitude
			valid = false;
			break;
		}

		horizontalAcceleration = std::sqrt(acceleration * acceleration - verticalAcceleration * verticalAcceleration);

		//Integrate
		time += timestep;
		currVel += horizontalAcceleration * timestep;

		//Printout
		//std::cout << verticalAcceleration << "\t" << horizontalAcceleration << "\t" << currVel << "\n";
	}

	double totalEngineMass = 0.0;
	for (int i = 0; i < numEngines; i++) {
		totalEngineMass += engineList[i].mass * (double)engineNum[i];
	}

	simulationResult res;
	res.finalMass = currMass - totalEngineMass;
	res.finalMass -= fuelUsed[0]/10.0; //LF tanks
	res.finalMass -= fuelUsed[1]/8.0; //LFOX tanks
	res.finalMass -= fuelUsed[2]/3.0; //Xenon tanks
	res.finalMass -= fuelUsed[3]/7.5; //Monopropellant tanks

	if (!valid) {
		res.finalMass = -1.0;
	} else {
		for (int i = 0; i < 4; i++) {
			res.fuelUsed[i] = fuelUsed[i];
		}
	}
	return res;
}

simulationResult optimizeLanding(double initialMass, std::vector<engine> engineList, std::vector<int> minEngineCounts) {
	simulationResult bestResult;
	bestResult.finalMass = -1.0;

	simulationResult bestRandomResult;
	simulationResult randomResult;
	simulationResult oldResult;
	simulationResult newResult;

	std::vector<int> randomEngineNum;
	std::vector<int> bestRandomEngineNum;

	std::vector<int> oldEngineNum;
	std::vector<int> newEngineNum;

	std::vector<double> key;
	int nKeys;
	int newAppendingNumEngines;

	int numEngines = engineList.size();

	for (int run = 0; run < 10; run++) {
		//Generate random runs
		bestRandomResult.finalMass = -1.0;
		for (int randRun = 0; randRun < 20000; randRun++) {
			randomEngineNum.clear();
			for (int i = 0; i < numEngines; i++) {
				newAppendingNumEngines = (int)(uniform(generator) * 50.0);
				if (uniform(generator) < 0.8) {
					newAppendingNumEngines = minEngineCounts[i];
				}
				if (newAppendingNumEngines < minEngineCounts[i]) {
					newAppendingNumEngines = minEngineCounts[i];
				}
				randomEngineNum.push_back(newAppendingNumEngines);
			}

			randomResult = simulateLandingAscent(initialMass, engineList, randomEngineNum);
			//std::cout << randomResult.finalMass << "\n";
			if (randomResult.finalMass < 0) {
				continue;
			}
			if (randomResult.finalMass < bestRandomResult.finalMass) {
				continue;
			}

			if (randomResult.finalMass > bestResult.finalMass) {
				std::cout << "Found good random result:\n" <<
				             "  - Payload: " << randomResult.finalMass << "kg (" << 100.0*randomResult.finalMass/initialMass << "% payload fraction)\n" <<
				             "  - Engine Numbers: {\n" <<
				             "    ";
				for (int i = 0; i < numEngines; i++) {
					std::cout << randomEngineNum[i] << ", ";
				}
				std::cout << "\n" <<
				             "  }\n";
				std::cout << "  - Fuel types used:\n" <<
				             "     - LF: " << randomResult.fuelUsed[0] << "kg\n" <<
				             "     - LFOX: " << randomResult.fuelUsed[1] << "kg\n" <<
				             "     - XE: " << randomResult.fuelUsed[2] << "kg\n" <<
				             "     - MP: " << randomResult.fuelUsed[3] << "kg\n\n";
			}

			bestRandomResult = randomResult;
			bestRandomEngineNum = randomEngineNum;
			//std::cout << randomResult.finalMass/initialMass << "\n";
		}

		if (bestRandomResult.finalMass < 0) {
			continue; //No random run could get a positive payload fraction. ABORT!
		}
		if (bestRandomResult.finalMass > bestResult.finalMass) {
			bestResult = bestRandomResult;
		}

		oldEngineNum = bestRandomEngineNum;
		oldResult = bestRandomResult;
		//Optimize random runs to make them better
		for (int step = 0; step < 20000; step++) {
			newEngineNum = oldEngineNum;

			for (int i = 0; i < numEngines; i++) {
				newEngineNum[i] += (int)(5.0*normal(generator));
				if (newEngineNum[i] < minEngineCounts[i]) {
					newEngineNum[i] = minEngineCounts[i];
				}
			}

			newResult = simulateLandingAscent(initialMass, engineList, newEngineNum);
			if ((newResult.finalMass > oldResult.finalMass) || (uniform(generator) < 0.25 - (double)step/10000.0)) {
				oldResult = newResult;
				oldEngineNum = newEngineNum;

				//std::cout << newResult.finalMass;

				if (newResult.finalMass > bestResult.finalMass) {
					std::cout << "Found good result:\n" <<
			    	         "  - Payload: " << newResult.finalMass << "kg (" << 100.0*newResult.finalMass/initialMass << "% payload fraction)\n" <<
			    	         "  - Engine Numbers: {\n" <<
			    	         "    ";
					for (int i = 0; i < numEngines; i++) {
						std::cout << newEngineNum[i] << ", ";
					}
					std::cout << "\n" <<
					             "  }\n" <<
				                 "  - Fuel types used:\n" <<
				                 "     - LF: " << newResult.fuelUsed[0] << "kg\n" <<
				                 "     - LFOX: " << newResult.fuelUsed[1] << "kg\n" <<
				                 "     - XE: " << newResult.fuelUsed[2] << "kg\n" <<
				                 "     - MP: " << newResult.fuelUsed[3] << "kg\n\n";
				}
				//std::cout << newResult.finalMass << "\n";
			}
		}

		if (oldResult.finalMass > bestResult.finalMass) {
			bestResult = oldResult;
		}
	}
	return bestResult;
}

int main() {
	const double initialMass = 1000000.0;
	//engine nerva = {800.0, 0, 60000.0, 3000.0};
	//std::cout << engine{800.0, 0, 60000.0, 3000.0}.fuelConsumption << "\n";

	//simulationResult landingResult = simulateLandingAscent(initialMass, std::vector<engine> {engine{800.0, 0, 60000.0, 3000.0}}, std::vector<int> {50}, std::vector<std::vector<double> > {{0.0, 1.0}});
	//std::cout << landingResult.finalMass << "\n";

	simulationResult optimizerResult = optimizeLanding(
		initialMass, 
		std::vector<engine> { //Put the allowed engines here, currently NERV, RAPIER, Wolfhound, Rhino, Dart, Dawn
			engine{800.0, 0, 60000.0, 3000.0}, 
			engine{305.0, 1, 180000.0, 2000.0}, 
			engine{380.0, 1, 375000.0, 3300.0},
			engine{340.0, 1, 180.0e+3, 1.0e+3},
			engine{4200.0, 2, 2000.0, 250.0}
		},
		std::vector<int> { //Minimum engine counts
			0,
			0,
			0,
			0,
			0
		}
	);

	return 0;
}
#include <iostream>
#include <cmath>
#include <vector>
#include <random>

//Physical constants
const double pi = 3.14159;
const double g = 9.81;

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
simulationResult simulateLandingAscent(double initialMass, std::vector<engine> engineList, std::vector<int> engineNum, std::vector<std::vector<double> > code) {
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

	double nextSegmentTime = code[0][0];
	int nextCodeSegment = 0;
	std::vector<double> currCodeSegment(code[0].size(), 0.0);

	double numEnginesFiring;

	bool valid = true;

	while (currVel > rotationalVelocity) {
		//Read code
		if (time >= nextSegmentTime) {
			if (nextCodeSegment < code.size()) {
				currCodeSegment = code[nextCodeSegment];
				nextCodeSegment += 1;
				nextSegmentTime += currCodeSegment[0];
			}
		}

		//Calculate thrust and fuel consumption
		thrust = 0.0;
		for (int i = 0; i < numEngines; i++) {
			numEnginesFiring = currCodeSegment[i + 1] * (double)engineNum[i];

			thrust += engineList[i].thrust * numEnginesFiring; //Thrust of each engine times number of engines
			engineGroupFuel = engineList[i].fuelConsumption * numEnginesFiring;
			
			fuelUsed[engineList[i].fueltype] += engineGroupFuel * timestep;
			currMass -= engineGroupFuel * timestep;
		}

		//Calculate acceleration
		acceleration = thrust/currMass;
		orbitalFraction = currVel/landingOrbitalVelocity;
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
		//Read code
		if (time >= nextSegmentTime) {
			if (nextCodeSegment < code.size()) {
				currCodeSegment = code[nextCodeSegment];
				nextCodeSegment += 1;
				nextSegmentTime += currCodeSegment[0];
			}
		}

		//Calculate thrust and fuel consumption
		thrust = 0.0;
		for (int i = 0; i < numEngines; i++) {
			numEnginesFiring = currCodeSegment[i + 1] * (double)engineNum[i];

			thrust += engineList[i].thrust * numEnginesFiring; //Thrust of each engine times number of engines
			engineGroupFuel = engineList[i].fuelConsumption * numEnginesFiring;
			
			fuelUsed[engineList[i].fueltype] += engineGroupFuel * timestep;
			currMass -= engineGroupFuel * timestep;
		}

		//Calculate acceleration
		acceleration = thrust/currMass;
		orbitalFraction = currVel/landingOrbitalVelocity;
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
	std::vector<std::vector<double> > randomCode;
	std::vector<int> bestRandomEngineNum;
	std::vector<std::vector<double> > bestRandomCode;

	std::vector<int> oldEngineNum;
	std::vector<std::vector<double> > oldCode;
	std::vector<int> newEngineNum;
	std::vector<std::vector<double> > newCode;

	std::vector<double> key;
	int nKeys;
	int newAppendingNumEngines;

	int numEngines = engineList.size();

	for (int run = 0; run < 10; run++) {
		//Generate random runs
		bestRandomResult.finalMass = -1.0;
		for (int randRun = 0; randRun < 1000; randRun++) {
			randomEngineNum.clear();
			for (int i = 0; i < numEngines; i++) {
				newAppendingNumEngines = (int)(uniform(generator) * 50.0);
				if (uniform(generator) < 0.5) {
					newAppendingNumEngines = minEngineCounts[i];
				}
				if (newAppendingNumEngines < minEngineCounts[i]) {
					newAppendingNumEngines = minEngineCounts[i];
				}
				randomEngineNum.push_back(newAppendingNumEngines);
			}

			randomCode.clear();
			nKeys = 1 + (int)(uniform(generator)*3.0);
			for (int i = 0; i < nKeys; i++) {
				key.clear();
				key.push_back(uniform(generator) * 300.0);
				for (int j = 0; j < numEngines; j++) {
					key.push_back(1.0);
				}
				randomCode.push_back(key);
			}
			randomCode[0][0] = 0.0;

			randomResult = simulateLandingAscent(initialMass, engineList, randomEngineNum, randomCode);
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
				             "  }\n" <<
				             "  - Code: {\n";
				for (int i = 0; i < nKeys; i++) {
					std::cout << "    {T:" << randomCode[i][0] << "s, ";
					for (int j = 1; j <= numEngines; j++) {
						std::cout << randomCode[i][j] << ", ";
					}
					std::cout << "}\n";
				}
				std::cout << "  }\n" <<
				             "  - Fuel types used:\n" <<
				             "     - LF: " << randomResult.fuelUsed[0] << "kg\n" <<
				             "     - LFOX: " << randomResult.fuelUsed[1] << "kg\n" <<
				             "     - XE: " << randomResult.fuelUsed[2] << "kg\n" <<
				             "     - MP: " << randomResult.fuelUsed[3] << "kg\n\n";
			}

			bestRandomResult = randomResult;
			bestRandomCode = randomCode;
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
		oldCode = bestRandomCode;
		oldResult = bestRandomResult;
		//Optimize random runs to make them better
		for (int step = 0; step < 10000; step++) {
			newEngineNum = oldEngineNum;
			newCode = oldCode;

			for (int i = 0; i < numEngines; i++) {
				newEngineNum[i] += (int)(5.0*normal(generator));
				if (newEngineNum[i] < minEngineCounts[i]) {
					newEngineNum[i] = minEngineCounts[i];
				}
			}
			for (int i = 0; i < oldCode.size(); i++) {
				newCode[i][0] += 50.0*normal(generator);
				if (newCode[i][0] < 0.0) {
					newCode[i][0] = 0.0;
				}
				for (int j = 1; j <= numEngines; j++) {
					newCode[i][j] += 0.1*normal(generator);
					if (newCode[i][j] < 0.0) {
						newCode[i][j] = 0.0;
					}
					if (newCode[i][j] > 1.0) {
						newCode[i][j] = 1.0;
					}
				}
			}
			newCode[0][0] = 0.0;

			newResult = simulateLandingAscent(initialMass, engineList, newEngineNum, newCode);
			if ((newResult.finalMass > oldResult.finalMass) || (uniform(generator) < 0.25 - (double)step/10000.0)) {
				oldCode = newCode;
				oldResult = newResult;
				oldEngineNum = newEngineNum;

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
					             "  - Code: {\n";
					for (int i = 0; i < newCode.size(); i++) {
						std::cout << "    {T:" << newCode[i][0] << "s, ";
						for (int j = 1; j <= numEngines; j++) {
							std::cout << newCode[i][j] << ", ";
						}
						std::cout << "}\n";
					}
					std::cout << "  }\n" <<
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
	const double initialMass = 466650;
	//engine nerva = {800.0, 0, 60000.0, 3000.0};
	//std::cout << engine{800.0, 0, 60000.0, 3000.0}.fuelConsumption << "\n";

	//simulationResult landingResult = simulateLandingAscent(initialMass, std::vector<engine> {engine{800.0, 0, 60000.0, 3000.0}}, std::vector<int> {50}, std::vector<std::vector<double> > {{0.0, 1.0}});
	//std::cout << landingResult.finalMass << "\n";

	simulationResult optimizerResult = optimizeLanding(
		initialMass, 
		std::vector<engine> { //Put the allowed engines here, currently NERV, RAPIER, Wolfhound, Dart, Dawn
			engine{800.0, 0, 60000.0, 3000.0}, 
			engine{305.0, 1, 180000.0, 2000.0}, 
			engine{380.0, 1, 375000.0, 3300.0},
			engine{340.0, 1, 180000.0, 1000.0},
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
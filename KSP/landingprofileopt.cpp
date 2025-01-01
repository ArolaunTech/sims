#include <iostream>
#include <vector>
#include <random>
#include <cmath>

//Planet settings
const double planetRadius = 250000.0;      //m
const double planetG = 2.7;                //m s^-2
const double landingHeight = 6800.0;       //m
const double rotationalPeriod = 1210000.0; //s

//Flight settings
const double initMass = 5450.0;     //kg
const double initXe = 450.0;        //kg
const double initLfOx = 30.0;       //kg
const double ionThrust = 14000.0;   //N
const double chemThrust = 180000.0; //N
const double ionIsp = 4200.0;       //s
const double chemIsp = 305.0;       //s

//Simulation settings
const double timestep = 0.5; //s

//Random generator settings
std::mt19937 generator(time(NULL));
std::uniform_real_distribution<double> uniform(0.0, 1.0);
std::normal_distribution<double> normal(0.0, 1.0);
std::uniform_int_distribution<int> numCodeSegments(4, 10);

//Calculations
const double pi = 3.14159265;
const double g = 9.80665; //m s^-2

const double planetSTDGP = planetRadius * planetRadius * planetG;
const double orbitalSpeed = std::sqrt(
	planetSTDGP / (planetRadius + landingHeight)
);
const double rotationalSpeed = 2 * pi * (planetRadius + landingHeight) / 
	rotationalPeriod;
const double landingGrav = planetSTDGP / 
	std::pow(planetRadius + landingHeight, 2.0);

const double ionXe = ionThrust/(ionIsp * g);
const double chemLfOx = chemThrust/(chemIsp * g);

struct craftState {
	double mass;
	double xe;
	double lfox;
};

bool works(craftState state) {
	return (state.xe > 0) && (state.lfox > 0);
}

bool greater(craftState a, craftState b) { //True if a > b
	return a.xe > b.xe;
}

craftState simulateLanding(std::vector<std::vector<double> > code, bool print) {
	double mass = initMass;
	double xe = initXe;
	double lfox = initLfOx;

	double vel = orbitalSpeed;
	double dist = 0.0;
	double time = 0.0;

	double accel;
	double vaccel;
	double haccel;

	int codeIdx = 0;
	double changeTime = 0.0;

	double currCodeChemThrust = code[0][1];
	std::size_t numCodes = code.size();
	double dir = -1.0;

	while (vel <= orbitalSpeed) {
		//Update code
		if ((codeIdx < (numCodes - 1)) && 
			(time > changeTime + code[codeIdx+1][0])) {
			codeIdx++;
			changeTime += code[codeIdx][0];
			currCodeChemThrust = code[codeIdx][1];
		}
		if (lfox == 0.0) {
			currCodeChemThrust = 0.0;
		}

		//Calculate thrust
		accel = (ionThrust + currCodeChemThrust * chemThrust)/mass;
		vaccel = landingGrav * (1.0 - std::pow(vel/orbitalSpeed, 2.0));
		if (vaccel > accel) { //Cannot support craft
			xe = -1.0;
			lfox = -1.0;
			break;
		}
		haccel = std::sqrt(accel * accel - vaccel * vaccel);

		//Integrate
		vel += dir * haccel * timestep;
		dist += vel * timestep;
		time += timestep;
		xe -= ionXe * timestep;
		lfox -= chemLfOx * currCodeChemThrust * timestep;
		mass -= (ionXe + chemLfOx * currCodeChemThrust) * timestep;
		if (xe < 0.0) //Ran out of fuel
			break;
		if (lfox < 0.0)
			lfox = 0.0;

		if (vel <= rotationalSpeed) {
			dir = 1.0;
			vel = rotationalSpeed;
		}

		if (print)
			std::cout << 
				"Time: " << time << "s,\t"
				"Velocity: " << vel << "m/s,\t"
				"Xe: " << xe << "kg,\t"
				"LfOx: " << lfox << "kg,\t"
				"Throttle: " << currCodeChemThrust << "\n";
	}

	craftState out;
	out.mass = mass;
	out.lfox = lfox;
	out.xe = xe;
	return out;
}

craftState optimizeLanding() {
	std::vector<std::vector<double> > code;
	std::vector<std::vector<double> > bestCode;

	craftState state;
	state.xe = -1.0;
	craftState bestState;
	bestState.xe = -1.0;

	std::vector<double> key = {0.0, 0.0};

	int numCodes;

	for (int i = 0; i < 100000; i++) {
		if ((i == 0) || (uniform(generator) < 0.1)) {
			code.clear();

			numCodes = numCodeSegments(generator);
			for (int j = 0; j < numCodes; j++) {
				key[0] = uniform(generator) * 200;
				key[1] = uniform(generator) * 0.001;
				code.push_back(key);
			}
			code[0][0] = 0.0;
		} else {
			code = bestCode;
			numCodes = code.size();
			for (int j = 0; j < numCodes; j++) {
				code[j][0] += normal(generator) * 5.0;
				if (code[j][0] < 0.0) code[j][0] = 0.0;

				code[j][1] += normal(generator) * 0.0005;
				if (code[j][1] < 0.0) code[j][1] = 0.0;
				if (code[j][1] > 1.0) code[j][1] = 1.0;
			}
			code[0][0] = 0.0;
		}

		state = simulateLanding(code, false);

		if ((i == 0) || (works(state) && greater(state, bestState))) {
			bestState = state;
			bestCode = code;
		}
	}

	//Logging
	simulateLanding(bestCode, true);

	std::cout << "\nCode:\n";
	std::size_t bestSize = bestCode.size();
	double time = 0.0;

	for (std::size_t i = 0; i < bestSize; i++) {
		time += bestCode[i][0];
		std::cout << 
			"Time: " << time << "s,\t"
			"Throttle: " << bestCode[i][1] << "\n";
	}

	return bestState;
}

int main() {
	craftState finalState = optimizeLanding();
	return 0;
}
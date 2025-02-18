#include <iostream>
#include <vector>
#include <array>
#include <cmath>

#include "random.h"
#include "curves.h"

//Planet settings
const double planetRadius = 250000.0;      //m
const double planetG = 2.7;                //m s^-2
const double initHeight = 10800.0;         //m
const double landingHeight = 6750.0;       //m
const double minHeight = 4000.0;           //m
const double rotationalPeriod = 1210000.0; //s

//Flight settings
const double initMass = 5361.0;     //kg
const double initXe = 261.0;        //kg
const double initLfOx = 10.7;       //kg
const double ionThrust = 14000.0;   //N
const double chemThrust = 180000.0; //N
const double ionIsp = 4200.0;       //s
const double chemIsp = 305.0;       //s

//Simulation settings
const double timestep = 0.5; //s

//Calculations
const double planetSTDGP = planetRadius * planetRadius * planetG;
const double orbitalSpeed = std::sqrt(
	planetSTDGP / (planetRadius + initHeight)
);
const double rotationalSpeed = 2 * pi * (planetRadius + landingHeight) / 
	rotationalPeriod;
const double landingGrav = planetSTDGP / 
	std::pow(planetRadius + landingHeight, 2.0);

const double ionXe = ionThrust/(ionIsp * g0);
const double chemLfOx = chemThrust/(chemIsp * g0);

struct craftState {
	double mass;
	double xe;
	double lfox;
};

struct Code {
	std::vector<std::array<double, 2> > code;

	double startAlt;
	TaylorCurve vland;
	TaylorCurve vascent;
	double takeOffAlt;

	double mutability;
};

bool works(craftState state) {
	return (state.xe > 0) && (state.lfox > 0);
}

bool greater(craftState a, craftState b) { //True if a > b
	return a.xe > b.xe;
}

craftState simulateLanding(Code code, bool print) {
	double mass = initMass;
	double deorbit;
	if (code.vland.coeffs[0] < 0) {
		deorbit = std::sqrt(-code.vland.coeffs[0] * ionThrust * std::pow(planetRadius + code.startAlt, 1.5)/(mass * std::sqrt(planetSTDGP)));
	} else {
		deorbit = code.vland.coeffs[0]; //Temporary, please fix
	}
	mass *= std::exp(-(
		std::abs(
			std::sqrt(planetSTDGP * (2/(planetRadius + initHeight) - 2/(planetRadius*2 + code.startAlt + initHeight))) - 
			orbitalSpeed
		) + std::abs(
			std::sqrt(planetSTDGP/(planetRadius + code.startAlt)) - 
			std::sqrt(planetSTDGP * (2/(planetRadius + code.startAlt) - 2/(planetRadius*2 + code.startAlt + initHeight)))
		) + deorbit
	)/(ionIsp * g0));
	double xe = initXe + mass - initMass;
	double lfox = initLfOx;

	double vel = std::sqrt(planetSTDGP/(planetRadius + code.startAlt)) - deorbit;
	TaylorCurve currVCurve = code.vland;
	double vvel = currVCurve.evaluate(0);
	double dist = 0.0;
	double time = 0.0;

	double alt = code.startAlt;

	double accel;
	double vaccel;
	double haccel;

	int codeIdx = 0;
	double changeTime = 0.0;

	double currCodeChemThrust = code.code[0][1];
	std::size_t numCodes = code.code.size();
	double dir = -1.0;

	vel = rotationalSpeed - 0.1;

	while ((vel <= std::sqrt(planetSTDGP/(planetRadius + alt))) || (dir < 0)) {
		//Update code
		if ((codeIdx < (numCodes - 1)) && 
			(time > changeTime + code.code[codeIdx+1][0])) {
			codeIdx++;
			changeTime += code.code[codeIdx][0];
			currCodeChemThrust = code.code[codeIdx][1];
		}
		if (lfox == 0.0) {
			currCodeChemThrust = 0.0;
		}

		//Calculate thrust
		accel = (ionThrust + currCodeChemThrust * chemThrust)/mass;
		//Gravity + Centrifugal force
		vaccel = 
			currVCurve.derivative(time) + 
			planetSTDGP/std::pow(planetRadius + alt, 2.0) - 
			vel * vel/(planetRadius + alt);
		if (vaccel > accel) { //Cannot support craft
			xe = -1.0;
			lfox = -1.0;
			break;
		}
		//Pythagoras + Coriolis force
		haccel = 
			std::sqrt(accel * accel - vaccel * vaccel) - 
			dir * vel * vvel/(planetRadius + alt);

		//Integrate
		vel += dir * haccel * timestep;
		alt += vvel * timestep;
		dist += vel * timestep;
		time += timestep;
		xe -= ionXe * timestep;
		lfox -= chemLfOx * currCodeChemThrust * timestep;
		mass -= (ionXe + chemLfOx * currCodeChemThrust) * timestep;
		vvel = currVCurve.evaluate(time);
		if (xe < 0.0) //Ran out of fuel
			break;
		if (lfox < 0.0)
			lfox = 0.0;

		if (alt < landingHeight) {
			xe = -1.0;
			lfox = -1.0;
			break;
		}

		if (vel <= rotationalSpeed) {
			time = 0;
			dir = 1.0;
			vel = rotationalSpeed;
			currVCurve = code.vascent;

			vaccel = planetSTDGP /
				std::pow(planetRadius + landingHeight, 2.0);
			accel = ionThrust/mass;
			if (vaccel > accel) {
				xe = -1.0;
				lfox = -1.0;
				break;
			}
			accel -= vaccel;
			/*
			//Crude approximation for suicide burn
			if (vvel * vvel > 2 * accel * (alt - landingHeight)) {
				xe = -1.0;
				lfox = -1.0;
				break;
			}*/
			double takeOffHeight = code.takeOffAlt - landingHeight;
			double neededHeight = std::pow(currVCurve.evaluate(time), 2)/(2 * accel);
			if (neededHeight > takeOffHeight) {
				takeOffHeight = neededHeight;
				if (print) std::cout << 
					"Note: taking off at " << takeOffHeight + landingHeight << "m, "
					"not " << code.takeOffAlt << "m\n";
			}

			xe -= ionXe * 
				(
					/*std::sqrt(
						(2 * vaccel * (alt - landingHeight) + vvel*vvel)/
						(accel * (vaccel + accel))) + */
					std::sqrt(
						(2 * vaccel * takeOffHeight + std::pow(currVCurve.evaluate(time), 2))/
						(accel * (vaccel + accel)))
				);
			alt = takeOffHeight + landingHeight;
		}

		if (print)
			std::cout << 
				"Time: " << time << "s,\t"
				"Velocity: " << vel << "m/s,\t"
				"Xe: " << xe << "kg,\t"
				"LfOx: " << lfox << "kg,\t"
				"Throttle: " << currCodeChemThrust << ",\t"
				"Dist: " << dist << "m,\t"
				"Alt: " << alt << "m\n";
	}

	craftState out;
	out.mass = mass;
	out.lfox = lfox;
	out.xe = xe;
	return out;
}

craftState optimizeLanding() {
	Code code;
	Code oldCode;
	Code bestCode;

	craftState state;
	state.xe = -1.0;
	craftState oldState;
	oldState.xe = -1.0;
	craftState bestState;
	bestState.xe = -1.0;

	std::array<double, 2> key;

	int numCodes;

	for (int i = 0; i < 1000000; i++) {
		if ((i == 0) || (uniform(generator) < 0.1)) {
			code.code.clear();

			numCodes = numCodeSegments(generator);
			for (int j = 0; j < numCodes; j++) {
				key[0] = uniform(generator) * 200;
				key[1] = uniform(generator) * 0.001;
				code.code.push_back(key);
			}
			code.code[0][0] = 0.0;

			code.startAlt = minHeight + (landingHeight-minHeight)*std::pow(uniform(generator), 1);
			//code.startAlt = initHeight;
			code.vland.initRandom(5, 1000, 5);
			code.vascent.initRandom(5, 500, 10);

			code.takeOffAlt = landingHeight;

			code.mutability = 0.05 * normal(generator);
		} else if (uniform(generator) < 0.8) {
			code = oldCode;
			numCodes = code.code.size();
			for (int j = 0; j < numCodes; j++) {
				code.code[j][0] += normal(generator) * 10.0;
				if (code.code[j][0] < 0.0) code.code[j][0] = 0.0;

				code.code[j][1] += normal(generator) * 0.0005;
				if (code.code[j][1] < 0.0) code.code[j][1] = 0.0;
				if (code.code[j][1] > 1.0) code.code[j][1] = 1.0;
			}
			code.code[0][0] = 0.0;
		} else if (uniform(generator) < 0.9) {
			code = oldCode;

			code.startAlt += 10 * uniform(generator) * normal(generator);
			if (code.startAlt < minHeight) code.startAlt = minHeight;

			code.vland.mutate(1, code.mutability);
			//code = oldCode;

			code.takeOffAlt += 10 * uniform(generator) * normal(generator);
			if (code.takeOffAlt < landingHeight) code.takeOffAlt = landingHeight;

			code.vascent.mutate(1, code.mutability);

			code.mutability += 0.001 * normal(generator);
		} else if (uniform(generator) < 0.5) {
			code = oldCode;

			code.startAlt = 2 * landingHeight - code.startAlt;
			if (code.startAlt < minHeight) continue;

			code.vland.multiply(-1);
		} else {
			code = oldCode;
			code.vascent.multiply(-1);
		}

		state = simulateLanding(code, false);

		if (uniform(generator) < 0.0001) {
			oldState = bestState;
			oldCode = bestCode;
		}

		if (
			(i == 0) || (
				works(state) && (
					greater(state, oldState) || 
					uniform(generator) < std::exp(0.5 * (state.xe - oldState.xe)/(std::pow(1-bestState.xe/200, 0.5)*(1-i/1e8)))
		))) {
			oldState = state;
			oldCode = code;
		}

		if ((i == 0) || (greater(state, bestState))) {
			bestState = state;
			bestCode = code;
			std::cout << 
				"\n" << i << " " << bestCode.mutability << " Record: " << bestState.xe << " (" << 
				0.01 * round(100 * ionIsp * g0 * std::log(1 + bestState.xe/(initMass - initLfOx - initXe))) << "m/s)\n";
			std::cout << "Start Alt: " << bestCode.startAlt << "m\nVvel 1:\n";
			bestCode.vland.print();
			std::cout << "\nTakeoff Alt: " << bestCode.takeOffAlt << "m\nVvel 2:\n";
			bestCode.vascent.print();
			std::cout << "\nCode:\n";
			std::size_t bestSize = bestCode.code.size();
			double time = 0.0;

			for (std::size_t i = 0; i < bestSize; i++) {
				time += bestCode.code[i][0];
				std::cout << 
					"Time: " << time << "s,\t"
					"Throttle: " << bestCode.code[i][1] << "\n";
			}
		}

		if (i % 100000 == 0) std::cout << i << "\n";
	}

	//Logging
	simulateLanding(bestCode, true);

	std::cout << "Start Alt: " << bestCode.startAlt << "m\n";
	bestCode.vland.print();
	std::cout << "Takeoff Alt: " << bestCode.takeOffAlt << "m\n";
	bestCode.vascent.print();
	std::cout << "\nCode:\n";
	std::size_t bestSize = bestCode.code.size();
	double time = 0.0;

	for (std::size_t i = 0; i < bestSize; i++) {
		time += bestCode.code[i][0];
		std::cout << 
			"Time: " << time << "s,\t"
			"Throttle: " << bestCode.code[i][1] << "\n";
	}

	return bestState;
}

int main() {
	craftState finalState = optimizeLanding();
	return 0;
}
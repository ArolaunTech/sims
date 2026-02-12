#include <iostream>
#include <cmath>
#include <vector>
#include <array>

#include "random.h"
#include "curves.h"
#include "consts.h"

//Planet consts
const double planetRadius = 250000.0; //m
const double planetG = 2.7;           //m/s/s
const double rotationalPeriod = 1210000.0;

//Craft consts
const double initMass = 5568.5;
const double initXe = 442.5;
const double initLFOX = 36.0;

const double ionThrust = 14000.0;
const double ionIsp = 4200.0;
const double rocketThrust = 180000.0;
const double rocketIsp = 305.0;

//Flight consts
const double orbitHeight = 10000.0;
const double landingHeight = 6750.0;
const double minHeight = 3000.0;
const double maxHeight = 20000.0;

//Sim consts
const double timestep = 0.5;

//Optimizer consts
const int minSegments = 8;
const int maxSegments = 13;

//Calculations
const double planetSTDGP = planetRadius * planetRadius * planetG;
const double ionFuel = ionThrust/(ionIsp * g0);
const double rocketFuel = rocketThrust/(rocketIsp * g0);
const double rotVel = 2 * pi * (planetRadius + landingHeight)/rotationalPeriod;
const double orbitVel = std::sqrt(planetSTDGP/(planetRadius + orbitHeight));
const double landingGrav = planetSTDGP/std::pow(planetRadius + landingHeight, 2);

//Structs
struct Code {
	FourierCurve xland; //Gives position with time
	FourierCurve yland;

	FourierCurve xascent;
	FourierCurve yascent;
};

struct craftState {
	double mass;
	double xe;
	double lfox;
};

//Functions
craftState calcFailedState(double mass, double time, double lfox, double state) {
	craftState out;
	out.mass = mass;
	out.lfox = lfox;
	out.xe = -1e6 + 1e5 * state + time;
	return out;
}

//Simulate
craftState simulate(Code code, bool print) {
	FourierCurve dxdt = code.xland.derivativeFunction();
	FourierCurve dydt = code.yland.derivativeFunction();
	FourierCurve axt = dxdt.derivativeFunction();
	FourierCurve ayt = dydt.derivativeFunction();

	double mass = initMass;
	double lfox = initLFOX;
	double xe = initXe;

	double x = code.xland.evaluate(0);
	double y = code.yland.evaluate(0);
	double vx, vy, vh, vv, ax, ay;
	vh = orbitVel;

	double accel, thrust;
	double currIonThrust, currChemThrust;

	double grav, gx, gy;

	double dist;

	double time = 0;
	while (vh > rotVel) {
		//Check distance
		dist = std::sqrt(x * x + y * y);
		if (dist > planetRadius + maxHeight) return calcFailedState(mass, time, lfox, 0);
		if (dist < planetRadius + minHeight) return calcFailedState(mass, time, lfox, 0.25);

		//Accel
		x = code.xland.evaluate(time);
		y = code.yland.evaluate(time);
		vx = dxdt.evaluate(time);
		vy = dydt.evaluate(time);
		ax = axt.evaluate(time);
		ay = ayt.evaluate(time);

		vh = (y*vx - x*vy)/dist;
		vv = (x*vx + y*vy)/dist;
		if (vh > orbitVel) return calcFailedState(mass, time, lfox, 0.375);

		//Gravity
		grav = planetSTDGP/(dist * dist);
		gx = -grav * x/dist;
		gy = -grav * y/dist;

		//Thrust
		ax -= gx;
		ay -= gy;
		accel = std::sqrt(ax * ax + ay * ay);
		thrust = mass * accel;
		if (thrust > ionThrust + rocketThrust) return calcFailedState(mass, time, lfox, 0.5);

		//Throttles
		currIonThrust = thrust;
		currChemThrust = 0;
		if (thrust > ionThrust) {
			currIonThrust = ionThrust;
			currChemThrust = thrust - ionThrust;
			if (lfox == 0) return calcFailedState(mass, time, 0, 1.75-vh/orbitVel);
		}

		//Integrate
		time += timestep;
		xe -= currIonThrust * timestep/(ionIsp * g0);
		lfox -= currChemThrust * timestep/(rocketIsp * g0);
		mass -= currIonThrust * timestep/(ionIsp * g0) + currChemThrust * timestep/(rocketIsp * g0);
		if (lfox < 0) lfox = 0;
		if (xe < 0) return calcFailedState(mass, time, lfox, 1);

		if (print) std::cout << 
			"Time: " << time << "s\t"
			"Mass: " << mass << "kg\t"
			"Xe: " << xe << "kg\t"
			"LFOX: " << lfox << "kg\t"
			"Pos: " << x << ", " << y << "m\t"
			"Dist: " << dist << "m\t"
			"Vel: " << vx << ", " << vy << "m/s\t"
			"Vrot: " << vh << ", " << vv << "m/s\t"
			"Accel: " << accel << "m/s/s\t"
			"Thrust: " << thrust << "N\n";
	}
	if ((dist < planetRadius + landingHeight))
		return calcFailedState(mass, time, 0, 2 + 0.01*dist/(planetRadius + landingHeight)-0.1 * std::abs(vv)/orbitVel);

	//Suicide burn
	if (vv > 0) vv *= -1;

	double suicideBurnMass = mass;
	double suicideBurnEnergy = -planetSTDGP/dist + 0.5 * vv * vv; // J/kg

	double veltime = std::abs(vv/(ionThrust/mass - landingGrav));
	double disttime = std::sqrt(2 * (dist-planetRadius-landingHeight)/std::abs(ionThrust/mass - landingGrav));
	if (veltime > disttime) disttime = veltime;

	mass -= xe;
	double vel = 0;
	double h = planetRadius + landingHeight;
	double energy = -planetSTDGP/h;
	while (energy < suicideBurnEnergy) {
		grav = planetSTDGP/(h * h);
		vel += (ionThrust/mass - grav) * timestep;
		h += vel * timestep;
		mass += ionFuel * timestep;

		energy = -planetSTDGP/h + 0.5 * vel * vel;
		//if (print) std::cout << 
		//	energy << " " << 
		//	h << " " <<
		//	vel << " " <<
		//	mass << " " << 
		//	suicideBurnEnergy << "\n";
		
		if (vel < 0) 
			return calcFailedState(
				suicideBurnMass, 0, lfox, 2.2 - 1e-7 * disttime
			);
		if (mass > suicideBurnMass) 
			return calcFailedState(
				suicideBurnMass, 0, lfox, 2.2 - 1e-7 * disttime
			);
		if (h > dist)
			return calcFailedState(
				suicideBurnMass, 0, lfox, 2.2 - 1e-7 * disttime
			);
	}

	double low = suicideBurnMass - xe;
	double high = suicideBurnMass;
	double middle;
	bool dec;
	while (high - low > 1e-3) {
		middle = (low + high)/2;

		mass = middle;
		vel = 0;
		h = planetRadius + landingHeight;
		energy = -planetSTDGP/h;
		dec = false;
		while (energy < suicideBurnEnergy) {
			grav = planetSTDGP/(h * h);
			vel += (ionThrust/mass - grav) * timestep;
			h += vel * timestep;
			mass += ionFuel * timestep;
	
			energy = -planetSTDGP/h + 0.5 * vel * vel;
			//if (print) std::cout << energy << " " << suicideBurnEnergy << "\n";
			if ((vel < 0) || (mass > suicideBurnMass) || (h > dist)) {
				dec = true;
				break;
			}
		}
		if (dec) {
			high = middle;
		} else {
			low = middle;
		}
		if (print) std::cout << low << " " << high << "\n";
	}
	if (print) std::cout << 
		"\nSuicide burn penalty: " << suicideBurnMass - middle << "kg\n"
		"Start burning at " << h-planetRadius << "m\n\n";

	xe -= suicideBurnMass - middle;
	mass = middle;

	x = 0;
	y = planetRadius + landingHeight;
	vx = rotVel;
	vy = 0;
	time = 0;

	dxdt = code.xascent.derivativeFunction();
	dydt = code.yascent.derivativeFunction();
	axt = dxdt.derivativeFunction();
	ayt = dydt.derivativeFunction();

	while (vh < orbitVel) {
		//Check distance
		dist = std::sqrt(x * x + y * y);
		if (dist > planetRadius + maxHeight) return calcFailedState(mass, time, lfox, 2.5+0.00001*middle);
		if (dist < planetRadius + minHeight) return calcFailedState(mass, time, lfox, 2.75+0.00001*middle);

		//Accel
		x = code.xascent.evaluate(time);
		y = code.yascent.evaluate(time);
		vx = dxdt.evaluate(time);
		vy = dydt.evaluate(time);
		ax = axt.evaluate(time);
		ay = ayt.evaluate(time);

		vh = (y*vx - x*vy)/dist;
		vv = (x*vx + y*vy)/dist;
		if (vh < rotVel) return calcFailedState(mass, time, lfox, 2.8+0.00001*middle);

		//Gravity
		grav = planetSTDGP/(dist * dist);
		gx = -grav * x/dist;
		gy = -grav * y/dist;

		//Thrust
		ax -= gx;
		ay -= gy;
		accel = std::sqrt(ax * ax + ay * ay);
		thrust = mass * accel;
		if (thrust > ionThrust + rocketThrust) return calcFailedState(mass, time, lfox, 3+0.00001*middle);

		//Throttles
		currIonThrust = thrust;
		currChemThrust = 0;
		if (thrust > ionThrust) {
			currIonThrust = ionThrust;
			currChemThrust = thrust - ionThrust;
			if (lfox == 0) return calcFailedState(mass, time, 0, 3+vh/orbitVel+0.00001*middle-std::abs(vv)/orbitVel);
		}

		//Integrate
		time += timestep;
		xe -= currIonThrust * timestep/(ionIsp * g0);
		lfox -= currChemThrust * timestep/(rocketIsp * g0);
		mass -= currIonThrust * timestep/(ionIsp * g0) + currChemThrust * timestep/(rocketIsp * g0);
		if (lfox < 0) lfox = 0;
		if (xe < 0) return calcFailedState(mass, time, lfox, 4+0.00001*middle+vh/orbitVel-std::abs(vv)/orbitVel);

		if (print) std::cout << 
			"Time: " << time << "s\t"
			"Mass: " << mass << "kg\t"
			"Xe: " << xe << "kg\t"
			"LFOX: " << lfox << "kg\t"
			"Pos: " << x << ", " << y << "m\t"
			"Dist: " << dist << "m\t"
			"Vel: " << vx << ", " << vy << "m/s\t"
			"Vrot: " << vh << ", " << vv << "m/s\t"
			"Accel: " << accel << "m/s/s\t"
			"Thrust: " << thrust << "N\n";
	}

	craftState out;
	out.mass = mass;
	out.lfox = lfox;
	out.xe = xe;
	return out;
}

craftState optimize() {
	/*Code code;
	Code oldCode;
	Code bestCode;
	craftState state;
	craftState oldState;
	craftState bestState;

	for (int i = 0; i < 100000000; i++) {
		if ((i == 0) || (uniform(generator) < 0.1)) {
			code.xland.initRandom(10, 700, 260000);
			code.yland.initRandom(10, 700, 260000);
			code.xland.coeffs[0] = 0;
			code.yland.coeffs[0] = planetRadius + orbitHeight;
			code.xland.coeffs[1] = orbitVel;
			code.yland.coeffs[1] = 0;

			code.xascent.initRandom(10, 700, 260000);
			code.yascent.initRandom(10, 700, 260000);
			code.xascent.coeffs[0] = 0;
			code.yascent.coeffs[0] = planetRadius + landingHeight;
			code.xascent.coeffs[1] = rotVel;
			code.yascent.coeffs[1] = 0;
		} else if (uniform(generator) < 0.2) {
			code.xland.initRandom(10, 700, 260000);
			code.yland.initRandom(10, 700, 260000);
			code.xland.coeffs[0] = 0;
			code.yland.coeffs[0] = planetRadius + orbitHeight;
			code.xland.coeffs[1] = orbitVel;
			code.yland.coeffs[1] = 0;
		} else if (uniform(generator) < 0.25) {
			code.xascent.initRandom(10, 700, 260000);
			code.yascent.initRandom(10, 700, 260000);
			code.xascent.coeffs[0] = 0;
			code.yascent.coeffs[0] = planetRadius + landingHeight;
			code.xascent.coeffs[1] = rotVel;
			code.yascent.coeffs[1] = 0;
		} else if (uniform(generator) < 0.5) {
			code = oldCode;
			code.xland.mutate(1, 0.001);
			code.yland.mutate(1, 0.001);
			code.xland.coeffs[0] = 0;
			code.yland.coeffs[0] = planetRadius + orbitHeight;
			code.xland.coeffs[1] = orbitVel;
			code.yland.coeffs[1] = 0;
		} else {
			code = oldCode;
			code.xascent.mutate(1, 0.001);
			code.yascent.mutate(1, 0.001);
			code.xascent.coeffs[0] = 0;
			code.yascent.coeffs[0] = planetRadius + landingHeight;
			code.xascent.coeffs[1] = rotVel;
			code.yascent.coeffs[1] = 0;
		}

		state = simulate(code, false);

		if (
			(i == 0) || 
			(state.xe > oldState.xe) || 
			(uniform(generator) < std::exp(0.5 * (state.xe - oldState.xe)/(std::pow(1-bestState.xe/200, 0.5)*(1-i/1e8))))
		) {
			oldState = state;
			oldCode = code;
		}

		if ((i == 0) || (state.xe > bestState.xe)) {
			bestState = state;
			bestCode = code;
			std::cout << 
				bestState.mass << " " << 
				bestState.xe << " " << 
				i << "\n";
		}

		if (i % 2000 == 0) {
			oldState = bestState;
			oldCode = bestCode;
		} 

		if (i % 100000 == 0) simulate(bestCode, true);
	}

	simulate(bestCode, true);

	return bestState;*/
	const int numCodes = 100;
	const int numCoeffs = 10;
	std::array<Code, numCodes> codes;
	std::array<Code, numCodes> vcodes;
	std::array<Code, numCodes> bestCodes;
	std::array<craftState, numCodes> bestStates;

	Code code;
	Code bestCode;

	craftState state;
	craftState bestState;
	bestState.xe = -1e9;

	for (int i = 0; i < numCodes; i++) {
		code.xland.initRandom(numCoeffs, 200000.0, 0.05, 1000.0);
		code.yland.initRandom(numCoeffs, 200000.0, 0.05, 1000.0);
		code.xland.setInitialPositionDerivative(0, orbitVel);
		code.yland.setInitialPositionDerivative(planetRadius + orbitHeight, 0);

		code.xascent.initRandom(numCoeffs, 200000.0, 0.05, 1000.0);
		code.yascent.initRandom(numCoeffs, 200000.0, 0.05, 1000.0);
		code.xascent.setInitialPositionDerivative(0, rotVel);
		code.yascent.setInitialPositionDerivative(planetRadius + landingHeight, 0);

		state = simulate(code, false);
		if (state.xe > bestState.xe) {
			bestCode = code;
			bestState = state;
		}

		codes[i] = code;
		bestCodes[i] = code;
		bestStates[i] = state;

		code.xland.initRandom(numCoeffs, 1000.0, 0.1, 1000.0);
		code.yland.initRandom(numCoeffs, 1000.0, 0.1, 1000.0);
		code.xascent.initRandom(numCoeffs, 1000.0, 0.1, 1000.0);
		code.yascent.initRandom(numCoeffs, 1000.0, 0.1, 1000.0);

		vcodes[i] = code;

		std::cout << "A " << bestState.xe << "\n";
	}

	for (int i = 0; i < 100000; i++) {
		for (int j = 0; j < numCodes; j++) {
			FourierCurve dir;
			dir.initRandom(numCoeffs, 100.0, 0.1, 1000.0);

			vcodes[j].xland = 
				0.95 * vcodes[j].xland + 
				dir + 
				0.5 * uniform(generator) * (bestCodes[j].xland - codes[j].xland) + 
				0.5 * uniform(generator) * (bestCode.xland - codes[j].xland);
			vcodes[j].yland = 
				0.95 * vcodes[j].yland + 
				dir + 
				0.5 * uniform(generator) * (bestCodes[j].yland - codes[j].yland) + 
				0.5 * uniform(generator) * (bestCode.yland - codes[j].yland);
			vcodes[j].xascent = 
				0.95 * vcodes[j].xascent + 
				dir + 
				0.5 * uniform(generator) * (bestCodes[j].xascent - codes[j].xascent) + 
				0.5 * uniform(generator) * (bestCode.xascent - codes[j].xascent);
			vcodes[j].yascent = 
				0.95 * vcodes[j].yascent + 
				dir + 
				0.5 * uniform(generator) * (bestCodes[j].yascent - codes[j].yascent) + 
				0.5 * uniform(generator) * (bestCode.yascent - codes[j].yascent);

			codes[j].xland = codes[j].xland + vcodes[j].xland;
			codes[j].yland = codes[j].yland + vcodes[j].yland;
			codes[j].xascent = codes[j].xascent + vcodes[j].xascent;
			codes[j].yascent = codes[j].yascent + vcodes[j].yascent;

			if (uniform(generator) < 0.01) {
				codes[j].xland.initRandom(numCoeffs, 250000.0, 0.1, 1000.0);
				codes[j].yland.initRandom(numCoeffs, 250000.0, 0.1, 1000.0);

				codes[j].xascent.initRandom(numCoeffs, 250000.0, 0.1, 1000.0);
				codes[j].yascent.initRandom(numCoeffs, 250000.0, 0.1, 1000.0);
			}

			if (uniform(generator) < 0.1) {
				codes[j] = bestCode;
				codes[j].xland.mutate(1, 0.1);
				codes[j].yland.mutate(1, 0.1);
				codes[j].xascent.mutate(1, 0.1);
				codes[j].yascent.mutate(1, 0.1);
			}

			codes[j].xland.setInitialPositionDerivative(0, orbitVel);
			codes[j].yland.setInitialPositionDerivative(planetRadius + orbitHeight, 0);

			codes[j].xascent.setInitialPositionDerivative(0, rotVel);
			codes[j].yascent.setInitialPositionDerivative(planetRadius + landingHeight, 0);

			state = simulate(codes[j], false);

			if (state.xe > bestStates[j].xe) {
				bestCodes[j] = codes[j];
				bestStates[j] = state;
			}

			if (state.xe > bestState.xe) {
				bestCode = codes[j];
				bestState = state;
				std::cout << i << " " << 
					bestState.xe << "\n";
			}
		}
		if (i % 1000 == 0) {
			std::cout << i << " " << 
				bestState.xe << "\n";
			bestCode.xland.print();
			bestCode.yland.print();
			bestCode.xascent.print();
			bestCode.yascent.print();

			simulate(bestCode, true);
		}
	}

	return state;
}

int main() {
	/*
	FourierCurve a;
	a.initRandom(10, 250000, 0.1, 1000);
	a.print();
	std::cout << a.evaluate(0) << " " << a.derivative(0) << "\n";
	a.setInitialPositionDerivative(0, orbitVel);
	a.print();
	std::cout << a.evaluate(0) << " " << a.derivative(0) << "\n";*/

	craftState final = optimize();
}
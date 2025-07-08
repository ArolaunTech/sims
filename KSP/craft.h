#include <map>
#include <functional>
#include <cmath>

#include "consts.h"
#include "engines.h"

#ifndef CRAFT_H
#define CRAFT_H

struct Craft {
	std::map<Engine, double, std::greater<Engine> > engines;
	std::map<FuelType, double> fuel_used;
	double mass;

	Craft() {
		fuel_used = std::map<FuelType, double> {
			{FuelType::FUEL_LFOX, 0},
			{FuelType::FUEL_LF, 0},
			{FuelType::FUEL_XENON, 0},
			{FuelType::FUEL_MP, 0}
		};
	}

	void use_dv(double dv) {
		Engine best_engine = engines.begin()->first;

		double isp = best_engine.isp;
		FuelType fuel = best_engine.fuel;

		double fuel_used_engine = mass * (1 - std::exp(-dv / (isp * g0)));

		fuel_used[fuel] += fuel_used_engine;
		mass -= fuel_used_engine;
	}

	double payload_fraction() const {
		double payload = mass;
		double total_mass = mass;

		for (auto const& iterator : fuel_used) {
			payload -= iterator.second / (fuel_ratios.at(iterator.first) - 1);
			total_mass += iterator.second;
		}

		for (auto const& iterator : engines) {
			payload -= iterator.first.mass * iterator.second;
		}

		return payload/total_mass;
	}
};

#endif
#include <map>
#include <string>

#ifndef ENGINES_H
#define ENGINES_H

enum FuelType {
	FUEL_LFOX,
	FUEL_LF,
	FUEL_XENON,
	FUEL_MP
};

const std::map<FuelType, double> fuel_ratios = {
	{FuelType::FUEL_LFOX, 9},
	{FuelType::FUEL_LF, 11},
	{FuelType::FUEL_XENON, 4},
	{FuelType::FUEL_MP, 8.5}
};

struct Engine {
	std::string name;
	double mass;
	double thrust;
	double isp;
	FuelType fuel;

	Engine(std::string newName, double newMass, double newThrust, double newIsp, FuelType newFuelType) {
		name = newName;
		mass = newMass;
		thrust = newThrust;
		isp = newIsp;
		fuel = newFuelType;
	}

	double effective_isp() const {
		return isp * (1 + 1 / (fuel_ratios.at(fuel) - 1));
	}
};

bool operator>(Engine const& lhs, Engine const& rhs) {
	double isp1 = lhs.effective_isp();
	double isp2 = rhs.effective_isp();

	if (isp1 > isp2) return true;
	if (isp2 > isp1) return false;
	return lhs.name > rhs.name;
}

bool operator<(Engine const& lhs, Engine const& rhs) {
	double isp1 = lhs.effective_isp();
	double isp2 = rhs.effective_isp();

	if (isp1 > isp2) return false;
	if (isp2 > isp1) return true;
	return lhs.name < rhs.name;
}

const std::map<std::string, Engine> engines = {
    {"Dawn",      Engine("Dawn",      250,   2000,    4200, FuelType::FUEL_XENON)},
    {"Nerv",      Engine("Nerv",      3000,  60000,   800,  FuelType::FUEL_LF)},
    {"Wolfhound", Engine("Wolfhound", 3300,  375000,  380,  FuelType::FUEL_LFOX)},
    {"Cheetah",   Engine("Cheetah",   1000,  125000,  355,  FuelType::FUEL_LFOX)},
    {"Poodle",    Engine("Poodle",    1750,  250000,  350,  FuelType::FUEL_LFOX)},
    {"Terrier",   Engine("Terrier",   500,   60000,   345,  FuelType::FUEL_LFOX)},
    {"Dart",      Engine("Dart",      1000,  180000,  340,  FuelType::FUEL_LFOX)},
    {"Rhino",     Engine("Rhino",     9000,  2000000, 340,  FuelType::FUEL_LFOX)},
    {"Skiff",     Engine("Skiff",     1600,  300000,  330,  FuelType::FUEL_LFOX)},
    {"Spark",     Engine("Spark",     130,   20000,   320,  FuelType::FUEL_LFOX)},
    {"Swivel",    Engine("Swivel",    1500,  215000,  320,  FuelType::FUEL_LFOX)},
    {"Skipper",   Engine("Skipper",   3000,  650000,  320,  FuelType::FUEL_LFOX)},
    {"Ant",       Engine("Ant",       20,    2000,    315,  FuelType::FUEL_LFOX)},
    {"Vector",    Engine("Vector",    4000,  1000000, 315,  FuelType::FUEL_LFOX)},
    {"Mammoth",   Engine("Mammoth",   15000, 4000000, 315,  FuelType::FUEL_LFOX)},
    {"Reliant",   Engine("Reliant",   1250,  240000,  310,  FuelType::FUEL_LFOX)},
    {"Mainsail",  Engine("Mainsail",  6000,  1500000, 310,  FuelType::FUEL_LFOX)},
    {"Cub",       Engine("Cub",       180,   32000,   310,  FuelType::FUEL_LFOX)},
    {"Bobcat",    Engine("Bobcat",    2000,  400000,  310,  FuelType::FUEL_LFOX)},
    {"Mastodon",  Engine("Mastodon",  5000,  1350000, 305,  FuelType::FUEL_LFOX)},
    {"Thud",      Engine("Thud",      900,   120000,  305,  FuelType::FUEL_LFOX)},
    {"Rapier",    Engine("Rapier",    2000,  180000,  305,  FuelType::FUEL_LFOX)},
    {"Kodiak",    Engine("Kodiak",    1250,  260000,  300,  FuelType::FUEL_LFOX)},
    {"Spider",    Engine("Spider",    20,    2000,    290,  FuelType::FUEL_LFOX)},
    {"Twitch",    Engine("Twitch",    80,    16000,   290,  FuelType::FUEL_LFOX)},
    {"Puff",      Engine("Puff",      90,    20000,   250,  FuelType::FUEL_MP)}
};

#endif
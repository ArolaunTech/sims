#include <string>
#include <filesystem>
#include <vector>
#include <limits>
#include <utility>
#include <unordered_map>

#include "consts.h"
#include "configreader.h"
#include "engines2.h"

#ifndef PARTLOADER_H
#define PARTLOADER_H

struct PartCategories {
	std::vector<std::pair<std::string, JetEngine> > jets;
};

PartCategories loadParts(std::string path) {
	std::string colloquial = readfile(path + "/colloquial.txt") + "\n";
	std::string word;
	std::vector<std::string> words;
	std::unordered_map<std::string, std::vector<std::string> > colloquialmap;

	bool second = false;

	for (std::size_t i = 0; i < colloquial.size(); i++) {
		switch(colloquial[i]) {
		case '=':
			words.push_back(word);
			word = "";
			second = true;
			break;
		case ',':
			if (second) {
				words.push_back(word);
				word = "";
			} else {
				word += ",";
			}
			break;
		case '\n':
			words.push_back(word);
			word = "";
			second = false;

			if (words.size() > 1) {
				colloquialmap[words[0]] = std::vector<std::string>{};

				for (std::size_t j = 1; j < words.size(); j++) {
					colloquialmap[words[0]].push_back(words[j]);
				}
			}

			words.clear();
			break;
		default:
			word += colloquial[i];
		}
	}

	PartCategories out;

	int id = 0;

	for (const auto& entry : std::filesystem::directory_iterator(path)) {
		if (std::filesystem::path(entry.path()).extension() != ".cfg") continue;

		// Read
		std::string config = readfile(entry.path());

		// Parse
		int end;
		ConfigTree tree = parsetree(config, 0, "root", end);

		// Get nodes
		const ConfigTree* part = nullptr;
		for (const auto& child : tree.children) {
			part = &child.second;
		}

		std::string partname = getvalue(*part, "name")[0];
		double partmass = 1000 * std::stod(getvalue(*part, "mass")[0]);

		// Create parts
		for (const auto& child : part->children) {
			if (child.first != "MODULE") continue;

			const ConfigTree* childnode = &child.second;

			std::string childname = getvalue(*childnode, "name")[0];
			std::string enginetype;

			for (std::size_t i = 0; i < childnode->words.size() - 2; i++) {
				if (childnode->words[i] == "EngineType" && childnode->words[i+1] == "=") {
					enginetype = childnode->words[i + 2];
				}
			}

			if (childname == "ModuleEngines" || childname == "ModuleEnginesFX") {
				if (enginetype == "Turbine") {
					// Jet Engine
					JetEngine partstruct;

					// id
					partstruct.id = id;

					// mass
					partstruct.mass = partmass;

					// Thrust
					partstruct.thrust = 1000 * std::stod(getvalue(*childnode, "maxThrust")[0]);

					// flowMultCap
					double flowmultcap = std::numeric_limits<double>::infinity();
					std::vector<std::string> flowmultval = getvalue(*childnode, "flowMultCap");

					if (flowmultval.size() > 0) flowmultcap = std::stod(flowmultval[0]);
					partstruct.flowmultcap = flowmultcap;

					// Isp
					for (const auto& child2 : childnode->children) {
						if (child2.first == "atmosphereCurve") {
							partstruct.isp = parsespline(child2.second).ys[0];
						}
					}

					// curves
					for (const auto& child2 : childnode->children) {
						if (child2.first == "velCurve") {
							partstruct.velCurve = parsespline(child2.second);
						}

						if (child2.first == "atmCurve") {
							partstruct.atmCurve = parsespline(child2.second);
						}
					}

					// Add jet engine
					for (const std::string& colloquialname : colloquialmap[partname]) {
						out.jets.push_back(std::make_pair(colloquialname, partstruct));
					}
				}
			}
		}

		//logtree(*part);

		id++;
	}

	return out;
}

#endif
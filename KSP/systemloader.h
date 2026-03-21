#include <filesystem>

#include "configreader.h"
#include "system.h"
#include "consts.h"

#ifndef SYSTEMLOADER_H
#define SYSTEMLOADER_H

System loadSystem(std::string filepath) {
	System out;

	for (const auto& entry : std::filesystem::directory_iterator(filepath + "/Configs")) {
		if (std::filesystem::path(entry.path()).extension() != ".cfg") continue;
		
		// Read file
		std::string config = readfile(entry.path());

		// Parse tree
		int end;
		ConfigTree tree = parsetree(config, 0, "root", end);
		const ConfigTree* curr = nullptr;

		for (const auto& child : tree.children) {
			curr = &child.second;
		}

		for (const auto& child : curr->children) {
			if (child.first == "Body") {
				curr = &child.second;
				break;
			}
		}

		// Get relevant nodes
		const ConfigTree* properties = nullptr;
		const ConfigTree* atmosphere = nullptr;
		
		for (const auto& child : curr->children) {
			if (child.first == "Properties") {
				properties = &child.second;
			}

			if (child.first == "Atmosphere") {
				atmosphere = &child.second;
			}
		}

		// Create body
		Body body;
		for (std::size_t i = 0; i < curr->words.size() - 2; i++) {
			if (curr->words[i] == "name" && curr->words[i + 1] == "=") {
				body.name = curr->words[i + 2];
			}
		}

		for (std::size_t i = 0; i < properties->words.size() - 2; i++) {
			if (properties->words[i] == "radius" && properties->words[i + 1] == "=") {
				body.radius = std::stod(properties->words[i + 2]);
			}

			if (properties->words[i] == "gravParameter" && properties->words[i + 1] == "=") {
				body.GM = std::stod(properties->words[i + 2]);
			}

			if (properties->words[i] == "rotationPeriod" && properties->words[i + 1] == "=") {
				body.rotperiod = std::stod(properties->words[i + 2]);
			}
		}

		body.atmospheric = false;
		body.atmosphere.oxygen = false;
		body.atmosphere.normalizedpressurecurve = false;
		body.atmosphere.normalizedtemperaturecurve = false;
		if (atmosphere != nullptr) {
			body.atmospheric = true;

			for (std::size_t i = 0; i < atmosphere->words.size() - 2; i++) {
				if (
					atmosphere->words[i] == "oxygen" && 
					atmosphere->words[i + 1] == "=" && 
					(atmosphere->words[i + 2] == "true" || atmosphere->words[i + 2] == "True")
				) {
					body.atmosphere.oxygen = true;
				}

				if (atmosphere->words[i] == "adiabaticIndex" && atmosphere->words[i + 1] == "=") {
					body.atmosphere.adiabaticindex = std::stod(atmosphere->words[i + 2]);
				}

				if (atmosphere->words[i] == "atmosphereMolarMass" && atmosphere->words[i + 1] == "=") {
					body.atmosphere.molarmass = std::stod(atmosphere->words[i + 2]);
				}

				if (
					(
						atmosphere->words[i] == "atmosphereDepth" || 
						atmosphere->words[i] == "altitude" || 
						atmosphere->words[i] == "maxAltitude"
					) && atmosphere->words[i + 1] == "="
				) {
					body.atmosphere.maxalt = std::stod(atmosphere->words[i + 2]);
				}

				if (
					atmosphere->words[i] == "pressureCurveIsNormalized" && 
					atmosphere->words[i + 1] == "=" && 
					(atmosphere->words[i + 2] == "true" || atmosphere->words[i + 2] == "True")
				) {
					body.atmosphere.normalizedpressurecurve = true;
				}

				if (
					atmosphere->words[i] == "temperatureCurveIsNormalized" && 
					atmosphere->words[i + 1] == "=" && 
					(atmosphere->words[i + 2] == "true" || atmosphere->words[i + 2] == "True")
				) {
					body.atmosphere.normalizedtemperaturecurve = true;
				}
			}

			for (const auto& child : atmosphere->children) {
				if (child.first == "pressureCurve") {
					body.atmosphere.pressurecurve = parsespline(child.second);
				}

				if (child.first == "temperatureCurve") {
					body.atmosphere.temperaturecurve = parsespline(child.second);
				}

				if (child.first == "temperatureSunMultCurve") {
					body.atmosphere.tempsunmultcurve = parsespline(child.second);
				}

				if (child.first == "temperatureLatitudeBiasCurve") {
					body.atmosphere.templatitudebiascurve = parsespline(child.second);
				}

				if (child.first == "temperatureLatitudeSunMultCurve") {
					body.atmosphere.templatitudesunmultcurve = parsespline(child.second);
				}

				if (child.first == "temperatureAxialSunBiasCurve") {
					body.atmosphere.tempaxialsunbiascurve = parsespline(child.second);
				}

				if (child.first == "temperatureAxialSunMultCurve") {
					body.atmosphere.tempaxialsunmultcurve = parsespline(child.second);
				}

				if (child.first == "temperatureEccentricityBiasCurve") {
					body.atmosphere.tempeccentricitybiascurve = parsespline(child.second);
				}
			}
		}

		// Add to system
		out.bodies[body.name] = body;

		//logtree(*curr);
	}

	return out;
}

#endif
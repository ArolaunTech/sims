#include <string>
#include <filesystem>
#include <vector>

#include "consts.h"
#include "configreader.h"
#include "engines.h"

#ifndef PARTLOADER_H
#define PARTLOADER_H

struct PartCategories {
	std::vector<Engine> engines;
};

PartCategories loadParts(std::string path) {
	PartCategories out;

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

		// Determine part type

		// Create part

		logtree(*part);
	}

	return out;
}

#endif
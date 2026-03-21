#include <filesystem>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <utility>

#include "curves.h"

#ifndef CONFIGREADER_H
#define CONFIGREADER_H

/*===== Quick, dirty cfg file reader that is able to decode planet cfg files. May be improved later. =====*/

struct ConfigTree {
	std::string name;

	std::vector<std::string> words;

	std::vector<std::pair<std::string, ConfigTree> > children;
};

void logtree(const ConfigTree& tree, int tabs = 0) {
	for (int i = 0; i < tabs; i++) std::cout << "\t";
	std::cout << tree.name << " {\n";

	for (int i = 0; i < tree.words.size(); i++) {
		for (int j = 0; j <= tabs; j++) std::cout << "\t";

		std::cout << tree.words[i] << "\n";
	}

	for (const auto& childpair : tree.children) {
		logtree(childpair.second, tabs + 1);
	}

	for (int i = 0; i < tabs; i++) std::cout << "\t";
	std::cout << "}\n";
}

std::string readfile(std::string path) {
	std::ifstream reader(path);

	std::string out((std::istreambuf_iterator<char>(reader)), std::istreambuf_iterator<char>());

	return out;
}

HermiteSpline parsespline(const ConfigTree& node) {
	HermiteSpline out;

	for (std::size_t i = 0; i < node.words.size() - 5; i++) {
		if (node.words[i] == "key" && node.words[i + 1] == "=") {
			out.xs.push_back(std::stod(node.words[i + 2]));
			out.ys.push_back(std::stod(node.words[i + 3]));
			out.ldys.push_back(std::stod(node.words[i + 4]));
			out.rdys.push_back(std::stod(node.words[i + 5]));
		}
	}

	return out;
}

ConfigTree parsetree(const std::string& config, int index, std::string name, int& end) {
	ConfigTree out;

	out.name = name;

	// If a pattern is a prefix of another then this is ambiguous
	const std::array<std::string, 6> patterns = {
		" ",
		"\n",
		"\t",
		"//",
		"{",
		"="
	};

	const std::array<bool, 6> addtowords = {
		false,
		false,
		false,
		true,
		true,
		true
	};

	std::vector<std::string> words;
	std::string word = "";

	int i = index;
	int childend;

	while (i < config.size()) {
		//std::cout << i << "\n";
		//for (int temp = 0; temp < words.size(); temp++) {
		//	std::cout << words[temp] << " ";
		//}
		//std::cout << "\n";
		//for (const auto& child : out.children) {
		//	logtree(child.second);
		//}
		//std::cout << "\n";

		if (config[i] == '}') {
			break;
		}

		bool found = false;
		for (int j = 0; j < patterns.size(); j++) {
			bool valid = true;
			for (int k = 0; k < patterns[j].size(); k++) {
				if (i + k >= config.size()) {
					valid = false;
					break;
				}

				if (config[i + k] != patterns[j][k]) {
					valid = false;
					break;
				}
			}

			if (valid) {
				if (word.size() > 0) words.push_back(word);
				word = "";

				if (addtowords[j]) words.push_back(patterns[j]);
				i += patterns[j].size();
				found = true;

				break;
			}
		}

		if (!found) {
			word += config[i];
			i++;
		} else if (words.size() > 0) {
			if (words[words.size() - 1] == "{") {
				words.pop_back();

				while (words.size() > 0) {
					bool pattern = false;

					for (int k = 0; k < patterns.size(); k++) {
						if (patterns[k] == words[words.size() - 1]) {
							pattern = true;
							break;
						}
					}

					if (!pattern) break;

					words.pop_back();
				}

				if (words.size() == 0) words.push_back(" ");

				out.children.push_back(std::make_pair(words[words.size() - 1], parsetree(config, i, words[words.size() - 1], childend)));
				i = childend;

				words.pop_back();
			} else if (words[words.size() - 1] == "//") {
				words.pop_back();

				while (i < config.size()) {
					if (config[i] == '\n') break;

					i++;
				}
			}
		}
	}
	if (word.size() > 0) words.push_back(word);

	out.words = words;

	end = i + 1;

	return out;
}

#endif
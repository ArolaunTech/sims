#include <iostream>
#include <iomanip>

#ifndef DEBUG_H
#define DEBUG_H

void logDoubleVector(std::vector<double> arr) {
	std::cout << "std::vector<double> {";

	//std::cout << std::fixed << std::setprecision(2);

	std::size_t numElements = arr.size();
	for (std::size_t i = 0; i < numElements; i++) {
		std::cout << i << ": " << arr[i] << ", ";
	}

	std::cout << "}";
}

#endif
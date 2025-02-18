#include <random>

#ifndef RANDOM
#define RANDOM

std::mt19937 generator(time(NULL));
std::uniform_real_distribution<double> uniform(0.0, 1.0);
std::normal_distribution<double> normal(0.0, 1.0);
std::uniform_int_distribution<int> randomint(0, 1000000000);

std::uniform_int_distribution<int> numCodeSegments(7,9);

int randint(int a, int b) {
	return a + (randomint(generator)) % (b - a + 1);
}

#endif
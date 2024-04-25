#include "Helpers.hpp"

#include <random>

std::mt19937 rng;

struct __SetupRand_ { public: __SetupRand_() { randSeed(time(NULL)); } }; __SetupRand_ __sr_;

void randSeed(unsigned int seed)
{
	rng.seed(seed);
}
int randI(int firstInt, int secondInt) ///inclusive
{
	if (firstInt > secondInt)
		std::swap(firstInt, secondInt);
	std::uniform_int_distribution<> dist(firstInt, secondInt);
	return dist(rng);
}
float randF(double firstFloat, double secondFloat)
{
	if (firstFloat > secondFloat)
		std::swap(firstFloat, secondFloat);
	std::uniform_real_distribution<> dist(firstFloat, secondFloat);
	return dist(rng);
}
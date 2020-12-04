#include "definitions.hpp"

std::random_device rd;
std::mt19937_64 rng(rd());
std::uniform_int_distribution<int> binary(0,1);
std::uniform_real_distribution<double> unif(0.0,1.0);

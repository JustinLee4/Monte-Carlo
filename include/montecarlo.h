#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include "atom.h"


#include <iostream>
#include <vector>
#include <random>

bool iterate_singly(std::vector<Water> &input_vec, std::vector<int> const indices, std::mt19937& gen);

std::tuple<std::vector<Water>,size_t,size_t> randomize_states(std::vector<Water> &input_vec);

std::tuple<double, bool> metropolis(double new_energy, double old_energy);

void zero_states(std::vector<Water> &input_vec);


#endif
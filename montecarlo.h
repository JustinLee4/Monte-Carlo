#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include "water.h"


#include <iostream>
#include <vector>


std::tuple<std::vector<Water>,size_t,size_t> randomize_states(std::vector<Water> &input_vec);

std::tuple<double, bool> metropolis(double new_energy, double old_energy);



#endif
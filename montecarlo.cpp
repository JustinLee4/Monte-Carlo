#include "montecarlo.h"

#include <random>
#include <algorithm>
#include <numeric>
#include <tuple>
#include <cmath>


std::tuple<std::vector<Water>,size_t,size_t> randomize_states(std::vector<Water> &input_vec) {


    static std::random_device rd;
    static std::mt19937 gen(rd());

    size_t N = input_vec.size();

    std::uniform_int_distribution<size_t> count_distrib(0, N);
    size_t k = count_distrib(gen);
    // size_t k = 1;
    for(int i = 0; i < N; i++){
        input_vec[i].setValue(0);
    }

    // std::cout << "N = " << N << ", k = " << k << std::endl;

    std::vector<size_t> indices(N);
    std::iota(indices.begin(), indices.end(), 0);

    std::shuffle(indices.begin(), indices.end(), gen);

    for(int i = 0; i < k; i++) {
        input_vec[indices[i]].setValue(true);
    }
    // for(int i = 0; i < N; i++) {
    //     input_vec[indices[i]].setValue(true);
    // }
    std::tuple<std::vector<Water>,size_t,size_t> ans = {input_vec, N, k};
    return ans;

}

std::tuple<double, bool> metropolis(double new_energy, double old_energy) {

    //dowser outputs energy in kcal/mol

    double kbT = .6;
    static std::random_device rd;
    static std::mt19937 gen(rd());

    
    double exponent = -(new_energy - old_energy) / (kbT);
    std::uniform_real_distribution<double> dis (0.0, 1.0);
    double rand = dis(gen);

    if(new_energy <= old_energy) {
        return {new_energy, true};
    } else if (pow(M_E, exponent) >= rand)
    {
        return {new_energy, true};
    } else {
        return {old_energy, false};
    }
}






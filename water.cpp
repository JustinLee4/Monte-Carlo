#include "water.h"

std::string Water::toString() {
    std::stringstream ss;

            // 2. Write the formatted data into the stream
            ss << "  Value (0/1): " << (value ? 1 : 0) << "\n";
            ss << "  Position (x, y, z): [";
            
            // Print the array elements
            ss << position[0] << ", " 
            << position[1] << ", " 
            << position[2] << "]" << "\n";
            
            ss << "energy = " << energy << "\n";
            // 3. Extract and return the final string from the stringstream
            return ss.str();
}

std::array<double, 3> Water::getCoords() const {
    return position;
}

void Water::setValue(bool new_value) {
    value = new_value;
}

bool Water::getValue() const{
    return value;
}
double Water::getEnergy() const{
    return energy;
}

void Water::setOverlap(bool new_isOverlapping){
    isOverlapping = new_isOverlapping;
}

bool Water::getOverlap() const{
    return isOverlapping;
}

void Water::reset_neighbors() {
    nearest_neighbors = {};
}

void Water::add_neighbors(Water* newNeighbor) {
    nearest_neighbors.push_back(newNeighbor);
}

double Water::constructive_interaction(){

    double ans = 0;

    for(int i = 0; i < nearest_neighbors.size(); i++) {


        //from https://docs.lammps.org/Howto_tip3p.html
        //coulomb interaction (constructive)
        // V =  k q^2 / r
        // k = 332.1 Å·kcal/(mol·e²)
        // q(oxygen) = -0.834 e
        // q(hydrogen) = 0.417 e



        //lennard-jones interaction (repulsive)
        // V = A / r^12 - B / r^6
        // A = 4(eps)(sigma^12)
        // B = 4(eps)(sigma^6)
        // eps = 0.1521 kcal/mol
        // sigma = 3.1507 Å, currently set to 3 Å for short range cutoff

        
        //for every water in first solvation shell, simply add 5 A
        
    }

    return ans;
}

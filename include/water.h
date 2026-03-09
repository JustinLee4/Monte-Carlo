#ifndef WATER
#define WATER

#include <iostream>
#include <sstream>
#include <array>
#include <vector>

class Water {
    bool value;
    double energy;
    bool isOverlapping;
    std::vector<Water*> nearest_neighbors = {};
    std::array<double, 3> position;
    
    public:
    Water(bool init_value, double init_energy, std::array<double,3> init_position):
        value(init_value),
        energy(init_energy),
        isOverlapping(false),
        nearest_neighbors({})
        {
            if (init_position.size() == 3) {
                std::copy(init_position.begin(), init_position.end(), position.begin());
            } else {
            // Handle error or default initialization if necessary
            std::cerr << "Warning: Position list must have exactly 3 elements." << std::endl;
            }
        }
        
    // print binary value and coordinates as a string
    std::string toString();

    //get coordinate information as an array of doubles
    std::array<double, 3> getCoords() const;

    void setValue(bool new_value);

    bool getValue() const;

    double getEnergy() const;

    void setOverlap(bool new_isOverlapping);

    bool getOverlap() const;

    void reset_neighbors();

    void add_neighbors(Water* newNeighbor);

    //return constructive interaction energy in kCal/mol
    double constructive_interaction();

    ~Water() {
    }
};

#endif
#ifndef WATER
#define WATER

#include <iostream>
#include <sstream>
#include <array>
#include <vector>
#include <unordered_map>
#include <random>


class ActiveList {
private:
    std::vector<int> active_indices;
    // Hash map maps arbitrary water_index -> position in active_indices array
    std::unordered_map<int, int> position_in_list; 

public:
    // We no longer need to know the total water count in advance!
    ActiveList() = default;

    // Optional: If you know the approximate cluster size, reserving space prevents memory reallocations
    void reserve(size_t expected_size) {
        active_indices.reserve(expected_size);
        position_in_list.reserve(expected_size);
    }

    void add(int index) {
        // If the index is already in the map, do nothing
        if (position_in_list.find(index) != position_in_list.end()) {
            return; 
        }
        
        active_indices.push_back(index);
        position_in_list[index] = active_indices.size() - 1;
    }

    void remove(int index) {
        // Find the index in the map
        auto it = position_in_list.find(index);
        if (it == position_in_list.end()) {
            return; // Not in the list, nothing to remove
        }

        int pos = it->second;
        int last_index = active_indices.back();
        
        // Swap the element we want to remove with the very last element in the vector
        active_indices[pos] = last_index;
        position_in_list[last_index] = pos; // Update the moved element's new position in the map

        // Pop the back of the vector (O(1) removal)
        active_indices.pop_back();
        
        // Completely erase the old index from the map
        position_in_list.erase(it);
    }

    int get_random(std::mt19937& gen) const {
        if (active_indices.empty()) return -1; // Simulation stuck / no valid moves
        std::uniform_int_distribution<int> distrib(0, active_indices.size() - 1);
        return active_indices[distrib(gen)];
    }

    bool is_empty() const { return active_indices.empty(); }
    
    size_t size() const { return active_indices.size(); }
};

class Atom {
    std::string resname;
    std::string atomname;
    double resnumber;
    double radius;
    std::array<double, 3> position;
    double b_factor;
    
    public:
    Atom(std::string init_resname, std::string init_atomname, std::array<double,3> init_position, double init_b_factor = 0.0):
        resname(init_resname),
        atomname(init_atomname),
        radius(0.0),
        b_factor(init_b_factor)
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

    void set_resnumber(double resnumber);

    double get_resnumber() const;

    void set_radius(double new_radius);

    double get_radius();

    std::string get_resname() const;

    std::string get_atomname() const;

    double get_bfactor() const;

    
    virtual ~Atom() {
    }
};

class Water: public Atom {
    public:
    bool value = false;
    int isOverlapping = 0;
    std::vector<int> nearby_water = {};


    Water(std::array<double, 3> init_position, double init_b_factor = 0.0) 
        : Atom("HOH", "O", init_position, init_b_factor) 
    {
        set_radius(1.25); 
    }

    void set_value(bool newval);

    bool get_value() const;

    void setOverlap(int new_isOverlapping);

    int getOverlap() const;

    void addOverlap();

    void subtractOverlap();

    void clear_Overlap();

    void add_neighbor(int newNeighbor);

    void add_overlap_with_neighbors(std::vector<Water>& input_vec, ActiveList& active_list);

    void subtract_overlap_with_neighbors(std::vector<Water>& input_vec, ActiveList& active_list);

};


#endif
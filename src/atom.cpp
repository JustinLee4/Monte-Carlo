#include "atom.h"

#include "common.h"

std::string Atom::toString() {
    std::stringstream ss;

            //Write the formatted data into the stream

            ss << "  Resname : " << resname;
            ss << "  Atomname : " << atomname;
            ss << "  Position (x, y, z): [";
            
            // Print the array elements
            ss << position[0] << ", " 
            << position[1] << ", " 
            << position[2] << "]" << "\n";
            
            return ss.str();
}

std::array<double, 3> Atom::getCoords() const {
    return position;
}

void Atom::set_resnumber(double new_resnumber) {
    resnumber = new_resnumber;
}

double Atom::get_resnumber() const {
    return resnumber;
}

void Atom::set_radius(double new_radius) {
    radius = new_radius;
}

double Atom::get_radius() {
    return radius;
}

std::string Atom::get_resname() const{
    return resname;
}

std::string Atom::get_atomname() const {
    return atomname;
}

double Atom::get_bfactor() const {
    return b_factor;
}

void Water::set_value(bool newval){
    value = newval;
}

bool Water::get_value() const {
    return value;
}

void Water::setOverlap(int new_isOverlapping) {
    isOverlapping = new_isOverlapping;
}

int Water::getOverlap() const {
    return isOverlapping;
}

void Water::addOverlap() {
    isOverlapping += 1;
}

void Water::subtractOverlap() {
    isOverlapping -= 1;
}

void Water::clear_Overlap() {
    isOverlapping = 0;
}

void Water::add_neighbor(int newNeighbor){
    nearby_water.push_back(newNeighbor);
}

void Water::add_overlap_with_neighbors(std::vector<Water>& input_vec, ActiveList& active_list){
    for (int neighbor_idx : nearby_water) {
        input_vec[neighbor_idx].addOverlap();
        active_list.remove(neighbor_idx);
    }
}

void Water::subtract_overlap_with_neighbors(std::vector<Water>& input_vec, ActiveList& active_list) {
    for (int neighbor_idx : nearby_water) {
        input_vec[neighbor_idx].subtractOverlap();
        if (input_vec[neighbor_idx].getOverlap() == 0) {
            active_list.add(neighbor_idx);
        }
    } 
}

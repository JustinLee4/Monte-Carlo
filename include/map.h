#ifndef MAP
#define MAP

#include <vector>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include <iostream>
#include <array>

#include "atom.h"


struct GridKey {
    int x, y, z;

    // We must tell C++ how to check if two keys are equal
    bool operator==(const GridKey& other) const {
        return x == other.x && y == other.y && z == other.z;
    }
};

namespace std {
    template <>
    struct hash<GridKey> {
        size_t operator()(const GridKey& k) const {
            return (
                (std::hash<int>()(k.x) ^
                 (std::hash<int>()(k.y) << 1)) >> 1
            ) ^ (std::hash<int>()(k.z) << 1);
        }
    };
}

GridKey getGridKey(const Atom& input_Atom, double gridCellSize);

GridKey getGridKey_pos(const std::array<double, 3> pos, double gridCellSize);

template<typename T>
std::unordered_map<GridKey, std::vector<int>>
buildSpatialGrid(std::vector<T> objects, double gridCellSize) {

    std::unordered_map<GridKey, std::vector<int>> grid;

    for (int i = 0; i < objects.size(); ++i) {
        GridKey key = getGridKey(objects[i], gridCellSize);
        grid[key].push_back(i);
    }
    return grid;
}


template<typename T>
std::unordered_map<int, std::vector<int>> buildClusterMap(const std::vector<T>& objects) {

    std::unordered_map<int, std::vector<int>> cluster_map;

    for (int i = 0; i < objects.size(); ++i) {
        int res_id = static_cast<int>(objects[i].get_resnumber()); 
        
        cluster_map[res_id].push_back(i);
    }
    
    return cluster_map;
}

void printSpatialGrid(const std::unordered_map<GridKey, std::vector<int>>& grid);

void getOverlap_cluster(const std::unordered_map<GridKey, std::vector<int>>& grid, std::vector<Water>& Watervector, Water& target, double gridCellSize, double diameter);

void testGrid();



#endif
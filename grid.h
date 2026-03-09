#ifndef GRID
#define GRID

#include <vector>
#include <unordered_map>
#include <tuple>
#include <cmath>
#include <iostream>
#include <array>

#include "water.h"


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
            // A common way to combine hashes: use xor (^) and bit-shifting.
            // This creates a good "hash" value from the three integers.
            return (
                (std::hash<int>()(k.x) ^
                 (std::hash<int>()(k.y) << 1)) >> 1
            ) ^ (std::hash<int>()(k.z) << 1);
        }
    };
}

GridKey getGridKey(const Water& input_water, double gridCellSize);

std::unordered_map<GridKey, std::vector<int>>
buildSpatialGrid(std::vector<Water> objects, double gridCellSize);

void printSpatialGrid(const std::unordered_map<GridKey, std::vector<int>>& grid);

void getOverlap_cluster(const std::unordered_map<GridKey, std::vector<int>>& grid, std::vector<Water>& watervector, int target, double gridCellSize);

void testGrid();



#endif
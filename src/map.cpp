#include "map.h"

#include "AtomicRadii.h"
#include "common.h"


GridKey getGridKey(const Atom& input_water, double gridCellSize) {
    std::array<double, 3> pos = input_water.getCoords();
    return {
        (int)std::floor(pos[0] / gridCellSize),
        (int)std::floor(pos[1] / gridCellSize),
        (int)std::floor(pos[2] / gridCellSize)
    };
}

GridKey getGridKey_pos(const std::array<double, 3> pos, double gridCellSize) {
    return {
        (int)std::floor(pos[0] / gridCellSize),
        (int)std::floor(pos[1] / gridCellSize),
        (int)std::floor(pos[2] / gridCellSize)
    };
}

void printSpatialGrid(const std::unordered_map<GridKey, std::vector<int>>& grid) {
    
    std::cout << "--- Spatial Grid Contents ---" << std::endl;
    // grid.size() tells you how many unique cells (keys) were created.
    std::cout << "Total unique cells (keys): " << grid.size() << std::endl;
    std::cout << "-------------------------------" << std::endl;

    // We loop through the map. 'pair' will be one {key, value} item.
    for (const auto& pair : grid) {
        // 'pair.first' is the GridKey
        const GridKey& key = pair.first;
        
        // 'pair.second' is the std::vector<int>
        const std::vector<int>& objectIndices = pair.second;
        
        // Print the Key
        std::cout << "Key (Cell): (" << key.x << ", " << key.y << ", " << key.z << ")" << std::endl;
        
        // Print the Value (the list of indices)
        std::cout << "  Indices : { ";
        for (int index : objectIndices) {
            std::cout << index << " ";
        }
        std::cout << "}" << std::endl;
    }
    std::cout << "-------------------------------" << std::endl;
}


//find overlaps and also update nearest neighbors
void getOverlap_cluster(const std::unordered_map<GridKey, std::vector<int>>& grid, std::vector<Water>& Watervector, Water& target, double gridCellSize, double diameter) { 


    // 1. Get the center coordinates of the target particle
    GridKey centerKey = getGridKey_pos(target.getCoords(), gridCellSize);

    // 2. Loop through X, Y, Z offsets (-1, 0, 1)
    for (int dx = -1; dx <= 1; dx++) {
        for (int dy = -1; dy <= 1; dy++) {
            for (int dz = -1; dz <= 1; dz++) {
                
                // 3. Construct the neighbor key
                // Note: Assuming GridKey is a struct with members x, y, z
                GridKey neighborKey;
                neighborKey.x = centerKey.x + dx;
                neighborKey.y = centerKey.y + dy;
                neighborKey.z = centerKey.z + dz;

                // 4. Try to find this neighbor box in the map
                auto it = grid.find(neighborKey);

                // 5. If the box exists (contains particles)
                if (it != grid.end()) {

                    const std::vector<int>& indices = it->second;
                    
                    for (int neighborIndex : indices) {  

                        std::array<double,3> origCoords = target.getCoords();
                        std::array<double,3> newCoords = Watervector[neighborIndex].getCoords();

                        double distance = 0;
                        double dX = origCoords[0] - newCoords[0];
                        double dY = origCoords[1] - newCoords[1];
                        double dZ = origCoords[2] - newCoords[2];

                        distance = (dX*dX) + (dY*dY) + (dZ*dZ);

                        if ((distance <= pow(diameter+0.01,2)) && (distance > 0.0001)){
                            target.add_neighbor(neighborIndex);
                        }
                    }
                }
            } 
        }
    }
}


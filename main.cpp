
#include <iostream>
#include <cstdio>

#include "pdbtovector.h"
#include "montecarlo.h"
#include "grid.h"
#include "water.h"

bool linux = true;
int total_reps = 100000;

int main(int argc, char* argv[]){

    double lowest_energy = INFINITY;
    double total_energy = 0;
    std::vector<Water> lowest_config;
    bool accepted;
    size_t N,k;
    std::string output_file;
    std::string monte_carlo_output_file;
    std::string input_file;
    std::string input_file_protein;
    std::string energy_file;
    std::ofstream monte_carlo_output;
    double hash_spacing = 10.0;


    // ------------------------------------------

    //setup printing energy file        

    if(argc == 6) {
        output_file = argv[1];      //file with lowest energy configuration
        monte_carlo_output_file = argv[2];      // file with iteration data
        input_file = argv[3];       //water molecules to test
        input_file_protein = argv[4];       //structure file
        energy_file = argv[5];      //dowser waters energy file

        monte_carlo_output.open(monte_carlo_output_file);

    } else {
        std::cout << "incorrect number of inputs\n" << "input should be \n"
        << "./main output_waters.pdb monte_carlo_output.csv input_waters.pdb input_structure.pdb dowser_energy.dow\n";
        return 0;
    }

    //------------------------------------------------------
    //calculate energy for all gridpoints (without water-water interaction)
    if(linux == true) {
    //start by reformatting the input protein file
        std::cout << "-> Entering reformatPDB" << std::endl;
    // int ret1 = std::system("$DOWSER/bin/linux/reformatPDB -pdbin 8OM1_no_ligands.pdb -pdbout reformat.pdb");
    //create placewat.dow w/ energy
        std::cout << "-> Entering placeWat" << std::endl;
    // int ret2 = std::system("$DOWSER/bin/linux/placeWat reformat.pdb clean_gridpoints.pdb both >placewat.dow");
    }

    //--------------------------------------------------------

    std::cout << "-> Entering pdbtovector_Waters" << std::endl;

    //create vector of waters
    std::vector<Water> watervector = pdbtovector_Waters(energy_file);
    // std::vector<Water> watervector = pdbtovector_Waters("placewat.dow");
    // std::vector<Water> watervector = pdbtovector_Waters("placewat_grid1.dow");

    // if(linux == true) {
    //     std::vector<Water> watervector = pdbtovector_Waters("placewat.dow");
    // }
    std::cout << "* total waters = " << watervector.size() << std::endl;

        
    std::cout << "-> Entering buildSpatialGrid" << std::endl;

    //create adjacency map
    std::unordered_map<GridKey, std::vector<int>> map = buildSpatialGrid(watervector, hash_spacing);


    // testGrid();

    // printSpatialGrid(map);



    //--------------------------------------------------------
    //This is the start of the loop
    

    std::cout << "-> Entering loop" << std::endl;

    //from here we are randomizing
    for(int z = 0; z < total_reps; z++) {

        std::cout << "\rIteration: " << z+1 << " of " << total_reps << std::flush;


        
        std::tuple<std::vector<Water>,size_t,size_t> newtuple = randomize_states(watervector);
        N = std::get<1>(newtuple);
        k = std::get<2>(newtuple);
        // std::cout << "Original waters: " << N << "\n";
        // std::cout << "Remaining waters: " << k << "\n";
        // printing all water molecules
        // for (int i = 0; i < watervector.size(); i++){
        //     std::string water = watervector[i].toString();
        //     std::cout << water;
        // }

        //find overlaps in selected waters

        //cleanup last iteration
        for(int i = 0; i < watervector.size(); i++){
            watervector[i].setOverlap(false);
            watervector[i].reset_neighbors();
        }

        for(int i = 0; i < watervector.size(); i++) {
            if(watervector[i].getValue() == 1) {
                if(watervector[i].getOverlap() == false){
                    getOverlap_cluster(map, watervector, i, hash_spacing);
                }
            }
        }


        //printing only water molecules that have been selected to a new pdb file
        // vectortopdb(watervector, "PDB_with_energies_and_overlaps.txt", N, k, true);
        // vectortopdb(watervector, output_file, N, k, false);



        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% printing combined PDB files
        //to print vector as pdb - not necessary if we find energies first

        //appending pdb files

        // std::string firstFile = "8OM1_no_ligands.pdb";
        // std::string secondFile = output_file;
        // std::string resultFile = "combined.pdb";
        // int insertAfterLine = 4250;

        // std::cout << "Attempting to merge files..." << std::endl;

        // if (append_pdb_files(firstFile, secondFile, resultFile, insertAfterLine)) {
        //     std::cout << "Success! Content inserted after line " << insertAfterLine << "." << std::endl;
        // } else {
        //     std::cout << "Operation failed." << std::endl;
        // }

        //TODO calculate energies using dowser
        //  $DOWSER/bin/linux/reformatPDB -pdbin combined.pdb -pdbout reformatted.pdb
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        //calculate energy
        total_energy = 0;
        // std::cout << total_energy << " ";
        for (int i = 0; i< watervector.size(); i++){
            if (watervector[i].getValue() == 1) {
                if(watervector[i].getOverlap()){
                    total_energy += 10;
                } else {
                    total_energy += watervector[i].getEnergy();
                    total_energy += watervector[i].constructive_interaction();
                    if(watervector[i].constructive_interaction()!=0){
                        std::cout << "UH OH" << std::endl;
                    }
                }
                // total_energy += watervector[i].getEnergy();

            }
        }
                // std::cout << total_energy << "\n ";

        //metropolis
        auto[temp, accepted] = metropolis(total_energy, lowest_energy);
        lowest_energy = temp;
        
        if(accepted) {
            lowest_config = std::get<0>(newtuple);
        }

        monte_carlo_output << z << ", " << N << ", " << k << ", " << total_energy << ", " << lowest_energy << "\n";
    
    }
    std::cout << std::endl;
    std::cout << "Completed " << total_reps << " iterations" << std::endl;

    vectortopdb(lowest_config,output_file,N,k,false);

    monte_carlo_output.close();

    return 0;
}

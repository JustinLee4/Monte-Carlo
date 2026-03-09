#include "atom.h"
#include "pdbtovector.h"
#include "montecarlo.h"
#include "map.h"


#include <iostream>
#include <cstdio>


bool linux = true;
int total_reps = 100000;

int main(int argc, char* argv[]){

    double lowest_energy = INFINITY;
    double total_energy = 0;
    std::vector<Water> lowest_config;
    bool accepted;
    size_t N,k;
    std::string monte_carlo_output_file;
    std::ofstream monte_carlo_output;
    double hash_spacing = 10.0;

    // input handling

    std::string input_energy_file = "";
    std::string input_cluster_file = "";
    std::string input_protein_file = "";
    std::string output_file = "";

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if ((arg == "-p" || arg == "--pdb") && i + 1 < argc) {
            input_protein_file = argv[++i]; 
        } 
        else if ((arg == "-e" || arg == "--energy") && i + 1 < argc) {
            input_energy_file = argv[++i];
        } 
        else if ((arg == "-c" || arg == "--clusters") && i + 1 < argc) {
            input_cluster_file = argv[++i];
        }
        else if ((arg == "-o" || arg == "--out") && i + 1 < argc) {
            output_file = std::string("results/") + argv[++i];
        }         
        else {
            std::cerr << "Error: Unknown or incomplete argument '" << arg << "'" << std::endl;
            std::cerr << "Usage: " << argv[0] << " -p <pdb> -e <energy> -c <clusters> -o <out>" << std::endl;
            return 1;
        }
    }

    if ((input_protein_file.empty() || input_energy_file.empty() || input_cluster_file.empty() || output_file.empty())) {
        std::cerr << "Error: Missing required arguments" << std::endl;
        std::cerr << "Usage: " << argv[0] << " -p <pdb> -e <energy> -c <clusters> -o <out>" << std::endl;
        return 1;
    }

    std::vector<int> input_clusters; 
    bool is_exclude_mode = false;   

    while (true) {
        std::cout << "\nEnter clusters to process.\n"
                << "Positive numbers include only those clusters.\n"
                << "Negative numbers exclude those clusters.\n"
                << "Input: ";
                
        std::string input;
        std::getline(std::cin, input);

        if (input.empty()) {
            continue;
        }

        std::stringstream ss(input);
        int num;
        std::vector<int> parsed_nums;
        bool has_positive = false;
        bool has_negative = false;

        // Read integers
        while (ss >> num) {
            parsed_nums.push_back(num);
            if (num > 0) has_positive = true;
            if (num < 0) has_negative = true;
        }

        // Validation 1: Mixed inputs
        if (has_positive && has_negative) {
            std::cout << "Error: Cannot mix positive and negative cluster numbers. Please try again.\n";
            continue;
        }
        
        // Validation 2: No valid numbers
        if (parsed_nums.empty()) {
            std::cout << "Invalid input. Please enter valid numbers.\n";
            continue;
        }

        // Save the rules!
        is_exclude_mode = has_negative;
        for (int parsed_num : parsed_nums) {
            // Store the absolute value so we just have a clean list of target IDs
            input_clusters.push_back(std::abs(parsed_num)); 
        }
        
        break; // Successfully got the rules, exit the loop
    }


    std::cout << "--- Files ---" << std::endl;
    std::cout << "Input PDB:  " << input_protein_file << std::endl;
    std::cout << "Input energy file: " << input_energy_file << std::endl;
    std::cout << "Input cluster file: " << input_cluster_file << std::endl;
    std::cout << "Output:     " << output_file << std::endl;

    std::string response;
    while (true) {
        std::cout << "\nProceed? (y/n): ";
        std::cin >> response;
        if (response == "y" || response == "Y") {
            break;
        } else if (response == "n" || response == "N") {
            std::cout << "Job cancelled by user." << std::endl;
            return 0;
        } else {
            std::cout << "Invalid input. Please enter 'y' or 'n'." << std::endl;
        }
    }


    //------------------------------------------------------
    //calculate energy for all gridpoints (without water-water interaction)
    // if(linux == true) {
    // //start by reformatting the input protein file
    //     std::cout << "-> Entering reformatPDB" << std::endl;
    // // int ret1 = std::system("$DOWSER/bin/linux/reformatPDB -pdbin 8OM1_no_ligands.pdb -pdbout reformat.pdb");
    // //create placewat.dow w/ energy
    //     std::cout << "-> Entering placeWat" << std::endl;
    // // int ret2 = std::system("$DOWSER/bin/linux/placeWat reformat.pdb clean_gridpoints.pdb both >placewat.dow");
    // }

    //--------------------------------------------------------

    std::cout << "-> Entering pdbtovector_Waters" << std::endl;

    //create vector of waters from energy file

    std::tuple<std::vector<Water>, double, double, double, double, double, double> energy_tuple = pdbtovector_Waters(input_energy_file);
    std::vector<Water> watervector_energy = std::get<0>(energy_tuple);

    std::tuple<std::vector<Water>, double, double, double, double, double, double> cluster_tuple = pdbtovector_Waters(input_cluster_file);
    std::vector<Water> watervector_cluster = std::get<0>(cluster_tuple);

    if(watervector_cluster.size() != watervector_energy.size()) {
        std::cout << "Error: Cluster file and Energy file not the same length" << std::endl;
        std::cout << "Waters in cluster file: " << watervector_cluster.size() << std::endl;
        std::cout << "Waters in energy file: " << watervector_energy.size() << std::endl;
        return 0;
    }

    std::cout << "* total waters = " << watervector_energy.size() << std::endl;

        
    std::cout << "-> Entering buildSpatialGrid" << std::endl;

    //create hashmap based on distance - for use in overlaps
    std::unordered_map<GridKey, std::vector<int>> distance_map = buildSpatialGrid(watervector_energy, hash_spacing);
    //create hashmap based on clusters - to loop through
    std::unordered_map<int, std::vector<int>> cluster_map = buildClusterMap(watervector_cluster);

    //--------------------------------------------------------
    std::vector<int> clusters;

    for (const auto& [cluster_id, indices] : cluster_map) {
        
        // Check if this specific cluster_id is in the user's list
        bool is_in_user_list = (std::find(input_clusters.begin(), input_clusters.end(), cluster_id) != input_clusters.end());
        
        if (is_exclude_mode) {
            // If they used negatives, process it ONLY if it is NOT in their list
            if (!is_in_user_list) {
                clusters.push_back(cluster_id);
            }
        } else {
            // If they used positives, process it ONLY if it IS in their list
            if (is_in_user_list) {
                clusters.push_back(cluster_id);
            }
        }
    }
    
    for (int cluster_id : clusters) {
        if (cluster_map.find(cluster_id) != cluster_map.end()) {        
            std::vector<int> current_cluster = {};
            int check = 0;

            std::cout << "-> Entering cluster " << cluster_id << std::endl;
            for (int index : cluster_map[cluster_id]) {
                check++;
                current_cluster.push_back(index);  
            }
            std::cout<< "there are " << check << " waters" << std::endl;

            //Montecarlo goes here

        }
    }
    std::cout << "-> Entering loop" << std::endl;

    //from here we are randomizing
    // for(int z = 0; z < total_reps; z++) {

    //     std::cout << "\rIteration: " << z+1 << " of " << total_reps << std::flush;


        
    //     std::tuple<std::vector<Water>,size_t,size_t> newtuple = randomize_states(watervector);
    //     N = std::get<1>(newtuple);
    //     k = std::get<2>(newtuple);
    //     // std::cout << "Original waters: " << N << "\n";
    //     // std::cout << "Remaining waters: " << k << "\n";
    //     // printing all water molecules
    //     // for (int i = 0; i < watervector.size(); i++){
    //     //     std::string water = watervector[i].toString();
    //     //     std::cout << water;
    //     // }

    //     //find overlaps in selected waters

    //     //cleanup last iteration
    //     for(int i = 0; i < watervector.size(); i++){
    //         watervector[i].setOverlap(false);
    //         watervector[i].reset_neighbors();
    //     }

    //     for(int i = 0; i < watervector.size(); i++) {
    //         if(watervector[i].getValue() == 1) {
    //             if(watervector[i].getOverlap() == false){
    //                 getOverlap_cluster(map, watervector, i, hash_spacing);
    //             }
    //         }
    //     }

    //     //calculate energy of configuration
    //     total_energy = 0;
    //     // std::cout << total_energy << " ";
    //     for (int i = 0; i< watervector.size(); i++){
    //         if (watervector[i].getValue() == 1) {
    //             if(watervector[i].getOverlap()){
    //                 total_energy += 10;
    //             } else {
    //                 total_energy += watervector[i].getEnergy();
    //                 total_energy += watervector[i].constructive_interaction();
    //                 if(watervector[i].constructive_interaction()!=0){
    //                     std::cout << "UH OH" << std::endl;
    //                 }
    //             }
    //             // total_energy += watervector[i].getEnergy();

    //         }
    //     }
    //             // std::cout << total_energy << "\n ";

    //     //metropolis
    //     auto[temp, accepted] = metropolis(total_energy, lowest_energy);
    //     lowest_energy = temp;
        
    //     if(accepted) {
    //         lowest_config = std::get<0>(newtuple);
    //     }

    //     monte_carlo_output << z << ", " << N << ", " << k << ", " << total_energy << ", " << lowest_energy << "\n";
    
    // }
    // std::cout << std::endl;
    // std::cout << "Completed " << total_reps << " iterations" << std::endl;

    // vectortopdb(lowest_config,output_file,N,k,false);

    // monte_carlo_output.close();

    return 1;
}

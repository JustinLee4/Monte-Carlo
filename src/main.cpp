#include "atom.h"
#include "pdbtovector.h"
#include "montecarlo.h"
#include "map.h"

#include <random>


#include <iostream>
#include <cstdio>


bool linux = true;
int total_reps = pow(10, 6);

int main(int argc, char* argv[]){

    double lowest_energy = INFINITY;
    double total_energy = 0;
    std::vector<std::tuple<int, std::vector<Water>>> lowest_config_cluster;
    bool accepted;
    std::ofstream monte_carlo_log_output;
    double hash_spacing = 5.0;
    std::random_device rd;
    std::mt19937 gen(rd());

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
        // else if ((arg == "-c" || arg == "--clusters") && i + 1 < argc) {
        //     input_cluster_file = argv[++i];
        // }
        else if ((arg == "-o" || arg == "--out") && i + 1 < argc) {
            output_file = std::string("results/") + argv[++i];
        }
        else if ((arg == "-reps" ) && i + 1 < argc) {
            total_reps = std::stod(argv[++i]);    
        }     
        else {
            std::cerr << "Error: Unknown or incomplete argument '" << arg << "'" << std::endl;
            std::cerr << "Usage: " << argv[0] << " -p <pdb> -e <energy> -c <clusters> -o <out> -reps <total reps>" << std::endl;
            return 1;
        }
    }

    if ((input_protein_file.empty() || input_energy_file.empty() || /*input_cluster_file.empty() ||*/ output_file.empty())) {
        std::cerr << "Error: Missing required arguments" << std::endl;
        std::cerr << "Usage: " << argv[0] << " -p <pdb> -e <energy> -c <clusters> -o <out> -reps <total reps>" << std::endl;
        return 1;
    }

    std::vector<int> input_clusters; 
    bool is_exclude_mode = false;   

    while (true) {
        std::cout << "\nEnter clusters to process. (space separated) \n"
                << "Positive numbers include only those clusters.\n"
                << "Negative numbers exclude those clusters.\n"
                << "Input: ";
                
        std::string input;
        std::getline(std::cin, input);

        if (input.empty()) {
            is_exclude_mode = true;
            std::cout << "all clusters" << std::endl;
            break;
        }

        std::stringstream ss(input);
        int num;
        std::vector<int> parsed_nums;
        bool has_positive = false;
        bool has_negative = false;

        while (ss >> num) {
            parsed_nums.push_back(num);
            if (num > 0) has_positive = true;
            if (num < 0) has_negative = true;
        }

        if (has_positive && has_negative) {
            std::cout << "Error: Cannot mix positive and negative cluster numbers. Please try again.\n";
            continue;
        }
        
        if (parsed_nums.empty()) {
            std::cout << "Invalid input. Please enter valid numbers.\n";
            continue;
        }

        is_exclude_mode = has_negative;
        for (int parsed_num : parsed_nums) {
            input_clusters.push_back(std::abs(parsed_num)); 
        }
        
        break;
    }


    std::cout << "--- Files ---" << std::endl;
    std::cout << "Input PDB:  " << input_protein_file << std::endl;
    std::cout << "Input energy file: " << input_energy_file << std::endl;
    std::cout << "Input cluster file: " << input_cluster_file << std::endl;
    std::cout << "Output:     " << output_file << std::endl;
    std::cout << std::scientific;
    std::cout <<  "Total steps:  " << static_cast<double>(total_reps) << std::endl;
    std::cout << std::defaultfloat;

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

    auto start_time = std::chrono::high_resolution_clock::now();

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

    // std::tuple<std::vector<Water>, double, double, double, double, double, double> cluster_tuple = pdbtovector_Waters(input_cluster_file);
    // std::vector<Water> watervector_cluster = std::get<0>(cluster_tuple);

    // if(watervector_cluster.size() != watervector_energy.size()) {
    //     std::cout << "Error: Cluster file and Energy file not the same length" << std::endl;
    //     std::cout << "Waters in cluster file: " << watervector_cluster.size() << std::endl;
    //     std::cout << "Waters in energy file: " << watervector_energy.size() << std::endl;
    //     return 0;
    // }

    std::cout << "* total waters = " << watervector_energy.size() << std::endl;

        
    std::cout << "-> Entering buildSpatialGrid" << std::endl;

    //create hashmap based on distance - for use in overlaps
    std::unordered_map<GridKey, std::vector<int>> distance_map = buildSpatialGrid(watervector_energy, hash_spacing);
    //create hashmap based on clusters - to loop through
    std::unordered_map<int, std::vector<int>> cluster_map = buildClusterMap(watervector_energy);


    //use overlap hashmap to populate neighbors for all gridpoints
    for (Water& water : watervector_energy) {
        getOverlap_cluster(distance_map, watervector_energy, water, hash_spacing, 2.5);
    }




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
    std::sort(clusters.begin(), clusters.end());

    std::filesystem::path root_path(output_file);

    std::filesystem::path dir_path = root_path.parent_path();

    if (!dir_path.empty() && !std::filesystem::exists(dir_path)) {
        try {
            std::filesystem::create_directories(dir_path);
        } catch (const std::filesystem::filesystem_error& e) {
            std::cerr << "Error creating directories: " << e.what() << '\n';
            return 0; 
        }
    }
    
    for (int cluster_id : clusters) {
        if (cluster_map.find(cluster_id) != cluster_map.end()) { 
            monte_carlo_log_output.open(output_file + "_cluster_" + std::to_string(cluster_id) + ".log");
            lowest_energy = INFINITY;
            double current_energy = INFINITY;       
            std::vector<int> current_cluster = {};
            int check = 0;


            std::cout << "-> Entering cluster " << cluster_id << std::endl;
            for (int index : cluster_map[cluster_id]) {
                check++;
                current_cluster.push_back(index);  
            }
            std::cout<< "* there are " << check << " waters in this cluster" << std::endl;
            std::vector<Water> temp_lowest_vector(check, Water({0.0, 0.0, 0.0}));

            ActiveList active_list;
            active_list.reserve(current_cluster.size());

            for (auto& water : watervector_energy) {
                water.set_value(false); 
                water.clear_Overlap();
            }

            for (int index : current_cluster) {
                active_list.add(index);
            }
            int rep_count = 0;

            while(rep_count < total_reps){
                //montecarlo randomization
                int changed_index = iterate_singly(watervector_energy, active_list, gen);

                rep_count++;

                //calculate energy of configuration
                total_energy = 0;
                // std::cout << total_energy << " ";
                for (int j : current_cluster) {
                    if (watervector_energy[j].get_value() == 1) {

                            total_energy += watervector_energy[j].get_bfactor();
                            // total_energy += watervector_energy[i].constructive_interaction();
                            // if(watervector[i].constructive_interaction()!=0){
                            //     std::cout << "UH OH" << std::endl;
                            // }
                    }
                }


                //metropolis
                auto[temp, accepted] = metropolis(total_energy, current_energy);
                
                if(accepted) {
                    current_energy = total_energy;
                } else {
                    bool current_state = watervector_energy[changed_index].get_value();
                    
                    if (current_state == true) {
                        watervector_energy[changed_index].set_value(false);
                        watervector_energy[changed_index].subtract_overlap_with_neighbors(watervector_energy, active_list);
                    } else {
                        watervector_energy[changed_index].set_value(true);
                        watervector_energy[changed_index].add_overlap_with_neighbors(watervector_energy, active_list);
                    }
                }
                if (current_energy < lowest_energy) {
                    lowest_energy = current_energy;
                    std::copy(watervector_energy.begin() + current_cluster[0], watervector_energy.begin() + current_cluster[0] + current_cluster.size(), temp_lowest_vector.begin());

                }

                monte_carlo_log_output << rep_count << ", " << total_energy << ", " << current_energy << ", " << lowest_energy << "\n";

                if (rep_count % 1000 == 0 || rep_count == total_reps) {
                    double percent = ((double)(rep_count) / total_reps) * 100.0;
                    std::cout << "\r\033[K  Progress:  " << std::fixed << std::setprecision(1) << percent << "%" << std::flush;                
                }
            }
            std::cout << std::endl;

            lowest_config_cluster.push_back({cluster_id, temp_lowest_vector});
            monte_carlo_log_output.close();


            std::cout << "energy of lowest configuration for cluster " << cluster_id << " is " << lowest_energy << std::endl;
        }
        
    }
    
    std::string final_output_name = output_file + "_all_clusters.pdb";
    std::ofstream combined_pdb_file(final_output_name);

    if (!combined_pdb_file.is_open()) {
        std::cerr << "Error: Could not open " << final_output_name << " for writing!\n";
        return 0; 
    }

    for(int i = 0; i < lowest_config_cluster.size(); i++) {
        
        int current_cluster_id = std::get<0>(lowest_config_cluster[i]);
        auto current_vector = std::get<1>(lowest_config_cluster[i]);
        
        vectortopdb(current_vector, combined_pdb_file, current_cluster_id);
    }

    combined_pdb_file.close();

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

        //--------- timer ----------
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end_time - start_time;

    long long total_seconds = static_cast<long long>(elapsed.count());
    long long hours = total_seconds / 3600;
    long long minutes = (total_seconds % 3600) / 60;
    long long seconds = total_seconds % 60;

    std::cout << "Time: " 
          << std::setfill('0') << std::setw(2) << hours << ":"
          << std::setfill('0') << std::setw(2) << minutes << ":"
          << std::setfill('0') << std::setw(2) << seconds 
          << std::endl;


    return 1;
}

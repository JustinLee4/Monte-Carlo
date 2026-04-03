#ifndef PDB_TO_VECTOR
#define PDB_TO_VECTOR

#include "atom.h"
#include "AtomicRadii.h"
#include "common.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>


std::array<double, 3> get_coords(const std::string& input);


std::tuple<std::string, std::string, double> get_data(const std::string& input);


//turn pdbfile to a vector of atoms
std::tuple<std::vector<Atom>, double, double, double, double, double, double> pdbtovector(std::string filename);

//turn pdbfile to a vector of Waters
std::tuple<std::vector<Water>, double, double, double, double, double, double> pdbtovector_Waters(std::string filename);


//---------------------------------------------------
//converting vector back into pdb file

// void vectortopdb(const std::vector<Atom> &atomvector, std::string output_filename);
template <typename T>
void vectortopdb(const std::vector<T> &atomvector, 
                 std::ofstream &out_file,
                 int cluster_number) {

    


    for(int i = 0; i < atomvector.size(); i++) {
        if(atomvector[i].get_value() == false) {
            continue;
        }

        std::array<double, 3> pos = atomvector[i].getCoords();
        int z = i + 1;
        if( z > 9999) {
            z = 9999;
        }

        out_file << std::fixed;
        out_file << std::setprecision(3);

        out_file << "HETATM"
         << std::setw(5) << std::right << z        // Col 7-11
         << " "                                        // Col 12
         << " O  "                                     // Col 13-16 (Atom Name)
         << " "                                        // Col 17
         << "HOH"                                      // Col 18-20 (ResName)
         << " "                                        // Col 21
         << "A"                                        // Col 22 (Chain)
         << std::setw(4) << std::right << cluster_number        // Col 23-26 (ResSeq)
         << "    "                                     // Col 27-30
         << std::setw(8) << std::fixed << std::right << pos[0] // X
         << std::setw(8) << std::fixed << std::right << pos[1] // Y
         << std::setw(8) << std::fixed << std::right << pos[2] // Z
         << std::setw(6) << "1.00"                    // Occ
         << std::fixed << std::setprecision(2) << std::setw(6) << atomvector[i].get_bfactor() // Temp
         << "          "                               // Spacing
         << " O" << "\n";                               // Element
    }   
    
}


//-----------------------------
//take a.pdb, append b.pdb, result c.pdb
bool append_pdb_files(const std::string& filepath_1, const std::string& filepath_2, const std::string& output_path, int targetLine);

#endif
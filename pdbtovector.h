#ifndef PDB_TO_VECTOR
#define PDB_TO_VECTOR

#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <sstream>
#include <vector>
#include <iterator>
#include <iomanip>

#include "water.h"

//helper functions for pdbtovector_Waters
//check for Water && check position of oxygen within a line
std::string find_water(std::string line, int line_case);
//get coords for Water class from HOH lines (HOH lines only)
std::array<double, 3> get_coords(std::string input);
//get energy for Water from HOH lines
double get_energy(std::string input);


//turn pdbfile to a vector of waters
std::vector<Water> pdbtovector_Waters(std::string filename);

void test();

//---------------------------------------------------
//converting vector back into pdb file

void vectortopdb(const std::vector<Water> &watervector, std::string output_filename, size_t N, size_t k, bool withEnergy_and_Overlaps);

//-----------------------------
//take a.pdb, append b.pdb, result c.pdb
bool append_pdb_files(const std::string& filepath_1, const std::string& filepath_2, const std::string& output_path, int targetLine);

#endif
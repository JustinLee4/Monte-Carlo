// #include <iostream>
// #include <fstream>
// #include <string>
// #include <filesystem>
// #include <sstream>
// #include <vector>
// #include <iterator>
#include "pdbtovector.h"

//helper functions for pdbtovector_Waters
//check for Water
std::string find_water(std::string line, int line_case)
{
    // std::cout << "entering reformat";
    std::vector<std::string> vec;
    //HETATM
    if(line_case == 1) {
        // std::cout << "entering line_case 1 \n";
        
        //choose line after HETATM
        std::string line1;
        line1 = line.substr(6);

        //partition line into vector by whitespaces
        std::string s;
        std::stringstream ss(line1);
        while(getline(ss, s, ' ')){
            std::size_t prev = 0, pos;
            while ((pos = s.find_first_of("-", prev)) != std::string::npos) {
                // std::cout << "pos = " << pos << std::endl;
                if(pos > prev){    
                    vec.push_back(s.substr(prev,pos-prev));

                }
                prev = pos+1;
                // std::cout << "minus found" << std::endl;
            }
            if (prev < s.length()) {
                if (prev == 0){
                    vec.push_back(s.substr(prev));
                } else {
                    vec.push_back(s.substr(prev-1));
                }
            }
        }

        ss.str(std::string());
        ss.clear();

        // for (std::size_t i = 0; i < vec.size(); i++){
        //         std::cout << vec[i] + " " ;
        // }
        // std::cout << "\n";

        //check vector for O & H
        if((vec[1] == "OW" || vec[1] == "O") && vec[2] == "HOH") {
            std::string output = "HETATM ";
            // std::cout << "found O or OW" << std::endl;
            for (std::size_t i = 0; i < vec.size(); i++){
                output = output + vec[i] + " ";
            }
            // std::cout << output << std::endl;
            return output;
        } else {
            // std::cout << "not O or OW" << std::endl;
            return "";
        }

    }
    //ATOM
    if(line_case == 2) {
        // std::cout << "entering line_case 2\n";
        //choose line after ATOM
        std::string line1;
        line1 = line.substr(4);
        // std::cout << line1 << "\n";

        //parse line into vector
        std::string s;
        std::stringstream ss(line1);
        while(getline(ss, s, ' ')){

            // if (s.length() != 0) {
            //     std::cout << s << std::endl;
            // }
            std::size_t prev = 0, pos;
            while ((pos = s.find_first_of("-", prev)) != std::string::npos) {
                // std::cout << "pos = " << pos << std::endl;
                if(pos > prev){    
                    vec.push_back(s.substr(prev,pos-prev));
                }
                prev = pos+1;
                // std::cout << "minus found" << std::endl;
            }
            if (prev < s.length()) {
                if (prev == 0){
                    vec.push_back(s.substr(prev));
                } else {
                    vec.push_back(s.substr(prev-1));
                }
            }
        }
        ss.str(std::string());
        ss.clear();

        // for (size_t i = 0; i < vec.size(); i++) {
        //     std::cout << "vec " << i << " = " << vec[i] << " ";
        // }
        // std::cout << "\n";
        

        //check vector for O & H
        if((vec[1] == "OW" || vec[1] == "O") && vec[2] == "HOH") {
            std::string output = "ATOM ";
            for (std::size_t i = 0; i < vec.size(); i++){
                output = output + vec[i] + " ";
            }
            return output;
        } else {
            return "";
        }
    }
    return "";
}

//get coords for Water class from HOH lines (HOH lines only)
std::array<double, 3> get_coords(std::string input) {
    std::array<double,3> output;
    std::string s;
    std::stringstream ss(input);
   
    // 1. Create the "start" iterator, initialized with the stringstream
    std::istream_iterator<std::string> start(ss);
    
    // 2. Create the "end" (default) iterator
    std::istream_iterator<std::string> end;
    
    // 3. Create the vector. This is now completely unambiguous.
    std::vector<std::string> tokens(start, end);
    // for(int i = 0; i < tokens.size(); i++){
    //     std::cout << tokens[i] << ", i= " << i << std::endl;
    // }
    for(int i = 0; i < 3; i++){
        output[i] = std::stod(tokens[i + 5]);
    }
    return output;
}

//get energy for Water class from HOH lines (HOH lines only)
double get_energy(std::string input) {
    
    double output;
    std::string s;
    std::stringstream ss(input);
   
    // 1. Create the "start" iterator, initialized with the stringstream
    std::istream_iterator<std::string> start(ss);
    
    // 2. Create the "end" (default) iterator
    std::istream_iterator<std::string> end;
    
    // 3. Create the vector. This is now completely unambiguous.
    std::vector<std::string> tokens(start, end);
    // for(int i = 0; i < tokens.size(); i++){
    //     std::cout << tokens[i] << ", i= " << i << std::endl;
    // }
   std::string last_token = tokens.back();
    // std::cout << last_token << "\n";
    return std::stod(last_token);
}

std::vector<Water> pdbtovector_Waters(std::string filename) {
    std::vector<Water> output;
    std::ifstream infile(filename);
    // std::ofstream out_file;
    // out_file.open("output.pdb");
    // out_file << "entered " << filename << std::endl;
    std::string line_string;
    int i = 0;
    while (std::getline(infile, line_string))
    {   
        i = i+1;
        // std::cout << i << std::endl;
        std::string answer;
        //Choosing HETATM and ATOM rows
        if(line_string.substr(0,6) == "HETATM") {
            //printing HETATM rows
            answer = find_water(line_string, 1);
            
            // std::cout << "ans = " << answer << std::endl;
            if(!answer.empty()) {
                // out_file << answer << "\n";
                std::array<double, 3> this_array = get_coords(answer);
                double this_energy = get_energy(answer);
                Water this_Water = Water(0, this_energy, this_array);
                output.push_back(this_Water);
            }  
        }
        if(line_string.substr(0,4) == "ATOM"){
            //Print ATOM rows
            answer = find_water(line_string, 2);

            if(!(answer.empty())) {
                // out_file << answer << "\n";
                std::array<double, 3> this_array = get_coords(answer);
                double this_energy = get_energy(answer);
                Water this_Water = Water(0, this_energy, this_array);
                output.push_back(this_Water);
            }
        }        
    }
    // out_file.close();
    return output;
}

void test() {
    ;
}

//-------------------------------------------

void vectortopdb(const std::vector<Water> &watervector, std::string output_filename, size_t N, size_t k, bool withEnergy_and_Overlaps) {
    std::vector<Water> output;
    std::ofstream out_file;
    out_file.open(output_filename);
    out_file << std::fixed;
    out_file << std::setprecision(3);

    // out_file << "REMARK  total number of initial water molecules = " << N << "\n"
            //  << "REMARK  total remaining water molecules = " << k << "\n";

    for(int i = 0; i < watervector.size(); i++) {

        bool value = watervector[i].getValue();
        if(value) {

            std::array<double, 3> pos = watervector[i].getCoords();
            
            if (withEnergy_and_Overlaps) {
                double energy = watervector[i].getEnergy();
                bool hasOverlap = watervector[i].getOverlap();

                out_file << "HETATM"                               // [ 1- 6] Record name
                        << std::setw(5) << std::right << i+1 // [ 7-11] Atom serial number
                        << " "
                        << std::setw(4) << std::left << "  O"      // [13-16] Atom name (left-justified)
                        << " "                                    // [17]    Alternate location
                        << std::setw(3) << std::right << "HOH"     // [18-20] Residue name
                        << " "                                    // [21]    (blank)
                        << " " //  << "A"                                     // [22]    Chain ID (using 'A' as default)
                        << std::setw(4) << std::right << i + 1     // [23-26] Residue seq num (assuming one res per obj)
                        << "    "                                  // [27-30] (insertion code + 3 blanks)
                        << std::setw(8) << std::right << pos[0]    // [31-38] X coordinate
                        << std::setw(8) << std::right << pos[1]    // [39-46] Y coordinate
                        << std::setw(8) << std::right << pos[2]    // [47-54] Z coordinate
                        << std::setw(6) << std::right << "1.00"    // [55-60] Occupancy (default 1.00)
                        //  << std::setw(6) << std::right << "0.00"    // [61-66] B-factor (default 0.00)
                        //  << "          "                             // [67-76] (blanks)
                        //  << std::setw(2) << std::right << "O"       // [77-78] Element symbol
                        //  << "  "                                   // [79-80] Charge
                        << " " << "Energy = " << energy
                        << " " << "Is Overlapping = " << hasOverlap
                        << "\n";                                  // End of line
            } else if (!withEnergy_and_Overlaps)
            {
                out_file << "HETATM"                               // [ 1- 6] Record name
                        << std::setw(5) << std::right << i+1 // [ 7-11] Atom serial number
                        << " "
                        << std::setw(4) << std::left << "  O"      // [13-16] Atom name (left-justified)
                        << " "                                    // [17]    Alternate location
                        << std::setw(3) << std::right << "HOH"     // [18-20] Residue name
                        << " "                                    // [21]    (blank)
                        << " " //  << "A"                                     // [22]    Chain ID (using 'A' as default)
                        << std::setw(4) << std::right << i + 1     // [23-26] Residue seq num (assuming one res per obj)
                        << "    "                                  // [27-30] (insertion code + 3 blanks)
                        << std::setw(8) << std::right << pos[0]    // [31-38] X coordinate
                        << std::setw(8) << std::right << pos[1]    // [39-46] Y coordinate
                        << std::setw(8) << std::right << pos[2]    // [47-54] Z coordinate
                        << std::setw(6) << std::right << "1.00"    // [55-60] Occupancy (default 1.00)
                        //  << std::setw(6) << std::right << "0.00"    // [61-66] B-factor (default 0.00)
                        //  << "          "                             // [67-76] (blanks)
                        //  << std::setw(2) << std::right << "O"       // [77-78] Element symbol
                        //  << "  "                                   // [79-80] Charge
                        << "\n";                                  // End of line
            }
            
        }
    }
    out_file.close();
}


bool append_pdb_files(const std::string& filepath_1, const std::string& filepath_2, const std::string& output_path, int targetLine) {
    std::ifstream file1(filepath_1);
    std::ifstream file2(filepath_2);
    std::ofstream outFile(output_path);

    // Check if files opened successfully
    if (!file1.is_open()) {
        std::cerr << "Error: Could not open input file: " << filepath_1 << std::endl;
        return false;
    }
    if (!file2.is_open()) {
        std::cerr << "Error: Could not open input file: " << filepath_2 << std::endl;
        return false;
    }
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not create output file: " << output_path << std::endl;
        return false;
    }

    std::string line;
    int currentLine = 0;

    // 1. Read File 1 up to the target line
    // We use a loop to read line by line so we can count them
    while (currentLine < targetLine && std::getline(file1, line)) {
        outFile << line << "\n";
        currentLine++;
    }

    // 2. Inject the entirety of File 2
    // rdbuf() is a highly efficient way to copy the entire buffer of a file
    outFile << file2.rdbuf();

    // Ensure File 2 ends with a newline before resuming File 1 (optional but recommended)
    // Un-comment the next line if File 2 might be missing a terminal newline
    // outFile << "\n"; 

    // 3. Write the remainder of File 1
    // rdbuf() will dump the rest of the file starting from where the previous loop left off
    outFile << file1.rdbuf();

    // Close streams (optional, destructors will handle this automatically)
    file1.close();
    file2.close();
    outFile.close();

    return true;
}

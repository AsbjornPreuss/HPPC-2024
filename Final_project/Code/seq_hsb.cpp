#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <fstream>
#include <cstdlib>


class spin_system {
    public:
    int flips = 100; // Number of flips the system will simulate.
    int n_spins = 25; // Number of spins in each direction.
    int n_dims = 2; // Number of dimensions the spins are placed in.
    int x_dist = 10;
    int y_dist = 10;
    int z_dist = 10;
    int nearest_neighbours = 1; // Number of nearest neighbour interactions to be calculated
    double J = 1; // Magnetization parameter, used in calculating energy.
    double H; // Total energy of the system;
    std::string filename = "seq_out.txt"; // Output file.
    std::vector<std::vector<double>> position; // Three by n_spins matrix, defining the spin's 3d position.
    std::vector<std::vector<double>> spin; // Three by n_spins matrix, defining the spin vector for each spin.

    spin_system(std::vector <std::string> argument){
    for (long unsigned int i = 1; i<argument.size() ; i += 2){
            std::string arg = argument.at(i);
            if(arg=="-h"){ // Write help
                std::cout << "Heisenberg_simulation -flips <number of flips performed> -n_spins <number of spins simulated>"
                          << " -ofile <filename> \n";
                exit(0);
                break;
            } else if(arg=="-flips"){
                flips = std::stoi(argument[i+1]);
            } else if(arg=="-n_spins"){
                n_spins = std::stoi(argument[i+1]);
            } else if(arg=="-ofile"){
                filename = argument[i+1];
            } else{
                std::cout << "---> error: the argument type is not recognized \n";
            }
        }
    }
};


// Function that generates rectangular positions for alle the spins in the system, 
void generate_positions_box(spin_system &sys){
    switch (sys.n_dims)
    {
    case 1:
        for (double i=0; i<sys.n_spins; i++)
            sys.position.push_back({i*sys.x_dist, 0, 0});
        break;
    case 2:
        for (double i=0; i<sqrt(sys.n_spins); i++)
        for (double j=0; j<sqrt(sys.n_spins); j++)
            sys.position.push_back({i*sys.x_dist, j*sys.y_dist, 0});
        break;
    case 3:
        for (double i=0; i<cbrt(sys.n_spins); i++)
        for (double j=0; j<cbrt(sys.n_spins); j++)
        for (double k=0; k<cbrt(sys.n_spins); k++)
            sys.position.push_back({i*sys.x_dist, j*sys.y_dist, k*sys.z_dist});
        break;
    default:
        break;
    }

};

// Function that generates random directions for all the spins in the system
void generate_spin_directions(spin_system &sys){
    
    for (int i = 0; i<sys.n_spins; i++){
        srand(i); // Seed is here to make it perform nicely when comparing to parallel
        double spin_azimuthal = (double) rand()/RAND_MAX;
        srand(i*rand()); // Seed is here to make it perform nicely when comparing to parallel
        double spin_polar = (double) rand()/RAND_MAX;
        
        sys.spin.push_back({sin(spin_azimuthal)*cos(spin_polar), 
                            sin(spin_azimuthal)*sin(spin_polar),
                            cos(spin_polar)});
    }   
};

// Function that calculates the energy of a single spin in 2d
double energy_calculation_2d(spin_system &sys, int spin){
    // Make list of nearest neighbour spin indexes
    // TODO: Make general for N-nearest neighbours
    int n_spins_row = sqrt(sys.n_spins);
    int spin_row = spin/n_spins_row; // Which row the spin is in
    int spin_col = spin%n_spins_row; // Which column the spin is in
    int spin_interactions[4] = {(spin_row - 1)%n_spins_row*n_spins_row + (spin_col)%n_spins_row,
                                 (spin_row)%n_spins_row*n_spins_row + (spin_col - 1)%n_spins_row, 
                                 (spin_row )%n_spins_row*n_spins_row + (spin_col + 1)%n_spins_row,
                                 (spin_row + 1)%n_spins_row*n_spins_row + (spin_col)%n_spins_row};
    double energy = 0;
    double dot_product;
    for (int i=0; i<4; i++){
        // Calculate the energy with the nearest neighbour with no corners
        if (spin_interactions[i]<0) spin_interactions[i] *= -1; // Should'nt be necessary, but modulo is not good with - 
        dot_product = sys.spin[spin][0]*sys.spin[spin_interactions[i]][0] 
                        + sys.spin[spin][1]* sys.spin[spin_interactions[i]][1]
                        + sys.spin[spin][2]* sys.spin[spin_interactions[i]][2];
        energy += - sys.J/2*dot_product;
    }
    return energy;
};

// Calculate the total energy of the system
void Calculate_h(spin_system& sys){
    for (int i=0; i<sys.n_spins; i++){
        sys.H += energy_calculation_2d(sys, i)*0.5; // Half the energy, because we calculate on all the spins
    }
};

// Write the spin configurations in the output file.
void Writeoutput(spin_system& sys, std::ofstream& file){  
    // Loop over all spins, and write out position and spin direction
    for (int i = 0; i<sys.n_spins; i++){
        file << sys.position[i][0] << " " << sys.position[i][1] << " "  << sys.position[i][2] << " "
            << sys.spin[i][0] << " " << sys.spin[i][1] << " "  << sys.spin[i][2] << " "
            << sys.H
            << '\n';
    }
};

void Simulate(spin_system& sys){
    for (int iteration=0; iteration<sys.flips; iteration++){
        // Choose a random spin site

        // Calculate if it lowers energy

        // If not, see if it should be randomised in direction

        // If yes, randomise

        // Change H to represent the total energy of the system

        
    }
}
//=============================================================================================
//=========================   MAIN FUNCTION   =================================================
//=============================================================================================
int main(int argc, char* argv[]){
    std::cout << "Hello Heisenberg!" << std::endl;
    spin_system sys({argv, argv+argc});
    generate_positions_box(sys);
    generate_spin_directions(sys);
    Calculate_h(sys);
    Simulate(sys);
    std::ofstream file(sys.filename); // open file
    Writeoutput(sys, file);



}
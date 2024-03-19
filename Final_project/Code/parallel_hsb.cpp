#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <mpi.h>

int mpi_size;
int mpi_rank;
int nproc_x = 2, nproc_y = 1,nproc_z=1;

bool verbose = true;
class spin_system {
    public:
    int flips = 100; // Number of flips the system will simulate.
    int n_spins = 25; // Number of spins in The system.
    int n_dims = 3; // Number of dimensions the spins are placed in.
    int n_spins_row; // Number of rows in the square/cube

    int x_offsets[6] = {-1,0,0,1,0,0}; // For calculating neighbor indices in energy calculation
    int y_offsets[6] = {0,-1,1,0,0,0};
    int z_offsets[6] = {0,0,0,0,-1,1};

    int x_dist = 10;
    int y_dist = 10;
    int z_dist = 10;

    int nearest_neighbours = 1; // Number of nearest neighbour interactions to be calculated
    int xlen, ylen, zlen;
    double J = 1; // Magnetization parameter, used in calculating energy.
    double H; // Total energy of the system;
    double B = 0; // Magnetic field in z direction
    double Temperature = 1; // Temperature of the system.
    std::string filename = "seq_out.txt"; // Output file.
    std::vector<std::vector<double>> position; // Three by n_spins matrix, defining the spin's 3d position.
    std::vector<std::vector<double>> spin; // Three by n_spins matrix, defining the spin vector for each spin.
    std::vector<std::vector<int>>    neighbours; // 2*n_dims by n_spins matrix, defining the neighbour indices of each cell, so they need only be calculated once.
    spin_system(std::vector <std::string> argument){
    for (long unsigned int i = 1; i<argument.size() ; i += 2){
            std::string arg = argument.at(i);
            if(arg=="-h"){ // Write help
                std::cout << "Heisenberg_simulation\n --flips <number of flips performed>\n --nspins <number of spins simulated>\n --ndims <number of dimensions to simulate> \n"
                          << " --ofile <filename>\n --magnet <strength of external magnetic field in z direction>\n"
                          << " --temp <temperature> \n";

                exit(0);
                break;
            } else if(arg=="--flips"){
                flips = std::stoi(argument[i+1]);
            } else if(arg=="--nspins"){
                n_spins = std::stoi(argument[i+1]);
            } else if(arg=="--ndims"){
                n_dims = std::stoi(argument[i+1]);
            } else if(arg=="--ofile"){
                filename = argument[i+1];
            } else if(arg=="--magnet"){
                B = std::stoi(argument[i+1]);
            } else if(arg=="--temp"){
                Temperature = std::stod(argument[i+1]);
            } else{
                std::cout << "---> error: the argument type is not recognized \n";
            }
        }
    n_spins_row = pow(float(n_spins),1./float(n_dims)); //Equal size in all dimensions

    xlen = x_dist*n_spins_row;
    ylen = y_dist*n_spins_row;
    zlen = z_dist*n_spins_row;
        
    }
};

class local_spins{
    public:
    int x_offsets[6] = {-1,0,0,1,0,0}; // For calculating neighbor indices in energy calculation
    int y_offsets[6] = {0,-1,1,0,0,0};
    int z_offsets[6] = {0,0,0,0,-1,1};

    int x_dist,y_dist,z_dist;

    int nearest_neighbours; // Number of nearest neighbour interactions to be calculated
    int xlen, ylen, zlen;
    int pad_xlen, pad_ylen, pad_zlen;
    int offset_x, offset_y, offset_z;
    int n_spins;
    double J; // Magnetization parameter, used in calculating energy.
    double H; // Total energy of the system;
    double B; // Magnetic field in z direction
    double Temperature; // Temperature of the system.
    std::string filename; // Output file.

    int no_in_padded_layer, no_in_layer;

    local_spins(spin_system &sys,
            int local_xlen, int local_ylen, int local_zlen,
            int offx, int offy, int offz){
        x_dist = sys.x_dist;
        y_dist = sys.y_dist;
        z_dist = sys.z_dist;
        
        nearest_neighbours = sys.nearest_neighbours;
        xlen = local_xlen;
        ylen = local_ylen;
        zlen = local_zlen;
        
        pad_xlen = xlen+2;
        pad_ylen = ylen+2;
        pad_zlen = zlen+2;
        
        n_spins = xlen*ylen*zlen;
        no_in_layer = xlen*ylen;
        no_in_padded_layer = (xlen+2)*(ylen+2);

        offset_x = offx;
        offset_y = offy;
        offset_z = offz;
        J = sys.J;
        H = sys.H;
        B = sys.B;
        Temperature = sys.Temperature;
        filename = sys.filename;
    };
    std::vector<std::vector<double>> position; // Three by n_spins matrix, defining the spin's 3d position.
    std::vector<std::vector<double>> spin; // Three by n_spins matrix, defining the spin vector for each spin.
    std::vector<std::vector<int>>    neighbours; // 2*n_dims by n_spins matrix, defining the neighbour indices of each cell, so they need only be calculated once.


    int index_to_padded_index(int index){
        int x = index%(no_in_layer)%xlen + 1;
        int y = (index%(no_in_layer))/xlen + 1;
        int z = index/(no_in_layer) + 1;

        return z*no_in_padded_layer + y*pad_xlen + x; 
    }

    void padded_index_to_padded_coordinates(int index, int& x, int& y, int& z){
        x = index % pad_xlen; // Which row the spin is in
        y = (index/pad_ylen)%pad_ylen; // Which column the spin is in
        z = index / no_in_padded_layer;
    }

};


// Function that generates rectangular positions for alle the spins in the system, 
void generate_positions_box(local_spins &sys){
        for (double i=0; i<sys.pad_xlen; i++)
        for (double j=0; j<sys.pad_ylen; j++)
        for (double k=0; k<sys.pad_zlen; k++)
            sys.position.push_back({i*sys.x_dist, j*sys.y_dist, k*sys.z_dist});
};

// Function that generates random directions for all the spins in the system
void generate_spin_directions(local_spins &sys){
    
    for (int i = 0; i<sys.pad_zlen*sys.no_in_padded_layer; i++){
        srand(i); // Seed is here to make it perform nicely when comparing to parallel
        double spin_azimuthal = (double) rand()/RAND_MAX * M_PI;
        srand(i*rand()); // Seed is here to make it perform nicely when comparing to parallel
        double spin_polar = (double) rand()/RAND_MAX * 2. * M_PI;
        
        sys.spin.push_back({sin(spin_azimuthal)*cos(spin_polar), 
                            sin(spin_azimuthal)*sin(spin_polar),
                            cos(spin_azimuthal)});
    }   
};

void generate_neighbours(local_spins &sys){
    
    for (int spin = 0; spin< sys.pad_zlen*sys.no_in_padded_layer; spin++){
        // Find position in square / cube
        int spin_x, spin_y, spin_z;
        sys.padded_index_to_padded_coordinates(spin, spin_x, spin_y, spin_z);

        // Find indices of neighbours
        std::vector<int> spin_interactions;
        for(int i = 0; i < 6; i++){
            spin_interactions.push_back((spin_x + sys.x_offsets[i])%sys.pad_xlen + 
                                        (spin_y + sys.y_offsets[i])%sys.pad_ylen * sys.pad_xlen + 
                                        (spin_z + sys.z_offsets[i])%sys.pad_zlen * sys.no_in_padded_layer);
        }
        sys.neighbours.push_back(spin_interactions);
    }
}
// Function that calculates the energy of a single spin in 2d
double energy_calculation_nd(local_spins &sys, int spin){
    double energy = 0;
    double dot_product;
    for (int i=0; i<6; i++){
        // Calculate the energy with the nearest neighbour with no corners
        dot_product = sys.spin[spin][0]*sys.spin[sys.neighbours[spin][i]][0] 
                        + sys.spin[spin][1]* sys.spin[sys.neighbours[spin][i]][1]
                        + sys.spin[spin][2]* sys.spin[sys.neighbours[spin][i]][2];
        energy -= sys.J/2*dot_product;
    }
    energy += sys.B*sys.spin[spin][2];
    return energy;
};

// Calculate the total energy of the system
void Calculate_h(local_spins& sys){
    sys.H = 0; // Set H to zero before recalculating it
    double mag_energy = 0;
    for (int i=0; i<sys.n_spins; i++){
        int pad_i = sys.index_to_padded_index(i);
        sys.H += energy_calculation_nd(sys, pad_i)*0.5; // Half the energy, because we calculate on all the spins
        mag_energy += sys.spin[pad_i][2]; 
    }
    sys.H += sys.B*mag_energy * 0.5; // Half of the magnetization energy is removed above
};

// Write the spin configurations in the output file.
void Writeoutput(local_spins& sys, std::ofstream& file){  
    // Loop over all spins, and write out position and spin direction
    // ONLY ONE RANK WRITES AT A TIME!!!!
    //
    int pad_i;
    for (int i = 0; i<sys.n_spins; i++){
        pad_i = sys.index_to_padded_index(i);
        file << sys.position[pad_i][0] << " " << sys.position[pad_i][1] << " "  << sys.position[pad_i][2] << " "
            << sys.spin[pad_i][0] << " " << sys.spin[pad_i][1] << " "  << sys.spin[pad_i][2] << " "
            << energy_calculation_nd(sys,pad_i)
            << std::endl;
    }
};

void exchange_ghost_cells(local_spins &local_sys,
                        MPI_Aint &sdispls, MPI_Aint &rdispls, 
                        MPI_Datatype &sendtypes, MPI_Datatype &recvtypes,
                        MPI_Comm cart_comm){
    int counts[6] = {1,1,1,1,1,1};
    // Make Proof of concept work
    std::vector<double> sx;
    for (uint64_t i=0; i<local_sys.spin.size(); i++){
        sx.push_back(local_sys.spin[i][0]);
    }
    // Send ghostcells
    std::cout << "MPI Error Code: "<< MPI_Neighbor_alltoallw (sx.data(), counts,  &sdispls, &sendtypes,
                            sx.data(), counts,  &rdispls, &recvtypes, cart_comm) << std::endl; 
};

void Simulate(spin_system& sys, local_spins& localsys,MPI_Aint &sdispls, MPI_Aint &rdispls, 
                        MPI_Datatype &sendtypes, MPI_Datatype &recvtypes,
                        MPI_Comm cart_comm){
    double old_energy, new_energy, spin_azimuthal, spin_polar, probability_of_change;
    std::vector<double> old_state(3);
    int not_flipped = 0;
    int flipped = 0;
    if(verbose) std::cout << old_state[0] << " " << old_state[1] << " " << old_state[2] << std::endl;
    
    auto begin = std::chrono::steady_clock::now();
    std::cout << "Temp: " <<sys.Temperature
                            << std::endl;

    int local_iterations = sys.flips/mpi_size;
    for (int iteration=0; iteration<local_iterations; iteration++){
        // Choose a random spin site
        srand(iteration);
        int rand_site = rand()%(localsys.n_spins);
        rand_site = localsys.index_to_padded_index(rand_site);
        
        // Calculate it's old energy
        old_energy = energy_calculation_nd(localsys, rand_site);

        // Store its old state. 
        old_state[0] = localsys.spin[rand_site][0];
        old_state[1] = localsys.spin[rand_site][1];
        old_state[2] = localsys.spin[rand_site][2];
      
        // Generate new state
        spin_azimuthal = (double) rand()/RAND_MAX * M_PI;
        srand(mpi_rank*iteration + iteration);
        spin_polar = (double) rand()/RAND_MAX * 2. * M_PI;
        localsys.spin[rand_site] = {sin(spin_azimuthal)*cos(spin_polar), 
                            sin(spin_azimuthal)*sin(spin_polar),
                            cos(spin_azimuthal)};
        flipped++; 

        // Calculate if it lowers energy
        new_energy = energy_calculation_nd(localsys, rand_site);
        if (new_energy > old_energy){
            // If not, see if it should be randomised in direction
            if (verbose) std::cout << "New Energy: " << new_energy << " Old energy: " << old_energy << std::endl;
            probability_of_change = exp(-(new_energy-old_energy)/localsys.Temperature); // FIgure out probability of change
            //std::cout << "Change prob: " << probability_of_change << " Exp factor :" << -(new_energy-old_energy)/(Boltzmann*sys.Temperature) << std::endl;
            srand((mpi_rank+1)*(iteration+1)*2);
            if (probability_of_change < (double) rand()/RAND_MAX){
                // If not, revert to old state
                localsys.spin[rand_site] = {
                    old_state[0], old_state[1], old_state[2]
                };
                not_flipped++;
                flipped--;
                new_energy = old_energy;
            }
        }
        
        // Change H to represent the total energy of the system. Gave wrong results. Unclear why. Currently commented out. H is just calculated at the end, as it is not used anywhere in the loop anyway.
        //sys.H = sys.H - old_energy + new_energy;

    
    } 
    exchange_ghost_cells(localsys,sdispls, rdispls, 
                        sendtypes, recvtypes,
                        cart_comm);
    Calculate_h(localsys);
    auto end = std::chrono::steady_clock::now();
    std::cout << "Elapsed Time: " << (end-begin).count() / 1000000000.0 << std::endl;
    std::cout << "Not flipped no. is " << not_flipped << std::endl;
    std::cout << "Flipped no. is " << flipped << std::endl;
    std::cout << "Total energy: " << localsys.H << std::endl;
    
}
//=============================================================================================
//=========================   MAIN FUNCTION   =================================================
//=============================================================================================
int main(int argc, char* argv[]){
    if (verbose) std::cout << "Hello Heisenberg!" << std::endl;

    //MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if (verbose) {
    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    std::cout << "Heisenberg running on " << processor_name
              << ", rank " << mpi_rank << " out of " << mpi_size << std::endl;

    }

    //Initialise and load config
    spin_system global_sys({argv, argv+argc});

    //Setup MPI
    int dims[3] = {nproc_x, nproc_y, nproc_z};
    int periods[3] = {1,1,1};
    int coords[3];
    MPI_Dims_create(mpi_size, 3, dims);
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 3, dims, periods,
            0, &cart_comm);
    MPI_Cart_coords(cart_comm, mpi_rank, 3, coords);

    int nleft, nright, nbottom, ntop, nfront, nback;
    MPI_Cart_shift(cart_comm, 0,1,&nleft,&nright);
    MPI_Cart_shift(cart_comm, 1,1,&nbottom,&ntop);
    MPI_Cart_shift(cart_comm, 2,1,&nfront,&nback);


    const long int offset_x = global_sys.xlen * coords[0] / nproc_x -1;
    const long int offset_y = global_sys.ylen * coords[1] / nproc_y -1;
    const long int offset_z = global_sys.zlen * coords[2] / nproc_z -1;

    const long int end_x = global_sys.xlen * (coords[0]+1) / nproc_x +1; 
    const long int end_y = global_sys.ylen * (coords[1]+1) / nproc_y +1; 
    const long int end_z = global_sys.zlen * (coords[2]+1) / nproc_z +1;


    local_spins local_sys(global_sys, 
            end_x-offset_x, end_y-offset_y, end_z-offset_z,
            offset_x, offset_y, offset_z);
    //========================================================================================
    //========================= START OF GHOST CELL COMMUNICATION SETUP ======================
    //========================================================================================

    // Define subarray types for ghost cell exchanges
    /* The following send blocks are defined as follows:
     * x_type sends to the nearest MPI block in the x direction
     * y_type sends to the nearest MPI block in the y direction
     * z_type sends to the nearest MPI block in the z direction
     * 
     * The blocks define the parts of data that will be sent in 
     * MPI_Neigbor_alltoallw.
     * 
     * HEAVILY INSPIRED BY https://github.com/essentialsofparallelcomputing/Chapter8/blob/master/GhostExchange/CartExchange3D_Neighbor/CartExchange.cc
     * LINE 110 AND FORWARD.
    */
    //MPI_Datatype Vector_type;
    //MPI_Type_contiguous(3,MPI_DOUBLE,&Vector_type);
    //MPI_Type_commit(&Vector_type);
    const int array_sizes[] = {local_sys.pad_xlen,local_sys.pad_ylen, local_sys.pad_zlen};
    // send subarrays
    int subarray_sizes_x[] = { 1,local_sys.ylen,local_sys.zlen};
    int subarray_x_start[] = {0,1,1};
    MPI_Datatype x_type;
    MPI_Type_create_subarray (3, array_sizes, subarray_sizes_x, subarray_x_start,
                            MPI_ORDER_C, MPI_DOUBLE, &x_type);
    MPI_Type_commit(&x_type);

    int subarray_sizes_y[] = {local_sys.xlen, 1, local_sys.zlen};
    int subarray_y_start[] = {1,0,1};
    MPI_Datatype y_type;
    MPI_Type_create_subarray (3, array_sizes, subarray_sizes_y, subarray_y_start,
                            MPI_ORDER_C, MPI_DOUBLE, &y_type);
    MPI_Type_commit(&y_type);

    int subarray_sizes_z[] = { local_sys.xlen, local_sys.ylen,1};
    int subarray_z_start[] = {1,1,0};
    MPI_Datatype z_type;
    MPI_Type_create_subarray (3, array_sizes, subarray_sizes_z, subarray_z_start,
                            MPI_ORDER_C, MPI_DOUBLE, &z_type);
    MPI_Type_commit(&z_type);
    int xyplane_mult = local_sys.pad_ylen*local_sys.pad_xlen*8; //8 because datatype is 3 doubles, 
    int xstride_mult = local_sys.pad_xlen*8;
    // Define displacements of send and receive in bottom top left right.
    MPI_Aint sdispls[6] = { 8,
                            local_sys.zlen*8 ,
                            xstride_mult,
                            local_sys.ylen*xstride_mult,
                            xyplane_mult,
                            local_sys.pad_zlen*xyplane_mult                          
                            };
    MPI_Aint rdispls[6] = {0,
                            (local_sys.xlen+1)*8,
                            0,
                            (local_sys.ylen+1)*xstride_mult,
                            0,
                            (local_sys.zlen+1)*xyplane_mult,
                            };
    
    MPI_Datatype sendtypes[6] = { x_type, x_type, y_type, y_type, z_type, z_type};
    MPI_Datatype recvtypes[6] = { x_type, x_type, y_type, y_type, z_type, z_type};

    //========================================================================================
    //========================= END OF GHOST CELL COMMUNICATION SETUP ========================
    //========================================================================================
    //Generate system TODO: Done in parallel
    generate_positions_box(local_sys);
    generate_spin_directions(local_sys);
    generate_neighbours(local_sys);
    //Magic TODO h as reduction
    Calculate_h(local_sys);
    Simulate(global_sys,local_sys,*sdispls, *rdispls, 
                        *sendtypes, *recvtypes,
                        cart_comm);
    
    if (mpi_rank == 0) std::cout << "Final energy: " << global_sys.H << std::endl;

    for (int i = 0; i < mpi_size; i++ ){
            if (mpi_rank == i){
                std::ofstream file(local_sys.filename); // open file
                Writeoutput(local_sys, file);
    }
    }
    MPI_Type_free(&x_type);
    MPI_Type_free(&y_type);
    MPI_Type_free(&z_type);
    //MPI_Type_free(&Vector_type);

    MPI_Finalize();
    
    return 0;
}

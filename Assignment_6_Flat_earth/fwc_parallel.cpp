#include <vector>
#include <iostream>
#include <H5Cpp.h>
#include <chrono>
#include <cmath>
#include <numeric>
#include <mpi.h>
bool verbose = false;
// Get the number of processes
int mpi_size;

// Get the rank of the process
int mpi_rank;

int nproc_lat, nproc_lon;

/** Representation of a flat world */
class World {
public:
    // The current world time of the world
    double time;
    // The size of the world in the latitude dimension and the global size
    uint64_t latitude, global_latitude;
    // The size of the world in the longitude dimension
    uint64_t longitude, global_longitude;
    // The offset for this rank in the latitude dimension
    long int offset_latitude;
    // The offset for this rank in the longitude dimension
    long int offset_longitude;
    // The temperature of each coordinate of the world.
    // NB: it is up to the calculation to interpret this vector in two dimension.
    std::vector<double> data;
    // The measure of the diffuse reflection of solar radiation at each world coordinate.
    // See: <https://en.wikipedia.org/wiki/Albedo>
    // NB: this vector has the same length as `data` and must be interpreted in two dimension as well.
    std::vector<double> albedo_data;

    /** Create a new flat world.
     *
     * @param latitude     The size of the world in the latitude dimension.
     * @param longitude    The size of the world in the longitude dimension.
     * @param temperature  The initial temperature (the whole world starts with the same temperature).
     * @param albedo_data  The measure of the diffuse reflection of solar radiation at each world coordinate.
     *                     This vector must have the size: `latitude * longitude`.
     */
    World(uint64_t latitude, uint64_t longitude, double temperature,
          std::vector<double> albedo_data) : latitude(latitude), longitude(longitude),
                                             data(latitude * longitude, temperature),
                                             albedo_data(std::move(albedo_data)) {}
};

double checksum(World &world) {
    // 
    double check = 0;
    double cs=0;
    // TODO: make sure checksum is computed globally
    // only loop *inside* data region -- not in ghostzones!
    for (uint64_t i = 1; i < world.latitude - 1; ++i)
    for (uint64_t j = 1; j < world.longitude - 1; ++j) {
        cs += world.data[i*world.longitude + j];
    }
    MPI_Reduce(&cs,&check,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    return check;
}

void stat(World &world, MPI_Comm comm) {
    // TODO: make sure stats are computed globally
    double mint_local = 1e99;
    double maxt_local = 0;
    double meant_local = 0;
    double mint = 1e99;
    double maxt = 0;
    double meant = 0;
    
    for (uint64_t i = 1; i < world.latitude - 1; ++i)
    for (uint64_t j = 1; j < world.longitude - 1; ++j) {
        mint_local = std::min(mint_local,world.data[i*world.longitude + j]);
        maxt_local = std::max(maxt_local,world.data[i*world.longitude + j]);
        meant_local += world.data[i*world.longitude + j];
    
    }
    MPI_Reduce(&mint_local, &mint,1, MPI_DOUBLE, MPI_MIN, 0, comm);
    MPI_Reduce(&maxt_local, &maxt,1, MPI_DOUBLE, MPI_MAX, 0, comm);
    MPI_Reduce(&meant_local, &meant,1, MPI_DOUBLE, MPI_SUM, 0, comm);
    meant = meant / (world.global_latitude * world.global_longitude);
    if (mpi_rank==0)
    std::cout <<   "min: " << mint
              << ", max: " << maxt
              << ", avg: " << meant << std::endl;
}

/** Exchange the ghost cells i.e. copy the second data row and column to the very last data row and column and vice versa.
 *
 * @param world  The world to fix the boundaries for.
 */
void exchange_ghost_cells(World &world,MPI_Aint &sdispls, MPI_Aint &rdispls, MPI_Datatype &sendtypes, MPI_Datatype &recvtypes, MPI_Comm cart_comm) {
    // TODO: figure out exchange of ghost cells bewteen ranks
    /*
    for (uint64_t i = 0; i < world.latitude; ++i) {
        world.data[i*world.longitude + 0] = world.data[i*world.longitude + world.longitude-2];
        world.data[i*world.longitude + world.longitude-1] = world.data[i*world.longitude + 1];
    }

    for (uint64_t j = 0; j < world.longitude; ++j) {
        world.data[0*world.longitude + j] = world.data[(world.latitude-2)*world.longitude + j];
        world.data[(world.latitude-1)*world.longitude + j] = world.data[1*world.longitude + j];
    }*/
    int counts[4] = {1,1,1,1};
    MPI_Neighbor_alltoallw (&world.data.data()[0], counts,  &sdispls, &sendtypes,
                            &world.data.data()[0], counts,  &rdispls, &recvtypes, cart_comm);
}

/** Warm the world based on the position of the sun.
 *
 * @param world      The world to warm.
 */
void radiation(World& world,MPI_Aint &sdispls, MPI_Aint &rdispls, MPI_Datatype &sendtypes, MPI_Datatype &recvtypes, MPI_Comm cart_comm) {
    double sun_angle = std::cos(world.time);
    double sun_intensity = 865.0;
    double sun_long = (std::sin(sun_angle) * (world.global_longitude / 2))
                      + world.global_longitude / 2.;
    double sun_lat = world.global_latitude / 2.;
    double sun_height = 100. + std::cos(sun_angle) * 100.;
    double sun_height_squared = sun_height * sun_height;
    
    for (uint64_t i = 1; i < world.latitude-1; ++i) {
        for (uint64_t j = 1; j < world.longitude-1; ++j) {
            // Euclidean distance between the sun and each earth coordinate
            double delta_lat  = sun_lat  - (i + world.offset_latitude);
            double delta_long = sun_long - (j + world.offset_longitude); 
            double dist = sqrt(delta_lat*delta_lat + 
                               delta_long*delta_long + 
                               sun_height_squared);
            world.data[i * world.longitude + j] += \
                (sun_intensity / dist) * (1. - world.albedo_data[i * world.longitude + j]);
        }
    }
    exchange_ghost_cells(world, sdispls, rdispls, sendtypes, recvtypes, cart_comm);
}

/** Heat radiated to space
 *
 * @param world  The world to update.
 */
void energy_emmision(World& world) {
    for (uint64_t i = 0; i < world.latitude * world.longitude; ++i) {
        world.data[i] *= 0.99;
    }
}

/** Heat diffusion
 *
 * @param world  The world to update.
 */
void diffuse(World& world,MPI_Aint &sdispls, MPI_Aint &rdispls, MPI_Datatype &sendtypes, MPI_Datatype &recvtypes, MPI_Comm cart_comm) {
    std::vector<double> tmp = world.data;
    for (uint64_t k = 0; k < 10; ++k) {
        for (uint64_t i = 1; i < world.latitude - 1; ++i) {
            for (uint64_t j = 1; j < world.longitude - 1; ++j) {
                // 5 point stencil
                double center = world.data[i * world.longitude + j];
                double left = world.data[(i - 1) * world.longitude + j];
                double right = world.data[(i + 1) * world.longitude + j];
                double up = world.data[i * world.longitude + (j - 1)];
                double down = world.data[i * world.longitude + (j + 1)];
                tmp[i * world.longitude + j] = (center + left + right + up + down) / 5.;
            }
        }
        std::swap(world.data, tmp);  // swap pointers for the two arrays
        exchange_ghost_cells(world, sdispls, rdispls, sendtypes, recvtypes, cart_comm); // update ghost zones
    }
}

/** One integration step at `world_time`
 *
 * @param world      The world to update.
 */
void integrate(World& world,MPI_Aint &sdispls, MPI_Aint &rdispls, MPI_Datatype &sendtypes, MPI_Datatype &recvtypes, MPI_Comm cart_comm) {
    radiation(world, sdispls, rdispls, sendtypes, recvtypes, cart_comm);
    energy_emmision(world);
    diffuse(world, sdispls, rdispls, sendtypes, recvtypes, cart_comm);
}

/** Read a world model from a HDF5 file
 *
 * @param filename The path to the HDF5 file.
 * @return         A new world based on the HDF5 file.
 */
World read_world_model(const std::string& filename) {
    //std::cout << filename << std::endl;
    H5::H5File file(filename, H5F_ACC_RDONLY);
    H5::DataSet dataset = file.openDataSet("world");
    H5::DataSpace dataspace = dataset.getSpace();

    if (dataspace.getSimpleExtentNdims() != 2) {
        throw std::invalid_argument("Error while reading the model: the number of dimension must be two.");
    }

    if (dataset.getTypeClass() != H5T_FLOAT or dataset.getFloatType().getSize() != 8) {
        throw std::invalid_argument("Error while reading the model: wrong data type, must be double.");
    }

    hsize_t dims[2];
    dataspace.getSimpleExtentDims(dims, NULL);
    std::vector<double> data_out(dims[0] * dims[1]);
    dataset.read(data_out.data(), H5::PredType::NATIVE_DOUBLE, dataspace, dataspace);
    std::cout << "World model loaded -- latitude: " << (unsigned long) (dims[0]) << ", longitude: "
              << (unsigned long) (dims[1]) << std::endl;
    return World(static_cast<uint64_t>(dims[0]), static_cast<uint64_t>(dims[1]), 293.15, std::move(data_out));
}

/** Write data to a hdf5 file
 *
 * @param group  The hdf5 group to write in
 * @param name   The name of the data
 * @param shape  The shape of the data
 * @param data   The data
 */
void write_hdf5_data(H5::Group& group, const std::string& name,
                     const std::vector <hsize_t>& shape, const std::vector<double>& data) {
    H5::DataSpace dataspace(static_cast<int>(shape.size()), &shape[0]);
    H5::DataSet dataset = group.createDataSet(name.c_str(), H5::PredType::NATIVE_DOUBLE, dataspace);
    dataset.write(&data[0], H5::PredType::NATIVE_DOUBLE);
}

/** Write a history of the world temperatures to a HDF5 file
 *s
 * @param world     world to write
 * @param filename  The output filename of the HDF5 file
 */
void write_hdf5(const World& world, const std::string& filename, uint64_t iteration) {

    static H5::H5File file(filename, H5F_ACC_TRUNC);

    H5::Group group(file.createGroup("/" + std::to_string(iteration)));
    write_hdf5_data(group, "world", {world.latitude, world.longitude}, world.data);
}

/** Simulation of a flat word climate
 *
 * @param num_of_iterations  Number of time steps to simulate
 * @param model_filename     The filename of the world model to use (HDF5 file)
 * @param output_filename    The filename of the written world history (HDF5 file)
 */
void simulate(uint64_t num_of_iterations, const std::string& model_filename, const std::string& output_filename) {

    // for simplicity, read in full model
    World global_world = read_world_model(model_filename); 

    // Define number of ranks in each direction
    nproc_lon = int(sqrt(mpi_size));
    nproc_lat = int(sqrt(mpi_size));
    if(mpi_rank==0&&verbose) std::cout << nproc_lon << " " << nproc_lat<<std::endl;
    
    // Set up cartesian topology 
    int dims[2] = {nproc_lat, nproc_lon};
    int periods[2] = {1, 1};
    int coords[2];
    MPI_Dims_create(mpi_size, 2, dims);
    MPI_Comm cart_comm;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &cart_comm);
    MPI_Cart_coords(cart_comm, mpi_rank, 2, coords);

    int nleft, nright, nbottom, ntop;
    MPI_Cart_shift(cart_comm, 1, 1, &nleft, &nright);
    MPI_Cart_shift(cart_comm, 0, 1, &nbottom, &ntop);
    
    // figure out size of domain for this rank
    //// Domain start point offset
    const long int offset_longitude = global_world.longitude*coords[1]/nproc_lon-1; // -1 because first cell is a ghostcell
    const long int offset_latitude  = global_world.latitude*coords[0]/nproc_lat-1;

    //// Domain end point
    const long int end_longitude = global_world.longitude*(coords[1]+1)/nproc_lon+1;
    const long int end_latitude  = global_world.latitude*(coords[0]+1)/nproc_lat+1;

    //// Domain size
    const uint64_t longitude = end_longitude-offset_longitude;//global_world.longitude/nproc_lon + 2; // one ghost cell on each end
    const uint64_t latitude  = end_latitude -offset_latitude ;//global_world.latitude/nproc_lat  + 2;
    if(verbose) std::cout << "Rank" << mpi_rank << ", " << longitude << " " << latitude << std::endl;
    
    
    // copy over albedo data to local world data
    std::vector<double> albedo(longitude*latitude);
    for (uint64_t i = 1; i < latitude-1; ++i)
    for (uint64_t j = 1; j < longitude-1; ++j) {
        uint64_t k_global = (i + offset_latitude) * global_world.longitude
                          + (j + offset_longitude);
        albedo[i * longitude + j] = global_world.albedo_data[k_global];
    }

    // create local world data
    World world = World(latitude, longitude, 293.15, albedo);
    world.global_latitude  = global_world.latitude;
    world.global_longitude = global_world.longitude;
    world.offset_latitude  = offset_latitude;
    world.offset_longitude = offset_longitude;

    // set up counters and loop for num_iterations of integration steps
    const double delta_time = world.global_longitude / 36.0;

    // Define subarray types for ghost cell exchanges
    const int array_sizes[] = {(latitude),(longitude)};
    // send subarrays
    int subarray_sizes_x[] = {latitude-2,1};
    int subarray_horiz_start[] = {1,0};
    MPI_Datatype horiz_type;
    MPI_Type_create_subarray (2, array_sizes, subarray_sizes_x, subarray_horiz_start,
                            MPI_ORDER_C, MPI_DOUBLE, &horiz_type);
    MPI_Type_commit(&horiz_type);

    int subarray_sizes_y[] = {1,longitude-2};
    int subarray_vert_start[] = {0,1};
    MPI_Datatype vert_type;
    MPI_Type_create_subarray (2, array_sizes, subarray_sizes_y, subarray_vert_start,
                            MPI_ORDER_C, MPI_DOUBLE, &vert_type);
    MPI_Type_commit(&vert_type);

    // Define displacements of send and receive in bottom top left right.
    MPI_Aint sdispls[4] = { longitude * 8,
                            ( latitude - 2 ) * longitude * 8,
                            8,
                            (longitude - 2 ) * 8};
    MPI_Aint rdispls[4] = {0,(latitude-1)*longitude*8,0,(longitude-1)*8};
    
    MPI_Datatype sendtypes[4] = {vert_type, vert_type, horiz_type, horiz_type};
    MPI_Datatype recvtypes[4] = {vert_type, vert_type, horiz_type, horiz_type};

    // Define gather_counts and gather_disp, to be used when gathering all data onto rank 0 for writing out. 
    int gather_counts[mpi_size];
    int gather_disp[mpi_size];
    int temp_disp = 0;
    
    for(int rank=0;rank<mpi_size;rank++){
        const int add =(global_world.latitude*(rank/nproc_lat+1)/nproc_lat-global_world.latitude*(rank/nproc_lat)/nproc_lat)*(global_world.longitude*(rank%nproc_lon+1)/nproc_lon-
global_world.longitude*(rank%nproc_lon)/nproc_lon);
        gather_counts[rank] = add;

        const int rank_disp = temp_disp;
        gather_disp[rank] = rank_disp;
        
        temp_disp += add;
        //if(mpi_rank==0) std::cout << rank << "  " <<gather_counts[rank] <<", " << gather_disp[rank] << std::endl;
    }
    
    // ==============================================================================
    // BEGIN SIMULATION
    // ==============================================================================
    

    auto begin = std::chrono::steady_clock::now();
    for (uint64_t iteration=0; iteration < num_of_iterations; ++iteration) {
        // Get current time
        world.time = iteration / delta_time;

        // Do integration
        integrate(world, *sdispls, *rdispls, *sendtypes, *recvtypes,cart_comm);
        
        // Gather all the data on rank 0 so it can be written to a file
        std::vector <double> world_WO_ghost; // World_WithOut_ghost cells
        for (int i = 0; i<(latitude)*(longitude); i++){
                if (i%longitude == 0 || i%longitude == longitude -1){continue;}
                if (i < longitude || i > (latitude-1)*longitude) continue;
                world_WO_ghost.push_back( world.data[i]);
        }
        if(verbose&&iteration==0)std::cout << "Rank" << mpi_rank << ", " <<world_WO_ghost.size();
        MPI_Gatherv(world_WO_ghost.data(), world_WO_ghost.size(), MPI_DOUBLE, global_world.data.data(), gather_counts,gather_disp, MPI_DOUBLE, 0,MPI_COMM_WORLD);
        
        
       if (!output_filename.empty()) {
            // Only rank zero writes water history to file
            if (mpi_rank == 0) {
                write_hdf5(global_world, output_filename, iteration);
                std::cout << iteration << " -- ";
                stat(world, cart_comm);
            }
        }
        // Wait for everyone to get here
        MPI_Barrier(cart_comm);
        
    }
    double check = checksum(world);
    if (mpi_rank==0){
    auto end = std::chrono::steady_clock::now();
    
    stat(world, cart_comm);
    std::cout << "checksum      : " << check << std::endl;
    std::cout << "elapsed time  : " << (end - begin).count() / 1000000000.0 << " sec" << std::endl;
    }
    MPI_Type_free(&horiz_type);
    MPI_Type_free(&vert_type);
}

/** Main function that parses the command line and start the simulation */
int main(int argc, char **argv) {

    // Initialize the MPI environment
    MPI_Init(&argc, &argv);

    // Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    // Get the rank of the process
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

    // Get the name of the processor
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message
    std::cout << "Flat World Climate running on " << processor_name
              << ", rank " << mpi_rank << " out of " << mpi_size << std::endl;

    uint64_t iterations=0;
    std::string model_filename;
    std::string output_filename;

    std::vector <std::string> argument({argv, argv+argc});

    for (long unsigned int i = 1; i<argument.size() ; i += 2){
        std::string arg = argument[i];
        if(arg=="-h"){ // Write help
            std::cout << "./fwc --iter <number of iterations>"
                      << " --model <input model>"
                      << " --out <name of output file>\n";
            exit(0);
        } else if (i == argument.size() - 1)
            throw std::invalid_argument("The last argument (" + arg +") must have a value");
        else if(arg=="--iter"){
            if ((iterations = std::stoi(argument[i+1])) < 0) 
                throw std::invalid_argument("iter most be positive (e.g. -iter 1000)");
        } else if(arg=="--model"){
            model_filename = argument[i+1];
        } else if(arg=="--out"){
            output_filename = argument[i+1];
        } else{
            std::cout << "---> error: the argument type is not recognized \n";
        }
    }
    if (model_filename.empty())
        throw std::invalid_argument("You must specify the model to simulate "
                                    "(e.g. --model models/small.hdf5)");
    if (iterations==0)
        throw std::invalid_argument("You must specify the number of iterations "
                                    "(e.g. --iter 10)");
    simulate(iterations, model_filename, output_filename);

    MPI_Finalize();

    return 0;
}

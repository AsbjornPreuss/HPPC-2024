#include <vector>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <numeric>
#include <cassert>
#include <array>
#include <algorithm>
#include <openacc.h>
using real_t = float;
constexpr size_t NX =4096, NY = 4096;
constexpr int number_of_gangs_in_run = 1;
constexpr int vector_lengths = 128;
using grid_t = std::array<std::array<real_t, NX>, NY>;

class Sim_Configuration {
public:
    int iter = 1000;  // Number of iterations
    double dt = 0.05;       // Size of the integration time step
    real_t g = 9.80665;     // Gravitational acceleration
    real_t dx = 1;          // Integration step size in the horizontal direction
    real_t dy = 1;          // Integration step size in the vertical direction
    int data_period = 100;  // how often to save coordinate to file
    std::string filename = "sw_output.data";   // name of the output file with history

    Sim_Configuration(std::vector <std::string> argument){
        for (long unsigned int i = 1; i<argument.size() ; i += 2){
            std::string arg = argument[i];
            if(arg=="-h"){ // Write help
                std::cout << "./par --iter <number of iterations> --dt <time step>"
                          << " --g <gravitational const> --dx <x grid size> --dy <y grid size>"
                          << "--fperiod <iterations between each save> --out <name of output file>\n";
                exit(0);
            } else if (i == argument.size() - 1)
                throw std::invalid_argument("The last argument (" + arg +") must have a value");
            else if(arg=="--iter"){
                if ((iter = std::stoi(argument[i+1])) < 0) 
                    throw std::invalid_argument("iter most be a positive integer (e.g. -iter 1000)");
            } else if(arg=="--dt"){
                if ((dt = std::stod(argument[i+1])) < 0) 
                    throw std::invalid_argument("dt most be a positive real number (e.g. -dt 0.05)");
            } else if(arg=="--g"){
                g = std::stod(argument[i+1]);
            } else if(arg=="--dx"){
                if ((dx = std::stod(argument[i+1])) < 0) 
                    throw std::invalid_argument("dx most be a positive real number (e.g. -dx 1)");
            } else if(arg=="--dy"){
                if ((dy = std::stod(argument[i+1])) < 0) 
                    throw std::invalid_argument("dy most be a positive real number (e.g. -dy 1)");
            } else if(arg=="--fperiod"){
                if ((data_period = std::stoi(argument[i+1])) < 0) 
                    throw std::invalid_argument("dy most be a positive integer (e.g. -fperiod 100)");
            } else if(arg=="--out"){
                filename = argument[i+1];
            } else{
                std::cout << "---> error: the argument type is not recognized \n";
            }
        }
    }
};

/** Representation of a water world including ghost lines, which is a "1-cell padding" of rows and columns
 *  around the world. These ghost lines is a technique to implement periodic boundary conditions. */
class Water {
public:
    grid_t u{}; // The speed in the horizontal direction.
    grid_t v{}; // The speed in the vertical direction.
    grid_t e{}; // The water elevation.
    Water() {
        for (size_t i = 1; i < NY - 1; ++i) 
        for (size_t j = 1; j < NX - 1; ++j) {
            real_t ii = 100.0 * (i - (NY - 2.0) / 2.0) / NY;
            real_t jj = 100.0 * (j - (NX - 2.0) / 2.0) / NX;
            e[i][j] = std::exp(-0.02 * (ii * ii + jj * jj));
        }
    }
};

/* Write a history of the water heights to an ASCII file
 *
 * @param water_history  Vector of the all water worlds to write
 * @param filename       The output filename of the ASCII file
*/
void to_file(const std::vector<grid_t> &water_history, const std::string &filename){
    std::ofstream file(filename);
    file.write((const char*)(water_history.data()), sizeof(grid_t)*water_history.size());
}

/** Exchange the horizontal ghost lines i.e. copy the second data row to the very last data row and vice versa.
 *
 * @param data   The data update, which could be the water elevation `e` or the speed in the horizontal direction `u`.
 * @param shape  The shape of data including the ghost lines.
 */
void exchange_horizontal_ghost_lines(grid_t& data0,grid_t& data1) {
    #pragma acc parallel loop gang num_gangs(number_of_gangs_in_run) vector_length(vector_lengths) present(data0,data1, NX, NY)
    for (uint64_t j = 0; j < NX; ++j) {
        data0[0][j]      = data0[NY-2][j]; 
        data0[NY-1][j]   = data0[1][j];
        data1[0][j]      = data1[NY-2][j]; 
        data1[NY-1][j]   = data1[1][j];
        
    }
}

/** Exchange the vertical ghost lines i.e. copy the second data column to the rightmost data column and vice versa.
 *
 * @param data   The data update, which could be the water elevation `e` or the speed in the vertical direction `v`.
 * @param shape  The shape of data including the ghost lines.
 */
void exchange_vertical_ghost_lines(grid_t& data0,grid_t& data1) {
    #pragma acc parallel loop gang num_gangs(number_of_gangs_in_run) vector_length(vector_lengths) present(data0,data1, NX, NY)
    for (uint64_t i = 0; i < NY; ++i) {
        data0[i][0] = data0[i][NX-2];
        data0[i][NX-1] = data0[i][1];
        data1[i][0] = data1[i][NX-2];
        data1[i][NX-1] = data1[i][1];
    }
}

/** One integration step
 *
 * @param w The water world to update.
 */
void integrate(Water &w, const real_t &dtdx, const real_t& dtdy, const real_t& g) {
    
    exchange_horizontal_ghost_lines(w.e,w.v);
    //exchange_horizontal_ghost_lines(w.v);
    exchange_vertical_ghost_lines(w.e,w.u);
    //exchange_vertical_ghost_lines(w.u);
    
    #pragma acc parallel loop gang vector collapse(2) num_gangs(number_of_gangs_in_run) vector_length(vector_lengths) present(w,dtdx,dtdy,g,NX,NY)
    for (uint64_t i = 0; i < NY - 1; ++i) {
    for (uint64_t j = 0; j < NX - 1; ++j){
        w.u[i][j] -= dtdx * g * (w.e[i][j+1] - w.e[i][j]);
        w.v[i][j] -= dtdy * g * (w.e[i + 1][j] - w.e[i][j]);
    }
    }
    #pragma acc parallel loop gang vector collapse(2) num_gangs(number_of_gangs_in_run) vector_length(vector_lengths) present(w,dtdx,dtdy,g,NX,NY)
    for (uint64_t i = 1; i < NY - 1; ++i) {
    for (uint64_t j = 1; j < NX - 1; ++j) {
        w.e[i][j] -= dtdx * (w.u[i][j] - w.u[i][j-1])
                   + dtdy * (w.v[i][j] - w.v[i-1][j]);
    }
    }
    
    
}

/** Simulation of shallow water
 *
 * @param num_of_iterations  The number of time steps to simulate
 * @param size               The size of the water world excluding ghost lines
 * @param output_filename    The filename of the written water world history (HDF5 file)
 */
void simulate(const Sim_Configuration config) {
    auto begin0 = 0.;
    Water water_world = Water();
    const real_t dt = config.dt; const real_t dx=config.dx; const real_t dy=config.dy; const real_t g=config.g;
    std::vector <grid_t> water_history;
    //std::array<grid_t, config.iter / config.data_period> water_history;
    auto begin = std::chrono::steady_clock::now();
    const real_t dtdx = dt/dx; const real_t dtdy = dt/dy;
    #pragma acc data copyin(water_world, dtdx,  dtdy, g, NX, NY)//water_world.e[0:NX][0:NY],water_world.u[0:NX][0:NY],water_world.v[0:NX][0:NY],
    {
    for (uint64_t t = 0; t < config.iter; ++t) {
        integrate(water_world, dtdx, dtdy, g);
        if (t % config.data_period == 0) {
            #pragma acc update host(water_world.e)
            water_history.push_back(water_world.e);
        }
    }
    }
    auto end = std::chrono::steady_clock::now();
    begin0 += (end-begin).count() / 1'000'000'000.;
    to_file(water_history, config.filename);

    std::cout << " 3 " << std::accumulate(water_world.e.front().begin(), water_world.e.back().end(), 0.0)<< " ";
    std::cout << (end - begin).count() / 1000000000.0 << " " << NX << " " << number_of_gangs_in_run << " " << vector_lengths << std::endl;
    
}

/** Main function that parses the command line and start the simulation */
int main(int argc, char **argv) {
    auto config = Sim_Configuration({argv, argv+argc});
    simulate(config);
    return 0;
}

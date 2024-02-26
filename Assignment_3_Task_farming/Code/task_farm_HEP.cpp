/*
  Assignment: Make an MPI task farm for analysing HEP data
  To "execute" a task, the worker computes the accuracy of a specific set of cuts.
  The resulting accuracy should be send back from the worker to the master.
*/

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <random>
#include <chrono>
#include <thread>
#include <array>
#include <vector>
#include <limits>
// To run an MPI program we always need to include the MPI headers
#include <mpi.h>

// Number of cuts to try out for each event channel.
// BEWARE! Generates n_cuts^8 permutations to analyse.
// If you run many workers, you may want to increase from 3.
const int n_cuts = 5;
const long n_settings = (long) pow(n_cuts,8);
const long NO_MORE_TASKS = n_settings+1;
const bool verbose = false;
// Class to hold the main data set together with a bit of statistics
class Data {
public:
    long nevents=0;
    std::string name[8] = { "averageInteractionsPerCrossing", "p_Rhad","p_Rhad1",
                            "p_TRTTrackOccupancy", "p_topoetcone40", "p_eTileGap3Cluster",
                            "p_phiModCalo", "p_etaModCalo" };
    std::vector<std::array<double, 8>> data; // event data
    std::vector<long> NvtxReco;              // counters; don't use them
    std::vector<long> p_nTracks;
    std::vector<long> p_truthType;           // authorative truth about a signal

    std::vector<bool> signal;                // True if p_truthType=2

    std::array<double,8> means_sig {0}, means_bckg {0}; // mean of signal and background for events
    std::array<double,8> flip; // flip sign if background larger than signal for type of event
};

// Routine to read events data from csv file and calculate a bit of statistics
Data read_data() {
    // name of data file
    std::string filename="mc_ggH_16_13TeV_Zee_EGAM1_calocells_16249871.csv";
    std::ifstream csvfile(filename); // open file

    std::string line;
    std::getline(csvfile, line); // skip the first line

    Data ds; // variable to hold all data of the file

    while (std::getline(csvfile,line)) {  // loop over lines until end of file
        if (line.empty()) continue;       // skip empty lines
        std::istringstream iss(line);
        std::string element;
        std::array<double,8> data;

        // read in one line of data in to class        
        std::getline(iss, element, ','); // line counter, skip it
        std::getline(iss, element, ','); // averageInteractionsPerCrossing
        data[0] = std::stod(element);
        std::getline(iss, element, ','); // NvtxReco
        ds.NvtxReco.push_back(std::stol(element));
        std::getline(iss, element, ','); // p_nTracks
        ds.p_nTracks.push_back(std::stol(element));
        // Load in a loop the 7 next data points: 
        // p_Rhad, p_Rhad1, p_TRTTrackOccupancy, p_topoetcone40, p_eTileGap3Cluster, p_phiModCalo, p_etaModCalo
        for(int i=1; i<8; i++) {
            std::getline(iss, element, ',');
            data[i] = std::stod(element);
        }
        std::getline(iss, element, ','); // p_truthType
        ds.p_truthType.push_back(std::stol(element));
        ds.data.push_back(data);
        ds.nevents++;
    }

    // Calculate means. Signal has p_truthType = 2
    ds.signal.resize(ds.nevents);
    long nsig=0, nbckg=0;
    for (long ev=0; ev<ds.nevents; ev++) {
        ds.signal[ev] = ds.p_truthType[ev] == 2;
        if (ds.signal[ev]) {
            for(int i=0; i<8; i++) ds.means_sig[i] += ds.data[ev][i];
            nsig++;
        } else {
            for(int i=0; i<8; i++) ds.means_bckg[i] += ds.data[ev][i];
            nbckg++;
        }
    }
    for(int i=0; i<8; i++) {
        ds.means_sig[i]  = ds.means_sig[i] / nsig;
        ds.means_bckg[i] = ds.means_bckg[i] / nbckg;
    }
    
    // check for flip and change sign of data and means if needed
    for(int i=0; i<8; i++) {
        ds.flip[i]= (ds.means_bckg[i] < ds.means_sig[i]) ? -1 : 1;
        for (long ev=0; ev<ds.nevents; ev++) ds.data[ev][i] *= ds.flip[i];
        ds.means_sig[i]  = ds.means_sig[i] * ds.flip[i];
        ds.means_bckg[i] = ds.means_bckg[i] * ds.flip[i];
    }

   return ds;
}

// call this function to complete the task. It calculates the accuracy of a given set of settings
double task_function(std::array<double,8>& setting, Data& ds) {
    // pred evalautes to true if cuts for events are satisfied for all cuts
    std::vector<bool> pred(ds.nevents,true);
    for (long ev=0; ev<ds.nevents; ev++)
        for (int i=0; i<8; i++)
            pred[ev] = pred[ev] and (ds.data[ev][i] < setting[i]);

    // accuracy is percentage of events that are predicted as true signal if and only if a true signal
    double acc=0;
    for (long ev=0; ev<ds.nevents; ev++) acc += pred[ev] == ds.signal[ev];

    return acc / ds.nevents;
}

void master (int nworkers, Data& ds) {
    std::array<std::array<double,8>,n_cuts> ranges; // ranges for cuts to explore

    // loop over different event channels and set up cuts
    for(int i=0; i<8; i++) {
        for (int j=0; j<n_cuts; j++)
            ranges[j][i] = ds.means_sig[i] + j * (ds.means_bckg[i] - ds.means_sig[i]) / n_cuts;
    }
    
    // generate list of all permutations of the cuts for each channel
    std::vector<std::array<double,8>> settings(n_settings);
    for (long k=0; k<n_settings; k++) {
        long div = 1;
        std::array<double,8> set;
        for (int i=0; i<8; i++) {
            long idx = (k / div) % n_cuts;
            set[i] = ranges[idx][i];
            div *= n_cuts;
        }
        settings[k] = set;
    }

    // results vector with the accuracy of each set of settings
    std::vector<double> accuracy(n_settings);

    auto tstart = std::chrono::high_resolution_clock::now(); // start time (nano-seconds)

//===============================================================================
    //                  OUR CODE STARTS HERE
    //===============================================================================
    if (verbose) std::cout << "Master  : Initialized\n";
    // Initialize vectors for keeping track of things
    std::vector<int> task_at_worker(nworkers);
    // Initialize variables used later.
    int tag = 1;
    int destination;
    //MPI_Request IRreq; // Only used for non-blocking mpi messages
    // Now we complete all the tasks.

    int task = 0;

    // Init tasks to all workers
    for (int worker = 0; worker < nworkers; worker++){
        if (verbose) std::cout <<"Master  : Working on task " << task <<"\n";
        tag = 1;
        destination = worker + 1;
        MPI_Send(&settings[task], 8, MPI_DOUBLE, destination, tag,
            MPI_COMM_WORLD); // Send the task to the worker queue
        task_at_worker[worker] = task;
    	task++;
    }

    // Receive tasks from any worker
    int finished_worker_rank;
    int finished_worker;
    double acc;
    MPI_Status status;
    while (task<n_settings){
        MPI_Recv(&acc, 1, MPI_DOUBLE,  MPI_ANY_SOURCE, MPI_ANY_TAG,
                    MPI_COMM_WORLD, &status);
    finished_worker_rank = status.MPI_SOURCE;
	finished_worker = finished_worker_rank-1;
	accuracy[task_at_worker[finished_worker]] = acc;
        // Give worker new task.
	task_at_worker[finished_worker] = task;
        MPI_Send(&settings[task], 8, MPI_DOUBLE, finished_worker_rank, tag,
            MPI_COMM_WORLD);
	task++;
    }
     
    // Now shut down the workers, once all tasks are done
    std::array<double,8> shut_down_data = {std::numeric_limits<double>::max(),0,0,0,0,0,0,0};
    for (int worker=0; worker < nworkers; worker++){
	finished_worker_rank = worker +1;
	
    MPI_Recv(&accuracy[task_at_worker[worker]], 1, MPI_DOUBLE, finished_worker_rank,
		       tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
    MPI_Send(&shut_down_data, 8, MPI_DOUBLE,  finished_worker_rank, tag,
                MPI_COMM_WORLD); 

    }
    //===============================================================================
    //                        OUR CODE ENDS HERE
    //===============================================================================
    


// THIS CODE SHOULD BE REPLACED BY TASK FARM
    // loop over all possible cuts and evaluate accuracy
    //for (long k=0; k<n_settings; k++)
    //    accuracy[k] = task_function(settings[k], ds);
    // THIS CODE SHOULD BE REPLACED BY TASK FARM
    // ================================================================

    auto tend = std::chrono::high_resolution_clock::now(); // end time (nano-seconds)
    // diagnostics
    // extract index and value for best accuracy
    double best_accuracy_score=0;
    long idx_best=0;
    for (long k=0; k<n_settings; k++)
        if (accuracy[k] > best_accuracy_score) {
            best_accuracy_score = accuracy[k];
            idx_best = k;
        }
    
    //std::cout << "Best accuracy obtained :" << best_accuracy_score << "\n";
    //std::cout << "Final cuts :\n";
    //for (int i=0; i<8; i++)
        //std::cout << std::setw(30) << ds.name[i] << " : " << settings[idx_best][i]*ds.flip[i] << "\n";
    
    //std::cout <<  "\n";
    //std::cout <<  "Number of settings:" << std::setw(9) << n_settings << "\n";
    //std::cout <<  "Elapsed time      :" << std::setw(9) << std::setprecision(4)
    //          << (tend - tstart).count()*1e-9 << "\n";
    //std::cout <<  "task time [mus]   :" << std::setw(9) << std::setprecision(4)
    //          << (tend - tstart).count()*1e-3 / n_settings << "\n";
    std::cout << best_accuracy_score << " ";
    for (int i=0; i<8; i++)
        std::cout << settings[idx_best][i]*ds.flip[i] << " ";
    std::cout << n_settings << " " << (tend-tstart).count()*1e-9 << " " << (tend - tstart).count()*1e-3 / n_settings << " " << nworkers +1 << "\n";
}

void worker (int rank, Data& ds) {
    //===============================================================================
    //                  OUR CODE STARTS HERE
    //===============================================================================
    std::array<double,8> setting;           // The task must be declared. It is given in the MPI_Recv function.
    int tag = 1;         // This worker will only take the process allocated to it, in the tag line.
    int master_rank = 0;    // The master rank is 0 as a default.
    double acc;
   
    while (true){
        //std::cout << "Worker " << rank << ": Initialized and waiting for a task\n";
        MPI_Recv(&setting, 8, MPI_DOUBLE, master_rank, tag,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (setting[0] == std::numeric_limits<double>::max()) break;
        //std::cout << "Worker " << rank << ": Received an MPI message to sleep for " << task << " milliseconds\n";
        acc = task_function(setting, ds);
        //std::cout << "Worker " << rank << ": Completed task. Sending result back to master\n";
        MPI_Send(&acc, 1, MPI_DOUBLE, master_rank, tag,
                        MPI_COMM_WORLD);
    }    
    if(verbose) std::cout << "Worker "<< rank << " Done " << "\n";
    
    //===============================================================================
    //                        OUR CODE ENDS HERE
    //===============================================================================
}

int main(int argc, char *argv[]) {
    int nrank, rank;

    MPI_Init(&argc, &argv);                // set up MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nrank); // get the total number of ranks
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // get the rank of this process

    // All ranks need to read the data
    Data ds = read_data();

    if (rank == 0)       // rank 0 is the master
        master(nrank-1, ds); // there is nrank-1 worker processes
    else                 // ranks in [1:nrank] are workers
        worker(rank, ds);

    MPI_Finalize();      // shutdown MPI
}

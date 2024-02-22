/*
  Assignment: Make an MPI task farm. A "task" is a randomly generated integer.
  To "execute" a task, the worker sleeps for the given number of milliseconds.
  The result of a task should be send back from the worker to the master. It
  contains the rank of the worker
*/

#include <iostream>
#include <random>
#include <chrono>
#include <thread>
#include <array>
#include <vector>

// To run an MPI program we always need to include the MPI headers
#include <mpi.h>

const int NTASKS=3000;  // number of tasks
const int RANDOM_SEED=12345;

void master (int nworkers) {
    // Initialize task, result and worker queue array
    std::array<int, NTASKS> tasks, result;
    // set up a random number generator
    std::random_device rd;
    //std::default_random_engine engine(rd());
    std::default_random_engine engine;
    engine.seed(RANDOM_SEED);
    // make a distribution of random integers in the interval [0:30]
    std::uniform_int_distribution<int> distribution(0, 30);

    for (int& t : tasks) {
        t = distribution(engine);   // set up some "tasks"   
    }

    //===============================================================================
    //                  OUR CODE STARTS HERE
    //===============================================================================
    std::cout << "Master  : Initialized\n";
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
        std::cout <<"Master  : Working on task " << task <<"\n";
        tag = 1;
        destination = worker + 1;
        MPI_Send(&tasks[task], 1, MPI_INT,  destination, tag,
            MPI_COMM_WORLD); // Send the task to the worker queue
        task_at_worker[worker] = task;
	task++;
    }

    // Receive tasks from any worker
    int finished_worker_rank;
    int finished_worker;
    while (task<NTASKS){
        MPI_Recv(&finished_worker_rank, 1, MPI_INT,  MPI_ANY_SOURCE, tag,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	finished_worker = finished_worker_rank-1;
	result[task_at_worker[finished_worker]] = finished_worker_rank;
        // Give worker new task.
	task_at_worker[finished_worker] = task;
        MPI_Send(&tasks[task], 1, MPI_INT, finished_worker_rank, tag,
            MPI_COMM_WORLD);
	task++;
    }
     
    // Now shut down the workers, once all tasks are done
    int shut_down_data = -1;
    for (int worker=0; worker < nworkers; worker++){
	finished_worker_rank = worker +1;
	MPI_Recv(&result[task_at_worker[worker]], 1, MPI_INT, finished_worker_rank,
		       tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	
        MPI_Send(&shut_down_data, 1, MPI_INT,  finished_worker_rank, tag,
                MPI_COMM_WORLD); 
    }
    //===============================================================================
    //                        OUR CODE ENDS HERE
    //===============================================================================

    // Print out a status on how many tasks were completed by each worker
    for (int worker=0; worker<nworkers; worker++) {
        int tasksdone = 0; int workdone = 0;
        for (int itask=0; itask<NTASKS; itask++)
        if (result[itask]==worker+1) {
            tasksdone++;
            workdone += tasks[itask];
        }
        std::cout << "Master  : Worker " << worker << " solved " << tasksdone << 
                    " tasks\n";    
    }
}

// call this function to complete the task. It sleeps for task milliseconds
void task_function(int task) {
    std::this_thread::sleep_for(std::chrono::milliseconds(task));
}

void worker (int rank) {
    //===============================================================================
    //                  OUR CODE STARTS HERE
    //===============================================================================
    int task = 1;           // The task must be declared. It is given in the MPI_Recv function.
    int tag = 1;         // This worker will only take the process allocated to it, in the tag line.
    int master_rank = 0;    // The master rank is 0 as a default.
    while (task >= 0){
        //std::cout << "Worker " << rank << ": Initialized and waiting for a task\n";
        MPI_Recv(&task, 1, MPI_INT, master_rank, tag,
                    MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //std::cout << "Worker " << rank << ": Received an MPI message to sleep for " << task << " milliseconds\n";
        task_function(task);
        //std::cout << "Worker " << rank << ": Completed task. Sending result back to master\n";
        MPI_Send(&rank, 1, MPI_INT, master_rank, tag,
                        MPI_COMM_WORLD);
        //std::cout << "Worker "<< rank << ": Sent result from task "<< task << " back\n";
    }
    //===============================================================================
    //                        OUR CODE ENDS HERE
    //===============================================================================
}

int main(int argc, char *argv[]) {
    int nrank, rank;
    enum {i_am_master, i_am_worker};
    MPI_Init(&argc, &argv);                // set up MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nrank); // get the total number of ranks
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // get the rank of this process

    if (rank == i_am_master){              // rank 0 is the master
        master(nrank-1);                   // There is nrank-1 workers
    }
    else{                                  // ranks in [1:nrank] are workers
        worker(rank);
    }

    MPI_Finalize();                        // shutdown MPI
}

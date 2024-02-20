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

// To run an MPI program we always need to include the MPI headers
#include <mpi.h>

const int NTASKS=5000;  // number of tasks
const int RANDOM_SEED=1234;

void master (int nworker) {
    std::array<int, NTASKS> task, result;

    // set up a random number generator
    std::random_device rd;
    //std::default_random_engine engine(rd());
    std::default_random_engine engine;
    engine.seed(RANDOM_SEED);
    // make a distribution of random integers in the interval [0:30]
    std::uniform_int_distribution<int> distribution(0, 30);

    for (int& t : task) {
        t = distribution(engine);   // set up some "tasks"
    }
    std::cout << "Master  : Initialized\n";
    // master(nrank-1); // there is nrank-1 worker processes
    MPI_Send(&task, 1, MPI_INT, 1, 1,
                MPI_COMM_WORLD);
    std::cout << "Master  : Task sent to worker "<< 1 <<"\n";
    MPI_Recv(&task, 1, MPI_INT, 1, 1,
                MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    std::cout << "Master  : Task finished and received\n";


    /*
    IMPLEMENT HERE THE CODE FOR THE MASTER
    ARRAY task contains tasks to be done. Send one element at a time to workers
    ARRAY result should at completion contain the ranks of the workers that did
    the corresponding tasks
    */

    // Print out a status on how many tasks were completed by each worker
    for (int worker=1; worker<=nworker; worker++) {
        int tasksdone = 0; int workdone = 0;
        for (int itask=0; itask<NTASKS; itask++)
        if (result[itask]==worker) {
            tasksdone++;
            workdone += task[itask];
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
    int task; // The task must be declared. It is given in the MPI_Recv function.
    int tag = rank; // This worker will only take the process allocated to it, in the tag line.
    int master_rank = 0; // The master rank is 0 as a default.
    std::cout << "Worker " << rank << ": Initialized and waiting for a task\n";
    MPI_Recv(&task, 1, MPI_INT, master_rank, tag,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    std::cout << "Worker " << rank << ": Received an MPI message to sleep for " << task << " milliseconds\n";
    task_function(task);
    std::cout << "Worker " << rank << ": Completed task. Sending result back to master\n";
    MPI_Send(&task, 1, MPI_INT, master_rank, tag,
                    MPI_COMM_WORLD);
    std::cout << "Worker "<< rank << ": Sent result from task "<< task << " back\n";
}

int main(int argc, char *argv[]) {
    int nrank, rank;
    enum {i_am_master, i_am_worker};
    MPI_Init(&argc, &argv);                // set up MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nrank); // get the total number of ranks
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // get the rank of this process

    if (rank == i_am_master){      // rank 0 is the master
        master(nrank);
    }
    else{                 // ranks in [1:nrank] are workers
        worker(rank);
    }

    MPI_Finalize();      // shutdown MPI
}

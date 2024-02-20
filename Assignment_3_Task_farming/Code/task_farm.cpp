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

const int NTASKS=100;  // number of tasks
const int RANDOM_SEED=12345;

void master (int nworkers) {
    std::cout << "Master  : Initialized\n";
    
    // Initialize task, result and worker queue array
    std::array<int, NTASKS> tasks, result;

    std::vector<int> worker_queue(nworkers);
    for (int i=0; i < nworkers; i++){
        worker_queue[i] = i +1;
        std::cout << "Worker queue[" << i << "] is " << worker_queue[i] << "\n";
    }
    // An iterator for removing the worker which has just had a message sent to it.
    std::vector<int>::iterator worker_position = worker_queue.begin();

    int tag = 1;
    int destination;
    int source;
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
    for (int task=0; task< NTASKS; task++) {
        std::cout <<"Working on task " << task <<"\n";
        // Send out tasks to all workers
        for (long unsigned int worker = 0; worker < worker_queue.size(); worker++){
            tag = 1;
            destination = worker + 1;
            MPI_Send(&tasks[task], 1, MPI_INT,  destination, tag,
                MPI_COMM_WORLD); // Send the task to the worker queue
            // Remove the worker after having sent a message to it
            worker_queue.erase(worker_position);
        }
        // Receive all tasks that have been sent out
        for (long unsigned int worker = 0; worker < worker_queue.size(); worker++){
            source = worker + 1;
            MPI_Recv(&result[task], 1, MPI_INT,  source, tag,
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            // Put the worker back in the worker queue.
            worker_queue.push_back(result[task]);
        }
    }
    // Now shut down the workers, once all tasks are done
    int shut_down_data = -1;
    MPI_Send(&shut_down_data, 1, MPI_INT,  1, tag,
            MPI_COMM_WORLD); 



    /*
    IMPLEMENT HERE THE CODE FOR THE MASTER
    ARRAY task contains tasks to be done. Send one element at a time to workers
    ARRAY result should at completion contain the ranks of the workers that did
    the corresponding tasks
    */

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
    int task = 1; // The task must be declared. It is given in the MPI_Recv function.
    int tag = rank; // This worker will only take the process allocated to it, in the tag line.
    int master_rank = 0; // The master rank is 0 as a default.
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
}

int main(int argc, char *argv[]) {
    int nrank, rank;
    enum {i_am_master, i_am_worker};
    MPI_Init(&argc, &argv);                // set up MPI
    MPI_Comm_size(MPI_COMM_WORLD, &nrank); // get the total number of ranks
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // get the rank of this process

    if (rank == i_am_master){      // rank 0 is the master
        master(nrank-1); // There is nrank-1 workers
    }
    else{                 // ranks in [1:nrank] are workers
        worker(rank);
    }

    MPI_Finalize();      // shutdown MPI
}

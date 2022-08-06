/*
 * Copyright: Davide Lomeo
 * 
 * Author: Davide Lomeo
 * 
 * Birth: 2021-03-12
 * 
 * Description: main file of the program that perform the domain decomposition using 
 *              parallelisation through MPI. This file reads in the user-defined 
 *              parameters from .txt file and controls the flow of the entire program.
 * 
 * Github Repository: https://github.com/davidelomeo/wave_equation_using_mpi
 * 
 * Revision: Last revised by Davide Lomeo on 2022-08-06
*/

#include <mpi.h>
#include "Decomposition.h"
#include "DataType.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <vector>

// (used Google's cpplint package, that apparenlty does not accept chrono as a header file)
#include <chrono>  // NOLINT

#include "Decomposition.cpp"
#include "DataType.cpp"

#define DO_TIMING

// Function to read the user-defined parameters for the problem from a textfile
void read_textfile(std::string filename, std::vector<double> &params) {

    // assigning ifstream values
    std::ifstream file;

    // opening the file using in-stream read-only format
    file.open(filename, std::ios_base::in);

    // adding the variables from the textfile to the the array
    if (file.is_open()) {
        double a;
        while (file >> a) {
            // adding the values to the vector, making sure that it is of type double
            params.push_back(static_cast<double>(a));
        }
    }
    file.close();
}
// ----------------------------------------------------------


int id, p;

int main(int argc, char *argv[]) {

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
	srand(time(NULL) + id * 10);

    // initialising the name of the textfile containing the parameters
    std::string filename = "parameters.txt";

    // initialising a vector to contain the parameters from textfile
    std::vector<double> params;

    // reading data from textfile
    read_textfile(filename, params);

    // from vector to array
    double *vars = &params[0];

// starting the timer here to avoid getting longer timings if the HD
// of the computer running the program is particularly slow.
#ifdef DO_TIMING
    // using a barrier to ensure that every processors gets here first
    // before starting to time the code
	MPI_Barrier(MPI_COMM_WORLD);
	auto start = std::chrono::high_resolution_clock::now();
#endif

    // step 1 - creating a decomposition object by loading the Decomposition
    //          class. This will handle the doamin decomposition
    auto *decomposition = new Decomposition<double>(id, p, vars);

    // step 2 - setting up the distribution of the domain according
    //          to the number of processors
    decomposition->domain_setup();

    // step 3 - creating 3 arrays for each processor with the size
    //          of their assigned portion of the domain
    decomposition->domain_per_node();

    // step 4 - outputting p matrices to show the domain at t = 0
    decomposition->grid_to_file();
    decomposition->out_cnt++;
    decomposition->t_out += decomposition->dt_out;

    // step 5 - computing the initial disturbance of the domain in each
    //          processor and output the result to file. Here the code
    //          iterates 'n_disturbances' times, according to how many
    //          initial disturbance have been setup by the user.
    // indices of first disturbance
    int r_splash_index = 11, x_splash_index = 12, y_splash_index = 13;
    int n_disturbances = vars[10];
    for (int i = 0; i < n_disturbances; i++) {
        decomposition->initial_disturbance(vars[r_splash_index], vars[x_splash_index], vars[y_splash_index]);
        r_splash_index += 3;
        x_splash_index += 3;
        y_splash_index += 3;
    }
    decomposition->grid_to_file();

    // step 6 - computing the iterations until reqching the user-defined
    //          t_max value
    decomposition->do_iteration();

#ifdef DO_TIMING
	MPI_Barrier(MPI_COMM_WORLD);
	auto finish = std::chrono::high_resolution_clock::now();

    // some print statements to help visualise the size of the matrix,
    // how many processors were used and the time taken in seconds
	if (id == 0) {
		std::chrono::duration<double> elapsed = finish - start;
		std::cout << std::setprecision(5);
		std::cout << "\nOutput:" << std::endl;
        std::cout << "Matrix size:    " << decomposition->imax << " x " << decomposition->jmax << std::endl;
        std::cout << "N. Processors:  " << p << std::endl;
        std::cout << "Time taken (s): " << elapsed.count() << "\n" << std::endl;
        std::cout << "The following are availebl ONLY if timing the central code:" << std::endl;
        std::cout << "Central iter. (s): " << decomposition->elapsed_iter.count() << std::endl;
        std::cout << "Grid to file (s): " << decomposition->elapsed_to_file.count() << "\n" << std::endl;
	}
#endif

    // deleting the decomposition and vars objects
    delete decomposition;

    MPI_Finalize();
    return 0;
}

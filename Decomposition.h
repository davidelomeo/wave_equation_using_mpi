/*
 * Copyright: Davide Lomeo
 * 
 * Author: Davide Lomeo
 * 
 * Birth: 2021-03-12
 * 
 * Description: Header file that contains the class that stores the variables and 
 *              calls the functions needed in the Decomposition.cpp source file.
 * 
 * Github Repository: https://github.com/davidelomeo/wave_equation_using_mpi
 * 
 * Revision: Last revised by Davide Lomeo on 2022-08-06
*/

#pragma once

#include "DataType.h"
#include <mpi.h>

#include <iostream>
#include <iomanip>
#include <string>

// (used Google's cpplint package, that apparenlty does not accept chrono as a header file)
#include <chrono>  // NOLINT

template <class T>
class Decomposition {
 public:
    // Constructor
    Decomposition(int &rank, int &proc, double *vars);

    // Destructor
    ~Decomposition();

    // ======================================== problem setup ========================================
    // Setting up the distribution of the domain across processors and handle odd row/col numbers
    void domain_setup(void);
    void distribution_check(int &reminder, int &loc_min, int &loc_max, int &dom, int &n_block);

    // Creating 3 arrays (current, old, new) for each of the chunks of the domain.
    void setup_array(double *&array_1d);
    void domain_per_node(void);

    // Setting an initial disturbance per node
    void initial_disturbance(double r_splash, double x_splash, double y_splash);


    // ======================================= domain evolution ======================================
    // Computing the evolution of the wave at each iteration, and apply boundary conditions.
    void do_iteration(void);
    void domain_evolution(void);

    // computing Neumann and Dirichlet boundary conditions in the case of fixed boundaries
    void neumann_boundaries(void);
    void dirichlet_boundaries(void);


    // ======================================== communications =======================================
    // Setting up communications between opposite boundaries in the case of a perodic boundaries problem
    void periodic_boundaries_comms(DataType<MPI_Datatype>& datatype, MPI_Request& requests);

    // Setting up communications between internal boundaries
    void internal_boundaries_comms(DataType<MPI_Datatype>& datatype, MPI_Request& requests);


    // ================================ printing and textfile reading ================================
    // Printing the arrays to file
    void grid_to_file(void);


    // ====================================== Public Variables =======================================
    // These variables have been made public to be accessible and modifiable from outside the current
    // class. This was done to allow to change the parameters of the problem from main.cpp, rather
    //  than only relying on the .txt file to read-in the input parameters.

    // rank and number of porcessors
    int id, p;

    // size of the matrix and physical size of the domain
    int imax, jmax, y_max, x_max;

    // time settings
    double t_max, t_out, dt_out;

    // speed of wave
    int c;

    // variables to store the iteration number
    int out_cnt = 0, it = 0;

    // flag needed to determine if using fixed or periodic boundary conditions and which boundary
    // type to use in the case of fixed boundaries.
    // Default -> Fixed boundaries with Neumann bounadry condition
    bool periodic_boundaries = false, fixed_bound_type = true;

    // number of disturbances 
    int initial_disturbances;

    // variables that will store the timings of the loop that solves the wave equation and the
    // section of the program that writes the solutions to file
    std::chrono::duration<double> elapsed_iter, elapsed_to_file;


    // ====================================== Private Variables ======================================
    // variables kept private to be unaccessible from outside of the current class. The purpose of
    // the variables being saved here is to make sure these are accessiblee from anywhere inside
    // the source code of this program.
    // Note that all the private variables to this class have a "_" at the end, to distinguish
    // them from variables accessible from outside the class. This was just a style decision

 private:
    // small spatial and temporal increments
    double dx_, dy_, dt_;
    double t_ = 0.0;

    // global size of the domain. These values are set to -1 for debugging purposes
    int imax_ = -1, jmax_ = -1;

    // i,j locations of processors on the domain
    int i_dom_, j_dom_;

    // variables that store the ids of the neighbours of each chunk of the domain
    int left_, right_, top_, bottom_;

    // row/ col length of the chunk of the domain, inlcusive of two ghost cells
    // setting the values to -1 for debugging purposes
    int block_imax_ = -1, block_jmax_ = -1;

    // chunk specific i,j locations without ghost cells
    int local_imin_, local_imax_, local_jmin_, local_jmax_;

    // Total number of neighbours. The value starts from the maximum possible neighbours
    // in this problem (4)
    int neighbours_ = 4;

    // creating null pointers that will be used to point to the boundaries of
    // each of the chunk of the domain
    double *send_left_ = nullptr, *send_right_ = nullptr, *send_top_ = nullptr, *send_bottom_ = nullptr;
    double *recv_in_right_ = nullptr, *recv_in_left_ = nullptr, *recv_in_bottom_ = nullptr, *recv_in_top_ = nullptr;

    // creating null pointers that will change to pointing to input arrays
    T *block_1d_ = nullptr, *new_block_1d_ = nullptr, *old_block_1d_ = nullptr;
};

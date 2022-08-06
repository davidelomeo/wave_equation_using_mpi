/*
 * Copyright: Davide Lomeo
 * 
 * Author: Davide Lomeo
 * 
 * Birth: 2021-03-12
 * 
 * Description: Source file that contains the core of the program. Here, the domain is
 * 				decomposed in n chunks according to the number of processors used, the 
 * 				wave equation is discretised, the communication between chunks of the
 * 				domain and the write to file occur. 
 * 
 * Github Repository: https://github.com/davidelomeo/wave_equation_using_mpi
 * 
 * Revision: Last revised by Davide Lomeo on 2022-08-06
*/

// Flag that determines if the program runs on a HPC system for the purpose of
// testing its efficeincy or if the program runs for the purpose of outputting
// many more matrices in order to obtain a smoother visual representation of the
// evolution of the domain. By default, the HPC flag will not write any file to
// disk,to prevent the pure timing of code to be polluted by a slow disk.
// #define HPC

// Flag that will enable the timings of blocks of the code
// #define DO_TIMING

// Flg that defines the most used mathematical constants
#define _USE_MATH_DEFINES

#include "Decomposition.h"
#include "DataType.h"
#include <mpi.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <utility>

// Integer array of dimensions of the domain (2 dimensions assumed for this problem).
// This is needed to split the domain using the MPI_Dims_create command. The variable
// was put here because MPI showed to have particular requirements about its location.
// I have tried to store the variable within the class or calling it from main trhough
// the constructor, but this did not work as expected and caused errors.
int dim[2];

// Constructor of the Decomposition object. Here, the used-defined parameters are assigned
// to variables within the Decompoistion class.
template <class T>
Decomposition<T>::Decomposition(int &rank, int &proc, double *vars) {

	// storing ids and processors numbers
	this->id = rank;
	this-> p = proc;

	// storing matrix imax and jmax read from file
	this->imax = vars[0];
	this->jmax = vars[1];

	// storing physical size of the domain read from file
	this->x_max = vars[2];
	this->y_max = vars[3];

	// storing time variables read from file
	this->t_max = vars[4];
	this->t_out = vars[5];
	this->dt_out = vars[6];

	// storing speed of wave read from file
	this->c = vars[7];

	// flag that changes the problem to periodic boundaries.
	// By default the problem uses fixed boundaries (false)
	if (vars[8] == 1)
		this->periodic_boundaries = true;

	// flag that changes the fixed Neumann boundary conditions to Dirichlet.
	// By default the problem uses Neumann boundary conditions (true)
	if (vars[9] == 0)
		this->fixed_bound_type = false;
	
	// storing the number of initial disturbances read from file
	this->initial_disturbances = vars[10];

    // copying the global imax/ jmax to private variables. These are originally
	// stored as public variables to allow to be read from outisde the class,
	// but the program only uses the private variable to prevent any changes
	// to these to be made after the Decomposition object has been created
    this->imax_ = this->imax;
    this->jmax_ = this->jmax;

    // setting spatial increments
    this->dx_ = this->x_max / (static_cast<double>(this->imax - 1));
	this->dy_ = this->y_max / (static_cast<double>(this->jmax - 1));

	// setting the time-step for the problem. This is striclty dependent on the
	// resolution due to the use of a discretisation in its explicit form
    this->dt_ = 0.1 * std::min(this->dx_, this->dy_) / this->c;
}
// ------------------------------------------------------------------------------------------------------------------------


// Destructor of the Decomposition object. The memory allocate to create the
// arrays is here freed.
template <class T>
Decomposition<T>::~Decomposition() {

	// deleting the created arrays as part of the object destruction process
    delete[] this->block_1d_;
    delete[] this->new_block_1d_;
	delete[] this->old_block_1d_;
}
// ------------------------------------------------------------------------------------------------------------------------


// Function that checks the size of the matrix and calculates the number of rows and cols
// to assign to each processor. This function virtually handles all the possible 2D domain
// size combinations, trying to distribute any extra row/ col evenly across chunks of the domain
template <class T>
void Decomposition<T>::distribution_check(int &reminder, int &loc_min, int &loc_max, int &dom, int &chunk) {

	// checking if the reminder between imax/ jmax and the size of the split domain is != 0
	// If this number is not 0, the function will distribute the additional rows/ cols across
	// the processors, starting from top to bottom
	if (reminder != 0) {
		// adding the extra row (col) to the bottom (right) end of the domain
		if (dom < reminder) {
			loc_min = (chunk +1) * dom;
			loc_max = loc_min +  chunk + 1;

		} else {
			// assigning the left-over row (col) to the current chunk of the domain when the remainder is 0
			loc_min = dom * chunk + reminder;
			loc_max = loc_min + chunk;
		}

	} else {
		// assigning the current chunk of the domain if the remainder is 0
		loc_min = chunk * dom;
		loc_max = chunk * (dom + 1 );
	}
}
// ------------------------------------------------------------------------------------------------------------------------


// Function to split the domain across n chunks (number of processors) allowing the lowest possible
// difference between m and n. If p is a odd number (e.g. 5, 7, 11), the domain will consist of a
// column of domain chunks (e.g. 5x1, 7x1, etc.)
template <class T>
void Decomposition<T>::domain_setup(void) {

	// Dims_create ustilise MPI built-in libraty to split the domain into even chunks according to
	// a cartesian grid of size n (here assumed 2D) split across p processors.
	MPI_Dims_create(p, 2, dim);

	// Computing the remiander to check if the input domain size can be evenly distributed across
	// the processors or if additional work is needed to allocate all the 'extra' rows/ cols.
	int i_remainder = this->imax_ % dim[0];
	int j_remainder = this->jmax_ % dim[1];

    // Assigning to each processor its location within the domain (previously split across n chunks).
	// This method follows a row-major ordering of the domain
    this->j_dom_ = id % dim[1];
    this->i_dom_ = id / dim[1];

    // Telling each processor its neighbouring cells, following a row-major ordering
    this->left_ = id - 1;
    this->right_ = id + 1;
    this->top_ = id - dim[1];
    this->bottom_ = id + dim[1];

	// assigning MPI_PROC_NULL (-2 if seen as integer) to the side of the chunks
	// of the domain that are global boundaries. The if statements also reduce of
	// 1 the maximum possible total number of neighbours (4) for each side of the
	// chunk of the domain that turns out to be a global boundary (i.e. no neighbours).
	// This information is used later on to detemrmine the number of requests needed in
	// the communications.
    if (this->j_dom_ == 0) {
        this->left_ = MPI_PROC_NULL;
		this->neighbours_ -= 1;
	}
    if (this->j_dom_ == dim[1] - 1) {
        this->right_ = MPI_PROC_NULL;
		this->neighbours_ -= 1;
	}
    if (this->i_dom_ == 0) {
        this->top_ = MPI_PROC_NULL;
		this->neighbours_ -= 1;
	}
    if (this->i_dom_ == dim[0] - 1) {
        this->bottom_ = MPI_PROC_NULL;
		this->neighbours_ -= 1;
	}

	// determining the m x n size of each of the chunks of the domain
    int chunk_imax = this->imax_ / dim[0];
	int chunk_jmax = this->jmax_ / dim[1];

    // Calling the function to assign to each processors its chunk of the domain. The function will
	// handle virtually all possible cases of even/ odd rows/ cols pairs. The function calls are two,
	// in order to handle row and col ordering separately.
	distribution_check(i_remainder, this->local_imin_, this->local_imax_, this->i_dom_, chunk_imax);
	distribution_check(j_remainder, this->local_jmin_, this->local_jmax_, this->j_dom_, chunk_jmax);

	// Assigning the imax and the jmax of the local chunk of the domain, adding a padding at the edjes
	// to facilitate the future communication between processors.
	this->block_imax_ = (this->local_imax_ - this->local_imin_) + 2;
	this->block_jmax_ = (this->local_jmax_ - this->local_jmin_) + 2;
}
// ------------------------------------------------------------------------------------------------------------------------


// Funciton that setup the values of a 1D array for each of the chunks of the domain
template <class T>
void Decomposition<T>::setup_array(double *&array_1d) {

	// iterating trhough the local imax/ jmax chunks of the domain
	for (int i = 0; i < this->block_imax_; i++)
		for (int j = 0; j < this->block_jmax_; j++)
			array_1d[i * this->block_jmax_ + j] = 0;
}
// ------------------------------------------------------------------------------------------------------------------------


// Function that creates 1D arrays to contain current, old and new values of the domain.
// These three arrays are needed to compute the evolution of the wave across the domain.
template <class T>
void Decomposition<T>::domain_per_node(void) {

	// creating arrays of double to store the values. As each node uses its own 'version'
	// of the program, the array is only of the size of the chunk of the domain that was
	// assigned to each of the nodes.
    this->block_1d_ = new T[this->block_imax_ * this->block_jmax_];
    this->new_block_1d_ = new T[this->block_imax_ * this->block_jmax_];
    this->old_block_1d_ = new T[this->block_imax_ * this->block_jmax_];

	// setting up the three arrays
	setup_array(this->block_1d_);
    setup_array(this->new_block_1d_);
	setup_array(this->old_block_1d_);
}
// ------------------------------------------------------------------------------------------------------------------------


// Function that computes the initial disturbance in the domain using user-defined values. This is also
// the only function the uses parameters that are external to the class. The choice was made to allow the
// user to have the freedom to have as many initial disturbances as desired, without interfering with the
// correct running of the code.
template <class T>
void Decomposition<T>::initial_disturbance(double r_splash, double x_splash, double y_splash) {
	for (int i = 1; i < this->block_imax_ - 1; i++)
		for (int j = 1; j < this->block_jmax_ - 1; j++) {

			// adding local is and js to asplit the intial disturbance across the
			// different chunks of the domain
			double x = this->dx_ * ((i-1) + this->local_imin_);
			double y = this->dy_ * ((j-1) + this->local_jmin_);

			double dist = sqrt(pow(x - x_splash, 2.0) + pow(y - y_splash, 2.0));

			if (dist < r_splash) {

				double h = 5.0*(cos(dist / r_splash * M_PI) + 1.0);

				this->block_1d_[i * this->block_jmax_ + j] = h;
				this->old_block_1d_[i * this->block_jmax_ + j] = h;
			}
		}
}
// ------------------------------------------------------------------------------------------------------------------------


// Function that performs the communication between opposite global boundaries in the case of a periodic domain.
// The function is capable of dealing with any number of processors distribution across the domain, including
// domain decompositions in only 1 processor that sends and receives to itself
template <class T>
void Decomposition<T>::periodic_boundaries_comms(DataType<MPI_Datatype>& datatype, MPI_Request& requests) {

	// initialising the communication's sources/ destinations to MPI_PROC_NULL to avoid the communications
	// going ahead if the processor does not deal with global boundaries
	int last_in_row = MPI_PROC_NULL, first_in_row  = MPI_PROC_NULL;
	int last_in_col  = MPI_PROC_NULL, first_in_col  = MPI_PROC_NULL;

	// assigning a real proc number to each of the sources/ destinations above, according to the distribution
	// of the domain. If the chunk of the domain deals with global boundaries, then it will send and receive
	// its boundaries to/ from the other end of the domain (which in thw case of 1 node used, will mean that
	// it will send to itself)
	if (this->left_ == MPI_PROC_NULL)
		last_in_row = id + (dim[1] - 1);
	if (this->right_ == MPI_PROC_NULL)
		first_in_row = id - (dim[1] - 1);
	if (this->top_ == MPI_PROC_NULL)
		last_in_col = id + (dim[0] * dim[1] - dim[1]);
	if (this->bottom_ == MPI_PROC_NULL)
		first_in_col = id - (dim[0] * dim[1] - dim[1]);

	// Series of non-blocking Irecv with custom datatypes
	MPI_Irecv(this->block_1d_, 1, datatype.receive_type_left, first_in_row, 1, MPI_COMM_WORLD, &requests);
	MPI_Irecv(this->block_1d_, 1, datatype.receive_type_right, last_in_row, 1, MPI_COMM_WORLD, &requests);
	MPI_Irecv(this->block_1d_, 1, datatype.receive_type_top, first_in_col, 1, MPI_COMM_WORLD, &requests);
	MPI_Irecv(this->block_1d_, 1, datatype.receive_type_bottom, last_in_col, 1, MPI_COMM_WORLD, &requests);

	// Series of non-blocking Isend with custom datatypes
	MPI_Isend(this->block_1d_, 1, datatype.send_type_left, last_in_row, 1, MPI_COMM_WORLD, &requests);
	MPI_Isend(this->block_1d_, 1, datatype.send_type_right, first_in_row, 1, MPI_COMM_WORLD, &requests);
	MPI_Isend(this->block_1d_, 1, datatype.send_type_top, last_in_col, 1, MPI_COMM_WORLD, &requests);
	MPI_Isend(this->block_1d_, 1, datatype.send_type_bottom, first_in_col, 1, MPI_COMM_WORLD, &requests);

	MPI_Barrier(MPI_COMM_WORLD);
}
// ------------------------------------------------------------------------------------------------------------------------


// Function that implements the Neumann boudary condition. With this condition, when the waves
// in the domain reach the edges, they get reflected and mantain their sign. This method also
// keeps the gradient of the displacement of the wave 0 a the boundary, such that the wave move
// upwards and the displacement comes down.
template <class T>
void Decomposition<T>::neumann_boundaries(void) {

	// applying condition to the top edge of the chunk of the domain
	if (this->local_imin_ == 0) {
		for (int j = 1; j < this->block_jmax_ - 1; j++)
			this->new_block_1d_[1 * this->block_jmax_ + j] = this->new_block_1d_[2 * this->block_jmax_ + j];
	}

	// applying condition to the left edge
	if (this->local_jmin_ == 0) {
		for (int i = 1; i < this->block_imax_ -1; i++)
			this->new_block_1d_[i * this->block_jmax_ + 1] = this->new_block_1d_[i * this->block_jmax_ + 2];
	}

	// applying condition to the bottom edge
	if (this->local_imax_ == this->imax_) {
		for (int j = 1; j < this->block_jmax_ -1; j++)
			this->new_block_1d_[(this->block_imax_ - 2) * this->block_jmax_ + j] = this->new_block_1d_[(this->block_imax_ - 3) * this->block_jmax_ + j];
	}

	// applying condition to the right edge
	if (this->local_jmax_ == this->jmax_) {
		for (int i = 1; i < this->block_imax_ - 1; i++)
			this->new_block_1d_[i * this->block_jmax_ + (this->block_jmax_ - 2)] = this->new_block_1d_[i * this->block_jmax_ + (this->block_jmax_ - 3)];
	}
}
// ------------------------------------------------------------------------------------------------------------------------


// Function to implement the Dirichlet boundary condition. This consist in setting a static
// displacement around the edges of the domain. The waves coming into the boundary will be
// reflected, but will change sign.
template <class T>
void Decomposition<T>::dirichlet_boundaries(void) {

	// applying condition to the top edge of the chunk of the domain
	if (this->local_imin_ == 0) {
		for (int j = 1; j < this->block_jmax_ - 1; j++)
			this->new_block_1d_[1 * this->block_jmax_ + j] = 0;
	}

	// applying condition to the left edge
	if (this->local_jmin_ == 0) {
		for (int i = 1; i < this->block_imax_ -1; i++)
			this->new_block_1d_[i * this->block_jmax_ + 1] = 0;
	}

	// applying condition to the bottom edge
	if (this->local_imax_ == this->imax_) {
		for (int j = 1; j < this->block_jmax_ -1; j++)
			this->new_block_1d_[(this->block_imax_ - 2) * this->block_jmax_ + j] = 0;
	}

	// applying condition to the right edge
	if (this->local_jmax_ == this->jmax_) {
		for (int i = 1; i < this->block_imax_ - 1; i++)
			this->new_block_1d_[i * this->block_jmax_ + (this->block_jmax_ - 2)] = 0;
	}
}
// ------------------------------------------------------------------------------------------------------------------------


// Function that implements the non-blocking communications between the processors in the
// internal boundaries of the domain. This function is used both in the case of periodic
// or fixed boundaries, as this are dealt with separately.
template <class T>
void Decomposition<T>::internal_boundaries_comms(DataType<MPI_Datatype>& datatype, MPI_Request& requests) {

	// non-blockin receive for each of the boundaries for the chunk of the domain that the processor deals with.
	// If one of the communications receive from MPI_NULL_PROC, the communication will automatically skip.
	MPI_Irecv(this->block_1d_, 1, datatype.receive_type_left, this->right_, 0, MPI_COMM_WORLD, &requests);
	MPI_Irecv(this->block_1d_, 1, datatype.receive_type_right, this->left_, 0, MPI_COMM_WORLD, &requests);
	MPI_Irecv(this->block_1d_, 1, datatype.receive_type_top, this->bottom_, 0, MPI_COMM_WORLD, &requests);
	MPI_Irecv(this->block_1d_, 1, datatype.receive_type_bottom, this->top_, 0, MPI_COMM_WORLD, &requests);

	// non-blockin send for each of the boundaries for the chunk of the domain that the processor deals with.
	// If one of the communications send to MPI_NULL_PROC, the communication will automatically skip.
	MPI_Isend(this->block_1d_, 1, datatype.send_type_left, this->left_, 0, MPI_COMM_WORLD, &requests);
	MPI_Isend(this->block_1d_, 1, datatype.send_type_right, this->right_, 0, MPI_COMM_WORLD, &requests);
	MPI_Isend(this->block_1d_, 1, datatype.send_type_top, this->top_, 0, MPI_COMM_WORLD, &requests);
	MPI_Isend(this->block_1d_, 1, datatype.send_type_bottom, this->bottom_, 0, MPI_COMM_WORLD, &requests);

	// barrier to make sure that all the processors continue the program only when ALL the communications are finished
	MPI_Barrier(MPI_COMM_WORLD);
}
// ------------------------------------------------------------------------------------------------------------------------


// Funciton that computes the movememnt of the wave across the domain. This gets called multiple times
// according to the timesteps and t_max set by the user. The function also calls the desired boundary
// condition when dealing with fixed boundaries.
template <class T>
void Decomposition<T>::domain_evolution() {

	// setting a deafult range of initial/ ending values for each chunk of the domain (accounting for ghost cells).
	// This will remain unchanged in the case of periodic boundaries, as the boundaries will be computed in the comms
	int start_ival = 1, start_jval = 1, end_ival = this->block_imax_ - 1, end_jval = this->block_jmax_ - 1;

	// The start and end values set above will increase and decrease of 1 respectively in the
	// presence of global boundaries. This will then allow the boundary conditions to
	// be applied after the iteration. This will not be run in the case of periodic boundaries
	if (!this->periodic_boundaries) {
		if (this->local_imin_ == 0)
			start_ival = 2;
		if (this->local_jmin_ == 0)
			start_jval = 2;
		if (this->local_imax_ == this->imax_)
			end_ival = this->block_imax_ - 2;
		if (this->local_jmax_ == this->jmax_)
			end_jval = this->block_jmax_ - 2;
	}

// Section to time the core of the program, that after the write to file, is the part of the code
// that takes the longest
#ifdef DO_TIMING
    // using a barrier to ensure that every processors gets here first
    // before starting to time the code
	MPI_Barrier(MPI_COMM_WORLD);
	auto start_iter_time = std::chrono::high_resolution_clock::now();
#endif

	// Solving the Wave Equation using the explicit scheme.
	for (int i = start_ival; i < end_ival; i++)
		for (int j = start_jval; j < end_jval; j++)
			this->new_block_1d_[i * this->block_jmax_ + j] = pow(this->dt_ * this->c, 2.0) * (
				(this->block_1d_[(i + 1) * this->block_jmax_ + j] - 2.0 * this->block_1d_[i * this->block_jmax_ + j] + this->block_1d_[(i - 1) * this->block_jmax_ + j]) / pow(this->dx_, 2.0) + (
					this->block_1d_[i * this->block_jmax_ + (j + 1)] - 2.0 * this->block_1d_[i * this->block_jmax_ + j] + this->block_1d_[i * this->block_jmax_ + (j - 1)]) / pow(this->dy_, 2.0)
					) + 2.0 * this->block_1d_[i * this->block_jmax_ + j] - this->old_block_1d_[i * this->block_jmax_ + j];

#ifdef DO_TIMING
	MPI_Barrier(MPI_COMM_WORLD);
	auto finish_iter_time = std::chrono::high_resolution_clock::now();
	if (id == 0)
		this->elapsed_iter += finish_iter_time - start_iter_time;
#endif

	// Checking if the chunk of the domain contains global boundaries, and calls the function
	// to compute the boundary conditions requested by the user (default is true - Neumann condition).
	// These are only called if the boundaries are not periodic
	if (!this->periodic_boundaries) {
		if (this->local_imin_ == 0 || this->local_jmin_ == 0 || this->local_jmax_ == this->jmax_ || this->local_imax_ == this->imax_) {
			if (this->fixed_bound_type)
				neumann_boundaries();
			else
				dirichlet_boundaries();
		}
	}

	// incrementing total time elpased by dt
	this->t_ += this->dt_;

	// stl swap function to swap the pointers of the arrays
	std::swap(this->old_block_1d_, this->new_block_1d_);
	std::swap(this->old_block_1d_, this->block_1d_);
}
// ------------------------------------------------------------------------------------------------------------------------


// Function that controls the evolution of the domain by calling the relevant functions
// according to the t_max set by the user. The function also calls the write-to-file
// function to allow each chunk to save its part of the domain to file at user-defined timesteps.
template <class T>
void Decomposition<T>::do_iteration(void) {

	// Creating a datatype object for the chunk of the domain that each processors deal with.
	// Creating the object faciliates the memory handling and allows for better code maintanability
	auto *datatype = new DataType<MPI_Datatype>(this->block_imax_, this->block_jmax_, this->block_1d_);

	// creating a requests object to support non-blocking communications. This is always equal to
	// 8 in the case of periodic boundaries and equal to the number of neighbouring processors
	// multiplied by 2 in the case of non periodic boundaries.
	MPI_Request *requests;

	if (this->periodic_boundaries)
		requests = new MPI_Request[8];
	else
		requests = new MPI_Request[this->neighbours_ * 2];

	// iteration that performns the central part of the program. This will continue as long as the
	// total elapsed time (t) is less than the maximum time set by the user
	while (this->t_ < this->t_max) {

		// calling the communication specific to periodic boundaries
		if (this->periodic_boundaries)
			periodic_boundaries_comms(*datatype, *requests);

		// calling the internal communications between neighbouring processors. The function is only called
		// if the number of processors is different from 1
		if (p-1 != 0)
			internal_boundaries_comms(*datatype, *requests);

		// calling the function to compute of the movement of the wave across the domain
		domain_evolution();

// Section to time the write to file section of the progrma. This was found to be the most
// time consuming task of the whole code
#ifdef DO_TIMING
                // using a barrier to ensure that every processors gets here first
                // before starting to time the code
                MPI_Barrier(MPI_COMM_WORLD);
                auto start_to_file_time = std::chrono::high_resolution_clock::now();
#endif

// if the HPC flag is active, then the program will not output the images, in order to prevent
// any costrained dictated by the speed of the local HPC hard drives.
#ifndef HPC
           	// calling the function to print the chunk of the domain
                if (this->t_out <= this->t_) {
                        if (this->out_cnt % 50 == 0)
                                std::cout << "Output: " << out_cnt << "\tt: " << t_ << "\titeration: " << it << std::endl;
                        grid_to_file();
                        this->out_cnt++;
                        this->t_out += this->dt_out;
                        }
#endif
                this->it++;

#ifdef DO_TIMING
        MPI_Barrier(MPI_COMM_WORLD);
        auto finish_to_file_time = std::chrono::high_resolution_clock::now();
        if (id == 0)
                this->elapsed_to_file += finish_to_file_time - start_to_file_time;
#endif
      	}

	// writing a .dat file that stores the i, j of the domain, the number of processors
	// the number of iterations, the type of boundaries and the number of initial disturbances.
	// These values are needed in the post-processing script.
	if (id == 0) {
		std::stringstream fname;
		std::fstream f1;
		fname << "./array_output/domain_parameters.dat";
		f1.open(fname.str().c_str(), std::ios_base::out);
		f1 << dim[0] << "\t" << dim[1] << "\t" << p << "\t" << out_cnt 
		   << "\t" << imax << "\t" << jmax << "\t" << periodic_boundaries 
		   << "\t" << fixed_bound_type << "\t" << initial_disturbances;
		f1.close();
	}

	// deleting datatype and requests objects to free previously allocated memory
	delete datatype;
	delete[] requests;
}
// ------------------------------------------------------------------------------------------------------------------------


// Function to save the output of the domain in a folder. The files will be named
// to include the n of iteration needed to obtain that output, followed by which
// processor has outputted that particular chunk of the domain.
template <class T>
void Decomposition<T>::grid_to_file() {

	std::stringstream fname;
	std::fstream f1;
	fname << "./array_output/out_iter_" << this->out_cnt << "_proc_" << id << ".dat";
	f1.open(fname.str().c_str(), std::ios_base::out);

	// using the following code to print the domain without the ghost cells. This is necessary
	// to facilitate the re-construcition of the domain in the post-processing code.
	for (int i = 1; i < this->block_imax_ - 1; i++) {
		for (int j = 1; j < this->block_jmax_ - 1; j++)
			f1 << this->block_1d_[i * this->block_jmax_ + j] << "\t";
		f1 << std::endl;
	}

	// // using the following code to print the domain with the ghost cells
	// for (int i = 0; i < this->block_imax_; i++) {
	// 	for (int j = 0; j < this->block_jmax_; j++)
	// 		f1 << this->block_1d_[i * this->block_jmax_ + j] << "\t";
	// 	f1 << std::endl;
	// }

	f1.close();
}
// ------------------------------------------------------------------------------------------------------------------------

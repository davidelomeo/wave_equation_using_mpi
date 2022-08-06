/*
 * Copyright: Davide Lomeo
 * 
 * Author: Davide Lomeo
 * 
 * Birth: 2021-03-15
 * 
 * Description: Source file that contains all the function needed to create custom
 *				MPI datatypes for each of the boundary of a domain of any size
 * 
 * Github Repository: https://github.com/davidelomeo/wave_equation_using_mpi
 * 
 * Revision: Last revised by Davide Lomeo on 2022-08-06
*/

#include "DataType.h"
#include <mpi.h>

#include <iomanip>
#include <iostream>

// costructor - storing external variables inside the class
template<class S>
DataType<S>::DataType(int row, int col, double *array) {

	// row/ cols lenght
    this->row_ = row;
    this->col_ = col;

	// storing the input array
    this->array_1d_ = array;

    // allocating memory to store the length on non-contiguous data, the type
    // of each block and the displacement in bytes between blocks.
    this->blocklengths_ = new int[this->row_];
    this->typelist_ = new S[this->row_];
    this->displacements_ = new MPI_Aint[this->row_];

    // function calls for each sending/ receiving boundary. These were placed here
    // as they only need to be called once, and never from outside the class. This
    // lead to the constructor to automatically generate all 8 datatypes for each
    // chunk of data (processor).
	create_left_type_toSend();
    create_right_type_toSend();
    create_top_type_toSend();
    create_bottom_type_toSend();
    create_left_type_toReceive();
    create_right_type_toReceive();
    create_top_type_toReceive();
    create_bottom_type_toReceive();
}
// ----------------------------------------------------------------------------------------------------------------


// Destructor used to free up allocated memory and datatypes after usage
template<class S>
DataType<S>::~DataType() {

    delete[] this->blocklengths_;
    delete[] this->displacements_;
    delete[] this->typelist_;
	MPI_Type_free(&send_type_left);
	MPI_Type_free(&send_type_right);
	MPI_Type_free(&send_type_top);
	MPI_Type_free(&send_type_bottom);
	MPI_Type_free(&receive_type_left);
	MPI_Type_free(&receive_type_right);
	MPI_Type_free(&receive_type_top);
	MPI_Type_free(&receive_type_bottom);
}
// ----------------------------------------------------------------------------------------------------------------


// Function that implements the left data type for the sending
template<class S>
void DataType<S>::create_left_type_toSend(void) {

    for (int i = 0; i < this->row_; i++) {

        // assigning length of one to each block of memory to send and
        // defining its type.
        this->blocklengths_[i] = 1;
        this->typelist_[i] = MPI_DOUBLE;

        // getting the address of each element of the left boundary that will be
        // sent and placing it in the array of addresses
        MPI_Get_address(&array_1d_[i  * this->col_ + 1], &this->temp_add_);
        this->displacements_[i] = this->temp_add_;
    }

    // getting the address of the beginning of the data
    MPI_Get_address(this->array_1d_, &this->add_start_);

    // creating a list of displacements
	for (int i = 0; i < this->row_; i++)
		this->displacements_[i] = this->displacements_[i] - this->add_start_;

    // generating the datatype and commiting it
    MPI_Type_create_struct(this->row_, this->blocklengths_, this->displacements_, this->typelist_, &this->send_type_left);
	MPI_Type_commit(&this->send_type_left);
}
// ----------------------------------------------------------------------------------------------------------------


// Function that implements the right data type for the sending
template<class S>
void DataType<S>::create_right_type_toSend(void) {

    for (int i = 0; i < this->row_; i++) {

        // assigning length of one to each block of memory to send and
        // defining its type.
        this->blocklengths_[i] = 1;
        this->typelist_[i] = MPI_DOUBLE;

        // getting the address of each element of the right boundary that will be
        // sent and placing it in the array of addresses
        MPI_Get_address(&array_1d_[i * this->col_ + (this->col_  - 2)], &this->temp_add_);
        this->displacements_[i] = this->temp_add_;
    }

    // getting the address of the beginning of the data
    MPI_Get_address(this->array_1d_, &this->add_start_);

    // creating a list of displacements
	for (int i = 0; i < this->row_; i++)
		this->displacements_[i] = this->displacements_[i] - this->add_start_;

    // generating the datatype and commiting it
    MPI_Type_create_struct(this->row_, this->blocklengths_, this->displacements_, this->typelist_, &this->send_type_right);
	MPI_Type_commit(&this->send_type_right);
}
// ----------------------------------------------------------------------------------------------------------------


// Function that implements the top data type that will be sent
template<class S>
void DataType<S>::create_top_type_toSend(void) {

    // getting the address of each element of the top boundary that will be sent
    // and of the start of the chunk of data
    MPI_Get_address(&array_1d_[1 * this->col_], &this->address_);
    MPI_Get_address(this->array_1d_, &this->add_start_);

    // storing the address of the contiguous chunk of data
    this->address_ = this->address_ - this->add_start_;

    // generating the datatype and commiting it
    MPI_Type_create_struct(1, &this->col_, &this->address_, &this->type_, &this->send_type_top);
	MPI_Type_commit(&this->send_type_top);
}
// ----------------------------------------------------------------------------------------------------------------



// Function that implements the bottom data type that will be sent
template<class S>
void DataType<S>::create_bottom_type_toSend(void) {

    // getting the address of each element of the top boundary that will be sent
    // and of the start of the chunk of data
    MPI_Get_address(&array_1d_[(this->row_ - 2) * this->col_], &this->address_);
    MPI_Get_address(this->array_1d_, &this->add_start_);

    // storing the address of the contiguous chunk of data
    this->address_ = this->address_ - this->add_start_;

    // generating the datatype and commiting it
    MPI_Type_create_struct(1, &this->col_, &this->address_, &this->type_, &this->send_type_bottom);
	MPI_Type_commit(&this->send_type_bottom);
}
// ----------------------------------------------------------------------------------------------------------------


// Function that implements the left data type for the receving chunk of data
template<class S>
void DataType<S>::create_left_type_toReceive(void) {


    for (int i = 0; i < this->row_; i++) {

        // assigning length of one to each block of memory to send and
        // defining its type.
        this->blocklengths_[i] = 1;
        this->typelist_[i] = MPI_DOUBLE;

        // getting the address of each element of the right boundary that will be
        // receving the data and placing it in the array of addresses
        MPI_Get_address(&array_1d_[i  * this->col_ + (this->col_ - 1)], &this->temp_add_);
        this->displacements_[i] = this->temp_add_;
    }

    // getting the address of the beginning of the data
    MPI_Get_address(this->array_1d_, &this->add_start_);

    // creating a list of displacements
	for (int i = 0; i < this->row_; i++)
		this->displacements_[i] = this->displacements_[i] - this->add_start_;

    // generating the datatype and commiting it
    MPI_Type_create_struct(this->row_, this->blocklengths_, this->displacements_, this->typelist_, &this->receive_type_left);
	MPI_Type_commit(&this->receive_type_left);
}
// ----------------------------------------------------------------------------------------------------------------

// Function that implements the right data type for the receving chunk of data
template<class S>
void DataType<S>::create_right_type_toReceive(void) {

    for (int i = 0; i < this->row_; i++) {

        // assigning length of one to each block of memory to send and
        // defining its type.
        this->blocklengths_[i] = 1;
        this->typelist_[i] = MPI_DOUBLE;

        // getting the address of each element of the left boundary that will be
        // receving the data and placing it in the array of addresses
        MPI_Get_address(&array_1d_[i * this->col_ + 0], &this->temp_add_);
        this->displacements_[i] = this->temp_add_;
    }

    // getting the address of the beginning of the data
    MPI_Get_address(this->array_1d_, &this->add_start_);

    // creating a list of displacements
	for (int i = 0; i < this->row_; i++)
		this->displacements_[i] = this->displacements_[i] - this->add_start_;

    // generating the datatype and commiting it
    MPI_Type_create_struct(this->row_, this->blocklengths_, this->displacements_, this->typelist_, &this->receive_type_right);
	MPI_Type_commit(&this->receive_type_right);
}
// ----------------------------------------------------------------------------------------------------------------


// Function that implements the top data type that will be receiving the data
template<class S>
void DataType<S>::create_top_type_toReceive(void) {

    // getting the address of each element of the bottom boundary that will be
    // receving the data and placing it in the array of addresses
    MPI_Get_address(&array_1d_[(this->row_ - 1) * this->col_], &this->address_);
    MPI_Get_address(this->array_1d_, &this->add_start_);

    // getting the address of the beginning of the data
    this->address_ = this->address_ - this->add_start_;

    // generating the datatype and commiting it
    MPI_Type_create_struct(1, &this->col_, &this->address_, &this->type_, &this->receive_type_top);
	MPI_Type_commit(&this->receive_type_top);
}
// ----------------------------------------------------------------------------------------------------------------


// Function that implements the bottom data type that will be receiving the data
template<class S>
void DataType<S>::create_bottom_type_toReceive(void) {

    // getting the address of each element of the top boundary that will be
    // receving the data and placing it in the array of addresses
    MPI_Get_address(&array_1d_[0], &this->address_);
    MPI_Get_address(this->array_1d_, &this->add_start_);

    // getting the address of the beginning of the data
    this->address_ = this->address_ - this->add_start_;

    // generating the datatype and commiting it
    MPI_Type_create_struct(1, &this->col_, &this->address_, &this->type_, &this->receive_type_bottom);
	MPI_Type_commit(&this->receive_type_bottom);
}
// ----------------------------------------------------------------------------------------------------------------

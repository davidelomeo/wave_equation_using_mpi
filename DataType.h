/*
 * Copyright: Davide Lomeo
 * 
 * Author: Davide Lomeo
 * 
 * Birth: 2021-03-15
 * 
 * Description: header file that contains the class that store the variables and call
 *              the function neede to create the custom MPI datatypes.
 * 
 * Github Repository: https://github.com/davidelomeo/wave_equation_using_mpi
 * 
 * Revision: Last revised by Davide Lomeo on 2022-08-06
*/

#pragma once
#include <iostream>

template<class S>
class DataType {

 public:
    // Constructor
    DataType(int row, int col, double *array);

    // Destructor
    ~DataType();

    // creating a type for the sending communication for each boundary
    void create_left_type_toSend(void);
    void create_right_type_toSend(void);
    void create_top_type_toSend(void);
    void create_bottom_type_toSend(void);

    // creating a type for the receiving communication for each boundary
    void create_left_type_toReceive(void);
    void create_right_type_toReceive(void);
    void create_top_type_toReceive(void);
    void create_bottom_type_toReceive(void);

    // ====================================== Public Variables =======================================
    // The following datatypes are kept publin in the class to be accessible from outside the class

    // generating MPI_Datatype variables for each of the boubdaries
    S send_type_left, send_type_right, send_type_top, send_type_bottom;
    S receive_type_left, receive_type_right, receive_type_top, receive_type_bottom;


    // ====================================== Private Variables ======================================
    // The following varibales do not need to be accessbile from outside the class. For this reason
    // keeping them private will preserve them from accidental overriding.

 private:
    // values set to -1 for debugging purpouses
    int row_ = -1, col_ = -1;

    // null pointer that will be replaced with the values of the input (row/ col) array
    double *array_1d_ = nullptr;

    // integer that will store an array containing the width of each left/ right
    // boundary that will be sent
    int *blocklengths_;

    // variables that will store the datatype of each chunk of data tht will be sent
    S *typelist_, type_ = MPI_DOUBLE;

    // variable that stores the displacement in bytes between non-contiguous memory chunks
    MPI_Aint *displacements_;

    // variables that store the addresses of each chunk of memory that will be sent/ received
    MPI_Aint temp_add_, add_start_, address_;
};

# Solving the wave equation using MPI

The program implements a parallel solver for the wave equation using MPI.

---


## Code Structure
The program is split across six files, five of which form the main program and are written in C++, and one that is a stand-alone post-processing and visualisation program written in Python, along with a small program to compare the results of the post-processing merging and the serial output.

### Main program
The files that compose the main program are:

1. `main.cpp`- file that loads in a .txt file containing the user-defined parameters and controls the flow of the entire program. This is also the only file that will need to be used to compile the program.

2. `Decomposition.cpp/ Decomposition.h` - Source and Header files that perform the domain decomposition, the discretisation, the communications between chunks of the domain and save outputs to files.

3. `DataType.cpp/ DataType.h` - Source and Header files that build custom MPI datatypes in order to allow the communications between the chunks of the domain.

### Post-processing program
The post-processing program constist only of the single python script `post_processing.py`.

### Input parameters
The program has been designed to allow the user to change the parameters of the problem from a textfile without the need to re-compile the program each time the settings are changed.

The file main.cpp loads in the textifle `parameters.txt` located in the same direcory of main.cpp (if the user wants to use a different name, or specifiy a different path to the file, the filename on main will need to be changed accordingly) and stores the parameters in an array of doubles. The set up of the textfile is as follow:

- `grid imax/ jmax sizes` *[indexes 0 and 1]* **-> size of the grid**
- `physical domain sizes` *[indexes 2 and 3]* **-> physical size of the domain**
- `tmax, t0, dt_out` *[indexes 4, 5, 6]*   **-> maximum time, initial time, timesteps to output**
- `c` *[index 7]* **-> speed of the wave**
- `periodic boundaries flag` *[indexes 8]* **-> flag to determine if using periodic or fixed boundaries (fixed boundaries `0` [False] default)**
- `fixed boundaries type` *[indexes 9]* **-> flag to determine if using Dirichlet or Neumann boundary conditions (Neumann `1` [True] default)**

Additionally, the user can decide to have $n$ `initial disturbances` by changin the number at *[index 10]* (11th number on the textfile - `1` is default), followed by $n$ series of three numbers that define the radius `r_splash` of the wave, the x location of the centre `x_splash` and the y location of the centre `y_splash`.

***Example:***<br>
Here there is an example of `parameters.txt` file, were the user setups a 100x100 grid, using fixed Dirichlet boundary conditions and having 3 initial disturbances:

    100  100

    10.0  10.0

    8.0  0.0  0.04

    1

    1

    1

    3

    1.0  3.0  3.0

    2.0  6.0  6.0

    1.0  1.0  1.0


The user **MUST** respect the ordering provided here, otherwise the program will crash/ fail. The read-from-file script is clever enough to only read-in inputs different from spaces  (but will throw an error if it finds anything different from numbers). This means that the spacing between the paramters does not matter, as long as there is at least one space in between them (the textfile has been setup with this sort of layout for cleanliness).


### Post-Processing, Visualisation
The post-processing script was run in Jupyter Notebook due to its easy UI and visusalisation, but the script is also executable from terminal. The program reads in the `domain_parameters.txt` file that was saved by the main program in the folder `array_output`.

The file contains:

- domain row/ cols numbers
- number of nodes used
- number of iterations done (+ 1 -> the grid at t=0)
- grid imax/ jmax sizes
- boundary condition flags
- number of initial disturbance(s)

and these are used to merge together all the outputs from each of the chunks of the domain for each iteration. The resultant merged grid is then saved in the `post_proc_arrays` folder, located in the main diretory of the progrma. Here, the grids are stored with the name `post_proc_<iter>.dat`, where `<iter>` is the number at which iteration that particular grid belongs.

The program ultimately uses the generated .dat files to output the domain in 3D animated plots trhough the use of GIFs.

There are some example GIFs in the folder `images` of this repository. Three files are examples of Dirichlet, Neumann and periodic boundaries using the (default) initial disturbance provided. The other three are also examples of all the boundary conditions, but using three initial disturbances, and using odd and not equal row/ cols sizes (to demonstrate that the code behaves as intended with any given grid size).

***Side note***:
There is a flag named `VISUALISATION` at the beginning of the file that it is set to `True` as default, but if turned to `False` will make the program skip the visualisation. This was thought to be useful if wanting to only run the post-processing for testing purposes (more below).

### Testing
There is only a small python script called `testing.py` that compares the result of the serial code with the result of the post processing.

There are two folders dedicated for pre-saved outputs from serial code for both periodic and fixed boundaries (the serial code used the 3 initial conditions in the paramaters.txt file in this repository). 

The script consist in using np.allclose to compare the outputs of the post-processing to the serial version of the code. 

---
## Program Compiling instructions
### C++ script
The main program used the c++17 library and the code was compiled using mpic++. If using VScode on a Mac, the user is invited to copy the following tasks.json setup in the directory in order to compile the program using  MPI for c++:

```
{
    "version": "2.0.0",
    "tasks": [
        {
            "type": "cppbuild",
            "label": "C/C++: mpic++ build active file",
            "command": "mpic++",
            "args": [
                "-std=c++17",
                "-stdlib=libc++",
                "-g",
                "${file}",
                "-o",
                "${fileDirname}/${fileBasenameNoExtension}"
            ],
            "options": {
                "cwd": "${workspaceFolder}"
            },
            "problemMatcher": [
                "$gcc"
            ],
            "group": {
                "kind": "build",
                "isDefault": true
            },
            "detail": "Task generated by Debugger."
        }
    ]
}
$
```
<br>
The program only needs main.cpp to compile. If compiling from terminal, the command line can simply be (assuming to use the default c++ version):

`mpic++ -g main.cpp `

### Python script
The post-processing can be setup in two ways:

**Using Conda** -> the user will need to type the following in terminal

`conda env create -f environment.yml` 

and activate the environment by typing `conda activate post_proc`

**Using pip** -> the user will need to type the following in terminal:

`pip install -r requirements.txt`

#### GIFs
The package used to convert the 3D animations to GIFs is **ImageMagick**, which was installed on the computer following the instructions on the website: https://imagemagick.org/script/download.php

If using Mac, it is enough to run `brew install imagemagick` and `brew install ghostscript` (provided that brew package manager is installed)

For Windows there is an executable downloadable from the link provided above.

Upon installation of the package, it is not needed for it to be imported in the Python script, as the package becomes part of the methods used by Python to save images.

---
## Running the program
### Main Program
If the main program is compiled using VScode, then it can be run by simply typing in the VScode terminal:

`mpiexec -n <n> main` 

if instead the program is compiled from terminal (not VScode), then it can be run by typing:

`mpiexec -n <n> a.out` (a.out is the standard name of the executables - this can be changed when compiling the code) 

where `<n>` is the number of processors that the user wants to use.

***Side Note:*** <br>
It was noticed that in a very small number of occasions MPI gave out an error for no apparent reasons. The error was the following:

```
--------------------------------------------------------------------------
A system call failed during shared memory initialization that should
not have.  It is likely that your MPI job will now either abort or
experience performance degradation.

  Local host:  host.hostname.com
  System call: unlink(2) /...
  Error:       No such file or directory (errno 2)
--------------------------------------------------------------------------
```
searching on the internet, there was no clear guidance on what the error might be or how to fix it, but the easiest and most effective 'fix' was found to run the program as follow:

`mpiexec -np <n> -oversubscribe main` 

It is still unclear why thie error came up, as it has only been raised at random and very rarely. Nevertheless, even when the errorr occurred, the program still executed and produced outputs as normal.

### Post_processing Program
The post-processing program can be run by simply typing the following from terminal:

`python3 post_processing.py` 

The program assumes that the user wants to have GIFs as output of the animation.

Alternatively, the user can copy the script in a Jupyter Notebook, for example, and un-comment line 188 `HTML(animation.to_jshtml())`. This will show the animation within the Jupyer Notebook, allowing to play, pause and skip frames through the interactive bar within the animation.

---

## Program Outputs and Flags
Every time that the program runs, it generates a number of `.dat` files according to the number of processors used and iterations performed, as well as printing some key information to the terminal.

### Grids outputs
At every iteration, each chunk of the domain write its part of the grid to a .dat file and stores them in the `array_output` folder located in the main directory.

Each file is named after the id of the processor that generated it and at what iteration it belongs, as follow:

`out_iter<iter>_proc_<id>.dat`

this naming convention is key for the post-processing script to run and generate the expected combined grids.

### Image outputs
By default, the `post_processing.py` scripts produces a GIF of the evolution of the domain according to its input parameters (as described above)

The output GIFs are saved in the `image` folder, and thery are in the form of:

`<n>_IDs_<Boundary_type>_<fixed_boundary_type>_<grid_size>.gif`

where `<n>` is the number of initial disturbances, `<Boundary_type>` tells if the boundaries are fixed or periodic, `<fixed_boundary_type>` if Neuman or Dirichlet (when fixed).

The GIFs files can be dragged onto an open browser internet, or through right clicking and selecting 'open with <browser_internet>' and they will play automatically. (The sample GIFS can be visualised diretly through Github)


### Terminal outputs
By default, the program outputs the iteration that it is currently working on to the terminal. This was setup to occur every 50 iterations, and it is needed only for debugging purposes (i.e. determining the progress of the runs). 

Additionally, the program outputs the following to terminal (this is just an example):

```
Output:
Matrix size:    100 x 100
N. Processors:  6
Time taken (s): 0.48208
```

This was done to provide general information about the setup of the problem and the timings of the program.

There are also some additional lines at the end, that are always outputted, but only meaningful if the user has decided to time the central section of the code on its own (i.e. discretisation and the write to file):
```
The following are availebl ONLY if timing the central code:
Central iter. (s): 2.1234e-314
Grid to file (s): 4.9407e-324
```

The reasoning behind this was to determne the timings of the two most computationally expensive sections of the program, that account for almost all the total running time.

### Flags
The program contains a number of flags to allow the user to run/ skip certain part of the code.

1. `main.cpp` contains the flag `DO_TIMING`, that is needed if the user intends to time the whole code and print the results to screen (by default the flag is active). It is only necessary to comment this flag out to skip the timings.

2. `Decomposition.cpp` contains two flags:
    1. `DO_TIMING`, that by default is off, but the user can un-comment it in order to time the discretisation and the write-to-file sections of the code.
    2. `HPC`, that by default is off, is a flag that the user should consider un-commenting if running the code on an HPC system (more info below)

---
## Running on HPC
The program can be run on a HPC system in order to test its performance on a greater number of cores (generally limited for consumer machines) and possibly on larger grid sizes.

As mentioned above, one of the most computationally expensive tasks of the program was found to be the writing the grids to files. It was, therefore, considered sensible to provide to the user the choice to skip this task by adding the `HPC` flag in the `Decomposition.cpp` source file. This way, the timing focuses on the computational efficiency of the program alone.

***Additional Notes:***
It is worth mentioning that the decision to provide the choice to skip the write-to-file was made after attempting to run the code on the DUG system and getting very long timings for very small problems.

Upon discussion with the DUG team, I was adviced that despite the HPC systems have multiple, very fast nodes, their hard drives might result slower than the standar consumer computer. It was then that I decided to include the flag.

### DUG jobs
In order to do a performance analysis of the program, this was copied to my personal folder in the DUG system (I only copied across the .cpp and the .h files and the parameters.txt file). There, the program was compiled using the command `mpic++ -g main.cpp`, but to execute it it was necessary to create several `.job` files that contained a series of bash commands similar to those used in the terminal.

Several .job files were used to run the program with increasingly larger grid sizes. The files used were:

- 100x100.job
- 500x500.job
- 1000x1000.job
- 2000x2000.job

The names reflect the size of the grid used. It is important to mention that with increasing grid size, the t_max in `parameters.txt` had to be reduced by a factor of 8 at every doubling of the size, due to constraines in the maximum working hours allowed.

The content of each of the .job files is virtually the same, and look as follow:

```
#!/bin/bash
#rj nodes=1 queue=iclmsc priority=100 schema=np.schema features= logdir=logs_100x100
set -euo pipefail

echo "Starting MPI example job"
#add an openmpi library to the environment
module add openmpi/4.0.5-mlnx-gcc

echo "Using these nodes:"
echo ${SLURM_JOB_NODELIST}

#Running the binary
mpirun -np $np -oversubscribe a.out

```
The only difference between the files is the `logdir` job command, that indicates the folder in which to save the logs (output) of each job.

Additionally, every job loads in the `np.schema` file, that store the `$np` keys (number of processors to use), and which content is as follow:

```
1 np=1
2 np=2
3 np=4
4 np=8
5 np=16
6 np=32
7 np=64

```

Using a 'helper' `.schema` file allowed to run seven jobs simultaneously, assigning to each a different `np` key, but with the advantage of keeping the same id for the running 'family' of jobs, and only adding `.<n>` (where `<n>` is the job number as per `np.schema` file) at the end of each `*.o`log file 

An example is:
```
100x100.o51595936.1
100x100.o51595936.2
100x100.o51595936.3
...
```
This format indicates that each file was generated from the same job (100x100.job), using different number of processors (e.g. `*.o*.3` indicates that the logfile relates to the job that used `np=4` processors.

To run the `.job` files it was enough to type:

`rjs <jobname>.job`

The log of the jobs were stored in folders named:

- logs_100x100
- logs_500x500
- logs_1000x1000
- logs_2000x2000

In order to have a substantial set of timings, each job was run three times, such that each of the folders stated above containes 21 `*.o*` files.
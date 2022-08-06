"""
This program is needed to check if the output of the post-processing matches
the output of the serial code
"""

import numpy as np


# loading the number of iterations from the domain_parameters .dat file
domain_parameters = np.loadtxt('./array_output/domain_parameters.dat')
n_iter = int(domain_parameters[3])

# checking the flag (periodic/ non-periodic)
if int(domain_parameters[6]) == 1:
    # if periodic, then the code is checked against the result of the
    # periodic serial
    for i in range(1, n_iter):
        serial = np.loadtxt(
            './serial_periodic_output/out_iter_{0:d}_proc_0.dat'.format(i))
        parallel = np.loadtxt(
            './post_proc_arrays/post_proc_{0:d}.dat'.format(i))
        print(np.allclose(serial, parallel) is True)

else:
    # if non-periodic, then the code is checked against the result of the
    # non-periodic serial
    for i in range(n_iter):
        serial = np.loadtxt(
            './serial_nonperiodic_output/out_iter_{0:d}_proc_0.dat'.format(i))
        parallel = np.loadtxt(
            './post_proc_arrays/post_proc_{0:d}.dat'.format(i))
        print(np.allclose(serial, parallel) is True)

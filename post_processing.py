"""
This program performs the post-processing needed to merge the individual
outputs of the chunks of the domain into sigle files and visualises the result.
"""

# ================================== Imports ==================================
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# import only needed if the user wants to use Jupyter Notebooks
# to visualise the animation
# from IPython.display import HTML

# increasing the maximum memory available to matplotlib to create an animation
from matplotlib import rcParams as pars
pars['animation.embed_limit'] = 2**128

# allowing the animation to be converted to a gif
plt.rcParams['animation.convert_path'] = 'convert'
# =============================================================================


# =============================================================================
# flag to change to false if NOT wanting to output images
VISUALISATION = True
# =============================================================================


# ========================= Post-Processing Function ==========================
def post_processing(n_iter, n_procs, rows, cols):
    '''The function uses the number of column/ row and processors
    used in the domain decomposition to determine how to concatenate
    the chunks of the domain that each prcessor was assigned.
    The script loads the .dat files stored in 'array_output' folder
    and assign a processor number and iteration to it, following
    the same naming convenction used in the main program to save
    the files.
    At each iteration, the program concatenates the .dat files
    column-wise first, and row-wise after. This approach allows to
    deal with any input i,j number, including single processors

    Parameters.
    n_iter:   number of iterations
    n_procs:  number of processors used in the domain decomposition
    rows:     number of rows of the decomposition
    cols:     number of cols of the decomposition

    Output.
    A list of np.arrays, each defining the combination of all the
    chunks that the domain was split into at each iteration. The
    file is saved as 'post_proc_<iter>' in the 'post_proc_arrays'
    folder.
    '''

    # empty lists to store the merged chunk of the domain as columns
    # and the fully merged arrays
    domain_col = []
    full_domain = []

    # storing a string that will be used to identify the processor and
    # and iteration within the main loop
    file_path = './array_output/out_iter_{:d}_proc_{:d}.dat'

    # iterating n_iter times
    for i in range(n_iter):
        for col in range(cols):
            # setting empty array
            array_0 = np.empty([])

            for row in range(rows):
                # getting the id of the node of the current chunk of the domain
                proc_id = row * cols + col

                # loading the .dat files. If row == 0, then the current chunk
                # of the domain will act as pivot, and all the remaining chunks
                # will be appended column-wise
                if row == 0:
                    array_0 = np.loadtxt(file_path.format(i, proc_id))
                else:
                    array = np.loadtxt(file_path.format(i, proc_id))

                    # concatenating the current chunk with the pivot
                    array_0 = np.concatenate((array_0, array), axis=0)

            # appending the column to the empty list
            domain_col.append(array_0)

        if cols != 1:
            # concatenating arrays row-wise
            full_domain.append(
                np.concatenate(
                    (domain_col[-2], domain_col[-1]),
                    axis=1))
        else:
            # assignign the column as the full domain if the domain
            # length = 1 across
            full_domain = domain_col
            continue

    return full_domain
# =============================================================================


# ========================= Visualisation Function ============================
def frame_plotting(frame_number, Z, plot):
    '''Function that outputs each individual timeframe from the
    input animation.
    '''
    ax.clear()
    Z = np.loadtxt(
        './post_proc_arrays/post_proc_{0:d}.dat'.format(frame_number))
    plot = ax.plot_surface(X, Y, Z, **plot_args)
    ax.set_zlim(-4, 4)
    ax.set_xlabel('$x$', fontsize=14)
    ax.set_ylabel('$y$', fontsize=14)
    ax.set_zlabel('$z$', fontsize=14)
    title = ('Domain Evolution - ' +
             boundaries_flag + ' ' +
             fixed_bound_type + ' - iter {0:d}')
    ax.set_title(title.format(frame_number), fontsize=16)
    return plot,
# =============================================================================


# ============================= Post_processing ===============================
# loading problem parameters from .dat file
domain_parameters = np.loadtxt('./array_output/domain_parameters.dat')
rows = int(domain_parameters[0])
cols = int(domain_parameters[1])
n_procs = int(domain_parameters[2])
n_iter = int(domain_parameters[3])
imax = int(domain_parameters[4])
jmax = int(domain_parameters[5])
ini_dist = int(domain_parameters[8])

# if statements to identify the type of boundaries used in the domain.
# This will then form the tile of the plot.
boundaries_flag = ""
fixed_bound_type = ""
if int(domain_parameters[6]) == 1:
    boundaries_flag = "Periodic Boundaries"
else:
    boundaries_flag = "Fixed Boundaries"
    if int(domain_parameters[7]) == 1:
        fixed_bound_type = "Neumann"
    else:
        fixed_bound_type = "Dirichlet"

# printing out the values for debugging purposes
print("rows:", rows, " cols:", cols, " procs:", n_procs,
      " n_iter:", n_iter, " imax:", imax, " jmax:", jmax,
      " - ", boundaries_flag, fixed_bound_type, " - IDs:", ini_dist)

# storing the list obtained from the post-processing function
domain_evolution = post_processing(n_iter, n_procs, rows, cols)

# writing .dat files for each iteration for future plotting purposes
for i in range(n_iter):
    with open('./post_proc_arrays/post_proc_{0:d}.dat'.format(i), 'w') as data:
        np.savetxt(data, domain_evolution[i], delimiter=' ')
# =============================================================================


# ============================ Visualisation ==================================

if VISUALISATION:
    # setting up the mesh size and loading the domain at t=0
    X = np.linspace(0, 100, imax)
    Y = np.linspace(0, 100, jmax)
    X, Y = np.meshgrid(Y, X)
    Z = np.loadtxt('./post_proc_arrays/post_proc_0.dat')

    # setting up plot parameters
    plot_args = {'rstride': 1, 'cstride': 1, 'cmap': 'coolwarm',
                 'linewidth': 0., 'antialiased': True, 'color': 'w'}

    # creating figure and setting 3d projection
    fig = plt.figure(figsize=(12, 12))
    ax = fig.add_subplot(111, projection='3d')

    # plotting the figure at t=0
    plot = ax.plot_surface(X, Y, Z, **plot_args)
    grid_size = str(imax) + " x " + str(jmax)
    plt.suptitle(grid_size + " matrix", fontsize=14)
    plt.close()

    # compute the animation for n iterations with
    animation = animation.FuncAnimation(fig, frame_plotting,
                                        frames=np.arange(0, n_iter),
                                        fargs=(Z, plot), interval=100,
                                        blit=True)

    # outputting the animation - Uncomment this line if wanting to show the
    # animation with a control bar within the Jupyter Notebook
    # HTML(animation.to_jshtml())

    # saving the animation to a GIF
    filename = (str(ini_dist) + ' IDs ' + boundaries_flag +
                ' ' + fixed_bound_type + ' ' +
                grid_size + '.gif').replace(" ", "_")
    animation.save('./images/' + filename,
                   writer='imagemagick', fps=80)
# =============================================================================

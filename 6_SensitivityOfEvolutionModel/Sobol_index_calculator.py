"""
Sobol Sensitivity as implemented in SAlib - Squid fluids model
@author: Steve
"""

## Import the needed libraries and functions
# Set the required public imports
import numpy as np
import matplotlib.pyplot as plt
from SALib import ProblemSpec
#from SALib.analyze import sobol
#from SALib.sample.sobol import sample

#from SALib.sample.sobol import sample

bounds = [[-1,1],[-1,1],[-1,1],[-1,1]]

## Define the problem dictionary for SaLib and get a collection of samples
# problem = {
#   'num_vars': 4,
#   'names': [r'$a$',r'$b$',r'$c$',r'$d$'],
#   'groups': None,
#   'bounds': bounds
# }

sp = ProblemSpec({
        "names": ["x1", "x2", "x3", "x4"],
        "groups": None,
        "bounds": [[-1,1]] * 4,
        "outputs": ["Y"],
    })

params_list = [r'$S_0$',r'$A_0$',r'$\chi$',r'$\sigma$']

# Nparams = 2**10 # Number of samples to take
# param_values = sample(sp, N=Nparams, calc_second_order=True) # Get parameter sets
# np.savetxt("data_outputs/params.txt", param_values) # Save the parameter values to a file

# ## Solve the system across the sample set to get the values
param_values = np.loadtxt("6_SensitivityOfEvolutionModel/params.txt") # Load the parameter values
QOI = np.loadtxt("6_SensitivityOfEvolutionModel/output_values.txt") # Preallocate storange for the quantity of interest in each run

# Provide the results to the interface
sp.set_samples(param_values)
sp.set_results(QOI)

# ## Perform the sobol index calculations
sp.analyze_sobol( calc_second_order=True, print_to_console=True)
#sp.analyze_sobol( calc_second_order=False, print_to_console=True)
print(np.mean(QOI))
print(np.var(QOI))

if False:

    ## Plot the outputs
    bar_width=.75 # Bar width
    # Plot the individual variances for each parameter
    plt.figure()

    COLOR1 = [0, 0.4470, 0.7410]
    COLOR2 = [0.8500, 0.3250, 0.0980]

    plt.bar(2*np.arange(len(params_list))-1.5*bar_width,sp.analysis['S1'],bar_width,  color=COLOR1)
    plt.bar(2*np.arange(len(params_list))-.5*bar_width,sp.analysis['ST'], bar_width, color=COLOR2)
    plt.legend(['First Order','Total Order'])
    plt.xticks(2*np.arange(len(params_list))-bar_width, params_list)
    plt.ylabel('Sobol Index',fontsize=16)
    plt.show()

error_bars = True

if error_bars:

    bar_width = 0.4
    COLOR1 = [0, 0.4470, 0.7410]
    COLOR2 = [0.8500, 0.3250, 0.0980]

    # Create the plot
    fig, ax = plt.subplots()

    # Plot bars with error bars
    x = np.arange(len(params_list))  # Positions for bars
    ax.bar(x - bar_width/2, sp.analysis['S1'], bar_width, label='First Order', color=COLOR1)
    ax.bar(x + bar_width/2, sp.analysis['ST'], bar_width, label='Total Order', color=COLOR2)

    # Add error bars
    yerr_S1 = np.array(sp.analysis['S1_conf']).T  # Transpose for easier indexing
    yerr_ST = np.array(sp.analysis['ST_conf']).T

    ax.errorbar(x - bar_width/2, sp.analysis['S1'], yerr=yerr_S1, fmt='none', ecolor='black', capsize=7)
    ax.errorbar(x + bar_width/2, sp.analysis['ST'], yerr=yerr_ST, fmt='none', ecolor='black', capsize=7)

    ax.axhline(0.05, color='y', linestyle='--') 

    # Set labels and title
    plt.ylim(0, 0.9)
    plt.xticks(x, params_list)
    plt.xlabel('Parameters', fontsize=16)
    plt.ylabel('Sobol Index', fontsize=16)
    plt.title('Sobol Indices with Confidence Intervals', fontsize=18)

    # Add legend
    #plt.legend()

    # Rotate x-axis labels for better readability if many parameters
    # if len(params_list) > 5:
    #     plt.xticks(rotation=45)

    # Show the plot
    plt.tight_layout()
    plt.show()

# # Plot the total variances for each parameter
# plt.figure()
# plt.bar(np.arange(len(params_list))-.5*bar_width,sp.analysis['ST'], bar_width)
# plt.xticks(np.arange(len(params_list)), params_list)
# plt.ylabel('Total Sobol Index',fontsize=16)
# plt.show()

# #

# sobol_indices = sp.analysis['S2']

# # Assuming 'sobol_indices' is your matrix of 2nd-order Sobol indices
# plt.imshow(sobol_indices, cmap='YlGnBu', interpolation='nearest')
# plt.colorbar(label='Sobol Index')

# # Set labels and title
# plt.title('2nd-Order Sobol Indices')
# plt.xlabel('Variable 2')
# plt.ylabel('Variable 1')

# # Add annotations
# for i in range(len(sobol_indices)):
#     for j in range(len(sobol_indices[0])):
#         plt.text(j, i, round(sobol_indices[i, j], 3),
#                  ha="center", va="center",
#                  color="black")

# plt.show()

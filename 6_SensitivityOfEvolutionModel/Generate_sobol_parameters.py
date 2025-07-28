""": 
Sobol Sensitivity as implemented in SAlib - Squid fluids model
@author: Steve
"""

## Import the needed libraries and functions
import numpy as np
import matplotlib.pyplot as plt
from SALib import ProblemSpec
from SALib.analyze import sobol
from SALib.sample.sobol import sample

sp = ProblemSpec({
        "names": ["x1", "x2", "x3", "x4"],
        "groups": None,
        "bounds": [[-1,1]] * 4,
        "outputs": ["Y"],
    })

Nparams = 2**12 # Number of samples to take
param_values = sample(sp, N=Nparams, calc_second_order=True) # Get parameter sets
np.savetxt("6_SensitivityOfEvolutionModel/params.txt", param_values) # Save the parameter values to a file

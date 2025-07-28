To make use of threads in MatLab the code uses multiple languages.
To avoid complications with executing MatLab within python script the files are kept isolated.
The following order should be used to produce the data.

In order to calculate the Sobol indices for the model...
1. Run the script Generate_sobol_parameters.py
2. Run the script calculateOutputs.m
3. Run the script Sobol_index_calculator.py

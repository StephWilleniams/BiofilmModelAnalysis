
# Import modules
include("ModelAnalysis.jl") 
import .ModelAnalysis

using DifferentialEquations, Distributions, Plots, Random, Interpolations

# Random seed
Random.seed!(13245)

N_synthetics = 3 # Number of synthetic data sets to generate
μ = 1.0 # Mean of the noise in the measurements/synthetic data
σ = 0.0025 # Standard deviation of the noise in the measurements/synthetic data
dist = LogNormal(ModelAnalysis.μ_for_mu(μ, σ), ModelAnalysis.σ_for_sigma(μ, σ)) # Log-normal distribution
PL_range = 0.29 # Range of the profile likelihood optimisation
n_PLs = 101 # Number of points in the profile likelihood optimisation
param_range = LinRange(1-PL_range,1+PL_range,n_PLs) # Univariate parameter
initial_varying_params = [0.05,0.05,0.05] # Initial values of the varying parameters
plot_colors = ["#0072BD","#D95319","#77AC30"] # Colors for sample frequency Plot
plot_styles = [:solid,:dash,:dot] # Styles for sample frequency Plot

## Constants
# Non-dimensionalised scaling factors
ND_T = 24.0 # Non dimension time - hours
ND_C =  1.0 # Non dimension concentration - ug/ml
Vol  =  3.0 # Volume of the dimensional system - ml
Surf =  8.0 # Surface area of the dimensional system - cm^2
# Efficiencies
e23 = 0.2 # Bacteria growth efficiency
e4  = 0.5 # Free predator growth efficiency
e5  = 0.33 # Biofilm predator growth efficiency
# Growth rates
r2 = (1/e23) *0.21  *ND_T # Free bacteria growth rate 
r3 = (1/e23) *0.007 *ND_T # Bound bacteria growth rate 
r4 = (1/e4)  *0.12  *ND_T # Predator growth rate 
r5 = (1/e5)  *0.09  *ND_T # Predator death rate 
# Half saturation constants
H23 = 1.0*Vol # Half saturations of bacteria 
H4  = 1.0*Vol # Half saturation of free predator
H5  = 0.1*Surf # Half saturation of bound predator
# Attatchment parameters
chi_on_max           = (Surf/Vol)*0.1*ND_T    # Attatchment parameters
chi_on_min           = (Surf/Vol)*0.0005*ND_T # Attatchment parameters
interaction_strength = Surf*0.01        # Interaction strength
chi_off              = 0.005*ND_T       # Dettachment rate
# ODE parameters
tspan        = (0.0, 1.0) # Time span of solution
output_times = 0.0:0.01:1.0 # Times at which to output the solutions
N_times      = length(output_times) # Number of output times
ode_algo     = Tsit5() # Non-stiff solver
# Initial conditions
C0 = 1.0 # Initial concentration of Free bacteria
F0 = 1.0 # Initial concentration of Free predator
B0 = 0.01 # Initial concentration of Biofilm bacteria
S0 = 0.01*Vol # Initial concentration of Ciliate predator
A0 = 0.001*Surf # Initial concentration of Biofilm predator
u0 = [C0,F0,B0,S0,A0]
# Model parameters
parameters = [r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, chi_on_max, chi_on_min, interaction_strength, chi_off]

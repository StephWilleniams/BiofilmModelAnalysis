# Import modules
include("ModelAnalysis.jl") 
import .ModelAnalysis

using DifferentialEquations, Distributions, Plots, Random

# Random seed
Random.seed!(13245)

# Define constants
t_start = 0.0 # Start time of the schedule
t_end = 1.0 # End time of the schedule
dt = 0.01 # Sample schedule minimum time step
n_times = 51 # Number of time points in the schedule
N_synthetics = 3 # Number of synthetic data sets (experimental repeat analog) to use for the profile likelihood optimisation
N_organisms = 100 # Number of organisms (schedules) in the population for selection
N_generations = 3 # Number of generations to run the genetic algorithm
cloning_noise = 9*dt/N_generations # Cloning mutation noise level
PL_range = 0.30 # Range of the profile likelihood optimisation
n_PLs = 101 # Number of points in the profile likelihood optimisation
initial_varying_params = [0.05,0.05,0.05] # Initial values of the varying parameters

# Preallocate organisms
global organisms = Vector{ModelAnalysis.organism}(undef,N_organisms)

# Define the model
ode_algo = Tsit5()
ND_T = 24.0 # Non-dimension time - hours
ND_C = 1.0 # Non-dimension concentration - ug/ml
Vol = 3.0 # Volume of the dimensional system - ml
Surf = 8.0 # Surface area of the dimensional system - cm^2
e23 = 0.2 # Free bacteria growth efficiency
e4 = 0.5 # Bound bacteria growth efficiency
e5 = 0.33 # Predator growth efficiency
r2 = (1/e23)*0.21*ND_T # Free bacteria growth rate - cells/hr * ND_hrs
r3 = (1/e23)*0.007*ND_T # Bound bacteria growth rate - cells/hr * ND_hrs
r4 = (1/e4)*0.12*ND_T # Predator growth rate - cells/hr * ND_hrs
r5 = (1/e5)*0.09*ND_T # Predator death rate - 1/hrs * ND_hrs
H23 = Vol*1.0 # Half saturations of bacteria
H4  = Vol*1.0 # Half saturation of free predator
H5  = Surf*0.1 # Half saturation of bound predator
chi_on_max = (Surf/Vol)*0.05*ND_T # Attatchment parameters
chi_on_min = (Surf/Vol)*0.0005*ND_T # Attatchment parameters
interaction_strength = Surf*0.01 # Interaction strength
chi_off = 0.005*ND_T
theta_baseline = [r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, chi_on_max, chi_on_min, interaction_strength, chi_off] # Parameters
C0 = 1.0 # Initial concentration of Free bacteria
F0 = 1.0 # Initial concentration of Free predator
B0 = 0.01 # Initial concentration of Biofilm bacteria
S0 = 0.01*Vol # Initial concentration of Ciliate predator
A0 = 0.001*Surf # Initial concentration of Biofilm predator
u0 = [C0,F0,B0,S0,A0]

# Define synthetic noise
μ = 1.0 # Mean of the noise in the measurements/synthetic data
σ = 0.0025 # Standard deviation of the noise in the measurements/synthetic data
dist = LogNormal(ModelAnalysis.μ_for_mu(μ, σ), ModelAnalysis.σ_for_sigma(μ, σ)) # Log-normal distribution for the noise in the measurements/synthetic data

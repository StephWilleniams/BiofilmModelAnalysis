# Calculate the Sobol sensitivity indices for the BFP model
# Author: Stephen Williams 

##--------------------------------------------------##

using DifferentialEquations, Plots, GlobalSensitivity, QuasiMonteCarlo, LaTeXStrings, CurveFit, Statistics
MatlabBlue = "#0072BD"
MatlabRed = "#D95319"

function BFP_model(du, u, params, t)
    # Extract the different species
    C, F, B, S, A = u
    # Extract parameters
    chi_on_max, chi_on_min, interaction_strength, chioff = params

    # Re-fix parameters not of interest to the sensitivity
    e23 = 0.2
    e4 = 0.5
    e5 = 0.33
    r2 = (1/e23)*0.21*ND_T
    r3 = (1/e23)*0.007*ND_T
    r4 = (1/e23)*0.12*ND_T
    r5 = (1/e23)*0.09*ND_T
    H23 = 3.0*1.0/ND_C # Half saturations of bacteria
    H4  = 3.0*1.0/ND_C # Half saturation of free predator
    H5  = 8.0*0.1/ND_C # Half saturation of bound predator

    CHfunct = C/(C+H23)
    FHfunct = F/(F+H4)
    BHfunct = B/(B+H5)

    chion = (chi_on_max*interaction_strength + chi_on_min*B)/(interaction_strength + B)
    # Calculate the ODEs
    du[1] = - r2*CHfunct*F - r3*CHfunct*B
    du[2] = + e23*r2*CHfunct*F - r4*FHfunct*S - chion*F + chioff*B
    du[3] = + e23*r3*CHfunct*B - r5*BHfunct*A + chion*F - chioff*B
    du[4] = + e4*r4*FHfunct*S 
    du[5] = + e5*r5*BHfunct*A
end

# Define a function of one input set only for simplicity - used for the sensitivity analysis
function model(theta)
    # Extract info from the input theta
    u0 = [1.0, 1.0, 0.01, 3.0*0.01, 8.0*0.001] # Initial conditions
    params = theta
    # Fixed values needed for the solver to work
    output_times = LinRange(0,1,100) # Solver times
    tspan = (output_times[1],output_times[end]) # End-points for the solver time
    # Run the solver for the system
    algo = Tsit5() # Optional stiff solver, required for optimisation
    prob = ODEProblem(BFP_model, u0, tspan, params) # Define the ODE problem
    sol = solve(prob, algo, saveat=output_times) # Solve the ODE problem
    # Set the parameter of interest for our outputs to check sensitivity with respect to
    QOI = sol[3,end]
    # Return the quantity of interest
    return QOI

end

### Constants
# Non-dimensionalised scaling factors
ND_T = 24.0 # Non dimension time - hours
ND_C = 1.0 # Non dimension concentration - ug/ml
Vol = 3.0 # Volume of the dimensional system - ml
Surf = 8.0 # Surface area of the dimensional system - cm^2
# Attatchment parameters
chi_on_max = (Surf/Vol)*0.05*ND_T # Attatchment parameters
chi_on_min = (Surf/Vol)*0.0005*ND_T # Attatchment parameters
interaction_strength = Surf*0.01 # Interaction strength
chi_off = 0.005*ND_T # Detachment rate
# Set the true values of the parameters
theta_true = (chi_on_max,chi_on_min,interaction_strength,chi_off) # Parameters
theta_names = ["chi_max","chi_min","a","chi_off"]

### Define and sensitivity problem parameters
delta = 0.9
N_samples = 20000
lb = [(1-delta)*theta_true[i] for i in eachindex(theta_true)]
ub = [(1+delta)*theta_true[i] for i in eachindex(theta_true)]
A, B = QuasiMonteCarlo.generate_design_matrices(N_samples, lb, ub, QuasiMonteCarlo.SobolSample())
sob_res = gsa(model, Sobol(nboot=100), A, B, batch=false)

## Plot the results
plt = plot()
bar!(plt,2:2:8,sob_res.ST, yerror = sob_res.ST_Conf_Int, legend=false, lw = 1, bar_width=0.8, color="#D95319")
bar!(plt,1:2:7,sob_res.S1, yerror = sob_res.S1_Conf_Int, legend=false, lw = 1, bar_width=0.8, color="#0072BD")
display(plt)

# savefig(plt, "SobolIndices.png")

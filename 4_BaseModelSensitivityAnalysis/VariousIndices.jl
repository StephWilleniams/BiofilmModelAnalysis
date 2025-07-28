# Calculate various sensitivity indices for the BFP model
# Author: Stephen Williams 

##--------------------------------------------------##

using DifferentialEquations, Plots, GlobalSensitivity, QuasiMonteCarlo, LaTeXStrings, CurveFit, Statistics

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
    u0 = [ 1.0, 1.0, 0.01, 3.0*0.01, 8.0*0.001 ] # Fixed values
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
theta_true = (chi_on_max, chi_on_min, interaction_strength, chi_off) # Parameters
theta_names = ["chi_max", "chi_min", "a", "chi_off"]

### Define and sensitivity problem parameters
param_index = 1 # Index of the parameter to plot
delta = 0.9 # Range of parameter values checked (1-delta, 1+delta)
lb = [(1-delta)*theta_true[i] for i in eachindex(theta_true)]
ub = [(1+delta)*theta_true[i] for i in eachindex(theta_true)]
N_samples = 40000 # Number of samples to use in the sensitivity analysis

## Calculate the sensitivity indices
X = QuasiMonteCarlo.sample(N_samples, lb, ub, QuasiMonteCarlo.SobolSample())
Y = reshape(model.([X[:, i] for i in 1:N_samples]), 1, N_samples)
reg_mat = gsa(X, Y, RegressionGSA(true))
xpts = X[param_index,:]/theta_true[param_index]
c,m = linear_fit(xpts, Y[1,:])

## Output the resulting indices
println("mean: ", mean(Y[1,:]))
#println(std(Y[1,:]), " ", m*sqrt(N_samples)/std(Y[1,:]))
println(reg_mat.pearson[1,param_index])
println(reg_mat.partial_correlation[1,param_index])
println(reg_mat.pearson_rank[1,param_index])
println(reg_mat.partial_rank_correlation[1,param_index])
println("intercept: ", c, " slope: ", m)

## Plot the result with a linear fit
# plt = plot()
# scatter!(plt, xpts,Y[1,:], alpha=0.1, label = false)
# display(plt)
# plot!(plt, xpts , m.*xpts .+ c , linewidth = 5, label = false)

#savefig(plt, "4_BaseModelSensitivityAnalysis/BFP_sensitivity_$(param_index).png")

# Probability districution of the output
# plt = histogram(Y[1,:], bins=100, normalize=:pdf, label = false, xlabel = "Quantity of Interest", ylabel = "Frequency", title = "Distribution of the Quantity of Interest")
# display(plt)

#savefig(plt, "4_BaseModelSensitivityAnalysis/histogram.png")

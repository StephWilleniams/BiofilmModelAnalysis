# Plot the profile likelihood of the non-dimensionalised BFP model
# Author: Stephen Williams 

##--------------------------------------------------##

# Define constants and modules
include("preamble.jl")

## Script variables
index_of_interest = 14 # Index of the parameter of interest
varying_indices = filter!(x -> x != index_of_interest, [11,12,13,14]) # Indices of the varying parameters
N_sample_times = [11,21,51] # Number of uniformly sampled times

## Plot preallocate
plt = hline([-quantile(Chisq(1),0.95)/2],label=false,lw=6,lc = :gold, xlims=(param_range[1],param_range[end]), ylims=(-7.5,0.3)) 

## Loop through cases
for i = eachindex(N_sample_times)
    global current_N_sample_times = N_sample_times[i]
    global organism = ModelAnalysis.organism(LinRange(0.0,1.0,current_N_sample_times),1,0.0) # Sampling schedule
    global pop = ModelAnalysis.SolvePopulation(ModelAnalysis.BFP_model, u0, organism, parameters, ode_algo)
    global synthetic_data = ModelAnalysis.SyntheticPopulation_fixIVP(pop, dist, N_synthetics)
    global output, theta_minus, theta_plus, width = ModelAnalysis.PLB(index_of_interest, PL_range, n_PLs, ModelAnalysis.BFP_model, ode_algo, u0, dist, organism, parameters, initial_varying_params, varying_indices, synthetic_data) 
    global itp = interpolate(output .- quantile(Chisq(1),0.95)/2, BSpline(Cubic()))
    println("Width with $(current_N_sample_times) sample times: $(width), $(theta_minus), $(theta_plus)")
    plot!(plt, param_range, itp(1:n_PLs), color = plot_colors[i], ls = plot_styles[i], label=false, lw=5)
end

## Final plot
vline!(plt, [1.0], label=false, lw=6, lc = :black, alpha=0.5)
display(plt)

# savefig("2_PracticalIdentifiability/$(index_of_interest).png") 

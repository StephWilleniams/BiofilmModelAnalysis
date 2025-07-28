# Perform the genetic algorithm to find the best sampling for the BFP model
# Author: Stephen Williams 

##--------------------------------------------------##
 
# Define constants and modules
include("preamble.jl")

index_of_interest = 14 # Index of the parameter of interest, takes values 11 to 14
varying_indices = filter!(x -> x != index_of_interest, [11,12,13,14]) # Indices of the varying parameters
for i in 1:N_organisms
    global organisms[i] = ModelAnalysis.organism(ModelAnalysis.GenerateInitialSchedule(t_start,t_end,dt,n_times),i,0.0)
end
for gen = 1:N_generations
    for i in 1:N_organisms
        local synthetic_data = ModelAnalysis.SyntheticPopulation(ModelAnalysis.SolvePopulation(ModelAnalysis.BFP_model, u0, organisms[i], theta_baseline, ode_algo), dist, N_synthetics)
        local output, root1, root2, width = ModelAnalysis.PLB(index_of_interest, PL_range, n_PLs, ModelAnalysis.BFP_model, ode_algo, u0, dist, organisms[i], theta_baseline, initial_varying_params, varying_indices, synthetic_data)
        global organisms[i].fitness = width
        println("Generation: ", gen, " Organism: ", i, " Fitness: ", organisms[i].fitness)
    end
    global organisms = ModelAnalysis.OrderOrganisms(organisms)
    global organisms = ModelAnalysis.CloneOrganisms(organisms,cloning_noise,dt,t_start,t_end)   
end
println("Fitest $(index_of_interest)_$(n_times) organism: ", organisms[1].fitness)

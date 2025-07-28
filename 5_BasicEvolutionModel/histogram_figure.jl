# Plot the districution of chion values in the BFP model across 500 generations
# Author: Stephen Williams 

##--------------------------------------------------##

# Define the model
include("preamble.jl")

plt = plot()

global S0 = 0.0
global A0 = 2.0
global u0 = vcat(C0,F0,B0,S0,A0)

for gen = 1:N_generations
    local prob = ODEProblem(BFP_compartment_model, u0, (0.0,1.0), theta)
    local sol = solve(prob, Tsit5())
    local NofT = sol(output_times)
    local F0 =       NofT[ind_B,end]./sum(NofT[ind_B,end])
    local B0 = 0.001*NofT[ind_B,end]./sum(NofT[ind_B,end])
    global u0 = vcat(C0,F0,B0,S0,A0)
    if gen % 50 == 0
        plot!(plt, chionVals, NofT[ind_B,end]./sum(NofT[ind_B,end]*bw), lw=3, legend=false, color=colors[gen],xlims=(2,12))
    end
end

display(plt)


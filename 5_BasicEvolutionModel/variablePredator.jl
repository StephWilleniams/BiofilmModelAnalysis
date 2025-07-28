# Plot the mean chion values in the BFP model across multiple predator cases
# Author: Stephen Williams 

##--------------------------------------------------##

# Define the model
include("preamble.jl")

plt = plot()

for case = 0:N_predCases
    println("Case: ", case)
    global S0 = pred0max*case/N_predCases
    global u0 = vcat(C0,F0,B0,S0,A0)
    local avgChi = zeros(N_generations)
    for gen = 1:N_generations
        local prob = ODEProblem(BFP_compartment_model, u0, (0.0,1.0), theta)
        local sol = solve(prob, Tsit5())
        local NofT = sol(output_times)
        local F0 = NofT[ind_B,end]./sum(NofT[ind_B,end])
        local B0 = 0.001*NofT[ind_B,end]./sum(NofT[ind_B,end])
        global u0 = vcat(C0,F0,B0,S0,A0)
        local alpha = 1/sum(NofT[ind_B,end]*bw)
        avgChi[gen] = alpha * bw * sum(chionVals .* NofT[ind_B,end])
    end
    if case == 0
        plot!(plt, avgChi, label=false, lw=5, color = colors2[case+1], style = :dash)
    else
        plot!(plt, avgChi, label=false, lw=5, color = colors2[case+1], style = :solid)
    end
end

# savefig(plt, "5_BasicEvolutionModel/variable_S_pred.png")

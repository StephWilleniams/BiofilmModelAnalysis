# Plot a single generation of the evolution model
# Author: Stephen Williams 

##--------------------------------------------------##

# Define the model
include("preamble.jl")

prob = ODEProblem(BFP_compartment_model, u0, (0.0,1.0), theta)
sol = solve(prob, Tsit5())
NofT = sol(output_times)
plt = plot()
az, el = 50, 30

# Define the ranges for x and y
x = output_times
y = chionVals
X = [i for i in x, j in y]
Y = [j for i in x, j in y]
Z = NofT[ind_F,:]
plot!(plt,x, y, Z, st = :surface, color = colors[end], alpha = 0.3, camera = (az,el) ,cbar = false)

x = output_times
y = chionVals
X = [i for i in x, j in y]
Y = [j for i in x, j in y]
Z = NofT[ind_B,:]
plot!(plt,x, y, Z, st = :surface, color = colors[1], alpha = 0.25, camera = (az,el) ,cbar = false) 

display(plt)

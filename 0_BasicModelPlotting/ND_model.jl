# Plot the solution of the non-dimensionalised BFP model
# Author: Stephen Williams 

##--------------------------------------------------##

### Preamble
using DifferentialEquations, Plots

### Define the ODE system
function BFP_model(du, u, theta, t)
    C, F, B, S, A = u 
    r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, chi_on_max, chi_on_min, interaction_strength, chioff = theta 
    CHfunct = C/(C+H23) 
    FHfunct = F/(F+H4) 
    BHfunct = B/(B+H5) 
    chion = (chi_on_max*interaction_strength + chi_on_min*B)/(interaction_strength + B) 
    du[1] = - r2*CHfunct*F - r3*CHfunct*B 
    du[2] = + e23*r2*CHfunct*F - r4*FHfunct*S - chion*F + chioff*B 
    du[3] = + e23*r3*CHfunct*B - r5*BHfunct*A + chion*F - chioff*B
    du[4] = + e4*r4*FHfunct*S 
    du[5] = + e5*r5*BHfunct*A
end

### Constants
# Non-dimensionalised scaling factors
ND_T = 24.0 # Non dimension time - hours
ND_C =  1.0 # Non dimension concentration - ug/ml
Vol  =  3.0 # Volume of the dimensional system - ml
Surf =  8.0 # Surface area of the dimensional system - cm^2
# Efficiencies
e23 = 0.2 # Free bacteria growth efficiency
e4  = 0.5 # Bound bacteria growth efficiency
e5  = 0.33 # Predator growth efficiency
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
ode_algo     = Tsit5() 
# Initial conditions
C0 = 1.0 # Initial Carbon source
F0 = 1.0 # Initial Planktonic bacteria
B0 = 0.01 # Initial Biofilm bacteria
S0 = 0.01*Vol # Initial Planktonic predator
A0 = 0.001*Surf # Initial Biofilm predator
u0 = [C0,F0,B0,S0,A0]
# Model parameters
parameters = (r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, chi_on_max, chi_on_min, interaction_strength, chi_off)

### Define and solve the ODE system
prob = ODEProblem(BFP_model, u0, tspan, parameters) # Define the ODE problem
sol = solve(prob, ode_algo) # Solve the ODE problem
NofT = sol(output_times) # Extract the solution at the output times

### Plot the solution
plt = plot()
[plot!(plt,output_times, NofT[i,:]./u0[1], ylims=(0,1.1), label=false, lw=5) for i in 1:5] # Plot the solution
display(plt)

savefig("0_BasicModelPlotting/1B.png") 

##--------------------------------------------------##


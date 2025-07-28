
using DifferentialEquations, Plots

function BFP_compartment_model(du, u, theta, t)
    N_compartments = 2000
    chionmin = 0.0
    chionmax = 0.2
    chionVals = 24*LinRange(chionmin,chionmax,N_compartments)*(8.0/3.0)
    # Extract parameters
    r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, a = theta

    Ftot = sum(u[2:1+N_compartments])
    Btot = sum(u[2+N_compartments:2+2*N_compartments-1])
    
    Hc = (u[1]/(u[1]+H23))
    Hp = Ftot/(Ftot + H5)
    Hb = Btot/(Btot + H4)

    du[1] = -r2*Hc*Ftot - r3*Hc*Btot
    for i in 1:N_compartments
        du[1+i]                = e23*r2*Hc*u[i+1]                - r4*Hp*(u[1+i]/Ftot)*u[end-1]              - a*chionVals[i]*u[i+1]/(a+Btot)# + chi_off*u[1+N_compartments+i] 
        du[1+N_compartments+i] = e23*r3*Hc*u[1+N_compartments+i] - r5*Hb*(u[1+N_compartments+i]/Btot)*u[end] + a*chionVals[i]*u[i+1]/(a+Btot)# - chi_off*u[1+N_compartments+i]  
    end
    du[end-1] = e4*r4*Hp*u[end-1]
    du[end]   = e5*r5*Hb*u[end]
end

# Model constants
ND_C = 1.0
ND_T = 24.0
e23 = 0.2 # Free bacteria growth efficiency
e4 = 0.5 # Bound bacteria growth efficiency
e5 = 0.33 # Predator growth efficiency
r2 = (1/e23)*0.21*ND_T # Free bacteria growth rate - cells/hr * ND_hrs
r3 = (1/e23)*0.007*ND_T # Bound bacteria growth rate - cells/hr * ND_hrs
r4 = (1/e4)*0.12*ND_T # Predator growth rate - cells/hr * ND_hrs
r5 = (1/e5)*0.09*ND_T # Predator death rate - 1/hrs * ND_hrs
H23 = 3.0*1.0/ND_C # Half saturations of bacteria
H4  = 3.0*1.0/ND_C # Half saturation of free predator
H5  = 8.0*0.1/ND_C # Half saturation of bound predator
#chi_off = 0.005*ND_T # Detachment rate
a = 0.01*8.0 # interaction strength
theta = [r2, r3, r4, r5, e23, e4, e5, H23, H4, H5, a] # Parameters

# Other constants
N_compartments = 2000
chionmin = 0.0
chionmax = 0.2
chionVals = LinRange(chionmin,chionmax,N_compartments)*ND_T*8.0/3.0
bw = chionVals[2] - chionVals[1]
output_times = 0:0.1:1.0
N_generations = 500
E_of_generations = zeros(N_generations)

# Preallocation
N_bins = 3+2*N_compartments
ind_F = 2:1+N_compartments
ind_B = 2+N_compartments:N_bins-2
C0 = 1.0 # Initial concentration of Free bacteria
sig = 45
F0 = zeros(N_compartments)
for i in 1:N_compartments
    F0[i] = exp(-0.5*((i-0.5*N_compartments)/sig)^2)/sig*sqrt(2*pi)
end
F0 = F0./sum(F0)
B0 = 0.001*F0./sum(F0) # Initial concentration of Biofilm bacteria
pred0max = 2.0
S0 = 0.0/ND_C # Initial concentration of Ciliate predator
A0 = 0.0/ND_C # Initial concentration\ of Biofilm predator
u0 = vcat(C0,F0,B0,S0,A0)
colors = palette(cgrad([:green, :blue]), N_generations); 
N_predCases = 10
colors2 = palette(cgrad([:grey, :purple1], N_predCases+1)); 
#colors2 = palette(cgrad([:grey, :goldenrod1], N_predCases+1)); 

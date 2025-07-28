# Structural Identifiability 
## Author: Stephen Williams

These equations can be used alongside the MapleCloud SIAN toolbox: https://maple.cloud/app/6509768948056064/ to determine structural identifiability.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

## State variables measured: C, F, B, S, A

dx1/dt = - r2*x1/(x1+H23)*x2 - r3*x1/(x1+H23)*x3,
dx2/dt = e23*r2*x1/(x1+H23)*x2 - r4*x2/(x2+H4)*x4 - ((chi_on_max*a + chi_on_min*x3)/(a+x3))*x2 + chioff*x3,
dx3/dt = e23*r3*x1/(x1+H23)*x3 - r5*x3/(x3+H5)*x5 + ((chi_on_max*a + chi_on_min*x3)/(a+x3))*x2 - chioff*x3,
dx4/dt = r4*x2/(x2+H4)*x4,
dx5/dt = r5*x3/(x3+H5)*x5,
y1=x1,
y2=x2,
y3=x3,
y4=x4,
y5=x5

RESULTS -- Fully Gloabally Identifiable.

## State variables measured: B
# Measuring chi_on_max, a, chi_on_min, chioff

dx1/dt = - r2*x1/(x1+H23)*x2 - r3*x1/(x1+H23)*x3,
dx2/dt = e23*r2*x1/(x1+H23)*x2 - r4*x2/(x2+H4)*x4 - ((chi_on_max*a + chi_on_min*x3)/(a+x3))*x2 + chioff*x3,
dx3/dt = e23*r3*x1/(x1+H23)*x3 - r5*x3/(x3+H5)*x5 + ((chi_on_max*a + chi_on_min*x3)/(a+x3))*x2 - chioff*x3,
dx4/dt = r4*x2/(x2+H4)*x4,
dx5/dt = r5*x3/(x3+H5)*x5,
y3=x3

RESULTS -- (chi_on_max, a, chi_on_min, chioff) globally identifiable. 

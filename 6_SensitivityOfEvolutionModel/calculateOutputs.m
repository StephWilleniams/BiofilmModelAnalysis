%% Expected run time (12 cores) -- 

clear workspace
close gcf

% Display the time in the default format
disp('Job start: ')
tic()
disp(datetime())

% Parallelization information
if isempty(gcp('nocreate'))
    pool = parpool('Threads');
end
params = load('params.txt');
N_runs = length(params(:,1));
NWORKERS = pool.NumWorkers;
jobsPthread = floor(N_runs/NWORKERS);
output_store = zeros(NWORKERS,jobsPthread);

% Input set ups
N_compartments = 2000;
chionmin = 0.0;
chionmax = 0.2;
chionVals = linspace(chionmin,chionmax,N_compartments)*(8.0/3.0)*24.0;
S0 = 1.0;
A0 = 1.0;
init_mean = N_compartments/2;
sig = 30;
baseline = [S0, A0, init_mean, sig];
delta = [1.0,1.0,0.3,1.0];
job_parameters = baseline .* (1 + delta.*params);

% Main loop
parfor node_index = 1:NWORKERS
    job_list = job_parameters;
    for job_index = 1:jobsPthread
        job_number = (node_index-1)*jobsPthread + job_index
        output_store(node_index,job_index) = model(job_list(job_number,:),chionVals);
    end
end

%%

% Final results not parellizable
missing_jobs = N_runs - jobsPthread*NWORKERS;
if missing_jobs>0
    final_outputs = zeros(missing_jobs,1);
    for i = 1:missing_jobs
        final_outputs(i) = model(job_parameters(N_runs-missing_jobs+i,:),chionVals);
    end
end

% Format outputs

%%

outputs = [];
for ii = 1:NWORKERS
    outputs = [outputs;output_store(ii,:)'];
end
outputs = [outputs;final_outputs];

elapsed_time = toc; % Get elapsed time
disp(['Total execution time: ', num2str(elapsed_time), ' seconds.']);
disp("Execution complete.")

%%

write_output = true;
if write_output

    disp(sum(isnan(outputs)));

    ind = find(isnan(outputs));
    outputs(ind) = 1.0;

    filename = 'output_values.txt';
    writematrix(outputs, filename);
end


%%

function dydt = BFP_compartment_model(~,y,chionVals)

ND_C = 1.0; % Non-dimensional concentration
ND_T = 24.0; % Non-dimensional time

e23 = 0.2; % Free bacteria growth efficiency
e4 = 0.5; % Bound bacteria growth efficiency
e5 = 0.33; % Predator growth efficiency

r2 = (1/e23)*0.21*ND_T; % Free bacteria growth rate - cells/hr * ND_hrs
r3 = (1/e23)*0.007*ND_T; % Bound bacteria growth rate - cells/hr * ND_hrs
r4 = (1/e4)*0.12*ND_T; % Predator growth rate - cells/hr * ND_hrs
r5 = (1/e5)*0.09*ND_T; % Predator death rate - 1/hrs * ND_hrs

H23 = 3.0*1.0/ND_C; % Half saturations of bacteria
H4 = 3.0*1.0/ND_C; % Half saturation of free predator
H5 = 8.0*0.1/ND_C; % Half saturation of bound predator

a = 8.0*0.01; 

N_compartments = 2000; % Number of compartments for Biofilm
N_bins = 3+2*N_compartments; % Total number of bins
ind_F = 2:1+N_compartments; % Indices for bound bacteria
ind_B = 2+N_compartments:N_bins-2; % Indices for bound bacteria

dydt = zeros(N_bins,1);

C = y(1);
F = y(2:1+N_compartments);
B = y(2+N_compartments:1+2*N_compartments);
S = y(end-1);
A = y(end);

sumF = sum(F);
sumB = sum(B);

HF = sumF/(sumF + H4);
HB = sumB/(sumB + H5);

dydt(1) = -r2*(C/(C+H23))*sumF - r3*(C/(C+H23))*sumB;

dydt(ind_F) = e23*r2*(C/(C+H23))*F ...
            - r4*HF*(F/sumF)*S ...
            - a*chionVals'.*F/(a + sumB);

dydt(ind_B) = e23*r3*(C/(C+H23))*B ...
            - r5*HB*(B/sumB)*A ...
            + a*chionVals'.*F/(a + sumB);

dydt(end-1) = e4*r4*HF*S;
dydt(end)   = e5*r5*HB*A;

end

function returner = model(params,chionVals)

QOI = zeros(2,1);

S0_input = params(1);
A0_input = params(2);
init_mean_input = params(3);
sig_input = params(4);

ND_C = 1.0; % Non-dimensional concentration

N_generations = 50; % Number of generations
N_compartments = 2000; % Number of compartments for Biofilm
N_bins = 3+2*N_compartments; % Total number of bins
ind_B = 2+N_compartments:N_bins-2; % Indices for bound bacteria

C0 = 1.0/ND_C;
F0 = exp(-0.5*(((1:N_compartments)-init_mean_input)/sig_input).^2)/(sig_input*sqrt(2*pi));
F0 = F0./sum(F0);
B0 = 0.001*F0./sum(F0);
S0 = S0_input/ND_C;
A0 = A0_input/ND_C;

y0 = [C0,F0,B0,S0,A0];

QOI(1) = (1/sum(F0)) * chionVals * F0';

for tt = 1:N_generations
    sol = ode45(@(t,y) BFP_compartment_model(t,y,chionVals),[0,1.0],y0);
    output = deval(sol,1.0);
    F0 = output(ind_B)'./sum(output(ind_B));
    B0 = 0.001*F0;
    y0 = [C0,F0,B0,S0,A0];
end

QOI(2) = (1/sum(F0)) * chionVals * F0';

returner = QOI(2)/QOI(1);
end
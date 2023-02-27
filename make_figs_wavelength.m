%% without growth
close all
clear all
addpath("SSB_scripts/")

%% Load parameters                            
params.domain_length = 60;                
params.time_step = 0.01;
params.number_gridpoints = 128;
params.random_seed = 1;
params.end_time = 2e2;
params.diffusion_constants = [1 1 0 30];
params.initial_conditions = [0.02 0.26 0.02 0.006];
params.reaction_parameters(1) = 0.1;     
params.reaction_parameters(2) = 6.7;    
params.reaction_parameters(3) = 1;   
params.reaction_parameters(4) = 1;
params.reaction_parameters(5) = 10;
params.reaction_parameters(6) = 1;
params.reaction_parameters(7) = 40;
params.reaction_parameters(8) = 40;
params.reaction_parameters(9) = 0.01;
params.reaction_parameters(10) = 8;
params.reaction_parameters(11) = 2;

%Define ODE
reaction_function=@(x,k)[ k(8).*x(4) - k(7).*x(1).*x(2)-k(1).*x(1)+ (k(4)./(1+((x(3)./k(9)).^k(10)))); ...
    k(8).*x(4) - k(7).*x(1).*x(2)-k(2).*x(2)+ (k(5)./(1+((x(3)./k(9)).^k(11)))); ...
    k(6).*x(1)- k(3).*x(3) ; ...
    k(7).*x(1).*x(2)-k(8).*x(4)];


%% Change k1
k1_big = params;
k1_big.reaction_parameters(1) = 0.1;
output_k1_big = simulate_without_growth(k1_big,reaction_function);

k1_small = params;
k1_small.reaction_parameters(1) = 0.02;
output_k1_small = simulate_without_growth(k1_small,reaction_function);

%% Change k2
k2_big = params;
k2_big.reaction_parameters(2) = 5;
output_k2_big = simulate_without_growth(k2_big,reaction_function);

k2_small = params;
k2_small.reaction_parameters(2) = 1;
output_k2_small = simulate_without_growth(k2_small,reaction_function);

%% Change k3
k3_big = params;
k3_big.reaction_parameters(3) = 4;
output_k3_big = simulate_without_growth(k3_big,reaction_function);

k3_small = params;
k3_small.reaction_parameters(3) = 0.8;
output_k3_small = simulate_without_growth(k3_small,reaction_function);

%% Change k4
k4_big = params;
k4_big.reaction_parameters(4) = 1;
output_k4_big = simulate_without_growth(k4_big,reaction_function);

k4_small = params;
k4_small.reaction_parameters(4) = 0.2;
output_k4_small = simulate_without_growth(k4_small,reaction_function);

%% Change k5
k5_big = params;
k5_big.reaction_parameters(5) = 50;
output_k5_big = simulate_without_growth(k5_big,reaction_function);

k5_small = params;
k5_small.reaction_parameters(5) = 10;
output_k5_small = simulate_without_growth(k5_small,reaction_function);

%% Change k6
k6_big = params;
k6_big.reaction_parameters(6) = 2;
output_k6_big = simulate_without_growth(k6_big,reaction_function);

k6_small = params;
k6_small.reaction_parameters(6) = 0.4;
output_k6_small = simulate_without_growth(k6_small,reaction_function);

%% Change k7
k7_big = params;
k7_big.reaction_parameters(7) = 200;
output_k7_big = simulate_without_growth(k7_big,reaction_function);

k7_small = params;
k7_small.reaction_parameters(7) = 40;
output_k7_small = simulate_without_growth(k7_small,reaction_function);

%% Change k8
k8_big = params;
k8_big.reaction_parameters(8) = 40;
output_k8_big = simulate_without_growth(k8_big,reaction_function);

k8_small = params;
k8_small.reaction_parameters(8) = 8;
output_k8_small = simulate_without_growth(k8_small,reaction_function);

%% Change k9
k9_big = params;
k9_big.reaction_parameters(9) = 0.025;
k9_big.end_time = 3e2;
output_k9_big = simulate_without_growth(k9_big,reaction_function);

k9_small = params;
k9_small.reaction_parameters(9) = 0.005;
k9_small.end_time = 3e2;
output_k9_small = simulate_without_growth(k9_small,reaction_function);

%% Change D_GDF5
D_GDF5_big = params;
D_GDF5_big.diffusion_constants = [1 1 0 30];
output_D_GDF5_big = simulate_without_growth(D_GDF5_big,reaction_function);

D_GDF5_small = params;
D_GDF5_small.diffusion_constants = [0.2 1 0 30];
output_D_GDF5_small = simulate_without_growth(D_GDF5_small,reaction_function);

%% Change D_NOG
D_NOG_big = params;
D_NOG_big.diffusion_constants = [1 1 0 30];
output_D_NOG_big = simulate_without_growth(D_NOG_big,reaction_function);

D_NOG_small = params;
D_NOG_small.diffusion_constants = [1 0.2 0 30];
output_D_NOG_small = simulate_without_growth(D_NOG_small,reaction_function);

%% Change D_C
D_C_big = params;
D_C_big.diffusion_constants = [1 1 0 150];
output_D_C_big = simulate_without_growth(D_C_big,reaction_function);

D_C_small = params;
D_C_small.diffusion_constants = [1 1 0 30];
output_D_C_small = simulate_without_growth(D_C_small,reaction_function);

%% Change D
D_big = params;
D_big.diffusion_constants = 5*[1 1 0 30];
output_D_big = simulate_without_growth(D_big,reaction_function);

D_small = params;
D_small.diffusion_constants = [1 1 0 30];
output_D_small = simulate_without_growth(D_small,reaction_function);

%% SAVE DATA

save data/wavelength.mat






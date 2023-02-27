%% without growth

close all
clear all
addpath("SSB_scripts/")

%% Load parameters                            
params.domain_length = 45;                
params.time_step = 0.02;
params.number_gridpoints = 256;
params.random_seed = 1;
params.end_time = 1e3;
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
reaction_function=@(x,k)[ k(8).*x(4) - k(7).*x(1).*x(2)-k(1).*x(1)+ (1./(1+((x(3)./k(9)).^k(10)))); ...
    k(8).*x(4) - k(7).*x(1).*x(2)-k(2).*x(2)+ (k(5)./(1+((x(3)./k(9)).^k(11)))); ...
    k(6).*x(1)- k(3).*x(3) ; ...
    k(7).*x(1).*x(2)-k(8).*x(4)];

%% Simulate
wildtype = simulate_without_growth(params,reaction_function);

%% repeat with different seed
params2 = params;
params2.random_seed = 2;
wildtype2 = simulate_without_growth(params2,reaction_function);


%% Simulations with different parameter set
params.diffusion_constants = [30 1 0 30];
params.initial_conditions = [0.003 2.4 0.003 0.03];
params.reaction_parameters(1) = 0.1;     
params.reaction_parameters(2) = 0.1;    
params.reaction_parameters(3) = 1;   
params.reaction_parameters(4) = 1;
params.reaction_parameters(5) = 2.08;
params.reaction_parameters(6) = 1;
params.reaction_parameters(7) = 46;
params.reaction_parameters(8) = 10;
params.reaction_parameters(9) = 0.001;
params.reaction_parameters(10) = 8;
params.reaction_parameters(11) = 2;
params.domain_length = 200;

params.time_step = 0.01;
params.end_time = 1e3;

no_shuttling = simulate_without_growth(params,reaction_function);


%% Simulations with equal diffusivities
params.diffusion_constants = [1 1 0 1];
params.initial_conditions = [0.02 0.87 0.027 0.02];

params.reaction_parameters(1) = 6.93;     
params.reaction_parameters(2) = 10;    
params.reaction_parameters(3) = 1;   
params.reaction_parameters(4) = 1;
params.reaction_parameters(5) = 73;
params.reaction_parameters(6) = 1;
params.reaction_parameters(7) = 100;
params.reaction_parameters(8) = 63;
params.reaction_parameters(9) = 0.0077;
params.reaction_parameters(10) = 2;
params.reaction_parameters(11) = 8;
params.reaction_parameters(12) = 0.015;
params.random_seed = 1;

params.domain_length = 15;
params.number_gridpoints = 256;
params.time_step = 0.005;
params.end_time = 1e3;

close all
reaction_function=@(x,k)[ k(8).*x(4) - k(7).*x(1).*x(2)-k(1).*x(1)+ (1./(1+((x(3)./k(9)).^k(10)))); ...
    k(8).*x(4) - k(7).*x(1).*x(2)-k(2).*x(2)+ (k(5)./(1+((x(3)./k(12)).^k(11)))); ...
    k(6).*x(1)- k(3).*x(3) ; ...
    k(7).*x(1).*x(2)-k(8).*x(4)];

equal_diffusion = simulate_without_growth(params,reaction_function);


%% save

save data/no_growth.mat




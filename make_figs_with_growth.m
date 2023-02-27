%% Clear and close things (optional)
clear all; close all;
addpath("SSB_scripts/")
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define geometry
params.total_growth_rate = 0.02;        
params.fraction_of_growth_at_tip = 0.8;                         
params.left_boundary_position = 0;              
params.length_patterning_zone = 8;           
params.domain_length = 45;             
params.initial_digit_length = 9;           
params.time_to_initalize_pattern = 150;         
params.time_step = 0.02;

params.number_gridpoints = 128;
params.random_seed = 1;
params.committed_species = [1 2 3 4];

% define ODE
reaction_function=@(x,k)[ k(8).*x(4) - k(7).*x(1).*x(2)-k(1).*x(1)+ (1./(1+((x(3)./k(9)).^k(10)))); ...
    k(8).*x(4) - k(7).*x(1).*x(2)-k(2).*x(2)+ (k(5)./(1+((x(3)./k(9)).^k(11)))); ...
    k(6).*x(1)- k(3).*x(3) ; ...
    k(7).*x(1).*x(2)-k(8).*x(4)];

% set parameters 
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

% set dynamics outside of the digit
params.boundary.reflection = 1; 

% We mainly choose reflective boundary conditions at the digit ends
% Note, other possibilities can be specified via boundary.function and
% boundary.parameters


% Reference param set
p_ref = params;
p_ref.reaction_parameters(2) = 3;   
p_ref.total_growth_rate = 0.1;
p_ref.time_to_initalize_pattern = 50;
p_ref.length_patterning_zone = 4;
p_ref.filenm = "ref";


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamics with PFR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p_PFR = p_ref;


p_PFR.boundary.function =@(x,k)[    k(1) - k(2) * x(1); ...
                    k(3) - k(4) * x(2); ...
                    k(5) - k(6) * x(3); ...
                    k(7) - k(8) * x(4); ...
    ];

%GDF5
p_PFR.boundary.parameters(1) = 0.1;
p_PFR.boundary.parameters(2) = 1;

%NOG
p_PFR.boundary.parameters(3) = 0;
p_PFR.boundary.parameters(4) = 1;

%PSMAD
p_PFR.boundary.parameters(5) = 0.1;
p_PFR.boundary.parameters(6) = 1;

%C
p_PFR.boundary.parameters(7) = 0;
p_PFR.boundary.parameters(8) = 1;

p_PFR.boundary.reflection = 0;

p_PFR.filenm = "PFR";
p_PFR.length_patterning_zone = 10;
p_PFR.initial_digit_length = 5;
p_PFR.total_growth_rate = 0.5;
p_PFR.domain_length = 50;

PFR_zoom = simulate_with_growth(p_PFR,reaction_function,[85 105 105.4 106 111]);
PFR_heatmap = simulate_with_growth(p_PFR,reaction_function,[60 87.6 111 139]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vary joint number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Increasing joint number: reference
p_reference_twojoints = p_ref; 
p_reference_twojoints.reaction_parameters(2) = 2;
p_reference_twojoints.initial_digit_length = 10;
p_reference_twojoints.time_to_initalize_pattern = 100;
p_reference_twojoints.filenm = "extra_joint_ref";
extra_joint_ref = simulate_with_growth(p_reference_twojoints,reaction_function,[]);

%% Increasing joint number: wavelength
p_wavelength_twojoints = p_reference_twojoints; 
p_wavelength_twojoints.reaction_parameters(2) = 4;

%we also need to alter initial condition settings, to ensure we begin with
%a single joint to mimic the reference
p_wavelength_twojoints.initial_digit_length = 7;
p_wavelength_twojoints.time_to_initalize_pattern = 80;

p_wavelength_twojoints.filenm = "shorter_wavelength";
extra_joint = simulate_with_growth(p_wavelength_twojoints,reaction_function,[]);

%% Increasing joint number: longer growth duration
p_longer_growth = p_reference_twojoints; 
p_longer_growth.domain_length = 60;
p_longer_growth.filenm = "p_longer_growth";
longer_growth = simulate_with_growth(p_longer_growth,reaction_function,[]);

%% Increasing joint number: faster growth, same duration
p_faster_growth = p_reference_twojoints; 

% growth duration is determined implicitly whenever digit grows to end of
% domain. Therefore, specify domain length. Then compute growth rate such
% that the time it takes to reach the end of the domain is the same
p_faster_growth.domain_length = 60;
p_faster_growth.total_growth_rate = p_reference_twojoints.total_growth_rate * (p_faster_growth.domain_length - p_faster_growth.initial_digit_length)/(p_reference_twojoints.domain_length - p_ref.initial_digit_length);
p_faster_growth.filenm = "p_faster_growth";
faster_growth = simulate_with_growth(p_faster_growth,reaction_function,[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MUTANTS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GDF5(-/-)
p_gdf5_mutant3 = p_ref;
p_gdf5_mutant3.reaction_parameters(6) = 0;
p_gdf5_mutant3.time_step = 0.002;
p_gdf5_mutant3.filenm = "gdf5_mutant3";
p_gdf5_mutant3.domain_length = 40;
gdf5_3 = simulate_with_growth(p_gdf5_mutant3,reaction_function,[]);

%% reference
p_reference3 = p_ref;
p_reference3.filenm = "reference3";
p_reference3.domain_length = 40;
reference3 = simulate_with_growth(p_reference3,reaction_function,[]);

%% NOG(-/-) (3 joints)
p_nog_mutant3 = p_ref;
p_nog_mutant3.filenm = "nog_mutant3";
p_nog_mutant3.reaction_parameters(5) = 0;
p_nog_mutant3.domain_length = 40;
nog3 = simulate_with_growth(p_nog_mutant3,reaction_function,[]);

%% GDF5 bead: control (3 joints)
p_gdf5_bead_control3 = p_ref;
p_gdf5_bead_control3.filenm = "gdf5_bead_control3";
p_gdf5_bead_control3.fraction_of_growth_at_tip = 1.0;
bead_params.species = 1;
bead_params.amplitude = 0;
bead_params.location = 17;
bead_params.width = 2;
p_gdf5_bead_control3.initial_digit_length = 8;
p_gdf5_bead_control3.domain_length = 33;
p_gdf5_bead_control3.time_to_initalize_pattern = 250;

bead_control3 = simulate_with_growth_and_bead(p_gdf5_bead_control3,reaction_function,[], bead_params);

%% GDF5 bead: treatment (3 joints)
p_gdf5_bead3= p_ref;
p_gdf5_bead3.filenm = "gdf5_bead3";

p_gdf5_bead3.fraction_of_growth_at_tip = 1.0;

bead_params.species = 1;
bead_params.amplitude = .1;
bead_params.location = 17;
bead_params.width = 2;

p_gdf5_bead3.initial_digit_length = 8;
p_gdf5_bead3.domain_length = 33;
p_gdf5_bead3.time_to_initalize_pattern = 250;

bead3 = simulate_with_growth_and_bead(p_gdf5_bead3,reaction_function,[], bead_params);



%% Save to file
save data/with_growth.mat




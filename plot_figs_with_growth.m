%% Clear and close things (optional)
clear all
close all
addpath("SSB_scripts/")
load data/with_growth.mat
load data/colormaps.mat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normalization = [.1 0.5 0.1 0.01 .02 0.4 0.1];


%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Fig S3 (dynamic, with PFR)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GDF5, reference
colAxis = [0 2];
cMap = cMap_GDF5;
heatmap_growth(PFR_heatmap,normalization, 5, "pfr_gdf5_1.svg", 1,colAxis, 0, cMap)
heatmap_growth(PFR_heatmap,normalization, 5, "pfr_gdf5_2.svg", 2,colAxis, 0, cMap)
heatmap_growth(PFR_heatmap,normalization, 5, "pfr_gdf5_3.svg", 3,colAxis, 0, cMap)
heatmap_growth(PFR_heatmap,normalization, 5, "pfr_gdf5_4.svg", 4,colAxis, 0, cMap)

%% PSMAD, reference

colAxis = [0.1 0.9];
cMap = cMap_pSmad;

heatmap_growth(PFR_heatmap,normalization, 3, "pfr_psmad_1.svg", 1,colAxis, 0, cMap)
heatmap_growth(PFR_heatmap,normalization, 3, "pfr_psmad_2.svg", 2,colAxis, 0, cMap)
heatmap_growth(PFR_heatmap,normalization, 3, "pfr_psmad_3.svg", 3,colAxis, 0, cMap)
heatmap_growth(PFR_heatmap,normalization, 3, "pfr_psmad_4.svg", 4,colAxis, 0, cMap)

%% Zoom
normalization = [.1 0.5 0.06 0.01 .02 0.4 0.1];

plotLineWithGrowth(PFR_zoom,normalization,[5 3], "pfr_1_zoom.pdf",1,[0 2])
plotLineWithGrowth(PFR_zoom,normalization,[5 3], "pfr_2_zoom.pdf",2,[0 2])
plotLineWithGrowth(PFR_zoom,normalization,[5 3], "pfr_3_zoom.pdf",3,[0 2])
plotLineWithGrowth(PFR_zoom,normalization,[5 3], "pfr_4_zoom.pdf",4,[0 2])
plotLineWithGrowth(PFR_zoom,normalization,[5 3], "pfr_5_zoom.pdf",5,[0 2])




%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Fig 4 (mutants) [3 JOINTS]
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colAxis = [0 1.2];
cMap = cMap_GDF5;
heatmap_growth(reference3,normalization, 5, "wildtype_heatmap3.svg", 1,colAxis, 0, cMap)
heatmap_growth(gdf5_3,normalization, 5, "gdf5_mutant_heatmap3.svg", 1,colAxis, 0, cMap)
heatmap_growth(nog3,normalization, 5, "nog_mutant_heatmap3.svg", 1,colAxis, 0, cMap)

bead_pos = bead_params.location * bead3.params.number_gridpoints / bead3.params.domain_length;
heatmap_growth(bead_control3,normalization, 5, "bead_control_heatmap3.svg", 1,colAxis, bead_pos, cMap)
heatmap_growth(bead3,normalization, 5, "bead_heatmap3.svg", 1,colAxis, bead_pos, cMap)

%% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Fig 5 (changes to joint number)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

colAxis = [0 2];

cMap = cMap_GDF5;


rel_length = extra_joint.params.domain_length / extra_joint_ref.params.domain_length;
heatmap_growth_vary_length(extra_joint_ref,normalization, 5, "reference.svg", 1,colAxis, 0, cMap,1)

rel_length = extra_joint.params.domain_length / extra_joint_ref.params.domain_length;
heatmap_growth_vary_length(extra_joint,normalization, 5, "wavelength.svg", 1,0.4*colAxis, 0, cMap,rel_length)

rel_length = faster_growth.params.domain_length / extra_joint_ref.params.domain_length;
heatmap_growth_vary_length(faster_growth,normalization, 5, "faster_growth.svg", 1,colAxis, 0, cMap,rel_length)

rel_length = longer_growth.params.domain_length / extra_joint_ref.params.domain_length;
heatmap_growth_vary_length(longer_growth,normalization, 5, "longer_growth.svg", 1,colAxis, 0, cMap,rel_length)




%% plot from scratch
clear all
load data/no_growth.mat

%% Adjust normalization
normalization = [.1 0.5 0.04 0.01 0.4 0.008 0.01];
plotind = [3 5 6];
filenm = "tmp.pdf";
plotLine(equal_diffusion,normalization,plotind, filenm)


%% wildtype
normalization = [.1 0.5 0.1 0.01 .02 0.4 0.1];
plotLine(wildtype,normalization,[5], "gdf_RNA.pdf")
plotLine(wildtype,normalization,[5 6], "nog_RNA.pdf")
plotLine(wildtype,normalization,[5 3], "pSmad.pdf")
plotLine(wildtype,normalization,[5 1], "free_GDF5.pdf")

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DIFFERENT MODELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% reference parameter set (different seed)
normalization = [.1 0.5 0.1 0.01 .02 0.4 0.1];
plotLine(wildtype2,normalization,[5 6], "wildtype2_nog_RNA.pdf")
plotLine(wildtype2,normalization,[5 3], "wildtype2_pSmad.pdf")

%% no shuttling
normalization = [0.03 5 0.02 0.05 .005 0.3 0.05];
plotLine(no_shuttling,normalization,[5 6], "no_shuttling_nog_RNA.pdf")
plotLine(no_shuttling,normalization,[5 3], "no_shuttling_pSmad.pdf")

%% equal diffusion
normalization = [.1 0.5 0.04 0.01 0.4 0.008 0.01];
plotLine(equal_diffusion,normalization,[5 6], "equal_diffusion_nog_RNA.pdf")
plotLine(equal_diffusion,normalization,[5 3], "equal_diffusion_pSmad.pdf")






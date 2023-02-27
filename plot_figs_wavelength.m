%% plot from scratch
clear all
addpath("SSB_scripts/")
load data/wavelength.mat
load data/colormaps.mat


%% wildtype
normalization = [.1 0.5 0.1 0.01 .02 0.4 0.1];
colAxis = [0 1.];

cMap = cMap_GDF5;

figure()

subplot(10,3,1)
heatmap_static_grid(output_k1_big,normalization, 5, colAxis, cMap, "k1")
subplot(10,3,4)
heatmap_static_grid(output_k1_small,normalization, 5, colAxis, cMap,"")

subplot(10,3,2)
heatmap_static_grid(output_k2_big,normalization, 5, colAxis, cMap, "k2")
subplot(10,3,5)
heatmap_static_grid(output_k2_small,normalization, 5, colAxis, cMap,"")

subplot(10,3,3)
heatmap_static_grid(output_k3_big,normalization, 5, colAxis, cMap, "k3")
subplot(10,3,6)
heatmap_static_grid(output_k3_small,normalization, 5, colAxis, cMap,"")

subplot(10,3,7)
heatmap_static_grid(output_k4_big,normalization, 5, colAxis, cMap, "k4")
subplot(10,3,10)
heatmap_static_grid(output_k4_small,normalization, 5, colAxis, cMap,"")

subplot(10,3,8)
heatmap_static_grid(output_k5_big,normalization, 5, colAxis, cMap, "k5")
subplot(10,3,11)
heatmap_static_grid(output_k5_small,normalization, 5, colAxis, cMap,"")

subplot(10,3,9)
heatmap_static_grid(output_k6_big,normalization, 5, colAxis, cMap, "k6")
subplot(10,3,12)
heatmap_static_grid(output_k6_small,normalization, 5, colAxis, cMap,"")

subplot(10,3,13)
heatmap_static_grid(output_k7_big,normalization, 5, colAxis, cMap, "k7")
subplot(10,3,16)
heatmap_static_grid(output_k7_small,normalization, 5, colAxis, cMap,"")

subplot(10,3,14)
heatmap_static_grid(output_k8_big,normalization, 5, colAxis, cMap, "k8")
subplot(10,3,17)
heatmap_static_grid(output_k8_small,normalization, 5, colAxis, cMap,"")

subplot(10,3,15)
heatmap_static_grid(output_k9_big,normalization, 5, colAxis, cMap, "k9")
subplot(10,3,18)
heatmap_static_grid(output_k9_small,normalization, 5, colAxis, cMap,"")

subplot(10,3,19)
heatmap_static_grid(output_D_GDF5_big,normalization, 5, colAxis, cMap, "D_{GDF5}")
subplot(10,3,22)
heatmap_static_grid(output_D_GDF5_small,normalization, 5, colAxis, cMap,"")

subplot(10,3,20)
heatmap_static_grid(output_D_NOG_big,normalization, 5, colAxis, cMap, "D_{NOG}")
subplot(10,3,23)
heatmap_static_grid(output_D_NOG_small,normalization, 5, colAxis, cMap,"")

subplot(10,3,21)
heatmap_static_grid(output_D_C_big,normalization, 5, colAxis, cMap, "D_{C}")
subplot(10,3,24)
heatmap_static_grid(output_D_C_small,normalization, 5, colAxis, cMap,"")

subplot(10,3,26)
heatmap_static_grid(output_D_big,normalization, 5, colAxis, cMap, "D")
subplot(10,3,29)
heatmap_static_grid(output_D_small,normalization, 5, colAxis, cMap,"")

saveas(gcf,strcat('Fig/grid_wavelength.svg'));

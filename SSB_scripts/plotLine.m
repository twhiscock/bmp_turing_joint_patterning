function plotLine(output,normalization,plotind, filenm)
close all
m = output.m;
params = output.params;
names = ["GDF5 protein" "NOG protein" "pSMAD" "Complex" "Gdf5 RNA" "Nog RNA" "total GDF5"];
colors = {'#77AC30','#4DBEEE','#FF00FF','#EDB120','#00FF00','#0072BD', '#A2142F'};
y(:,1) = m(:,1)/normalization(1);
y(:,2) = m(:,2)/normalization(2);
y(:,3) = m(:,3)/normalization(3);
y(:,4) = m(:,4)/normalization(4);
y(:,5) = (1/normalization(5))*(1./(1+((m(:,3)./params.reaction_parameters(9)).^params.reaction_parameters(10))));
y(:,6) = (1/normalization(6))*(1./(1+((m(:,3)./params.reaction_parameters(9)).^params.reaction_parameters(11))));
y(:,7) = (m(:,1) + m(:,4))/normalization(7);
x = params.domain_length/params.number_gridpoints*(1:params.number_gridpoints);


figure()
hold on 
for i=1:length(plotind)
    plot(x,y(:,plotind(i)),'LineWidth',3, 'Color',colors{plotind(i)})
end
hold off
% legend(names(plotind),'Location','bestoutside')
% xlabel('position')
% ylabel('normalized concentration')

pbaspect([2 1 1])
xlim([0 params.domain_length])
ylim([0 1])
xticks([])
yticks([])

saveas(gcf,strcat('Fig/',filenm));

end
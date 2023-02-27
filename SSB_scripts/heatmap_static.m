function heatmap_static(output,normalization, plotind, filenm,colAxis,cMap)

m = output.m(:,:);
params = output.params;

y(:,1) = m(:,1)/normalization(1);
y(:,2) = m(:,2)/normalization(2);
y(:,3) = m(:,3)/normalization(3);
y(:,4) = m(:,4)/normalization(4);
y(:,5) = (1/normalization(5))*(1./(1+((m(:,3)./params.reaction_parameters(9)).^params.reaction_parameters(10))));
y(:,6) = (1/normalization(6))*(1./(1+((m(:,3)./params.reaction_parameters(9)).^params.reaction_parameters(11))));
y(:,7) = (m(:,1) + m(:,4))/normalization(7);

figure()
imagesc(y(:,plotind)')
xlim([1 length(y)])
caxis(colAxis)
set(gcf,'color','w');
pbaspect([6 1 1])
colormap(cMap)

set(gca,'xtick',[])
set(gca,'ytick',[])
saveas(gcf,strcat('Fig/',filenm));

end
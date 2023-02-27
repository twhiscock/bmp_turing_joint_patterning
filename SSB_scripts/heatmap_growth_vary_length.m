function heatmap_growth_vary_length(output,normalization, plotind, filenm, t_index,colAxis, bead_pos,cMap, relative_length)

m = output.m(:,:,t_index);
ind = 1:output.l_digit(t_index);

params = output.params;


y(:,1) = m(ind,1)/normalization(1);
y(:,2) = m(ind,2)/normalization(2);
y(:,3) = m(ind,3)/normalization(3);
y(:,4) = m(ind,4)/normalization(4);
y(:,5) = (1/normalization(5))*(1./(1+((m(ind,3)./params.reaction_parameters(9)).^params.reaction_parameters(10))));
y(:,6) = (1/normalization(6))*(1./(1+((m(ind,3)./params.reaction_parameters(9)).^params.reaction_parameters(11))));
y(:,7) = (m(ind,1) + m(ind,4))/normalization(7);

figure()
imagesc(y(ind,plotind)')


xlim([1 max(output.l_digit)])


caxis(colAxis)
set(gcf,'color','w');
pbaspect([6*relative_length 1 1])
colormap(cMap)
if (bead_pos > 0)
    hold on
    scatter(bead_pos,0.5,100,'r')
    hold off
end
axis off
set(gca,'xtick',[])
set(gca,'ytick',[])
saveas(gcf,strcat('Fig/',filenm));

end
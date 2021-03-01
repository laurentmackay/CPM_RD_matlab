if plotting

tp__0=tic;

plotCellIm(panel1,reshape(RacRatio0,shape),cell_mask,i0,j0);
caxis(panel1,'auto');
colorbar(panel1);
title(panel1,'RacRatio0', 'Fontsize', 24);

plotCellIm(panel2,reshape(RacRatio,shape),cell_mask,i0,j0);
caxis(panel2,'auto');
colorbar(panel2);
title(panel2,'RacRatio', 'Fontsize', 24);

plotCellIm(panel3,reshape(RhoRatio,shape),cell_mask,i0,j0);
caxis(panel3,'auto');
colorbar(panel3);
title(panel3,'RhoRatio', 'Fontsize', 24);

plotCellIm(panel4,reshape(PaxRatio,shape),cell_mask,i0,j0);
caxis(panel4,'auto');
colorbar(panel4);
title(panel4,'PaxRatio', 'Fontsize', 24);

sgtitle(pic_fig,['t=' num2str(time) ', t_{plot}=' num2str(double(tic-tp__0)*1e-6), ', t_{sim}=' num2str(toc)], 'Fontsize', 10,'FontWeight','bold')

end
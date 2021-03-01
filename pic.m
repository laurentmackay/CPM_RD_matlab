if plotting

tp__0=tic;

plotCellIm(panel1,reshape(u(:, 1),shape),cell_mask,i0,j0);
caxis(panel1,'auto');
colorbar(panel1);
title(panel1,'Raci', 'Fontsize', 24);

plotCellIm(panel2,reshape(u(:, 2),shape),cell_mask,i0,j0);
caxis(panel2,'auto');
colorbar(panel2);
title(panel2,'Rac', 'Fontsize', 24);

plotCellIm(panel3,reshape(u(:, 3),shape),cell_mask,i0,j0);
caxis(panel3,'auto');
colorbar(panel3);
title(panel3,'Rhoi', 'Fontsize', 24);

plotCellIm(panel4,reshape(u(:, 4),shape),cell_mask,i0,j0);
caxis(panel4,'auto');
colorbar(panel4);
title(panel4,'Rho', 'Fontsize', 24);

plotCellIm(panel5,reshape(u(:, 5),shape),cell_mask,i0,j0);
caxis(panel5,'auto');
colorbar(panel5);
title(panel5,'Paxi', 'Fontsize', 24);

plotCellIm(panel6,reshape(u(:, 6),shape),cell_mask,i0,j0);
caxis(panel6,'auto');
colorbar(panel6);
title(panel6,'Paxs', 'Fontsize', 24);

plotCellIm(panel7,reshape(u(:, 7),shape),cell_mask,i0,j0);
caxis(panel7,'auto');
colorbar(panel7);
title(panel7,'Paxt', 'Fontsize', 24);

plotCellIm(panel8,reshape(u(:, 8),shape),cell_mask,i0,j0);
caxis(panel8,'auto');
colorbar(panel8);
title(panel8,'FAK', 'Fontsize', 24);

plotCellIm(panel9,reshape(u(:, 9),shape),cell_mask,i0,j0);
caxis(panel9,'auto');
colorbar(panel9);
title(panel9,'FP', 'Fontsize', 24);

plotCellIm(panel10,reshape(u(:,10),shape),cell_mask,i0,j0);
caxis(panel10,'auto');
colorbar(panel10);
title(panel10,'xFAK', 'Fontsize', 24);

plotCellIm(panel11,reshape(u(:,11),shape),cell_mask,i0,j0);
caxis(panel11,'auto');
colorbar(panel11);
title(panel11,'GPP', 'Fontsize', 24);

plotCellIm(panel12,reshape(u(:,12),shape),cell_mask,i0,j0);
caxis(panel12,'auto');
colorbar(panel12);
title(panel12,'xGPP', 'Fontsize', 24);

sgtitle(pic_fig,['t=' num2str(time) ', t_{plot}=' num2str(double(tic-tp__0)*1e-6), ', t_{sim}=' num2str(toc)], 'Fontsize', 10,'FontWeight','bold')

end
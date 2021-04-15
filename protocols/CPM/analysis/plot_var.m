


figure(1);
plotCellIm(gca,reshape(x(:,:,i_var),shape),cell_mask,i0,j0);
title([var ' (t=' num2str(time) ')'], 'Fontsize', 20);
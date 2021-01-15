if plotting
%     close all
    % to make a video all frames must be the same size setting its position just
    % stops some bugs
%     figure('Position', [200 75 1000 900])
    fs=14; %axis font size
    tp__0=tic;

%     imagesc(cell_mask,[0 1]);
    plotCellIm(panelA,double(cell_mask),cell_mask,i0,j0)
%     plotCellIm(panelA,alpha_chem(:,:,5),cell_mask,i0,j0)
%      plotCellIm(panelA,I_Ks,cell_mask,i0,j0)

    hold(panelA,'on')
    try
    plot(panelA, center(2,1:iter),center(1,1:iter),'r')
    catch e
        disp(e)
    end
    hold(panelA, 'off')
    ax = panelA;
    ax.FontSize = fs;
    colorbar(ax);
    % title('Cell', 'Fontsize', 24)
    % xlabel('X')
    % ylabel('Y')
    
    
%     subplot(2,2,2)
    plotCellIm(panelB,reshape(RhoRatio,shape),cell_mask,i0,j0)

    ax = panelB;
    caxis(ax,'auto')
    colorbar(ax);
    ax.FontSize = fs;
    title(panelB,'Rho', 'Fontsize', 24)
    % xlabel('X')
    % ylabel('Y')
    
    plotCellIm(panelC,reshape(RacRatio0,shape),cell_mask,i0,j0)
    ax = panelC;
    colorbar(ax);
    caxis(ax,'auto')
    ax.FontSize = fs;
    set(ax,'Color',[1 1 1]*1)
    title(ax, 'Rac', 'Fontsize', 24)
    axis(ax,'tight')
    
    
    plotCellIm(panelD,reshape(PaxRatio,shape),cell_mask,i0,j0)
    ax = panelD;
    caxis(ax,'auto')
    colorbar(ax)
    ax.FontSize = fs;
    set(ax,'Color',[1 1 1]*1)
    title(ax,'Pax', 'Fontsize', 24)
    axis tight
    
    

    
    
    title(panelA,['t=' num2str(time) ', t_{plot}=' num2str(double(tic-tp__0)*1e-6), ', t_{sim}=' num2str(toc)], 'Fontsize', 10)
    drawnow

end
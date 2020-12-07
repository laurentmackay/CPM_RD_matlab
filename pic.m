if plotting
%     close all
    % to make a video all frames must be the same size setting its position just
    % stops some bugs
%     figure('Position', [200 75 1000 900])
    fs=14; %axis font size
    
    subplot(2,2,1)
%     imagesc(cell_mask,[0 1]);
    plotCellIm(panelA,double(cell_mask),cell_mask,i0,j0)
%     plotCellIm(panelA,alpha_chem(:,:,5),cell_mask,i0,j0)
     plotCellIm(panelA,I_Ks,cell_mask,i0,j0)
    colorbar;
    title(['t=' num2str(time)], 'Fontsize', 24)
    hold on
    try
    plot(center(2,1:iter),center(1,1:iter),'r')
    catch e
        disp(e)
    end
    hold off
    ax = gca;
    ax.FontSize = fs;
    % title('Cell', 'Fontsize', 24)
    % xlabel('X')
    % ylabel('Y')
    
    
    subplot(2,2,2)
    plotCellIm(panelB,RhoRatio,cell_mask,i0,j0)
    caxis('auto')
    colorbar;
    ax = gca;
    ax.FontSize = fs;
    title('K', 'Fontsize', 24)
    % xlabel('X')
    % ylabel('Y')
    
    subplot(2,2,3)
    plotCellIm(panelC,RacRatio0,cell_mask,i0,j0)
    caxis('auto')
    colorbar
    ax = gca;
    ax.FontSize = fs;
    set(gca,'Color',[1 1 1]*1)
    title('Rac', 'Fontsize', 24)
    axis tight
    
    
    
     subplot(2,2,4)
    plotCellIm(panelD,PaxRatio,cell_mask,i0,j0)
    caxis('auto')
    colorbar
    ax = gca;
    ax.FontSize = fs;
    set(gca,'Color',[1 1 1]*1)
    title('Pax', 'Fontsize', 24)
    axis tight
    
    
    
    
    % xlabel('X')
    % ylabel('Y')
%     
%     subplot(2,2,4)
%     imagesc(PaxRatio,[0 0.4]);
%     colorbar
%     ax = gca;
%     ax.FontSize = fs;
%     title('Pax', 'Fontsize', 24)
    % xlabel('X')
    % ylabel('Y')
    
    %colormap jet
    %saveas(gcf,'CPM.png') %if you want a an image of a frame
    
    
%     Results(:,:,1,z)=cell_mask;
%     Results(:,:,2:(N_species+1),z)=x; %storing the results
    drawnow
    %adding videos the frame
    
    %frames=[frames getframe(gcf)];
    frame=getframe(gcf);
%     writeVideo(vid,frame);
end
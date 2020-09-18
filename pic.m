if usejava('desktop')
%     close all
    % to make a video all frames must be the same size setting its position just
    % stops some bugs
%     figure('Position', [200 75 1000 900])
    fs=14; %axis font size
    
    subplot(2,2,1)
    imagesc(cell_mask,[0 1]);
    colorbar;
    center(z,:)=com(cell_mask);
    hold on
    plot(center(1:z,2),center(1:z,1),'r')
    hold off
    ax = gca;
    ax.FontSize = fs;
    % title('Cell', 'Fontsize', 24)
    % xlabel('X')
    % ylabel('Y')
    
    
    subplot(2,2,2)
    plotCellIm(panelB,RhoRatio,cell_mask,i0,j0)
    colorbar;
    ax = gca;
    ax.FontSize = fs;
    title('Rho', 'Fontsize', 24)
    % xlabel('X')
    % ylabel('Y')
    
    subplot(2,2,3)
    plotCellIm(panelC,RacRatio,cell_mask,i0,j0)
    colorbar
    ax = gca;
    ax.FontSize = fs;
    set(gca,'Color',[1 1 1]*1)
    title('Rac', 'Fontsize', 24)
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
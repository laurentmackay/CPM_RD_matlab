if plotting
%     close all
    % to make a video all frames must be the same size setting its position just
    % stops some bugs
%     figure('Position', [200 75 1000 900])
    fs=14; %axis font size
  for i=1:length(h_vec)

        plotCellIm(h_vec{i},x(:,:,i)/scaling,cell_mask,i0,j0)
        caxis('auto')
        colorbar(h_vec{i})
        ax = gca;
        ax.FontSize = fs;
        set(gca,'Color',[1 1 1]*1)
    %     title('R_i', 'Fontsize', 24)
        axis tight
  end
    drawnow

end
%-----------This is for extracting video from a result mat file----------%
%NOTE: variables Results and center must be stored properly in the mat file

close all
clear
set(0,'DefaultFIgureVisible','on');

wksp=strcat('result1.mat');% change to the mat file you wanna extract from
load(wksp,'Results','center');

vid = VideoWriter(['results'],'MPEG-4');
open(vid);

nf=size(Results,4);%number of timepoints/frames

for i=1:nf
    close
    % to make a video all frames must be the same size setting its position just
    % stops some bugs
    figure('Position', [200 200 1000 900]) 
    fs=14; %axis font size 

    cell_mask=Results(:,:,1,i);
    RacRatio=Results(:,:,3,i)./sum(Results(:,:,[3 5 8],i),3);
    RhoRatio=Results(:,:,4,i)./sum(Results(:,:,[2 4],i),3);
    PaxRatio=Results(:,:,6,i)./sum(Results(:,:,[6 7 9],i),3);

    subplot(2,2,1)
    imagesc(cell_mask,[0 1]);
    % colorbar;
    hold on
    plot(center(1:i,2),center(1:i,1),'r')
    hold off
    ax = gca;
    ax.FontSize = fs;
    % title('Cell', 'Fontsize', 24)
    % xlabel('X')
    % ylabel('Y')


    subplot(2,2,2)
    imagesc(RhoRatio,[0.1 0.6]);
    colorbar;
    ax = gca;
    ax.FontSize = fs;
    title('Rho', 'Fontsize', 24)
    % xlabel('X')
    % ylabel('Y')

    subplot(2,2,3)
    imagesc(RacRatio,[0 0.3]);
    colorbar
    ax = gca;
    ax.FontSize = fs;
    title('Rac', 'Fontsize', 24)
    % xlabel('X')
    % ylabel('Y')

    subplot(2,2,4)
    imagesc(PaxRatio,[0 0.4]);
    colorbar
    ax = gca;
    ax.FontSize = fs;
    title('Pax', 'Fontsize', 24)
    % xlabel('X')
    % ylabel('Y')

    %colormap jet
    %saveas(gcf,'CPM.png') %if you want a an image of a frame 


    drawnow
    %adding videos the frame 

    %frames=[frames getframe(gcf)];
    frame=getframe(gcf);
    writeVideo(vid,frame);
end
close(vid)
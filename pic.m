% close all
% to make a video all frames must be the same size setting its position just
% stops some bugs
% figure('Position', [200 200 1000 900]) 

fs=14; %axis font size 
% figure(1);

if time==0
    rac_induced=[]
    rac_uninduced=[];
end
    rac_induced=[rac_induced [time; mean(RacRatio(induced_mask&cell_mask))]];
    rac_uninduced=[rac_uninduced [time; mean(RacRatio(~induced_mask&cell_mask))]];
    

htitle=suptitle(['t = ' num2str(time)]);
set(htitle,'FontSize',28)

plot(panelA,rac_induced(1,:),rac_induced(2,:),rac_uninduced(1,:),rac_uninduced(2,:))

% imagesc(panelB,RhoRatio);
imagesc(panelB,RhoRatio);
colorbar(panelB);
caxis(panelB,[0, 0.6])

% panelB.FontSize = fs;
title(panelB,'Rho', 'Fontsize', 24)
% xlabel('X')
% ylabel('Y')

imagesc(panelC,RacRatio);
% imagesc(panelC,x(:,:,2));
colorbar(panelC)
caxis(panelC,[0, 0.3])
% ax = gca;
% panelC.FontSize = fs;
title(panelC,'Rac', 'Fontsize', 24)
% xlabel('X')
% ylabel('Y')


imagesc(panelD,PaxRatio);
% imagesc(panelD,x(:,:,4));
colorbar(panelD)
caxis(panelD,[0, 0.4])

panelD.FontSize = fs;
title(panelD,'Pax', 'Fontsize', 24)
% xlabel('X')
% ylabel('Y')

%colormap jet
%saveas(gcf,'CPM.png') %if you want a an image of a frame 

drawnow

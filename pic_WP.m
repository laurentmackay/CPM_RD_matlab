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

% plot(panelA,rac_induced(1,:),rac_induced(2,:),rac_uninduced(1,:),rac_uninduced(2,:))
Rac_profile=zero2nan(mean(RacRatio,2));
Results=[Results Rac_profile];
N_thick=5;
N_thin=15;
N_lines=size(Results,2);
Times=[Times time];
i_thin=[];

if N_lines>1
i_thick=unique(round(linspace(2,N_lines,5)));

for i=1:length(i_thick)-1
    temp=round(linspace(i_thick(i),i_thick(i+1),floor(N_thin/(N_thick-1))+2));
    i_thin=[i_thin temp(2:end-1)];
end
i_thin=unique(i_thin);
else
    i_thick=[];
end
%%

yaxis=(0:N-1)*h;
hA=plot(panelA,yaxis,Results(:,[1 i_thin i_thick]));


set(hA(1),'LineWidth',3,'LineStyle','--','Color','k');
cstart=0.75;
cend=0;

[~,colours]=meshgrid([1 1 1],linspace(cstart,cend,length(i_thick)));
for i=1:length(i_thick)
set(hA(i+length(i_thin)+1),'LineWidth',2,'LineStyle','-','Color',colours(i,:),'Marker','.','MarkerSize',15);
end

[~,colours]=meshgrid([1 1 1],linspace(cstart,cend,length(i_thin)));
for i=1:length(i_thin)
%      alpha(hA(i+1),0.5)
    set(hA(i+1),'LineWidth',0.01,'LineStyle','-','Color',[colours(i,:) 0.5]);
end

%%

% set(hA(2:(length(i_thick)+1)),'LineWidth',1,'LineStyle','--','Color',colours);
% if length(hA)>1
% 
%     i=2:length(hA);
%     i_thick=unique(round(linspace(2,N_lines,5)));
%     i_thin=setdiff(i,i_thick);
%     set(hA(i_thin),'LineWidth',0.01,'Color',[1 1 1]*0.85);
%     set(hA(i_thick),'LineWidth',1);
% end

% if time==0
% 
% else
%     hold(panelA,'on')
%     hA=plot(panelA,(0:N-1)*h,zero2nan(mean(RacRatio,2)),'-k');
%     set(hA,'LineWidth',0.5);
%     hold(panelA,'off')
% end
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

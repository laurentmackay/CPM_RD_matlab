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
    

% plot(panelA,rac_induced(1,:),rac_induced(2,:),rac_uninduced(1,:),rac_uninduced(2,:))
Rac_profile=zero2nan(mean(RacRatio,2));
Results=[Results Rac_profile];
N_thick=8;
N_thin=15;
N_lines=size(Results,2);
Times=[Times time];
i_thin=[];


rmin=min(Results);
rmax=max(Results);

thresh=(rmin+rmax)/2;
above=Results>thresh;

pin=zeros(size(Results,2),1);
for i=1:size(above,2)
    ind=find(~above(2:end,i),1);
    
    
    pin(i)=ind+(thresh(i)-Results(ind+1,i))/(Results(ind+1,i)-Results(ind,i));
    
end



if N_lines>1
    
    [pmin,imin]=min(pin);
    [pmax,imax]=max(pin);
    thresh=pmin+(1/N_thick:1/N_thick:1)*(pmax-pmin);
    i_thick=zeros(1,N_thick);
    for i=1:N_thick
        ind=find(pin(imin:imax)>thresh(i),1);
        if ~isempty(ind)
            i_thick(i)=ind;
        else
            i_thick(i)=imax;
        end
    end
    
i_thick=unique(i_thick);

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
    set(hA(i+1),'LineWidth',0.01,'LineStyle','-','Color',[colours(i,:) 0.5]);
end
axis(panelA,'tight');
xlabel(panelA,'Y (\mum)');
ylabel(panelA,'Scaled Rac');
%%


% panelB.FontSize = fs;
plot(panelB,Times,pin,'-k')
ylabel(panelB,'Wave Front')
xlabel(panelB,'Time')

axis(panelB,[0 time, 0 round(max(pin))+1])



imagesc(panelD,(1:shape(1))*h,(1:shape(2))*h,RacRatio);
colorbar(panelD)
caxis(panelD,[0, 0.3])
% ax = gca;
% panelC.FontSize = fs;
title(panelD,['Rac t=' num2str(time,6) ' s'], 'Fontsize', 18)
xlabel(panelD,'X (\mum)')
ylabel(panelD,'Y (\mum)')


imagesc(panelC,(1:shape(1))*h,(1:shape(2))*h,RacRatio0);
% imagesc(panelD,x(:,:,4));
colorbar(panelC)
caxis(panelC,[0, 0.3])

panelC.FontSize = fs;
title(panelC,'Rac t=0 s', 'Fontsize', 18)
xlabel(panelC,'X (\mum)')
ylabel(panelC,'Y (\mum)')

gca=panelC;


%colormap jet
%saveas(gcf,'CPM.png') %if you want a an image of a frame 

drawnow

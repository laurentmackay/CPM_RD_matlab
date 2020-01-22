% close all
% figure('Position', [400 200 700 700])
fs=14;
cell_mask=x(:,:,1);
subplot(2,2,1)
imagesc(cell_mask,[0 1]);
colorbar;
center(z,:)=com(cell_mask);
hold on
plot(center(1:z,1),center(1:z,2),'r')
hold off
ax = gca;
ax.FontSize = fs;
title('Cell', 'Fontsize', 24)
xlabel('X')
ylabel('Y')


subplot(2,2,2)
%imagesc(x(:,:,4),[0 ceil(.75*totalRho/a)]);
imagesc(squeeze(RhoRatio),[0.1 0.6]);
colorbar;
ax = gca;
ax.FontSize = fs;
title('Rho', 'Fontsize', 24)
xlabel('X')
ylabel('Y')

subplot(2,2,3)
imagesc(squeeze(RacRatio),[0 0.3]);
%imagesc(x(:,:,5),[0 ceil(totalRac/a)]);
colorbar
ax = gca;
ax.FontSize = fs;
title('Rac', 'Fontsize', 24)
xlabel('X')
ylabel('Y')

subplot(2,2,4)
imagesc(squeeze(PaxRatio),[0 0.3]);
%imagesc(x(:,:,7),[0 ceil(.5*totalPax/a)]);
colorbar
ax = gca;
ax.FontSize = fs;
title('Pax', 'Fontsize', 24)
xlabel('X')
ylabel('Y')
set(gcf,'position',[100 100 1100 1000]) 


%colormap jet
% Results(:,:,:,z)=x;
drawnow
% frame=getframe(gcf);
% writeVideo(vid,frame);



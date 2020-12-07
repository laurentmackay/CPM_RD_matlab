file='final_B_5.mat';

n=5;
% m1=45;
% m2=65
% s2=200;
% v0=get_instant_velocity(file,1:m1,n);
% vf=get_instant_velocity(file,s2-m2:s2,n);

vtot = get_instant_velocity(file,[],n);
acf_tot = get_vel_acf_windowed(vtot,12,3); 

figure(1);clf();
plot(acf_tot)

 figure(2);clf();plot(getHalflife(acf_tot))
 

% 
% 
% s0=sqrt(sum(v0.^2,1));
% sf=sqrt(sum(vf.^2,1));
% 
% nbin=25;
% edges=linspace(0,max(max(sf),max(s0)),nbin);
% norm='Probability';
% 
% 
% [h0,e0]=histcounts(s0,edges,'Normalization',norm);
% [hf,ef]=histcounts(sf,edges,'Normalization',norm);
% b0=e0(1:end-1)+diff(e0)/2;
% bf=ef(1:end-1)+diff(ef)/2;
% % figure(1);clf();
% % plot(b0,h0,bf,hf);
% 
% 
% figure(2);clf();
% plot(s0);
% hold on
% plot(sf)
% hold off
% 
% 
% 
% figure(3);clf();
% [acf,lags]=get_vel_acf(v0)
% plot(lags,acf)
% hold on
% [acf,lags]=get_vel_acf(vf)
% plot(lags,acf)
% hold off
file='../_chem_Rx_Pax_Kathy/results/random_walk_B_sweep/final_B_2.25_copy1.mat';
file='../_chem_Rx_Pax_Kathy/results/random_walk_B_sweep/final_B_6.2586_copy2.mat';

% file='results/final_B_5.mat';
n=3;
% m1=45;
% m2=65
% s2=200;
% v0=get_instant_velocity(file,1:m1,n);
% vf=get_instant_velocity(file,s2-m2:s2,n);

vtot = get_instant_velocity(file,[],n);

%%

win_len=10;
[acf_tot,~,t_acf] = get_vel_acf_windowed(vtot,win_len,1); 

% figure(4);

xapprox=cumsum(vtot,2)';%plot(xapprox(:,1),xapprox(:,2))

cpmstep=1.3;
cpmsteps=15;

cpm_time= n*cpmstep*cpmsteps;
t_vec=cpm_time*(1:size(vtot,2));
t_acf = cpm_time*t_acf;

% figure(3);clf();
% plot(t_vec,sqrt(sum(vtot.^2)));
% xlabel('MCS')
% ylabel('instantaneous velocity')


% v_com = vtot(1,:) + 1j *vtot(2,:);
 
% corrgram(v_com,v_com);
% acf_2 = acf_2(floor(end/2)+1:end,:);
% 


% figure(1);clf();
% subplot(2,1,1);
% plot(acf_tot);

% subplot(2,1,2);
% plot(acf_2)

 figure(2);clf();
 subplot(2,1,1);
 
plot(t_vec,sqrt(sum(vtot.^2)));
xlabel('MCS')
ylabel('instantaneous velocity')

subplot(2,1,2);
plot(t_acf(1:end),getHalflife(acf_tot(:,1:end)))
xlabel('time');
ylabel('ACF Halflife')
yline(1.5)
% subplot(2,1,2);

% plot(getHalflife(acf_2,0.5))

% 
% 
% s0=sqrt(sum(v0.^2,1));
% sf=sqrt(sum(vf.^2,1));
% % 
% % nbin=25;
% % edges=linspace(0,max(max(sf),max(s0)),nbin);
% % norm='Probability';
% % 
% % 
% [h0,e0]=histcounts(s0,edges,'Normalization',norm);
% [hf,ef]=histcounts(sf,edges,'Normalization',norm);
% b0=e0(1:end-1)+diff(e0)/2;
% bf=ef(1:end-1)+diff(ef)/2;
% % % figure(1);clf();
% plot(b0,h0,bf,hf);
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
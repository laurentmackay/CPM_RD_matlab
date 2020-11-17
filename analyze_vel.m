file='final_B_5.mat';

n=round(2.5*15);

v0=get_instant_velocity(file,1:1000,n);
vf=get_instant_velocity(file,3000:5000,n);

s0=sqrt(sum(v0.^2,1));
sf=sqrt(sum(vf.^2,1));

nbin=25;
edges=linspace(0,max(max(sf),max(s0)),nbin);
norm='Probability';


[h0,e0]=histcounts(s0,edges,'Normalization',norm);
[hf,ef]=histcounts(sf,edges,'Normalization',norm);
b0=e0(1:end-1)+diff(e0)/2;
bf=ef(1:end-1)+diff(ef)/2;
figure(1);clf();
plot(b0,h0,bf,hf);


figure(2);clf();
plot(s0);
hold on
plot(sf)
hold off



figure(3);clf();
plot(get_vel_acf(v0))
hold on
plot(get_vel_acf(vf))
hold off
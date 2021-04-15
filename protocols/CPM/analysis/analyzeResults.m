function [outputArg1,outputArg2] = analyzeResults(f,iters)

%%load in the results
load(['results/' f]);
iter=min(iter,length(Times));

if nargin==1
    iters=1:iter;
end

if iters(1)>1
T0=Times(iters(1)-1);
else
T0=0;
end
center=center(:,iters);
Times=Times(iters);
areas=areas(iters);
perims=perims(iters);

%% initial compuation
lags=1:iter-1;
dr_mean=zeros(size(lags));
dt_mean=zeros(size(lags));
drs={};
dts={};
for l=lags
    drs{l}=sqrt(sum((center(:,l+1:end)-center(:,1:end-l)).^2,1));
    dts{l}=Times(l+1:end)-Times(1:end-l);
    
    dr_mean(l)=mean(drs{l});
    dt_mean(l)=mean(dts{l});
end

dr=sqrt(sum((center(:,1)-center(:,2:end)).^2,1));
dt=Times(2:end)-T0;

figure(1);clf();


%% dispersion exponent
subplot(2,3,1);

% loglog([dts{:}],[drs{:}],'.k','MarkerSize',0.05/log(length(lags)))
% hold on
hplot=loglog([dt(1) dt(end)],[dr(1) dr(1)*dt(end)^0.5],...
    [dt(1) dt(end)],[dr(1) dr(1)*dt(end)],...
    dt,dr);
set(hplot,'LineWidth',2)
% set(hplot(1),'Color','r');
set(hplot(1),'Color','g','LineStyle',':');
set(hplot(2),'Color','g','LineStyle','-.');
set(hplot(3),'Color','g');
% hold off


xlim([dt(1) dt(end)]);
ylim([min(dr), max(dr)]);
xlabel('\Deltat')
ylabel('\Deltar (\mum)');
legend(hplot,{'Brownian', 'Ballistic', 'Random Walk'},'Location','Best')

%% displacement distributions
v_insta = [drs{1}]./[dts{1}];
v_insta2 = [drs{1:2}]./[dts{1:2}];
n0=3;

% sqrt(sum((center(:,1:n0:end-n0)-center(:,1+n0:n0:end)).^2,1))./diff(Times)
v_inst = @(i) [drs{1:n0:i}]./[dts{1:n0:i}];

subplot(2,3,2);
hold on
for i=1:1
    [counts,bin_centers]=histcounts(v_inst(i),50,'Normalization','Probability');
    bin_centers=bin_centers(1:end-1)+diff(bin_centers);
    plot(bin_centers,counts);
end
hold off


% [counts0,bin_centers0]=histcounts(v_insta,50,'Normalization','Probability');
% bin_centers0=bin_centers0(1:end-1)+diff(bin_centers0)/2;
% [counts,bin_centers]=histcounts([drs{1:2}]./[dts{1:2}],50,'Normalization','Probability');
% bin_centers=bin_centers(1:end-1)+diff(bin_centers);
% plot(bin_centers,counts,bin_centers0,counts0);


%% persistence time
subplot(2,3,3);

vxx=autocorr(v_inst(1));
plot(vxx)
axis tight


%% cell geometry
% perims=zeros(1,iter);
% areas=zeros(1,iter);
% for i=1:iter
%     cell_mask=Results(:,:,1,i);
%     perims(i)=perim(cell_mask);
%     areas(i)=nnz(cell_mask);
% end

% perims(1)=per;
% areas(1)=a;

subplot(2,3,4);

plot(Times,perims,[Times(1) Times(end)],[per per]);
title('Cell Perimeter')

subplot(2,3,5);
plot(Times,areas,[Times(1) Times(end)],[a a]);
title('Cell Area')
subplot(2,3,6);

% H0=lam_a*(a-areas).^2+lam_p*(per-perims).^2+J*perims;
plot(Times,Ham0(iters));
title('H_0')

%% plotting traj
figure(2); clf();
plot(center(1,:),center(2,:))
end


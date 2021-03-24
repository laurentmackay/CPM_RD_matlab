
t0=6000;
tf=4e4;

thresh=0.05;

figure(34);
hold(gca,'on')

exp_coeff={}    ;




B_vals = {};
set_experiment('lam_p_0_Inf')




results = ls_results()';
for result = results
    B_str = regexp(result.name,'B_([^\_]+)','tokens');
    B_vals{end+1}=str2num(B_str{1}{1});
end

B_vals_res = cell2mat(B_vals);
[B_vals, ~ ,inds_B] = unique(B_vals_res);
disp([ int2str(length(inds_B)) ' results total, ' int2str(length(B_vals)) ' unique B values'])



thresh_up=cell(length(B_vals),1);
thresh_down=cell(length(B_vals),1);
stpnt=cell(length(B_vals),1);

exp_coeff=cell(length(B_vals),1);


for i=1:length(inds_B)
    
    [t,y]=timecourse(results(i).name,"[min(min(x(inds))), max(max(x(inds)))]",...
        "i_rac = find(strcmp(chems,'Rac')); inds=cell_inds(1:A)+sz*(i_rac-1); ");
    delta=abs( y(:,1)-y(:,2));
    
    
    if max(delta)>0.05
        
        inds2=t>t0;
        i0=find(inds2,1);
        dt=mean(diff(t));
        lambda_estimate=([diff(delta(inds2));  Inf]/1e2)./delta(inds2);
        %         figure(55);semilogy(lambda_estimate);
        metric=abs(log(lambda_estimate(20))-log(lambda_estimate));
        i_thresh = i0 + find(metric<0.05,1,'last');
        
        
        
%         figure(56); semilogy(metric)
%         
        inds=t>t0 & t<t(i_thresh);
        
        log_delta=log(delta);
        
        
        lin_fit = fit(t(inds)',log(delta(inds)),'poly1');
        exp_coeff{inds_B(i)}{end+1}=lin_fit.p1;
        
        i_thresh = i0 + find(abs(log(lin_fit.p1)-log(lambda_estimate))>0.05,1);
        
%         figure(34);
%         semilogy(t,y)
%         semilogy(t(i_thresh),y(i_thresh,1), 'rx')
%         semilogy(t(i_thresh),y(i_thresh,2), 'rx')
        
        thresh_down{inds_B(i)}{end+1}=y(i_thresh,1);
        thresh_up{inds_B(i)}{end+1}=y(i_thresh,2);
        
        stpnt{inds_B(i)}{end+1}=mean(y(1,:));
        
        
%         figure(55);semilogy(t(inds2),lambda_estimate, t(i_thresh), lambda_estimate(i_thresh-i0), 'rx');
        
    end
end

figure(34);
hold(gca,'off')

%%
mean_thresh_up=cell(length(B_vals),1);
mean_thresh_down=cell(length(B_vals),1);
std_thresh_up=cell(length(B_vals),1);
std_thresh_down=cell(length(B_vals),1);
mean_stpnt=cell(length(B_vals),1);
mean_exp_coeff=cell(length(B_vals),1);

for i=1:size(thresh_up,1)
    
    if length(thresh_up{i})>0
        mean_thresh_up{i}=mean([thresh_up{i}{:}]);
        mean_thresh_down{i}=mean([thresh_down{i}{:}]);
        std_thresh_up{i}=std([thresh_up{i}{:}]);
        std_thresh_down{i}=std([thresh_down{i}{:}]);
        mean_stpnt{i}=mean([stpnt{i}{:}]);
        mean_exp_coeff{i}=mean([exp_coeff{i}{:}]);
    else
        mean_thresh_up{i}=NaN;
        mean_thresh_down{i}=NaN;
        std_thresh_up{i}=0;
        std_thresh_down{i}=0;
        mean_stpnt{i}=NaN;
        mean_exp_coeff{i}=NaN;
    end
    
end

mean_thresh_up=cell2mat(mean_thresh_up);
mean_thresh_down=cell2mat(mean_thresh_down);
mean_exp_coeff = cell2mat(mean_exp_coeff);

figure(35);
subplot(1,2,1);
plot(B_vals,cell2mat(mean_stpnt))


hold on
errorbar(B_vals, mean_thresh_up, cell2mat(std_thresh_up))
errorbar(B_vals, mean_thresh_down, cell2mat(std_thresh_down))
hold off

subplot(1,2,2);
plot(B_vals,mean_exp_coeff)

%%

figure(36); 
actual_time= 1.0e+05 *[ Inf
    1.3823
    0.9703
    0.7553
    0.6577
    0.5253
    0.4863
    0.4267
    0.4073
    0.3620
    0.3163
    0.3083
    0.2787
    0.2780
    0.2460
    0.2423
    0.2437
    0.2160
    0.1927
    0.1987]



plot(B_vals,actual_time,B_vals,9000-log(mean_thresh_up-mean_thresh_down)./mean_exp_coeff, B_vals,9000-mean(log(mean_thresh_up-mean_thresh_down),'omitnan')./mean_exp_coeff)
ylim([0, 2e5]);
% set(gca, 'Yscale','Log')



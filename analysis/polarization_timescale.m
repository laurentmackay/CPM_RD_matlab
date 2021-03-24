    
t0=4000;
tf=4e4;

thresh=0.05;

    figure(34);
    hold(gca,'on')

exp_coeff={}    


set_experiment('lam_p_0_1')

thresh_up={};
thresh_down={};


for j=1:3

[t,y]=timecourse(['final_B_2.4105_copy' int2str(j) '.mat'],"[min(min(x(inds))), max(max(x(inds)))]",...
        "i_rac = find(strcmp(chems,'Rac')); inds=cell_inds(1:A)+sz*(i_rac-1); ");
        delta=abs( y(:,1)-y(:,2));
        
        

        
        inds2=t>t0;
        i0=find(inds2,1);
        dt=mean(diff(t));
        lambda_estimate=([diff(delta(inds2));  Inf]/1e2)./delta(inds2);
        figure(55);semilogy(lambda_estimate);
       
        i_thresh = i0 + find(abs(log(lambda_estimate(1))-log(lambda_estimate))>0.05,1);
        

        
        figure(56); semilogy(abs(log(lambda_estimate(1))-log(lambda_estimate)))
        
        inds=t>t0 & t<t(i_thresh);
        
        log_delta=log(delta);
        
        
        lin_fit = fit(t(inds)',log(delta(inds)),'poly1');
        exp_coeff{end+1}=lin_fit.p1;
        
        i_thresh = i0 + find(abs(log(lin_fit.p1)-log(lambda_estimate))>0.05,1);
        
        figure(34);
        semilogy(t,delta)
        semilogy(t(i_thresh),delta(i_thresh), 'rx')
%         semilogy(t(i_thresh),y(i_thresh,2), 'rx')
        
%         semilogy(t,y)
%         semilogy(t(i_thresh),y(i_thresh,1), 'rx')
%         semilogy(t(i_thresh),y(i_thresh,2), 'rx')
        
        thresh_down{end+1}=y(i_thresh,1)
        thresh_up{end+1}=y(i_thresh,2)
        
        
        
end

    figure(34);
    hold(gca,'off')

    
    

    

    

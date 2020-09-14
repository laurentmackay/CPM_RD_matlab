function [x,alpha_rx,alpha_chem,time,PaxRatio,RhoRatio,K_is,K,RacRatio,...
            I_Ks,reaction,jump,ir0,id0,cell_inds,k_X,PIX,k_G,k_C,GIT,...
            Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,alpha,...
            PAKtot,rx_count,dt_diff] =  CPM_chem_func(x,alpha_rx,...
            alpha_chem,time,PaxRatio,RhoRatio,K_is,K,rx_speedup,RacRatio,...
            I_Ks,reaction,jump,ir0,id0,cell_inds,k_X,PIX,k_G,k_C,GIT,...
            Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,alpha,...
            PAKtot,nrx,A,pi,pT0,dt_diff,rx_count,P_diff,Rac_Square,...
            Pax_Square,Rho_Square)
N_instantaneous=50; % the number of complexing reactions iterated at a moment
sz=size(x,1)*size(x,2);
% rng(round(mod(cputime,1e-2)/1e-5))

%variables declared so that the c encoder is happy
vox=cell_inds(1);
neg=any(x<0);

mpT0=max(pT0);

for kk=1:nrx    
    %---------------------Activation/Inactivation------------------------%
    %total propensity
    a_total=sum(alpha_rx);
    
    tau  = (1/rx_speedup)*(1/a_total)*log(1/rand()); % time increment
    time=time+tau;
    reacted=false;
    neg=false;
    RN=rand();
    temp=0;
    
    %determine what reaction and where the reaction happens
    for rx=1:size(alpha_chem,3)
        i0=(rx-1)*sz;
        if alpha_rx(rx)>=RN*a_total-temp&&~reacted
            %find the lattice square where reaction happens
            ii=Alg2(alpha_chem(i0+cell_inds(1:A))/a_total,RN,temp/a_total);
            vox=cell_inds(ii(1));
            reacted = true;
            rx_count(vox)=rx_count(vox)+1;
        elseif ~reacted
            temp=temp+alpha_rx(rx);
        end
        if reacted
            break;
        end
    end

    if reacted==false
        error('Oh no! D: chemical reaction propensites did not sum correctly')
    end
    
    %Inactive rho to active rho
    if rx==1
        x(vox+(3-1)*sz) = x(vox+(3-1)*sz)+1;
        x(vox+(1-1)*sz) = x(vox+(1-1)*sz)-1;
    end

    %Inactive Rac to active Rac
    if rx==2
        x(vox+(2-1)*sz)=x(vox+(2-1)*sz)-1;
        x(vox+(4-1)*sz)=x(vox+(4-1)*sz)+1;
    end

    %Active rho to inactive rho
    if rx==3
        x(vox+(3-1)*sz) = x(vox+(3-1)*sz)-1;
        x(vox+(1-1)*sz) = x(vox+(1-1)*sz)+1;
    end

    %Active Rac to inactive Rac
    if rx==4
        x(vox+(2-1)*sz)=x(vox+(2-1)*sz)+1;
        x(vox+(4-1)*sz)=x(vox+(4-1)*sz)-1;
    end

    %Unphosphorylated Pax to phosphorylated Pax
    if rx==5
        x(vox+(5-1)*sz)=x(vox+(5-1)*sz)-1;
        x(vox+(6-1)*sz)=x(vox+(6-1)*sz)+1;
    end

    %Phosphorylated Pax to unphosphorylated
    if rx==6
        x(vox+(5-1)*sz)=x(vox+(5-1)*sz)+1;
        x(vox+(6-1)*sz)=x(vox+(6-1)*sz)-1;
    end


    %---------------------Complexing & Decomplexing----------------------%

    if rx==2||rx==4%Rac
        K_R=(1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(vox))*alpha*PAKtot*K_is(vox);
        for j=1:N_instantaneous
            if x(vox+(7-1)*sz)/(x(vox+(7-1)*sz)+x(vox+(4-1)*sz)*K_R)>rand()%decomplex
                x(vox+(4-1)*sz)=x(vox+(4-1)*sz)+1;
                x(vox+(7-1)*sz)=x(vox+(7-1)*sz)-1;
            elseif x(vox+(4-1)*sz)>0
                x(vox+(4-1)*sz)=x(vox+(4-1)*sz)-1;
                x(vox+(7-1)*sz)=x(vox+(7-1)*sz)+1;
            end
        end
    elseif rx==5||rx==6%Pax
        K_P=k_G*k_X*k_C*GIT*PIX*K_is(vox)*PAKtot*(1+alpha_R*RacRatio(vox));
        for j=1:N_instantaneous
            if x(vox+(8-1)*sz)/(x(vox+(8-1)*sz)+x(vox+(6-1)*sz)*K_P)>rand()%decomplex
                x(vox+(6-1)*sz)=x(vox+(6-1)*sz)+1;
                x(vox+(8-1)*sz)=x(vox+(8-1)*sz)-1;
            elseif x(vox+(6-1)*sz)>0
                x(vox+(8-1)*sz)=x(vox+(8-1)*sz)+1;
                x(vox+(6-1)*sz)=x(vox+(6-1)*sz)-1;
            end
        end
    end
    neg=x(vox+(rx-1)*sz)<0;

    %-----------recalculating value that would have changed--------------%

    if length(vox)>1
    [tmp,tmp2]=meshgrid(ir0,cell_inds(1:A));
        I_rx=tmp+tmp2;
    else
        I_rx=vox+ir0;
    end

    a_c_0=alpha_chem(I_rx);

    RacRatio(vox)=x(vox+(4-1)*sz)./Rac_Square;
    RhoRatio(vox)=x(vox+(3-1)*sz)./Rho_Square;
    PaxRatio(vox)=x(vox+(6-1)*sz)./Pax_Square;
        
    K_is(vox)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(vox)).*(1+alpha_R*RacRatio(vox))+k_G*k_X*GIT*PIX);
    K(vox)=alpha_R*RacRatio(vox).*K_is(vox).*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio(vox));  
    I_Ks(vox)=I_K*(1-K_is(vox).*(1+alpha_R*RacRatio(vox)));

    reaction(vox+(1-1)*sz) = I_rho*(L_R^m./(L_R^m +(RacRatio(vox)+gamma*K(vox)).^m));       
    reaction(vox+(2-1)*sz) = (I_R+I_Ks(vox)).*(L_rho^m./(L_rho^m+RhoRatio(vox).^m)); 
    reaction(vox+(5-1)*sz) = B_1*(K(vox).^m./(L_K^m+K(vox).^m));

    alpha_chem(I_rx) = reaction(I_rx).*x(I_rx); %chemical reaction
    alpha_rx=alpha_rx+sum(alpha_chem(I_rx)-a_c_0,1);
    a_total=sum(alpha_rx);
    
    if neg
        error('Oh no! D: (negtive numbers): ')
    end
    
    %-----------------------------Diffusion-------------------------------%
    dt_diff=dt_diff+tau;
    ind_diff=mpT0.*dt_diff>P_diff; %determine indices for bulk diffusions
    
    if nnz(ind_diff)>0
        
        %wow diffusing
        x=Alg3(x,dt_diff,jump',pT0,pi,cell_inds,A,ind_diff);
        
        dt_diff(ind_diff)=0;
        vox=cell_inds(1:A);
        
        %---------recalculating value that would have changed------------%

        if length(vox)>1
        [tmp,tmp2]=meshgrid(ir0,cell_inds(1:A));
            I_rx=tmp+tmp2;
        else
            I_rx=vox+ir0;
        end
        
        a_c_0=alpha_chem(I_rx);
        
        
        RacRatio(vox)=x(vox+(4-1)*sz)./Rac_Square;
        RhoRatio(vox)=x(vox+(3-1)*sz)./Rho_Square;
        PaxRatio(vox)=x(vox+(6-1)*sz)./Pax_Square;
        
        K_is(vox)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(vox)).*(1+alpha_R*RacRatio(vox))+k_G*k_X*GIT*PIX);
        K(vox)=alpha_R*RacRatio(vox).*K_is(vox).*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio(vox));         
        I_Ks(vox)=I_K*(1-K_is(vox).*(1+alpha_R*RacRatio(vox)));
        
        reaction(vox+(1-1)*sz) = I_rho*(L_R^m./(L_R^m +(RacRatio(vox)+gamma*K(vox)).^m));            
        reaction(vox+(2-1)*sz) = (I_R+I_Ks(vox)).*(L_rho^m./(L_rho^m+RhoRatio(vox).^m));            
        reaction(vox+(5-1)*sz) = B_1*(K(vox).^m./(L_K^m+K(vox).^m));
        
        alpha_chem(I_rx) = reaction(I_rx).*x(I_rx); %chemical reaction
        alpha_rx=alpha_rx+sum(alpha_chem(I_rx)-a_c_0,1);
        a_total=sum(alpha_rx);
        
    end
    
    
end
end
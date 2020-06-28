function [x,diffusing_species_sum,D,h,alpha_rx,num_diffuse,...
    alpha_chem,time,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
    RacRatio,RbarRatio,I_Ks,reaction,jump,ir0,id0,cell_inds,...
    k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
    alpha,PAKtot,rx_count,dt_diff] =  CPM_chem_func(x,diffusing_species_sum,D,h,alpha_rx,num_diffuse,...
    alpha_chem,time,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
    RacRatio,RbarRatio,I_Ks,reaction,jump,ir0,id0,cell_inds,...
    k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
    alpha,PAKtot,nrx,A,num_vox_diff,pi,pT0,dt_diff,rx_count,P_diff,rx_speedup)



N_instantaneous=50; % the number of steady reaction itterated at a moment in time
sz=size(x,1)*size(x,2);
% rng(round(mod(cputime,1e-2)/1e-5))
%variables declared so that the c encoder is happy
vox=cell_inds(1);
p=1;
i2=cell_inds(1);
rx=1;
vox0=id0+vox;
xi0=x(id0+vox);
i_update=[vox i2];
neg=any(x<0);

diff_err=0.01;
dt_max=sqrt(diff_err)*(h^2)/max(D);

a_total_new=sum(alpha_rx);


% dt_diff=-log(p_diff)/a_total_new;
% dt_diff=min(dt_diff,dt_max);

mpT0=max(pT0);

do_rx=true;



for kk=1:nrx
    %      alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
    
    
    
    %total propensity for rxn+diff
    a_total= a_total_new;
    
    tau  = (1/rx_speedup)*(1/a_total)*log(1/rand()); % time increment
    
    time=time+tau;
    
    neg=false;
    RN=rand();
    
    reacted=false;
    if do_rx
        temp=0;
        %determing where the reaction happens
        for rx=1:size(alpha_chem,3)
            i0=(rx-1)*sz;
            if temp+alpha_rx(rx)>=RN*a_total&&~reacted
                %                 ii=find(cumsum(alpha_chem(i0+cell_inds(1:A)))>=RN*a_total-temp,1);
                ii=Alg22(alpha_chem(i0+cell_inds(1:A))/a_total,RN,temp/a_total);
                vox=cell_inds(ii(1));
                vox0=id0+vox;
                xi0=x(vox0);
                reacted = true;
                %             if rx==2 || rx==1 || rx==5
                rx_count(vox)=rx_count(vox)+1;
                %             end
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
        %             numReac=numReac+1;%ellie
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
        %     I moved the following 2 lines into my dependence tree due to
        %     complexing (should uncomment if replicating this algorthim)
        %         dxi=xi0-x(id0+i(1));
        %         diffusing_species_sum = diffusing_species_sum - (diffuse_mask(:,i)*dxi);
        %         neg=x(i+(rx-1)*sz)<0;
        
        %-----------recalculating value that would have changed--------------
        
        %deal with complexing and decomplexing
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
        dxi=xi0-x(vox0);
        diffusing_species_sum = diffusing_species_sum - (diffuse_mask(:,vox)*dxi);
        neg=x(vox+(rx-1)*sz)<0;
        
        
if length(vox)>1
[tmp,tmp2]=meshgrid(ir0,vox);
    I_rx=tmp+tmp2;
else
    I_rx=vox+ir0;
end

a_c_0=alpha_chem(I_rx);

RacRatio(vox)=nan2zero(x(vox+(4-1)*sz)./(x(vox+(4-1)*sz)+x(vox+(2-1)*sz)+x(vox+(7-1)*sz)));
RbarRatio(vox)=nan2zero(x(vox+(7-1)*sz)./(x(vox+(4-1)*sz)+x(vox+(2-1)*sz)+x(vox+(7-1)*sz)));
RhoRatio(vox)=nan2zero(x(vox+(3-1)*sz)./(x(vox+(3-1)*sz)+x(vox+(1-1)*sz)));
PaxRatio(vox)=nan2zero(x(vox+(6-1)*sz)./(x(vox+(6-1)*sz)+x(vox+(5-1)*sz)+x(vox+(8-1)*sz)));

gamma=0.3;


if sum([nnz(isnan(RhoRatio)), nnz(isnan(RacRatio)), nnz(isnan(PaxRatio))])~=0
    disp("woah")
end


K_is(vox)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(vox)).*(1+alpha_R*RacRatio(vox))+k_G*k_X*GIT*PIX);
K(vox)=alpha_R*RacRatio(vox).*K_is(vox).*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio(vox));%RbarRatio(I)/gamma;         %changed from paper
I_Ks(vox)=I_K*(1-K_is(vox).*(1+alpha_R*RacRatio(vox)));

reaction(vox+(1-1)*sz) = I_rho*(L_R^m./(L_R^m +(RacRatio(vox)+gamma*K(vox)).^m));            %From inactive rho to active rho changed from model
reaction(vox+(2-1)*sz) = (I_R+I_Ks(vox)).*(L_rho^m./(L_rho^m+RhoRatio(vox).^m));                %From inactive Rac to active Rac
reaction(vox+(5-1)*sz) = B_1*(K(vox).^m./(L_K^m+K(vox).^m));





alpha_chem(I_rx) = reaction(I_rx).*x(I_rx); %chemical reaction
alpha_rx=alpha_rx+sum(alpha_chem(I_rx)-a_c_0,1);
        
        a_total_new=sum(alpha_rx);
        
    end
    
    if neg
        error('Oh no! D: (negtive numbers): ')
    end
    
    
    dt_diff=dt_diff+tau;
    ind_diff=mpT0.*dt_diff>P_diff;
    
    if nnz(ind_diff)>0

        x=Alg322(x,dt_diff,D,h,jump',diffuse_mask,pT0,pi,cell_inds,A,ind_diff);
        
        dt_diff(ind_diff)=0;
        i=cell_inds(1:A);

if length(vox)>1
[tmp,tmp2]=meshgrid(ir0,vox);
    I_rx=tmp+tmp2;
else
    I_rx=vox+ir0;
end

a_c_0=alpha_chem(I_rx);

RacRatio(vox)=nan2zero(x(vox+(4-1)*sz)./(x(vox+(4-1)*sz)+x(vox+(2-1)*sz)+x(vox+(7-1)*sz)));
RbarRatio(vox)=nan2zero(x(vox+(7-1)*sz)./(x(vox+(4-1)*sz)+x(vox+(2-1)*sz)+x(vox+(7-1)*sz)));
RhoRatio(vox)=nan2zero(x(vox+(3-1)*sz)./(x(vox+(3-1)*sz)+x(vox+(1-1)*sz)));
PaxRatio(vox)=nan2zero(x(vox+(6-1)*sz)./(x(vox+(6-1)*sz)+x(vox+(5-1)*sz)+x(vox+(8-1)*sz)));

gamma=0.3;


if sum([nnz(isnan(RhoRatio)), nnz(isnan(RacRatio)), nnz(isnan(PaxRatio))])~=0
    disp("woah")
end


K_is(vox)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(vox)).*(1+alpha_R*RacRatio(vox))+k_G*k_X*GIT*PIX);
K(vox)=alpha_R*RacRatio(vox).*K_is(vox).*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio(vox));%RbarRatio(I)/gamma;         %changed from paper
I_Ks(vox)=I_K*(1-K_is(vox).*(1+alpha_R*RacRatio(vox)));

reaction(vox+(1-1)*sz) = I_rho*(L_R^m./(L_R^m +(RacRatio(vox)+gamma*K(vox)).^m));            %From inactive rho to active rho changed from model
reaction(vox+(2-1)*sz) = (I_R+I_Ks(vox)).*(L_rho^m./(L_rho^m+RhoRatio(vox).^m));                %From inactive Rac to active Rac
reaction(vox+(5-1)*sz) = B_1*(K(vox).^m./(L_K^m+K(vox).^m));





alpha_chem(I_rx) = reaction(I_rx).*x(I_rx); %chemical reaction
alpha_rx=alpha_rx+sum(alpha_chem(I_rx)-a_c_0,1);
        
        a_total_new=sum(alpha_rx);
        
        
    end
    
    
end


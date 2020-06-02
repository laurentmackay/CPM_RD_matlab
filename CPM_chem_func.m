function [x,diffusing_species_sum,D,h,alpha_rx,num_diffuse,...
    alpha_chem,time,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
    RacRatio,RbarRatio,I_Ks,reaction,ij_diffuse,jump,ir0,id0,cell_inds,...
    k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
    alpha,PAKtot,numDiff,numReac] =  CPM_chem_func(x,diffusing_species_sum,D,h,alpha_rx,num_diffuse,...
    alpha_chem,time,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
    RacRatio,RbarRatio,I_Ks,reaction,ij_diffuse,jump,ir0,id0,cell_inds,...
    k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
    alpha,PAKtot,nrx,A,numDiff,numReac)


N_instantaneous=50; % the number of steady reaction itterated at a moment in time
sz=size(x,1)*size(x,2);

%variables declared so that the c encoder is happy
i=1;
p=1;
i2=1;
rx=1;
xi0=x(id0+i(1));
i_update=[i i2];

for kk=1:nrx
    %      alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
    
    %total propensity to diffuse for all species
    alpha_diff=sum(diffusing_species_sum).*D/(h*h);
    
    
    %total propensity for rxn+diff
    a_total=sum(alpha_diff)+sum(alpha_rx);
    
    tau  = (1/a_total)*log(1/rand()); % time increment
    time=time+tau;
    reacted=false;
    diffused=false;
    neg=false;
    RN=rand();
    
    if RN*a_total <= sum(alpha_diff(:))
        
        diffused=true;
        
        %-------------diffusion-----------
        p=find((RN*a_total<=cumsum(alpha_diff)),1); %which protein diffuses
        p=p(1); %for encoder
        i0=(p-1)*sz; % the vectorized location of the appropriate proteins
        
        %adding all previous reactions to the propensity
        if p==1
            temp=0;
        else
            temp=sum(alpha_diff(1:p-1));
        end
        
        for drx=1:length(num_diffuse) %iterate over possible diffusion directions
            %check if it will diffuse in this direction
            if temp+(D(p)/h^2)*diffusing_species_sum(drx,p)>RN*a_total&&~reacted
                
                %find the point that diffuses
                ii=find((D(p)/h^2)*cumsum(x(i0+ij_diffuse(drx,1:num_diffuse(drx)')))>RN*a_total-temp,1);
                i(1)=ij_diffuse(drx,ii);
                %find the point they diffuses to
                i2=jump(i,drx);
                i_update=[i i2];
                
                %carry out the diffusion reaction
                x(i0+i) = x(i0+i)-1;
                x(i0+i2) = x(i0+i2)+1;
                %update the sum of diffusing species
                %in case you have diffused to an edge
                diffusing_species_sum(:,p)=diffusing_species_sum(:,p)+(diffuse_mask(:,i2)-diffuse_mask(:,i));
                
                reacted=true;
                neg=any(x(i+i0)<0);
                
            elseif ~reacted
                %keep adding
                temp=temp+(D(p)/h^2)*diffusing_species_sum(drx,p);
            else
                break;
            end
        end
        
        if reacted==false %making sure the code works
            error('Oh no! D: diffusion propensites did not sum correctly')
        end
        numDiff=numDiff+1;%ellie
    else% ---------- Reaction Time -------------
        reacted=false;
        temp=sum(alpha_diff(:));
        %determing where the reaction happens
        for rx=1:size(alpha_chem,3)
            i0=(rx-1)*sz;
            if temp+alpha_rx(rx)>=RN*a_total&&~reacted
                ii=find(cumsum(alpha_chem(i0+cell_inds(1:A)))>=RN*a_total-temp,1);
                i=cell_inds(ii(1));
                i_update=i;
                xi0=x(id0+i(1));
                reacted = true;
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
        numReac=numReac+1;%ellie
        %Inactive rho to active rho
        if rx==1
            x(i+(3-1)*sz) = x(i+(3-1)*sz)+1;
            x(i+(1-1)*sz) = x(i+(1-1)*sz)-1;
        end
        
        %Inactive Rac to active Rac
        if rx==2
            x(i+(2-1)*sz)=x(i+(2-1)*sz)-1;
            x(i+(4-1)*sz)=x(i+(4-1)*sz)+1;
        end
        
        %Active rho to inactive rho
        if rx==3
            x(i+(3-1)*sz) = x(i+(3-1)*sz)-1;
            x(i+(1-1)*sz) = x(i+(1-1)*sz)+1;
        end
        
        %Active Rac to inactive Rac
        if rx==4
            x(i+(2-1)*sz)=x(i+(2-1)*sz)+1;
            x(i+(4-1)*sz)=x(i+(4-1)*sz)-1;
        end
        
        %Unphosphorylated Pax to phosphorylated Pax
        if rx==5
            x(i+(5-1)*sz)=x(i+(5-1)*sz)-1;
            x(i+(6-1)*sz)=x(i+(6-1)*sz)+1;
        end
        
        %Phosphorylated Pax to unphosphorylated
        if rx==6
            x(i+(5-1)*sz)=x(i+(5-1)*sz)+1;
            x(i+(6-1)*sz)=x(i+(6-1)*sz)-1;
        end
        % I moved the following 2 lines into my dependence tree due to
        % complexing (should uncomment if replicating this algorthim)
        %         dxi=xi0-x(id0+i(1));
        %         diffusing_species_sum = diffusing_species_sum - (diffuse_mask(:,i)*dxi);
        %         neg=x(i+(rx-1)*sz)<0;
    end
    
    %-----------recalculating value that would have changed--------------

    if ~diffused %this is a hack to prevent errors, but it is technically wrong...complexation should always be run
        %deal with complexing and decomplexing
        
        if rx==2||rx==4%Rac
            K_R=(1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(i))*alpha*PAKtot*K_is(i);
            for j=1:N_instantaneous
                if x(i+(7-1)*sz)/(x(i+(7-1)*sz)+x(i+(4-1)*sz)*K_R)>rand()%decomplex
                    x(i+(4-1)*sz)=x(i+(4-1)*sz)+1;
                    x(i+(7-1)*sz)=x(i+(7-1)*sz)-1;
                elseif x(i+(4-1)*sz)>0
                    x(i+(4-1)*sz)=x(i+(4-1)*sz)-1;
                    x(i+(7-1)*sz)=x(i+(7-1)*sz)+1;
                end
            end
        elseif rx==5||rx==6%Pax
            K_P=k_G*k_X*k_C*GIT*PIX*K_is(i)*PAKtot*(1+alpha_R*RacRatio(i));
            for j=1:N_instantaneous
                if x(i+(8-1)*sz)/(x(i+(8-1)*sz)+x(i+(6-1)*sz)*K_P)>rand()%decomplex
                    x(i+(6-1)*sz)=x(i+(6-1)*sz)+1;
                    x(i+(8-1)*sz)=x(i+(8-1)*sz)-1;
                elseif x(i+(6-1)*sz)>0
                    x(i+(8-1)*sz)=x(i+(8-1)*sz)+1;
                    x(i+(6-1)*sz)=x(i+(6-1)*sz)-1;
                end
            end
        end
        dxi=xi0-x(id0+i(1));
        diffusing_species_sum = diffusing_species_sum - (diffuse_mask(:,i)*dxi);
        neg=x(i+(rx-1)*sz)<0;
    end
    

    for j=1:length(i_update)
        j=i_update(j);
        i_chem=ir0+j;
        
        a_c_0=alpha_chem(i_chem);
        
        %update Ratios
        RacRatio(j)=nan2zero(x(j+(4-1)*sz)./(x(j+(4-1)*sz)+x(j+(2-1)*sz)+x(j+(7-1)*sz)));
        RbarRatio(j)=nan2zero(x(j+(7-1)*sz)./(x(j+(4-1)*sz)+x(j+(2-1)*sz)+x(j+(7-1)*sz)));
        RhoRatio(j)=nan2zero(x(j+(3-1)*sz)./(x(j+(3-1)*sz)+x(j+(1-1)*sz)));
        PaxRatio(j)=nan2zero(x(j+(6-1)*sz)./(x(j+(6-1)*sz)+x(j+(5-1)*sz)+x(j+(8-1)*sz)));
        
        if sum([nnz(isnan(RhoRatio)), nnz(isnan(RacRatio)), nnz(isnan(PaxRatio))])~=0
            disp("woah")
        end
        
        
        
        %update other parameters
        K_is(j)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(j)).*(1+alpha_R*RacRatio(j))+k_G*k_X*GIT*PIX);
        K(j)=alpha_R*RacRatio(j).*K_is(j).*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio(j));%RbarRatio(j)/gamma;         %changed from paper
        I_Ks(j)=I_K*(1-K_is(j).*(1+alpha_R*RacRatio(j)));
        
        reaction(j+(1-1)*sz) = I_rho*(L_R^m./(L_R^m +(RacRatio(j)+RbarRatio(j)).^m));            %From inactive rho to active rho changed from model
        reaction(j+(2-1)*sz) = (I_R+I_Ks(j)).*(L_rho^m./(L_rho^m+RhoRatio(j).^m));                %From inactive Rac to active Rac
        reaction(j+(5-1)*sz) = B_1*(K(j).^m./(L_K^m+K(j).^m));
        
        
        alpha_chem(i_chem) = reaction(i_chem).*x(i_chem); %chemical reaction
        alpha_rx=alpha_rx+alpha_chem(i_chem)-a_c_0;
        
        
    end
    
    
    if neg
        error('Oh no! D: (negtive numbers)')
    end
    
end

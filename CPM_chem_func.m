function [x,diffusing_species_sum,D,h,alpha_rx,num_diffuse,...
    alpha_chem,time,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
    RacRatio,RbarRatio,I_Ks,reaction,ij_diffuse,jump,ir0,id0,...
    k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
    alpha,PAKtot] =  CPM_chem_func(x,diffusing_species_sum,D,h,alpha_rx,num_diffuse,...
    alpha_chem,time,diffuse_mask,PaxRatio,RhoRatio,K_is,K,...
    RacRatio,RbarRatio,I_Ks,reaction,ij_diffuse,jump,ir0,id0,cell_inds,...
    k_X,PIX,k_G,k_C,GIT,Paxtot,alpha_R,gamma,I_K,I_rho,I_R,L_R,m,L_rho,B_1,L_K,...
    alpha,PAKtot,nrx,A)


N_instantaneous=50; % the number of steady reaction itterated at a moment in time 
sz=size(x,1)*size(x,2);

%variables declared so that the c encoder is happy 
i=1;
p=1;
i2=1;
rx=1;
xi0=x(id0+i(1));
I=[i i2];
for kk=1:nrx
    %total propensity to diffuse for all species
    alpha_diff=sum(diffusing_species_sum).*D/(h*h);
    
    
    %total propensity for rxn+diff
    a_total=sum(alpha_diff)+sum(alpha_rx(:));
    
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
                %find the point it diffuses to
                i2=jump(i,drx);
                I=[i i2];
                
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
        
    else% ---------- Reaction Time -------------
        reacted=false;      
        temp=sum(alpha_diff(:));
        %determing where the reaction happens 
        for rx=1:size(alpha_chem,3)
            i0=(rx-1)*sz;
            if temp+alpha_rx(rx)>=RN*a_total&&~reacted
                ii=find(cumsum(alpha_chem(i0+cell_inds(1:A)))>=RN*a_total-temp,1);
                i(1)=cell_inds(ii);
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
    %the reaction trees are seprate as 2 cells change when diffusion
    %happens
    
    if diffused %diff reaction tree
        ai0=alpha_chem(ir0+i(1));
        ai20=alpha_chem(ir0+i2(1));
        if p==2||p==4
            RacRatio(I)=x(I+(4-1)*sz)./(x(I+(4-1)*sz)+x(I+(2-1)*sz)+x(I+(7-1)*sz));
            RbarRatio(I)=x(I+(7-1)*sz)./(x(I+(4-1)*sz)+x(I+(2-1)*sz)+x(I+(7-1)*sz));
            if any(isnan(RacRatio(I)))
                RacRatio(isnan(RacRatio))=0;
                RbarRatio(isnan(RbarRatio))=0;
            end
            K_is(I)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(I)).*(1+alpha_R*RacRatio(I))+k_G*k_X*GIT*PIX);
            K(I)=RbarRatio(I)/gamma;         %changed from paper
            I_Ks(I)=I_K*(1-K_is(I).*(1+alpha_R*RacRatio(I)));
            reaction(I+(1-1)*sz) = I_rho*(L_R^m./(L_R^m +(RacRatio(I)+RbarRatio(I)).^m));            %From inactive rho to active rho changed from model
            reaction(I+(2-1)*sz) = (I_R+I_Ks(I)).*(L_rho^m./(L_rho^m+RhoRatio(I).^m));                %From inactive Rac to active Rac
            reaction(I+(5-1)*sz) = B_1*(K(I).^m./(L_K^m+K(I).^m));
        end
        if p==1||p==3
            RhoRatio(I)=x(I+(3-1)*sz)./(x(I+(3-1)*sz)+x(I+(1-1)*sz));
            if any(isnan(RhoRatio(I)))
                RhoRatio(isnan(RhoRatio))=0;
            end
            reaction(I+(2-1)*sz) = (I_R+I_Ks(I)).*(L_rho^m./(L_rho^m+RhoRatio(I).^m));                %From inactive Rac to active Rac
        end
        if p==5||p==6
            PaxRatio(I)=x(I+(6-1)*sz)./(x(I+(6-1)*sz)+x(I+(5-1)*sz)+x(I+(8-1)*sz));
            if any(isnan(PaxRatio(I)))
                PaxRatio(isnan(PaxRatio))=0;
            end
            K_is(I)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(I)).*(1+alpha_R*RacRatio(I))+k_G*k_X*GIT*PIX);
            I_Ks(I)=I_K*(1-K_is(I).*(1+alpha_R*RacRatio(I)));
            reaction(I+(2-1)*sz) = (I_R+I_Ks(I)).*(L_rho^m./(L_rho^m+RhoRatio(I).^m));                %From inactive Rac to active Rac
        end
        alpha_chem(ir0+i(1)) = reaction(i(1)+ir0).*x(i(1)+ir0);
        alpha_chem(ir0+i2(1)) = reaction(i2(1)+ir0).*x(i2(1)+ir0);
        alpha_rx=alpha_rx+(alpha_chem(ir0+i(1))-ai0)+(alpha_chem(ir0+i2(1))-ai20);
    else
        ai0=alpha_chem(ir0+i(1));
        %-----chem reaction tree
        
        if rx==2||rx==4
            K_R=(1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(i))*alpha*PAKtot*K_is(i);
            for j=1:N_instantaneous
                if x(i+(7-1)*sz)/(x(i+(7-1)*sz)+x(i+(4-1)*sz)*K_R)>rand()
                    x(i+(4-1)*sz)=x(i+(4-1)*sz)+1;
                    x(i+(7-1)*sz)=x(i+(7-1)*sz)-1;
                elseif (x(i+(7-1)*sz)+x(i+(4-1)*sz))>0
                    x(i+(4-1)*sz)=x(i+(4-1)*sz)-1;
                    x(i+(7-1)*sz)=x(i+(7-1)*sz)+1;
                end
            end
            RacRatio(i)=x(i+(4-1)*sz)/(x(i+(4-1)*sz)+x(i+(2-1)*sz)+x(i+(7-1)*sz));
            RbarRatio(i)=x(i+(7-1)*sz)/(x(i+(4-1)*sz)+x(i+(2-1)*sz)+x(i+(7-1)*sz));
            if any(isnan(RacRatio(i)))
                RacRatio(i)=0;
                RbarRatio(i)=0;
            end
            K_is(i)=1/((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(i))*(1+alpha_R*RacRatio(i))+k_G*k_X*GIT*PIX);
            K(i)=RbarRatio(i)/gamma;         %changed from paper
            I_Ks(i)=I_K*(1-K_is(i)*(1+alpha_R*RacRatio(i)));
            reaction(i+(1-1)*sz) = I_rho*(L_R^m/(L_R^m +(RacRatio(i)+RbarRatio(i))^m));            %From inactive rho to active rho changed from model
            reaction(i+(2-1)*sz) = (I_R+I_Ks(i))*(L_rho^m/(L_rho^m+RhoRatio(i)^m));                %From inactive Rac to active Rac
            reaction(i+(5-1)*sz) = B_1*(K(i)^m/(L_K^m+K(i)^m));
        end
        if rx==1||rx==3
            RhoRatio(i)=x(i+(3-1)*sz)/(x(i+(3-1)*sz)+x(i+(1-1)*sz));
            if any(isnan(RhoRatio(i)))
                RhoRatio(i)=0;
            end
            reaction(i+(2-1)*sz) = (I_R+I_Ks(i))*(L_rho^m/(L_rho^m+RhoRatio(i)^m));                %From inactive Rac to active Rac
        end
        if rx==5||rx==6
            K_P=k_G*k_X*k_C*GIT*PIX*K_is(i)*PAKtot*(1+alpha_R*RacRatio(i));
            for j=1:N_instantaneous
                if x(i+(8-1)*sz)/(x(i+(8-1)*sz)+x(i+(6-1)*sz)*K_P)>rand()
                    x(i+(6-1)*sz)=x(i+(6-1)*sz)+1;
                    x(i+(8-1)*sz)=x(i+(8-1)*sz)-1;
                elseif (x(i+(8-1)*sz)+x(i+(6-1)*sz))>0
                    x(i+(8-1)*sz)=x(i+(8-1)*sz)+1;
                    x(i+(6-1)*sz)=x(i+(6-1)*sz)-1;
                end
            end
            PaxRatio(i)=x(i+(6-1)*sz)/(x(i+(6-1)*sz)+x(i+(5-1)*sz)+x(i+(8-1)*sz));
            if any(isnan(PaxRatio(i)))
                PaxRatio(i)=0;
            end
            K_is(i)=1/((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(i))*(1+alpha_R*RacRatio(i))+k_G*k_X*GIT*PIX);
            I_Ks(i)=I_K*(1-K_is(i)*(1+alpha_R*RacRatio(i)));
            reaction(i+(2-1)*sz) = (I_R+I_Ks(i))*(L_rho^m/(L_rho^m+RhoRatio(i)^m));                %From inactive Rac to active Rac
        end
        dxi=xi0-x(id0+i(1));
        diffusing_species_sum = diffusing_species_sum - (diffuse_mask(:,i)*dxi);
        neg=x(i+(rx-1)*sz)<0;
        alpha_chem(ir0+i(1)) = reaction(i(1)+ir0).*x(i(1)+ir0); %chemical reaction
        alpha_rx=alpha_rx+(alpha_chem(ir0+i(1))-ai0);
    end
    
    if neg
        error('Oh no! D: (negtive numbers)')
    end
    
end

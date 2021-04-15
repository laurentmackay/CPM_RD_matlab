%total propensity to diffuse for all species
alpha_diff=sum(diffusing_species_sum).*D/(h*h);

%total propensity for rxn+diff
a_total=sum(alpha_diff)+sum(alpha_rx(:));

tau  = (1/a_total)*log(1/rand()); % time increment
time=time+tau;
reacted=false;
diffused=false;
RN=rand();
% if min(min(min(x)))<0
%     error('Oh no! D: (negtive numbers)')
% end
if RN*a_total <= sum(alpha_diff(:))
    
    diffused=true;
    
    diff_reactions=diff_reactions+1;
    
    %-------------diffusion-----------
    p=find((RN*a_total<=cumsum(alpha_diff)),1); %which protein diffuses
    i0=(p-1)*sz;
    if p==1
        temp=0;
    else
        temp=sum(alpha_diff(1:p-1));
    end

    for drx=1:length(ij_diffuse) %iterate over possible diffusion directions
        %check if it will diffuse in this direction
        if temp+(D(p)/h^2)*diffusing_species_sum(drx,p)>RN*a_total&&~reacted
            
            %find the point that diffuses
            ii=find((D(p)/h^2)*cumsum(x(i0+ij_diffuse{drx}))>RN*a_total-temp,1);
            i=ij_diffuse{drx}(ii);
            
            %find the point it diffuses to
            i2=jump(i,drx);
			
            %carry out the diffusion reaction
            x(i0+[i i2]) = x(i0+[i i2]) + [-1, 1];
            
            %update the sum of diffusing species
            %in case you have diffused to an edge 
            diffusing_species_sum(:,p)=diffusing_species_sum(:,p)+(diffuse_mask(:,i2)-diffuse_mask(:,i));
            reacted=true;

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
    
else% ---------- Reaction Time! -------------
    reacted=false;
%     cell_mask=find(x(:,:,(1-1)));
    
    temp=sum(alpha_diff(:));
    for rx=1:size(alpha_chem,3)
        i0=(rx-1)*sz;
        if temp+alpha_rx(rx)>=RN*a_total&&~reacted
            ii=find(cumsum(alpha_chem(i0+cell_inds))>=RN*a_total-temp,1);%find the spatial location in the linear reference frame of the cell
            i=cell_inds(ii);%convert it back to the grid frame of reference
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
    xi0=x(id0+i);
    
    if rx==1
        x(i+(4-1-1)*sz) = x(i+(4-1-1)*sz)+1;
        x(i+(2-1-1)*sz) = x(i+(2-1-1)*sz)-1;

        Rho_reactions = Rho_reactions + 1;
    end
    
    %Inactive Rac to active Rac
    if rx==2
        K_R=(1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(i))*alpha*PAKtot*K_is(i);
        x(i+(3-1-1)*sz)=x(i+(3-1-1)*sz)-1;
        x(i+(4-1-1)*sz)=x(i+(4-1-1)*sz)+1;
        
        %complexing reaction

        
        for j=1:N_instantaneous
            if x(i+(8-1-1)*sz)/(x(i+(8-1-1)*sz)+x(i+(5-1-1)*sz)*K_R)>rand()
                x(i+(5-1-1)*sz)=x(i+(5-1-1)*sz)+1;
                x(i+(8-1-1)*sz)=x(i+(8-1-1)*sz)-1;
            elseif (x(i+(8-1-1)*sz)+x(i+(5-1-1)*sz))>0
                x(i+(5-1-1)*sz)=x(i+(5-1-1)*sz)-1;
                x(i+(8-1-1)*sz)=x(i+(8-1-1)*sz)+1;
            end
        end
%         Rac_tot=(x(i+(5-1-1)*sz)+x(i+(8-1-1)*sz));
%         x(i+(5-1-1)*sz)=poissrnd(Rac_tot/(1+K_R));
%         x(i+(8-1-1)*sz)=Rac_tot-x(i+(5-1-1)*sz);
%         disp(['1target: ' num2str() ', reality:' num2str(x(i+(5-1-1)*sz))])

        
        Rac_reactions = Rac_reactions + 1;
    end
    %Active rho to inactive rho
    if rx==3
        x(i+(4-1-1)*sz) = x(i+(4-1-1)*sz)-1;
        x(i+(2-1-1)*sz) = x(i+(2-1-1)*sz)+1;
        

        
        Rho_reactions = Rho_reactions + 1;
    end
    
    %Active Rac to inactive Rac
    if rx==4
        K_R=(1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(i))*alpha*PAKtot*K_is(i);
        x(i+(3-1-1)*sz)=x(i+(3-1-1)*sz)+1;
        x(i+(5-1-1)*sz)=x(i+(5-1-1)*sz)-1;
        
        %complexing reaction

        
        for j=1:N_instantaneous
            if x(i+(8-1-1)*sz)/(x(i+(8-1-1)*sz)+x(i+(5-1-1)*sz)*K_R)>rand()
                x(i+(5-1-1)*sz)=x(i+(5-1-1)*sz)+1;
                x(i+(8-1-1)*sz)=x(i+(8-1-1)*sz)-1;
            elseif (x(i+(8-1-1)*sz)+x(i+(5-1-1)*sz))>0
                x(i+(5-1-1)*sz)=x(i+(5-1-1)*sz)-1;
                x(i+(8-1-1)*sz)=x(i+(8-1-1)*sz)+1;
            end
        end
%         Rac_tot=(x(i+(5-1-1)*sz)+x(i+(8-1-1)*sz));
%         x(i+(5-1-1)*sz)=poissrnd(Rac_tot/(1+K_R));
%         x(i+(8-1-1)*sz)=Rac_tot-x(i+(5-1-1)*sz);
        
%         disp(['2target: ' num2str(Rac_tot/(1+K_R)) ', reality:' num2str(x(i+(5-1-1)*sz))])
        
        
        Rac_reactions = Rac_reactions + 1;
    end
    
    %Unphosphorylated Pax to phosphorylated Pax
    if rx==5
        K_P=k_G*k_X*k_C*GIT*PIX*K_is(i)*PAKtot*(1+alpha_R*RacRatio(i));
        x(i+(6-1-1)*sz)=x(i+(6-1-1)*sz)-1;
        x(i+(7-1-1)*sz)=x(i+(7-1-1)*sz)+1;
        
        %complexing reaction
        Pax_tot=x(i+(9-1-1)*sz)+x(i+(7-1-1)*sz);
        
        for j=1:N_instantaneous
            if x(i+(9-1-1)*sz)/(x(i+(9-1-1)*sz)+x(i+(7-1-1)*sz)*K_P)>rand()
                x(i+(7-1-1)*sz)=x(i+(7-1-1)*sz)+1;
                x(i+(9-1-1)*sz)=x(i+(9-1-1)*sz)-1;
            elseif (x(i+(9-1-1)*sz)+x(i+(7-1-1)*sz))>0
                x(i+(9-1-1)*sz)=x(i+(9-1-1)*sz)+1;
                x(i+(7-1-1)*sz)=x(i+(7-1-1)*sz)-1;
            end
        end
        %complexing reaction
%         Pax_tot=x(i+(9-1-1)*sz)+x(i+(7-1-1)*sz);
%         x(i+(7-1-1)*sz)=poissrnd(Pax_tot/(1+K_P));
%         x(i+(9-1-1)*sz)=Pax_tot-x(i+(7-1-1)*sz);
%        disp(['3target: ' num2str(Pax_tot/(1+K_P)) ', reality:' num2str(x(i+(7-1-1)*sz))])
        
        Pax_reactions = Pax_reactions + 1;
    end
    
    %Phosphorylated Pax to unphosphorylated
    
    if rx==6
        K_P=k_G*k_X*k_C*GIT*PIX*K_is(i)*PAKtot*(1+alpha_R*RacRatio(i));
        x(i+(6-1-1)*sz)=x(i+(6-1-1)*sz)+1;
        x(i+(7-1-1)*sz)=x(i+(7-1-1)*sz)-1;

        for k=1:N_instantaneous
            if x(i+(9-1-1)*sz)/(x(i+(9-1-1)*sz)+x(i+(7-1-1)*sz)*K_P)>rand()
                x(i+(7-1-1)*sz)=x(i+(7-1-1)*sz)+1;
                x(i+(9-1-1)*sz)=x(i+(9-1-1)*sz)-1;
            elseif (x(i+(9-1-1)*sz)+x(i+(7-1-1)*sz))>0
                x(i+(9-1-1)*sz)=x(i+(9-1-1)*sz)+1;
                x(i+(7-1-1)*sz)=x(i+(7-1-1)*sz)-1;
            end
        end
        
        %complexing reaction
%         Pax_tot=x(i+(9-1-1)*sz)+x(i+(7-1-1)*sz);
%         x(i+(7-1-1)*sz)=poissrnd(Pax_tot/(1+K_P));
%         x(i+(9-1-1)*sz)=Pax_tot-x(i+(7-1-1)*sz);
        
%         disp(['4target: ' num2str(Pax_tot/(1+K_P)) ', reality:' num2str(x(i+(7-1-1)*sz))])
%         
        Pax_reactions = Pax_reactions + 1;
    end
            
        dxi=xi0-x(id0+i);
        diffusing_species_sum = diffusing_species_sum - (diffuse_mask(:,i).*dxi);
      
end


ai0=alpha_chem(ir0+i);
if diffused
    ai20=alpha_chem(ir0+i2);
    rx=p;%shourd be other way arpund
end
if rx==2||rx==4
    RacRatio(i)=x(i+(5-1-1)*sz)./(x(i+(5-1-1)*sz)+x(i+(3-1-1)*sz)+x(i+(8-1-1)*sz));
    RbarRatio(i)=x(i+(8-1-1)*sz)./(x(i+(5-1-1)*sz)+x(i+(3-1-1)*sz)+x(i+(8-1-1)*sz));
    if any(isnan(RacRatio(i)))
        RacRatio(i)=0;
        RbarRatio(i)=0;
    end
    K_is(i)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(i)).*(1+alpha_R*RacRatio(i))+k_G*k_X*GIT*PIX);
    K(i)=RbarRatio(i)/gamma;         %changed from paper
    I_Ks(i)=I_K*(1-K_is(i).*(1+alpha_R*RacRatio(i)));
    reaction(i+(1-1)*sz) = I_rho*(L_R^m./(L_R^m +(RacRatio(i)+RbarRatio(i)).^m));            %From inactive rho to active rho changed from model
    reaction(i+(2-1)*sz) = (I_R+I_Ks(i)).*(L_rho^m./(L_rho^m+RhoRatio(i).^m));                %From inactive Rac to active Rac
    reaction(i+(5-1)*sz) = B_1*(K(i).^m./(L_K^m+K(i).^m));    
end
if rx==1||rx==3
    RhoRatio(i)=x(i+(3+1-1-1)*sz)./(x(i+(3+1-1-1)*sz)+x(i+sz));
    if any(isnan(RhoRatio(i)))
        RhoRatio(i)=0;
    end
    reaction(i+(2-1)*sz) = (I_R+I_Ks(i)).*(L_rho^m./(L_rho^m+RhoRatio(i).^m));                %From inactive Rac to active Rac
end
if rx==5||rx==6
    PaxRatio(i)=x(i+(7-1-1)*sz)./(x(i+(7-1-1)*sz)+x(i+(6-1-1)*sz)+x(i+(9-1-1)*sz));
    if any(isnan(PaxRatio(i)))
        PaxRatio(i)=0;
    end
    K_is(i)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(i)).*(1+alpha_R*RacRatio(i))+k_G*k_X*GIT*PIX);
    I_Ks(i)=I_K*(1-K_is(i).*(1+alpha_R*RacRatio(i)));
    reaction(i+(2-1)*sz) = (I_R+I_Ks(i)).*(L_rho^m./(L_rho^m+RhoRatio(i).^m));                %From inactive Rac to active Rac
end


neg=false;
if ~diffused
    alpha_chem(ir0+i) = reaction(i+ir0).*x(i+(ir0+sz)); %chemical reaction
    alpha_rx=alpha_rx+(alpha_chem(ir0+i)-ai0);
	
    neg=x(i+(rx-1)*sz)<0;
else
    alpha_chem(ir0+i) = reaction(i+ir0).*x((ir0+sz)+i); %chemical reaction
    alpha_chem(ir0+i2) = reaction(i2+ir0).*x((ir0+sz)+i2); %chemical reaction
    
    alpha_rx=alpha_rx+(alpha_chem(ir0+i)-ai0)+(alpha_chem(ir0+i2)-ai20);
    neg=x(i+i0)<0||x(i2+i0)<0;
    
end




if neg
% if sum(x(:)<0)>0
    error('Oh no! D: (negtive numbers)')
end

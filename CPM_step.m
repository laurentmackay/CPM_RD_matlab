adj_empty = ~cell_mask(up) | ~cell_mask(down) | ~cell_mask(left) | ~cell_mask(right); %cells at the boundry 
adj_full = cell_mask(up) | cell_mask(down) | cell_mask(left) | cell_mask(right); %cells at the boundry 

%finding boundry points 
bndry_cell = cell_mask & adj_empty;
bndry_empty = ~cell_mask & adj_full;
bndry = find( bndry_cell | bndry_empty );

if any(cell_maskp~=cell_mask)
    error('not reseting')
end
ij = bndry(randi(length(bndry)));

r=randi(4);
cell_maskp(ij) = cell_mask(jump(sub2ind([sz,4],ij,r)));% make a new trial configuration

Per=perim(cell_maskp); % perimter
A=nnz(cell_maskp); % area
HA=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; % the hamiltonian after the possible change
dH=HA-H0;
Ncell_mask=squeeze(sum(sum(x))); %for a sanity check 

if getfield(bwconncomp(cell_maskp,4),'NumObjects')==1 %makes sure the cell stays connected 
    grow= cell_maskp(ij) & ~cell_mask(ij);
    shrink= ~cell_maskp(ij) & cell_mask(ij);
    if grow
        dH=dH+B_rho*(RhoRatio(jump(sub2ind([sz,4],ij,r)))-rho_eq)-B_R*(RacRatio(jump(sub2ind([sz,4],ij,r)))-R_eq);
        if rand<exp(-(dH+Hb)/T) %step raises energy boltzman prob forward
            cell_mask=cell_maskp; %changing cell shape
            
            for j=0:(N_species-1) %splitting the molecules with the new lattice
                x(ij+j*sz)=floor(x(jump(sub2ind([sz,4],ij,r))+j*sz)/2);
                x(jump(sub2ind([sz,4],ij,r))+j*sz)=ceil(x(jump(sub2ind([sz,4],ij,r))+j*sz)/2);
            end
            
            
            
            i=[ij jump(sub2ind([sz,4],ij,r))]; %places where molecule number has changed
            
            H0=HA; %changing the hamiltonn to the new one
            
            %similar idea now to the CPM_chem_func dependence tree
            
            RacRatio(i)=x(i+(4-1)*sz)./(x(i+(4-1)*sz)+x(i+(2-1)*sz)+x(i+(7-1)*sz));
            RbarRatio(i)=x(i+(7-1)*sz)./(x(i+(4-1)*sz)+x(i+(2-1)*sz)+x(i+(7-1)*sz));
            RhoRatio(i)=x(i+(3-1)*sz)./(x(i+(3-1)*sz)+x(i+(1-1)*sz));
            PaxRatio(i)=x(i+(6-1)*sz)./(x(i+(6-1)*sz)+x(i+(5-1)*sz)+x(i+(8-1)*sz));
            if any(isnan(PaxRatio(i)))
                PaxRatio(isnan(PaxRatio))=0;
            end
            if any(isnan(RhoRatio(i)))
                RhoRatio(isnan(RhoRatio))=0;
            end
            if any(isnan(RacRatio(i)))
                RacRatio(isnan(RacRatio))=0;
                RbarRatio(isnan(RbarRatio))=0;
            end
            
            %----reactions that vary lattice to lattice
            K_is(i)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(i)).*(1+alpha_R*RacRatio(i))+k_G*k_X*GIT*PIX);
            K(i)=RbarRatio(i)/gamma;         %changed from paper
            I_Ks(i)=I_K*(1-K_is(i).*(1+alpha_R*RacRatio(i)));
            reaction(i+(1-1)*sz) = I_rho*(L_R^m./(L_R^m +(RacRatio(i)+RbarRatio(i)).^m));            %From inactive rho to active rho changed from model
            reaction(i+(2-1)*sz) = (I_R+I_Ks(i)).*(L_rho^m./(L_rho^m+RhoRatio(i).^m));                %From inactive Rac to active Rac
            reaction(i+(5-1)*sz) = B_1*(K(i).^m./(L_K^m+K(i).^m));
            
            
            for ip=i
                ai=alpha_chem(ir0+ip(1));
                alpha_chem(ir0+ip) = reaction(ip+ir0).*x(ip+ir0);
                alpha_rx=alpha_rx+(alpha_chem(ir0+ip)-ai);
            end
            
        else
            cell_maskp=cell_mask;
        end
    elseif shrink
        dH=dH-B_rho*(RhoRatio(ij)-rho_eq)+B_R*(RacRatio(ij)-R_eq);
        
        if rand<exp(-(dH+Hb)/T)%step raises energy
            cell_mask=cell_maskp;  %changing cell shape 
            
            neighbors=[jump(ij,1) jump(ij,2) jump(ij,3) jump(ij,4)]; %finding places the molecules will go to
            neighbors=neighbors(find(cell_maskp(neighbors)));
            
            for j=0:(N_species-1) %dumping out molecules from the retracting site
                x(neighbors+j*sz)=x(neighbors +j*sz)+diff(round(linspace(0,x(ij+j*sz),length(neighbors)+1)));
                x(ij+j*sz)=0;
            end
            
            i=[ij neighbors]; %indices where molecule number changed 
            
            
            %Same as grow
            H0=HA;

            
            RacRatio(neighbors)=x(neighbors+(4-1)*sz)./(x(neighbors+(4-1)*sz)+x(neighbors+(2-1)*sz)+x(neighbors+(7-1)*sz));
            RbarRatio(neighbors)=x(neighbors+(7-1)*sz)./(x(neighbors+(4-1)*sz)+x(neighbors+(2-1)*sz)+x(neighbors+(7-1)*sz));
            RhoRatio(neighbors)=x(neighbors+(3-1)*sz)./(x(neighbors+(3-1)*sz)+x(neighbors+(1-1)*sz));
            PaxRatio(neighbors)=x(neighbors+(6-1)*sz)./(x(neighbors+(6-1)*sz)+x(neighbors+(5-1)*sz)+x(neighbors+(8-1)*sz));
            if any(isnan(PaxRatio(i)))
                PaxRatio(isnan(PaxRatio))=0;
            end
            if any(isnan(RhoRatio(i)))
                RhoRatio(isnan(RhoRatio))=0;
            end
            if any(isnan(RacRatio(i)))
                RacRatio(isnan(RacRatio))=0;
                RbarRatio(isnan(RbarRatio))=0;
            end
            RacRatio(ij)=0;
            RbarRatio(ij)=0;
            RhoRatio(ij)=0;
            PaxRatio(ij)=0;
            
            %----reactions that vary lattice ot lattice
            K_is(i)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(i)).*(1+alpha_R*RacRatio(i))+k_G*k_X*GIT*PIX);
            K(i)=RbarRatio(i)/gamma;         %changed from paper
            I_Ks(i)=I_K*(1-K_is(i).*(1+alpha_R*RacRatio(i)));
            reaction(i+(1-1)*sz) = I_rho*(L_R^m./(L_R^m +(RacRatio(i)+RbarRatio(i)).^m));            %From inactive rho to active rho changed from model
            reaction(i+(2-1)*sz) = (I_R+I_Ks(i)).*(L_rho^m./(L_rho^m+RhoRatio(i).^m));                %From inactive Rac to active Rac
            reaction(i+(5-1)*sz) = B_1*(K(i).^m./(L_K^m+K(i).^m));
            
            
            for ip=i
                ai=alpha_chem(ir0+ip(1));
                alpha_chem(ir0+ip) = reaction(ip+ir0).*x(ip+ir0);
                alpha_rx=alpha_rx+(alpha_chem(ir0+ip)-ai);
            end
            
        else
            %reaction doesn't happen
            cell_maskp=cell_mask;
        end
    else
        cell_maskp=cell_mask;
    end
else
    cell_maskp=cell_mask;
end

Ncell_maskp=squeeze(sum(sum(x)));
%sanity checks 
if any(Ncell_mask~=Ncell_maskp)
    errror('molecule loss')
end

if min(cell_mask(:))<0
    error('Oh no! D: (negtive numbers)')
end


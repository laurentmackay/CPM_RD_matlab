adj_empty = ~cell_mask(up) | ~cell_mask(down) | ~cell_mask(left) | ~cell_mask(right); %cells at the boundry 
adj_full = cell_mask(up) | cell_mask(down) | cell_mask(left) | cell_mask(right); %cells at the boundry 

%finding boundry points 
bndry_cell = cell_mask & adj_empty;
bndry_empty = ~cell_mask & adj_full;
bndry = find( bndry_cell | bndry_empty );

if any(cell_maskp~=cell_mask)
    error('not reseting')
end
vox_trial = bndry(randi(length(bndry)));

r=randi(size(jump,2));
vox_ref=jump(sub2ind([sz,4],vox_trial,r));
cell_maskp(vox_trial) = cell_mask(vox_ref);% make a new trial configuration

Per=perim(cell_maskp); % perimeter
A=nnz(cell_maskp); % area
HA=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; % the hamiltonian after the possible change
dH=HA-H0;
Ncell_mask=squeeze(sum(sum(x))); %for a sanity check 

%set the equilibrium point as the overall Ratio such that in the selected boundary point ij:
%higher RacRatio contributes to protrusion
%higher RhoRatio contributes to retraction
%vice versa
rho_eq=mean(RhoRatio(find(cell_mask)));
R_eq=mean(RacRatio(find(cell_mask)));

no_holes = getfield(bwconncomp(cell_maskp,4),'NumObjects')==1 && getfield(bwconncomp(~cell_maskp,4),'NumObjects')==1 ;%makes sure the cell stays connected and no hole   
reacted = 0;
if  no_holes
    %check if growing or shrinking
    grow= cell_maskp(vox_trial) & ~cell_mask(vox_trial);
    shrink= ~cell_maskp(vox_trial) & cell_mask(vox_trial);
    
    
    if grow 
        f=1;
        dH=dH+B_rho*(RhoRatio(vox_ref)-rho_eq)-B_R*(RacRatio(vox_ref)-R_eq);
    elseif shrink
        f=-1;
        dH=dH-B_rho*(RhoRatio(vox_trial)-rho_eq)+B_R*(RacRatio(vox_trial)-R_eq);
    end
        
        
        if (grow || shrink) && rand<exp(-(dH+Hb)/T) 
            reacted=1;
            cell_mask=cell_maskp; %changing cell shape
            
            for j=0:(N_species-1) %splitting the molecules with the new lattice
%                 P1=D(j+1)*0.5*cpmstep/(h^2)

            counts=x(cell_inds(1:A-1)+j*sz);
            if grow
                tmp=x(vox_ref+j*sz);
            else
                tmp=x(vox_trial+j*sz);
                counts(cell_inds(1:A-1)==vox_trial)=0;
            end

                
                Ntot=sum(counts);
                p=counts/sum(counts);
                
%                 if shrink
%                     p=1-p;
%                 end
                
                for k=1:tmp
                    l=Alg2(p,rand(),0);
                    x(cell_inds(l)+j*sz)=x(cell_inds(l)+j*sz)-f;
%                     counts(l)=counts(l)-1;
                end
                if grow
                    x(vox_trial+j*sz)=tmp;
                else
                    x(vox_trial+j*sz)=0;
                end
            end

            
            I=[vox_trial vox_ref]; %places where molecule number has changed
            H0=HA; %changing the hamiltonn to the new one
            
            %recalculate parameters
            Per=perim(cell_mask); 
            A=nnz(cell_mask);
            cell_inds(1:A)=find(cell_mask);

            vox=cell_inds(1:A);
            update_alpha_chem
            
            alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
            if grow 
                disp('grow');
                grow_count=grow_count+1;
            else
                disp('shrink');
                shrink_count=shrink_count+1;
            end
             
        end
end
        
    if ~reacted 
        %no move
        cell_maskp=cell_mask;
        Per=perim(cell_mask); % perimter
        A=nnz(cell_mask); % area
        cell_inds(1:A)=find(cell_mask);
    end



Ncell_maskp=squeeze(sum(sum(x)));
%sanity checks 
if any(Ncell_mask~=Ncell_maskp)
    error('molecule loss')
end

if min(cell_mask(:))<0
    error('Oh no! D: (negtive numbers)')
end
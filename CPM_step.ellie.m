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

Per=perim(cell_maskp); % perimeter
A=nnz(cell_maskp); % area
HA=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; % the hamiltonian after the possible change
dH=HA-H0;
Ncell_mask=squeeze(sum(sum(x))); %for a sanity check 

%set the equilibrium point as the overall Ratio such that in the selected boundary point ij:
%higher RacRatio contributes to protrusion
%higher RhoRatio contributes to retraction
%vice versa
rho_eq=sum(sum(x(:,:,4))/(sum(sum(sum(x(:,:,[3 4]))))));
R_eq=sum(sum(x(:,:,2))/(sum(sum(sum(x(:,:,[1 2]))))));

if getfield(bwconncomp(cell_maskp,4),'NumObjects')==1 && getfield(bwconncomp(~cell_maskp,4),'NumObjects')==1 %makes sure the cell stays connected and no hole    
    grow= cell_maskp(ij) & ~cell_mask(ij);
    shrink= ~cell_maskp(ij) & cell_mask(ij);
    if grow
        dH=dH+B_rho*(RhoRatio(jump(sub2ind([sz,4],ij,r)))-rho_eq)-B_R*(RacRatio(jump(sub2ind([sz,4],ij,r)))-R_eq);
        if rand<exp(-(dH+Hb)/T) 
            cell_mask=cell_maskp; %changing cell shape
            
            for j=0:(N_species-1) %splitting the molecules with the new lattice
                x(ij+j*sz)=floor(x(jump(sub2ind([sz,4],ij,r))+j*sz)/2);
                x(jump(sub2ind([sz,4],ij,r))+j*sz)=ceil(x(jump(sub2ind([sz,4],ij,r))+j*sz)/2);
            end
            %{
            %Before we found a way to do bulk diffusion, this somewhat works better
            for j=0:(N_species-1) %copy the ratios to the new square
                x(ij+j*sz)=x(jump(sub2ind([sz,4],ij,r))+j*sz);
            end
            %}
            
            I=[ij jump(sub2ind([sz,4],ij,r))]; %places where molecule number has changed
            H0=HA; %changing the hamiltonn to the new one
            
            %recalculate parameters
            Per=perim(cell_mask); 
            A=nnz(cell_mask);
            cell_inds(1:A)=find(cell_mask);
            
            if grow
                vox=cell_inds(1:A);
            else
                vox=[cell_inds(1:A); vox_trial];
            end

            update_alpha_chem
            
            alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
            
            
        else
            %no grow
            cell_maskp=cell_mask;
            Per=perim(cell_maskp); % perimeter
            A=nnz(cell_maskp); % area
        end
    elseif shrink
        dH=dH-B_rho*(RhoRatio(ij)-rho_eq)+B_R*(RacRatio(ij)-R_eq);
        
        if rand<exp(-(dH+Hb)/T)
            cell_mask=cell_maskp;  %changing cell shape 
            
            neighbors=[jump(ij,1) jump(ij,2) jump(ij,3) jump(ij,4)]; %finding places the molecules will go to
            neighbors=neighbors(find(cell_maskp(neighbors)));
            
            for j=0:(N_species-1) %dumping out molecules from the retracting site
                x(neighbors+j*sz)=x(neighbors +j*sz)+diff(round(linspace(0,x(ij+j*sz),length(neighbors)+1)));
                x(ij+j*sz)=0;
            end
            %{
            %Before we found a way to do bulk diffusion, this somewhat works better
            for j=0:(N_species-1) %completely remove the molecules
                x(ij+j*sz)=0;
            end
            %}
            
            I=[ij neighbors]; %indices where molecule number changed
            %Same as grow
            H0=HA;

            Per=perim(cell_maskp); % perimeter
            A=nnz(cell_mask);
            cell_inds(1:A)=find(cell_mask);
            
            if grow
                vox=cell_inds(1:A);
            else
                vox=[cell_inds(1:A); vox_trial];
            end

            update_alpha_chem
            
            alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
            
        else
            %no shrink
            cell_maskp=cell_mask;
            Per=perim(cell_maskp); % perimeter
            A=nnz(cell_maskp); % area
        end
    else
        %no move
        cell_maskp=cell_mask;
        Per=perim(cell_maskp); % perimter
        A=nnz(cell_maskp); % area
    end
else
    %doesn't protrude or retract due to connection issue
    cell_maskp=cell_mask;
    Per=perim(cell_maskp); % perimter
    A=nnz(cell_maskp); % area
end

Ncell_maskp=squeeze(sum(sum(x)));
%sanity checks 
if any(Ncell_mask~=Ncell_maskp)
    errror('molecule loss')
end

if min(cell_mask(:))<0
    error('Oh no! D: (negtive numbers)')
end
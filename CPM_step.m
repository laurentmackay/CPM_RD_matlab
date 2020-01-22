adj_empty = ~cell_mask(up) | ~cell_mask(down) | ~cell_mask(left) | ~cell_mask(right);
adj_full = cell_mask(up) | cell_mask(down) | cell_mask(left) | cell_mask(right);

bndry_cell = cell_mask & adj_empty;
bndry_empty = ~cell_mask & adj_full;
bndry = find( bndry_cell | bndry_empty );

if any(cell_maskp~=cell_mask)
    error('not reseting')
end
ij = bndry(randi(length(bndry),[1,N_update]))';


r=randi(4,[1,N_update]);
cell_maskp(ij) = cell_mask(jump(sub2ind([sz,4],ij,r)));% make a new trial configuration

Per=perim(cell_maskp(:,:,1)); % perimter

A=nnz(cell_maskp(:,:,1)); % area
HA=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; % the hamiltonian after the possible change
dH=HA-H0;

if getfield(bwconncomp(cell_maskp,4),'NumObjects')==1
    grow= cell_maskp(ij) & ~cell_mask(ij);
    shrink= ~cell_maskp(ij) & cell_mask(ij);
    rg=r(find(grow)); % for finding the neighbors of where the cell grew
    grow=ij(find(grow)); % finding the indecies where the cell growed
    shrink=ij(find(shrink)); % finding the indecies where the cell shrunk
    
    dH=dH+...
B_rho*sum(RhoRatio(jump(sub2ind([sz,4],grow,rg)))-rho_eq)-...
B_R*sum(RacRatio(jump(sub2ind([sz,4],grow,rg)))-R_eq);
    for j=1:8
        cell_maskp(grow+j*sz)=floor(cell_mask(jump(sub2ind([sz,4],grow,rg))+j*sz)/2);
        cell_maskp(jump(sub2ind([sz,4],grow,rg))+j*sz)=ceil(cell_mask(jump(sub2ind([sz,4],grow,rg))+j*sz)/2);
    end
    
    
    dH=dH-B_rho*sum(RhoRatio(shrink)-rho_eq)+B_R*sum(RacRatio(shrink)-R_eq);
    for shrinkp=shrink
        neighbors=[up(shrinkp) down(shrinkp) right(shrinkp) left(shrinkp)];
        neighbors=neighbors(find(cell_maskp(neighbors)));
        for j=1:8
            cell_maskp(neighbors+j*sz)=cell_maskp(neighbors +j*sz)+diff(round(linspace(0,cell_maskp(shrinkp+j*sz),length(neighbors)+1)));
            cell_maskp(shrinkp+j*sz)=0;
        end
    end
    
    
    % Ncell_mask=sum(cell_mask,[1 2]);
    Ncell_mask=squeeze(sum(sum(cell_mask)));
    % Ncell_maskp=sum(cell_maskp,[1 2]);
    Ncell_maskp=squeeze(sum(sum(cell_maskp)));
    if any(Ncell_mask(2:9)~=Ncell_maskp(2:9))
        dH=inf;
        errror('molecule loss')
    end
    
    if rand<ecell_maskp(-(dH+Hb)/T)%step raises energy
        cell_mask=cell_maskp;
        H0=HA;
        %work out a more specfic way latter !!!!!!!!!!!!!!!!!!!
        RhoRatio=cell_mask(:,:,(4-1))./(cell_mask(:,:,(4-1))+cell_mask(:,:,(2-1)));
        RhoRatio(isnan(RhoRatio))=0;
        RacRatio=cell_mask(:,:,(5-1))./(cell_mask(:,:,(5-1))+cell_mask(:,:,(3-1))+cell_mask(:,:,(8-1)));
        RacRatio(isnan(RacRatio))=0;
        Pacell_maskRatio=cell_mask(:,:,(7-1))./(cell_mask(:,:,(7-1))+cell_mask(:,:,(6-1))+cell_mask(:,:,(9-1)));
        Pacell_maskRatio(isnan(Pacell_maskRatio))=0;
        RbarRatio=cell_mask(:,:,(8-1))./(cell_mask(:,:,(5-1))+cell_mask(:,:,(3-1))+cell_mask(:,:,(8-1)));  % this is gamma*K
        RbarRatio(isnan(RbarRatio))=0;
        
        %----reactions that vary lattice ot lattice
        K_is=1./((1+k_cell_mask*PIcell_mask+k_G*k_cell_mask*k_C*GIT*PIcell_mask*Pacell_masktot*Pacell_maskRatio).*(1+alpha_R*RacRatio)+k_G*k_cell_mask*GIT*PIcell_mask);
        K=RbarRatio/gamma;
        I_Ks=I_K*(1-K_is.*(1+alpha_R*RacRatio));
        reaction(:,:,1) = I_rho*(L_R^m./(L_R^m +(RacRatio+RbarRatio).^m));            %From inactive rho to active rho changed from model
        reaction(:,:,2) = (I_R+I_Ks).*(L_rho^m./(L_rho^m+RhoRatio.^m));                %From inactive Rac to active Rac
        reaction(:,:,5) = B_1*(K.^m./(L_K^m+K.^m));
        alpha_chem=reaction(:,:,1:6).*cell_mask(:,:,2:7);
    else
        cell_maskp=cell_mask;
    end
    
else
    cell_maskp=cell_mask;
end
if min(cell_mask(:))<0
    error('Oh no! D: (negtive numbers)')
end

if time>=timecheck+picstep
%     pic
    z=z+1;
    center(z,:)=com(cell_mask(:,:,(1-1)));
    timecheck=timecheck+picstep;
    time;
    reactions;
    toc
end
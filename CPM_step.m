adj_empty = ~cell_mask(up) | ~cell_mask(down) | ~cell_mask(left) | ~cell_mask(right); %cells at the boundry
adj_full = cell_mask(up) | cell_mask(down) | cell_mask(left) | cell_mask(right); %cells at the boundry

%finding boundry points
bndry_cell = cell_mask & adj_empty;
bndry_empty = ~cell_mask & adj_full;
bndry = find( bndry_cell | bndry_empty );

bndry_up=cell_mask  & ~cell_mask(up);
bndry_down=cell_mask  & ~cell_mask(down);
bndry_l=cell_mask  & ~cell_mask(left);
bndry_r=cell_mask  & ~cell_mask(right);

bndry_ud= bndry_up | bndry_down;
bndry_lr= bndry_l | bndry_r;

bndrys=[bndry_up(:) bndry_down(:) bndry_l(:) bndry_r(:)];

% bndry_lr=cell_mask & (~cell_mask(left) | ~cell_mask(right));

if any(cell_maskp~=cell_mask)
    error('not reseting')
end

%         for j=0:(N_species-1) %sanity checks
%             %                 P1=D(j+1)*0.5*cpmstep/(h^2)
%             u=x(:,:,j+1);
% 
%             if any(u(:)<0)
%                 disp('we were given something negative :(')
%             end
%             if any(u(~cell_mask)~=0)
%                 error("we got a wild ass")
%             end
%             disp(['species ' num2str(j) 'cjecks out!'])
%         end

%set the equilibrium point as the overall Ratio such that in the selected boundary point ij:
%higher RacRatio contributes to protrusion
%higher RhoRatio contributes to retraction
%vice versa
rho_eq=mean(RhoRatio(find(cell_mask)));
R_eq=mean(RacRatio(find(cell_mask)));
Ncell_mask=squeeze(sum(sum(x))); %for a sanity check
A0=A;

no_holes=false;

while ~no_holes
    vox_trial = bndry(randi(length(bndry)));
    
    r=randi(size(jump,2));
    vox_ref=jump(sub2ind([sz,4],vox_trial,r));
    cell_maskp(vox_trial) = cell_mask(vox_ref);% make a new trial configuration
    
    Per=perim(cell_maskp); % perimeter
    A=nnz(cell_maskp); % area
    HA=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; % the hamiltonian after the possible change
    dH=HA-H0;
    no_holes = getfield(bwconncomp(cell_maskp,4),'NumObjects')==1 && getfield(bwconncomp(~cell_maskp,4),'NumObjects')==1 ;%makes sure the cell stays connected and no hole
    if ~no_holes
        cell_maskp(vox_trial)=cell_mask(vox_trial);%revert change
    end
end


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
%         if grow
            cm0=cell_mask;
%         else
%             cm0=cell_mask;
%         end
        cell_mask=cell_maskp; %changing cell shape
        
        if shrink
            bndry_up=cell_mask  & ~cell_mask(up);
            bndry_down=cell_mask  & ~cell_mask(down);
            bndry_l=cell_mask  & ~cell_mask(left);
            bndry_r=cell_mask  & ~cell_mask(right);

            bndry_ud= bndry_up | bndry_down;
            bndry_lr= bndry_l | bndry_r;

        end
        
            %recalculate parameters
        Per=perim(cell_mask);
        A=nnz(cell_mask);

        
        
        if grow
            inds=cell_inds(1:A-1);
            cell_inds(1:A)=find(cell_mask);
        else
            cell_inds(1:A)=find(cell_mask);
            inds=cell_inds(1:A);
        end
        
        
        
            if grow
                dist=sqrt(((i0(vox_ref)-i0(inds)).^2)+(j0(vox_ref)-j0(inds)).^2)*h;
            else
                dist=sqrt(((i0(vox_ref)-i0(inds)).^2)+(j0(vox_ref)-j0(inds)).^2)*h;
            end

            min_dist=3*h;
            transport_mask=D/min(D)*min_dist>dist;
%             if ~grow
            transport_mask(find(vox_trial==inds),:)=false;
%             end
            
        x0=x;
            
            
            for i=1:length(D)
 
                inds2=inds(transport_mask(:,i))+(i-1)*sz;
                i_trial=vox_trial+(i-1)*sz;
                if grow
                    us=x(vox_ref+(i-1)*sz);
                    Ts=sum(x(inds2));  
                    f=Ts/(Ts+us);
                    x(i_trial)=us;
                    inds2=[inds2; i_trial];
                else
                    
                    f=1+(x(i_trial)/sum(x(inds2)));
                    x(i_trial)=0;
                end
                
                x(inds2)=floor(x(inds2)*f)+[0; diff(floor(cumsum(rem(x(inds2)*f,1))+1e-5))]; %the 1e-5 is a fudge-factor to prevent underflow erros, they are typically of the order 1e-10 so the 1e-5 dominates
%                 x(inds2)=x(inds2)*f;
%                 if sum(x(inds+(i-1)*sz))-sum0~=0
%                     disp('check this out boss')
%                 end
            end
        
        
        I=[vox_trial vox_ref]; %places where molecule number has changed
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
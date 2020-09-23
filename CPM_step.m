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
        if grow
            cm0=cell_mask;
        else
            cm0=cell_maskp;
        end
        cell_mask=cell_maskp; %changing cell shape
        
        if shrink
            bndry_up=cell_mask  & ~cell_mask(up);
            bndry_down=cell_mask  & ~cell_mask(down);
            bndry_l=cell_mask  & ~cell_mask(left);
            bndry_r=cell_mask  & ~cell_mask(right);

            bndry_ud= bndry_up | bndry_down;
            bndry_lr= bndry_l | bndry_r;

        end
        if grow
            inds=cell_inds(1:A-1);
        else
            inds=cell_inds(1:A);
        end
        
        
        dist=sqrt(((i0(vox_trial)-i0(inds)).^2)+(j0(vox_trial)-j0(inds)).^2)*h;
            di=i0(vox_trial)-i0;
            dj=j0(vox_trial)-j0;
            dist2=sqrt((di.^2)+(dj.^2))*h;
            di=di./dist2;
            dj=dj./dist2;
            di(isnan(di))=1;
            dj(isnan(dj))=1;
            lambda=1;
            vmag=-lambda./dist2;
            if shrink
                vmag=-vmag;
            end
%             dvmag=exp(-dist2/lambda)/lambda;
            vi=vmag.*di;
            vj=vmag.*dj;
            vi(vox_trial)=0;
            vj(vox_trial)=0;
            
            ju=zeros(shape);
            jd=zeros(shape);
            jr=zeros(shape);
            jl=zeros(shape);
        for j=0:(N_species-1) %splitting the molecules with the new lattice
            %                 P1=D(j+1)*0.5*cpmstep/(h^2)
            u=x(:,:,j+1);
            if any(u(:)<0)
                disp('we were given something negative :(')
            end
            if grow
                tmp=u(vox_ref);
                
                %                 x(vox_trial+j*sz)=tmp;
            else
                tmp=u(vox_trial);
            end
            if grow
                samps=randi(A0,tmp,1);
                counts= histcounts(samps,(0:A0)+0.5)';
            else
                samps=randi(A0,tmp,1);
                ind_trial=find(vox_trial==cell_inds(1:A0),1);
                samps(samps==ind_trial)=A0;
                counts= histcounts(samps,(0:A0)+0.5)';

            end
            
            
            
            
            %                 pdiff=exp(-(dist2.^2/(4*corr_len^2*h^2)));
%             pdiff=(tanh((dist2-2*corr_len)/2)+1)/2;

            
            

            

            jud = @(dir,ex)  (u(dir(~ex))-u(~ex)).*(vj(dir(~ex))+vj(~ex))/2;
            jlr = @(dir,ex)  -(u(dir(~ex))-u(~ex)).*(vi(dir(~ex))+vi(~ex))/2;
%             jud = @(dir,ex)  (u(dir(~ex))-u(~ex)).*(vi(dir(~ex))+vi(~ex))/2;
%             jlr = @(dir,ex)  (u(dir(~ex))-u(~ex)).*(vj(dir(~ex))+vj(~ex))/2;

            ju(~bndry_up)=ju(~bndry_up)+jud(up,bndry_up);
            ju(~cm0)=0;
            jd(~bndry_down)=jd(~bndry_down)+jud(down,bndry_down);
            jd(~cm0)=0;           

            
            jr(~bndry_r)=jr(~bndry_r)+jlr(right,bndry_r);
            jr(~cm0)=0;
            jl(~bndry_l)=jl(~bndry_l)+jlr(left,bndry_l);
            jl(~cm0)=0;
%             sum(sum(jr+jl))
            
            jtot=ju+jd+jl+jr;
            jtot(~cm0)=0;
            udot=-jtot;
        
            if grow
                Pa0=u(vox_ref)*lambda;
            else
                Pa0=u(vox_trial)*lambda;
            end
            

            ind_plus=udot>0&isfinite(jtot);
            ind_minus=udot<0&isfinite(jtot);
            Pp0=sum(udot(ind_plus));
            Pm0=-sum(udot(ind_minus));

            Pa=Pa0/(Pp0+Pm0+Pa0);
            
            
            %                 jdiff=1-(pdiff).^1;
            udot=udot*(1-Pa);
            
            l=zeros(prod(shape),1);
            g=zeros(prod(shape),1);
            
            l(udot<0)=udot(udot<0);
            g(udot>0)=udot(udot>0);
            
            l=[floor(-l(1)); diff(floor(cumsum(-l)))];
            g=[floor(g(1)); diff(floor(cumsum(g)))];
            if sum(g-l)~=0
                disp('we done f-d up')
            end
            u=u+reshape(g-l,shape);
            
            if grow
                u(cell_inds(1:A0))=u(cell_inds(1:A0))-counts;
                u(vox_trial)=tmp;
            else
                u(cell_inds(1:A0))=u(cell_inds(1:A0))+counts;
                u(vox_trial)=0;
            end
            
            %                 diffw=sqrt(exp(-(dist.^2)/(4*N*h^2)));
            %                 p=jdiff(inds);
%             udot0=udot;
%             ind_freeze=true(size(u));
%             l=1;
%             g=l;
%             inc=true;
%             while inc
%                 
%                 ind_freeze([l jump(l,:)])=true;
%                 ind_freeze([g jump(g,:)])=true;
%                 
%                 ind_minus=udot<0&isfinite(jtot)&cm0;
%                 ind_plus=udot>0&isfinite(jtot)&cm0;
%                 ind_minus(vox_trial)=0;
%                 ind_plus(vox_trial)=0;
%                 
%                 
%                 p=u(ind_minus).*udot(ind_minus)/sum(udot(ind_minus));
%                 
%                 l=Alg2(p,rand(),0);
%                 it=find(ind_minus,l);
%                 l=it(l);
%                 
%                 ind_freeze([l jump(l,:)])=false;
%                 r2=rand();
%                 if r2<Pp
%                     
%                     pp=u(ind_plus).*udot(ind_plus)/sum(udot(ind_plus));
%                     g=Alg2(pp,rand(),0);
%                     it=find(ind_plus,g);
%                     g=it(g);
%                     ind_freeze([g jump(g,:)])=false;
%                     u(g)=u(g)+1;
%                     updates=unique([l jump(l,:) g jump(g,:)]);
%                     
%                 else
%                     if grow
%                         u(vox_trial)=u(vox_trial)+1;
%                     end
%                     %                     ll=l;
%                     updates=[l jump(l,:)];
%                 end
%                 if grow
%                     u(l)=u(l)-1;
%                 else
%                     if r2<Pm
%                         u(l)=u(l)-1;
%                     else
%                         u(vox_trial)=u(vox_trial)-1;
%                     end
%                     
%                 end
%                 jvec=[jud(up,ind_freeze) jud(down,ind_freeze) jlr(left,ind_freeze) jlr(right,ind_freeze)];
%                 jvec(bndrys(updates,:))=0;
%                 udot(updates)=-sum(jvec,2);
%                 
%                 if grow
%                     inc=u(vox_trial)~=tmp;
%                 else
%                     inc=u(vox_trial)~=0;
%                 end
%             end
            
            plotCellIm(panelC,x(:,:,j+1),cm0,i0,j0);
            caxis(panelC,'auto')
            title(panelC,'Old')
            
            title(panelB,'New')
            plotCellIm(panelB,u,cell_mask,i0,j0);
            caxis(panelB,'auto')
            %                  plotCellIm(panelA,sqrt(vi.^2+vj.^2),cm0,i0,j0);
            plotCellQuiver(panelA,-vj,vi,cm0,i0,j0);
            ut= ux.*vj+uy.*vi;
            title(panelA,'velocity')
            plotCellIm(panelD,udot0,cm0,i0,j0);
            colorbar(panelD)
            if sum(u(:))~=sum(sum(x(:,:,j+1)))
                error('we lost or gained some molecules')
            end
            x(:,:,j+1)=u;
            if any(u(:)<0)
                disp('we made something negative :(')
            end
            
%             sum(u(:))
%             
            
            %                 if shrink
            %                     p=1-p;
            %                 end
            
            %                 for k=1:tmp
            %                     l=Alg2(p,rand(),0);
            %                     x(cell_inds(l)+j*sz)=x(cell_inds(l)+j*sz)-f;
            % %                     counts(l)=counts(l)-1;
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
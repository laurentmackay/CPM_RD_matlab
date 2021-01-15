function [A,B_1,B_R,B_rho,D,GIT,Hb,I_K,I_R,I_rho,J,L_K,L_R,L_rho,N_dim,...
N_species,PAKtot,PIX,Pax_Square,Paxtot,Rac_Square,Rho_Square,T,T_integration,a,...
alpha_PAK,alpha_R,alpha_chem,bndry,bndrys,cell_inds,cell_mask,delta_P,delta_R,...
delta_rho,down,dt,e,grow_count,h,i0,i2,i_chem_0,ir0,j0,j2,jump,k_C,k_G,k_X,...
lam_a,lam_p,left,m,n,per,perim,right,shape,shrink_count,sz,time,up,v,x] = ...
FVM_integrator_2SBF_slim_fun(A,B_1,B_R,B_rho,D,GIT,Hb,I_K,I_R,I_rho,J,L_K,L_R,...
L_rho,N_dim,N_species,PAKtot,PIX,Pax_Square,Paxtot,Rac_Square,Rho_Square,T,...
T_integration,a,alpha_PAK,alpha_R,alpha_chem,bndry,bndrys,cell_inds,cell_mask,...
delta_P,delta_R,delta_rho,down,dt,e,grow_count,h,i0,i2,i_chem_0,ir0,j0,j2,jump,k_C,...
k_G,k_X,lam_a,lam_p,left,m,n,per,perim,right,shape,shrink_count,sz,time,up,v,x)



u = x(cell_inds(1:A) + i_chem_0);

interior=~bndrys&cell_mask(:);

vox=repmat((1:sz)',1,N_dim);
vox=vox(interior)';




row = find(interior);
N_ind = length(row);
i = [row'; row'];
j=[ vox(:)';  jump(row)';];

Delta=repmat([-1/h; +1/h],N_ind,1);
u_x = sparse(i,j,Delta,numel(interior),sz);



[i2,j2,v] = find(u_x);
i2=mod(i2-1,sz)+1;
u_xx = sparse(i2,j2,v/h,sz,sz);
u_xx=u_xx(cell_inds(1:A),cell_inds(1:A));



RacRatio0 = u(:,2) ./ Rac_Square;
RacRatio = (u(:,2) + u(:,7)) ./ Rac_Square;
RhoRatio = u(:,4) ./ Rho_Square;
PaxRatio = u(:,6) ./ Pax_Square;
K_is=1./((1+k_X.*PIX+k_G.*k_X.*k_C.*GIT.*PIX.*Paxtot.*PaxRatio).*(1+alpha_R.*RacRatio0)+k_G.*k_X.*GIT.*PIX);
K=alpha_R.*RacRatio0.*K_is.*(1+k_X.*PIX+k_G.*k_X.*k_C.*Paxtot.*GIT.*PIX.*PaxRatio);
I_Ks=I_K.*(1-K_is.*(1+alpha_R.*RacRatio0));
Q_R = (I_R+I_Ks).*(L_rho.^m./(L_rho.^m+RhoRatio.^m));
Q_rho = I_rho.*(L_R.^m./(L_R.^m +(RacRatio).^m));
Q_P = B_1.*(K.^n./(L_K.^n+K.^n));
f_Raci = -(Q_R.*u(:,1))+(delta_R.*u(:,2));
f_Rhoi = -(Q_rho.*u(:,3))+(delta_rho.*u(:,4));
f_Paxi = -(Q_P.*u(:,5))+(delta_P.*u(:,6));
subs__0 = 1./(((u(:,2).*alpha_R)./Rac_Square + 1).*(PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1) + GIT.*PIX.*k_G.*k_X);
subs__1 = 1./(Pax_Square.*Rac_Square + Pax_Square.*u(:,2).*alpha_R + PIX.*Pax_Square.*Rac_Square.*k_X + PIX.*Pax_Square.*u(:,2).*alpha_R.*k_X + GIT.*PIX.*Pax_Square.*Rac_Square.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*Rac_Square.*k_C.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*u(:,2).*alpha_R.*k_C.*k_G.*k_X).^2;
subs__2 = PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1;
subs__3 = (u(:,2).*alpha_R)./Rac_Square;
subs__4 = subs__3 + 1;
subs__5 = k_X.^2;
subs__6 = k_G.^2;
subs__7 = PIX.^2;
subs__8 = GIT.^2;
subs__9 = subs__0.^2;
subs__10 = GIT.*PAKtot.*PIX.*k_C.*k_G.*k_X.*subs__0.*subs__4;
subs__11 = PAKtot.*alpha_PAK.*subs__0.*subs__2;
affinity__1=subs__11;
affinity__2=subs__10;
J_gamma_1_2=subs__11 - PAKtot.*alpha_PAK.*subs__2.^2.*subs__3.*subs__9 + 1;
J_gamma_2_2=PAKtot.*u(:,6).*Pax_Square.^2.*Rac_Square.*alpha_R.*k_C.*subs__1.*subs__5.*subs__6.*subs__7.*subs__8;
J_gamma_1_6=PAKtot.*Paxtot.*Pax_Square.*u(:,2).*Rac_Square.^2.*alpha_PAK.*k_C.*subs__1.*subs__5.*subs__6.*subs__7.*subs__8;
J_gamma_2_6=subs__10 - (PAKtot.*u(:,6).*Paxtot.*k_C.^2.*subs__4.^2.*subs__5.*subs__6.*subs__7.*subs__8.*subs__9)./Pax_Square + 1;
subs2__0 = 1./(J_gamma_1_2 + J_gamma_2_6 + J_gamma_1_2.*J_gamma_2_6 - J_gamma_1_6.*J_gamma_2_2 + 1);
subs2__1 = J_gamma_2_2.*f_Raci;
subs2__2 = J_gamma_1_6.*f_Paxi;
subs2__3 = -subs2__2;
subs2__4 = -subs2__1;
subs2__5 = J_gamma_1_2.*f_Raci;
subs2__6 = J_gamma_1_2.*f_Paxi;

Rx = [f_Raci,...
-subs2__0.*(f_Raci + subs2__3 + J_gamma_2_6.*f_Raci),...
f_Rhoi,...
-f_Rhoi,...
f_Paxi,...
-subs2__0.*(f_Paxi + subs2__4 + subs2__6),...
-subs2__0.*(subs2__2 + subs2__5 + J_gamma_1_6.*subs2__4 + J_gamma_2_6.*subs2__5),...
-subs2__0.*(subs2__1 + J_gamma_2_6.*f_Paxi + J_gamma_2_2.*subs2__3 + J_gamma_2_6.*subs2__6)];
adj_empty = ~cell_mask(up) | ~cell_mask(down) | ~cell_mask(left) | ~cell_mask(right); 
adj_full = cell_mask(up) | cell_mask(down) | cell_mask(left) | cell_mask(right); 

bndry_cell = cell_mask & adj_empty;
bndry_empty = ~cell_mask & adj_full;
bndry_mask = bndry_cell | bndry_empty;
bndry = find( bndry_mask );

bndry_up=cell_mask  & ~cell_mask(up);
bndry_down=cell_mask  & ~cell_mask(down);
bndry_l=cell_mask  & ~cell_mask(left);
bndry_r=cell_mask  & ~cell_mask(right);

bndry_ud= bndry_up | bndry_down;
bndry_lr= bndry_l | bndry_r;

bndrys=[bndry_up(:) bndry_down(:) bndry_l(:) bndry_r(:)];is_discrete = all(mod(x(cell_inds(1:A)),1)==0);

if any(cell_maskp~=cell_mask)
    error('not reseting')
end


rho_eq=mean(RhoRatio(find(cell_mask)));
R_eq=mean(RacRatio(find(cell_mask)));
Ncell_mask=squeeze(sum(sum(x))); 
A0=A;

no_holes=false;

while ~no_holes
    vox_trial = bndry(randi(length(bndry)));
    
    r=randi(size(jump,2));
    vox_ref=jump(sub2ind([sz,4],vox_trial,r));
    cell_maskp(vox_trial) = cell_mask(vox_ref);
    
    Per=perim(cell_maskp); 
    A=nnz(cell_maskp); 
    HA=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; 
    dH=HA-H0;
    no_holes = getfield(bwconncomp(cell_maskp,4),'NumObjects')==1 && getfield(bwconncomp(~cell_maskp,4),'NumObjects')==1 ;
    if ~no_holes
        cell_maskp(vox_trial)=cell_mask(vox_trial);
    end
end


reacted = 0;
if  no_holes
    
    grow= cell_maskp(vox_trial) & ~cell_mask(vox_trial);
    shrink= ~cell_maskp(vox_trial) & cell_mask(vox_trial);
    
    
    if grow
        f=1;
        dH_chem=B_rho*(RhoRatio(vox_ref)-rho_eq)-B_R*(RacRatio(vox_ref)-R_eq);
        
    elseif shrink
        f=-1;
        dH_chem=-B_rho*(RhoRatio(vox_trial)-rho_eq)+B_R*(RacRatio(vox_trial)-R_eq);
       
    end
    
    
    if (grow || shrink) && rand<exp(-(dH+dH_chem+Hb)/T)
        reacted=1;
            cm0=cell_mask;
        cell_mask=cell_maskp; 
        
        if shrink
            bndry_up=cell_mask  & ~cell_mask(up);
            bndry_down=cell_mask  & ~cell_mask(down);
            bndry_l=cell_mask  & ~cell_mask(left);
            bndry_r=cell_mask  & ~cell_mask(right);

            bndry_ud= bndry_up | bndry_down;
            bndry_lr= bndry_l | bndry_r;

        end
        
            
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
                dist=max(abs(i0(vox_ref)-i0(inds)),abs(j0(vox_ref)-j0(inds)));
            else
                dist=max(abs(i0(vox_trial)-i0(inds)),abs(j0(vox_trial)-j0(inds)));
            end

            min_dist=5600;
            transport_mask=((D~=0).*D/min(D(D~=0))+(D==0).*prod(shape))*min_dist>dist;

            transport_mask(find(vox_trial==inds),:)=false;
            
        x0=x;
            inds2=inds+i_chem_0;
            i_trial=vox_trial+i_chem_0;
            
            
            if grow
                us=x(vox_ref+i_chem_0);
                Ts=sum(x(inds2));
                f=Ts./(Ts+us);
                x(i_trial)=us;
                inds2=[inds2; i_trial];
            else
                ut=x(i_trial);
                f=1+(ut./sum(x(inds2)));
                x(i_trial)=0;
            end

            if is_discrete
                x(inds2)=floor(x(inds2).*f)+[zeros(1,N_species); diff(floor(cumsum(rem(x(inds2).*f,1.0))+1e-5))]; 
            else
                x(inds2)=x(inds2).*f;
            end
            
        
        
        I=[vox_trial vox_ref]; 
        H0=HA; 
        
        
        Per=perim(cell_mask);
        A=nnz(cell_mask);
        cell_inds(1:A)=find(cell_mask);

        
        if grow
            vox=cell_inds(1:A);
        else
            vox=[cell_inds(1:A); vox_trial];
        end
        
       if length(vox)>1
[tmp,tmp2]=meshgrid(ir0,vox);
    I_rx=tmp+tmp2;
else
    I_rx=vox+ir0;
end
a_c_0=alpha_chem(I_rx);


RacRatio0(vox) = x(vox+1*sz) / Rac_Square;
RacRatio(vox) = (x(vox+1*sz) + x(vox+6*sz)) / Rac_Square;
RhoRatio(vox) = x(vox+3*sz) / Rho_Square;
PaxRatio(vox) = x(vox+5*sz) / Pax_Square;
K_is(vox)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(vox)).*(1+alpha_R*RacRatio0(vox))+k_G*k_X*GIT*PIX);
K(vox)=alpha_R*RacRatio0(vox).*K_is(vox).*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio(vox));
I_Ks(vox)=I_K*(1-K_is(vox).*(1+alpha_R*RacRatio0(vox)));
Q_R(vox) = (I_R+I_Ks(vox)).*(L_rho^m./(L_rho^m+RhoRatio(vox).^m));
Q_rho(vox) = I_rho*(L_R^m./(L_R^m +(RacRatio(vox)).^m));
Q_P(vox) = B_1*(K(vox).^n./(L_K^n+K(vox).^n));
alpha_chem(vox+0*sz)=(Q_R(vox)).*x(vox+0*sz);
alpha_chem(vox+1*sz)=(delta_R).*x(vox+1*sz);
alpha_chem(vox+2*sz)=(Q_rho(vox)).*x(vox+2*sz);
alpha_chem(vox+3*sz)=(delta_rho).*x(vox+3*sz);
alpha_chem(vox+4*sz)=(Q_P(vox)).*x(vox+4*sz);
alpha_chem(vox+5*sz)=(delta_P).*x(vox+5*sz);

alpha_rx=alpha_rx+sum(alpha_chem(I_rx)-a_c_0,1);alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
        if grow
            grow_count=grow_count+1;
        else
            shrink_count=shrink_count+1;
        end
        
    end
end

if ~reacted
    
    cell_maskp=cell_mask;
    Per=perim(cell_mask); 
    A=nnz(cell_mask); 
    cell_inds(1:A)=find(cell_mask);
else
end



Ncell_maskp=squeeze(sum(sum(x)));
if (is_discrete & any(Ncell_mask~=Ncell_maskp)) | (~is_discrete & any(abs(Ncell_mask-Ncell_maskp)>1e-5))
    error('molecule loss')
end

if min(cell_mask(:))<0
    error('Oh no! D: (negtive numbers)')
endu_prev = u;
t0=time;
eye = speye(A);
MAT_list = arrayfun(@(Di)( eye*3/(2*dt)-Di*u_xx),D,'UniformOutput',0);





while time-t0<T_integration
    

    
    Rx_prev=Rx;
    
   RacRatio0 = u(:,2) ./ Rac_Square;
RacRatio = (u(:,2) + u(:,7)) ./ Rac_Square;
RhoRatio = u(:,4) ./ Rho_Square;
PaxRatio = u(:,6) ./ Pax_Square;
K_is=1./((1+k_X.*PIX+k_G.*k_X.*k_C.*GIT.*PIX.*Paxtot.*PaxRatio).*(1+alpha_R.*RacRatio0)+k_G.*k_X.*GIT.*PIX);
K=alpha_R.*RacRatio0.*K_is.*(1+k_X.*PIX+k_G.*k_X.*k_C.*Paxtot.*GIT.*PIX.*PaxRatio);
I_Ks=I_K.*(1-K_is.*(1+alpha_R.*RacRatio0));
Q_R = (I_R+I_Ks).*(L_rho.^m./(L_rho.^m+RhoRatio.^m));
Q_rho = I_rho.*(L_R.^m./(L_R.^m +(RacRatio).^m));
Q_P = B_1.*(K.^n./(L_K.^n+K.^n));
f_Raci = -(Q_R.*u(:,1))+(delta_R.*u(:,2));
f_Rhoi = -(Q_rho.*u(:,3))+(delta_rho.*u(:,4));
f_Paxi = -(Q_P.*u(:,5))+(delta_P.*u(:,6));
subs__0 = 1./(((u(:,2).*alpha_R)./Rac_Square + 1).*(PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1) + GIT.*PIX.*k_G.*k_X);
subs__1 = 1./(Pax_Square.*Rac_Square + Pax_Square.*u(:,2).*alpha_R + PIX.*Pax_Square.*Rac_Square.*k_X + PIX.*Pax_Square.*u(:,2).*alpha_R.*k_X + GIT.*PIX.*Pax_Square.*Rac_Square.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*Rac_Square.*k_C.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*u(:,2).*alpha_R.*k_C.*k_G.*k_X).^2;
subs__2 = PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1;
subs__3 = (u(:,2).*alpha_R)./Rac_Square;
subs__4 = subs__3 + 1;
subs__5 = k_X.^2;
subs__6 = k_G.^2;
subs__7 = PIX.^2;
subs__8 = GIT.^2;
subs__9 = subs__0.^2;
subs__10 = GIT.*PAKtot.*PIX.*k_C.*k_G.*k_X.*subs__0.*subs__4;
subs__11 = PAKtot.*alpha_PAK.*subs__0.*subs__2;
affinity__1=subs__11;
affinity__2=subs__10;
J_gamma_1_2=subs__11 - PAKtot.*alpha_PAK.*subs__2.^2.*subs__3.*subs__9 + 1;
J_gamma_2_2=PAKtot.*u(:,6).*Pax_Square.^2.*Rac_Square.*alpha_R.*k_C.*subs__1.*subs__5.*subs__6.*subs__7.*subs__8;
J_gamma_1_6=PAKtot.*Paxtot.*Pax_Square.*u(:,2).*Rac_Square.^2.*alpha_PAK.*k_C.*subs__1.*subs__5.*subs__6.*subs__7.*subs__8;
J_gamma_2_6=subs__10 - (PAKtot.*u(:,6).*Paxtot.*k_C.^2.*subs__4.^2.*subs__5.*subs__6.*subs__7.*subs__8.*subs__9)./Pax_Square + 1;
subs2__0 = 1./(J_gamma_1_2 + J_gamma_2_6 + J_gamma_1_2.*J_gamma_2_6 - J_gamma_1_6.*J_gamma_2_2 + 1);
subs2__1 = J_gamma_2_2.*f_Raci;
subs2__2 = J_gamma_1_6.*f_Paxi;
subs2__3 = -subs2__2;
subs2__4 = -subs2__1;
subs2__5 = J_gamma_1_2.*f_Raci;
subs2__6 = J_gamma_1_2.*f_Paxi;

Rx = [f_Raci,...
-subs2__0.*(f_Raci + subs2__3 + J_gamma_2_6.*f_Raci),...
f_Rhoi,...
-f_Rhoi,...
f_Paxi,...
-subs2__0.*(f_Paxi + subs2__4 + subs2__6),...
-subs2__0.*(subs2__2 + subs2__5 + J_gamma_1_6.*subs2__4 + J_gamma_2_6.*subs2__5),...
-subs2__0.*(subs2__1 + J_gamma_2_6.*f_Paxi + J_gamma_2_2.*subs2__3 + J_gamma_2_6.*subs2__6)];
adj_empty = ~cell_mask(up) | ~cell_mask(down) | ~cell_mask(left) | ~cell_mask(right); 
adj_full = cell_mask(up) | cell_mask(down) | cell_mask(left) | cell_mask(right); 

bndry_cell = cell_mask & adj_empty;
bndry_empty = ~cell_mask & adj_full;
bndry_mask = bndry_cell | bndry_empty;
bndry = find( bndry_mask );

bndry_up=cell_mask  & ~cell_mask(up);
bndry_down=cell_mask  & ~cell_mask(down);
bndry_l=cell_mask  & ~cell_mask(left);
bndry_r=cell_mask  & ~cell_mask(right);

bndry_ud= bndry_up | bndry_down;
bndry_lr= bndry_l | bndry_r;

bndrys=[bndry_up(:) bndry_down(:) bndry_l(:) bndry_r(:)];is_discrete = all(mod(x(cell_inds(1:A)),1)==0);

if any(cell_maskp~=cell_mask)
    error('not reseting')
end


rho_eq=mean(RhoRatio(find(cell_mask)));
R_eq=mean(RacRatio(find(cell_mask)));
Ncell_mask=squeeze(sum(sum(x))); 
A0=A;

no_holes=false;

while ~no_holes
    vox_trial = bndry(randi(length(bndry)));
    
    r=randi(size(jump,2));
    vox_ref=jump(sub2ind([sz,4],vox_trial,r));
    cell_maskp(vox_trial) = cell_mask(vox_ref);
    
    Per=perim(cell_maskp); 
    A=nnz(cell_maskp); 
    HA=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; 
    dH=HA-H0;
    no_holes = getfield(bwconncomp(cell_maskp,4),'NumObjects')==1 && getfield(bwconncomp(~cell_maskp,4),'NumObjects')==1 ;
    if ~no_holes
        cell_maskp(vox_trial)=cell_mask(vox_trial);
    end
end


reacted = 0;
if  no_holes
    
    grow= cell_maskp(vox_trial) & ~cell_mask(vox_trial);
    shrink= ~cell_maskp(vox_trial) & cell_mask(vox_trial);
    
    
    if grow
        f=1;
        dH_chem=B_rho*(RhoRatio(vox_ref)-rho_eq)-B_R*(RacRatio(vox_ref)-R_eq);
        
    elseif shrink
        f=-1;
        dH_chem=-B_rho*(RhoRatio(vox_trial)-rho_eq)+B_R*(RacRatio(vox_trial)-R_eq);
       
    end
    
    
    if (grow || shrink) && rand<exp(-(dH+dH_chem+Hb)/T)
        reacted=1;
            cm0=cell_mask;
        cell_mask=cell_maskp; 
        
        if shrink
            bndry_up=cell_mask  & ~cell_mask(up);
            bndry_down=cell_mask  & ~cell_mask(down);
            bndry_l=cell_mask  & ~cell_mask(left);
            bndry_r=cell_mask  & ~cell_mask(right);

            bndry_ud= bndry_up | bndry_down;
            bndry_lr= bndry_l | bndry_r;

        end
        
            
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
                dist=max(abs(i0(vox_ref)-i0(inds)),abs(j0(vox_ref)-j0(inds)));
            else
                dist=max(abs(i0(vox_trial)-i0(inds)),abs(j0(vox_trial)-j0(inds)));
            end

            min_dist=5600;
            transport_mask=((D~=0).*D/min(D(D~=0))+(D==0).*prod(shape))*min_dist>dist;

            transport_mask(find(vox_trial==inds),:)=false;
            
        x0=x;
            inds2=inds+i_chem_0;
            i_trial=vox_trial+i_chem_0;
            
            
            if grow
                us=x(vox_ref+i_chem_0);
                Ts=sum(x(inds2));
                f=Ts./(Ts+us);
                x(i_trial)=us;
                inds2=[inds2; i_trial];
            else
                ut=x(i_trial);
                f=1+(ut./sum(x(inds2)));
                x(i_trial)=0;
            end

            if is_discrete
                x(inds2)=floor(x(inds2).*f)+[zeros(1,N_species); diff(floor(cumsum(rem(x(inds2).*f,1.0))+1e-5))]; 
            else
                x(inds2)=x(inds2).*f;
            end
            
        
        
        I=[vox_trial vox_ref]; 
        H0=HA; 
        
        
        Per=perim(cell_mask);
        A=nnz(cell_mask);
        cell_inds(1:A)=find(cell_mask);

        
        if grow
            vox=cell_inds(1:A);
        else
            vox=[cell_inds(1:A); vox_trial];
        end
        
       if length(vox)>1
[tmp,tmp2]=meshgrid(ir0,vox);
    I_rx=tmp+tmp2;
else
    I_rx=vox+ir0;
end
a_c_0=alpha_chem(I_rx);


RacRatio0(vox) = x(vox+1*sz) / Rac_Square;
RacRatio(vox) = (x(vox+1*sz) + x(vox+6*sz)) / Rac_Square;
RhoRatio(vox) = x(vox+3*sz) / Rho_Square;
PaxRatio(vox) = x(vox+5*sz) / Pax_Square;
K_is(vox)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(vox)).*(1+alpha_R*RacRatio0(vox))+k_G*k_X*GIT*PIX);
K(vox)=alpha_R*RacRatio0(vox).*K_is(vox).*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio(vox));
I_Ks(vox)=I_K*(1-K_is(vox).*(1+alpha_R*RacRatio0(vox)));
Q_R(vox) = (I_R+I_Ks(vox)).*(L_rho^m./(L_rho^m+RhoRatio(vox).^m));
Q_rho(vox) = I_rho*(L_R^m./(L_R^m +(RacRatio(vox)).^m));
Q_P(vox) = B_1*(K(vox).^n./(L_K^n+K(vox).^n));
alpha_chem(vox+0*sz)=(Q_R(vox)).*x(vox+0*sz);
alpha_chem(vox+1*sz)=(delta_R).*x(vox+1*sz);
alpha_chem(vox+2*sz)=(Q_rho(vox)).*x(vox+2*sz);
alpha_chem(vox+3*sz)=(delta_rho).*x(vox+3*sz);
alpha_chem(vox+4*sz)=(Q_P(vox)).*x(vox+4*sz);
alpha_chem(vox+5*sz)=(delta_P).*x(vox+5*sz);

alpha_rx=alpha_rx+sum(alpha_chem(I_rx)-a_c_0,1);alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
        if grow
            grow_count=grow_count+1;
        else
            shrink_count=shrink_count+1;
        end
        
    end
end

if ~reacted
    
    cell_maskp=cell_mask;
    Per=perim(cell_mask); 
    A=nnz(cell_mask); 
    cell_inds(1:A)=find(cell_mask);
else
end



Ncell_maskp=squeeze(sum(sum(x)));
if (is_discrete & any(Ncell_mask~=Ncell_maskp)) | (~is_discrete & any(abs(Ncell_mask-Ncell_maskp)>1e-5))
    error('molecule loss')
end

if min(cell_mask(:))<0
    error('Oh no! D: (negtive numbers)')
endb=(2*u-u_prev/2)/dt + 2*Rx-Rx_prev;
    u_prev=u;
    for i = 1:N_species
        u(:,i) = MAT_list{i}\b(:,i);
    end

    time=time+dt;


       

end
x(cell_inds(1:A) + i_chem_0) = u(:);
u = reshape(x,[sz ,size(x,3)]);RacRatio0 = u(:,2) ./ Rac_Square;
RacRatio = (u(:,2) + u(:,7)) ./ Rac_Square;
RhoRatio = u(:,4) ./ Rho_Square;
PaxRatio = u(:,6) ./ Pax_Square;
K_is=1./((1+k_X.*PIX+k_G.*k_X.*k_C.*GIT.*PIX.*Paxtot.*PaxRatio).*(1+alpha_R.*RacRatio0)+k_G.*k_X.*GIT.*PIX);
K=alpha_R.*RacRatio0.*K_is.*(1+k_X.*PIX+k_G.*k_X.*k_C.*Paxtot.*GIT.*PIX.*PaxRatio);
I_Ks=I_K.*(1-K_is.*(1+alpha_R.*RacRatio0));
Q_R = (I_R+I_Ks).*(L_rho.^m./(L_rho.^m+RhoRatio.^m));
Q_rho = I_rho.*(L_R.^m./(L_R.^m +(RacRatio).^m));
Q_P = B_1.*(K.^n./(L_K.^n+K.^n));
f_Raci = -(Q_R.*u(:,1))+(delta_R.*u(:,2));
f_Rhoi = -(Q_rho.*u(:,3))+(delta_rho.*u(:,4));
f_Paxi = -(Q_P.*u(:,5))+(delta_P.*u(:,6));
subs__0 = 1./(((u(:,2).*alpha_R)./Rac_Square + 1).*(PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1) + GIT.*PIX.*k_G.*k_X);
subs__1 = 1./(Pax_Square.*Rac_Square + Pax_Square.*u(:,2).*alpha_R + PIX.*Pax_Square.*Rac_Square.*k_X + PIX.*Pax_Square.*u(:,2).*alpha_R.*k_X + GIT.*PIX.*Pax_Square.*Rac_Square.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*Rac_Square.*k_C.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*u(:,2).*alpha_R.*k_C.*k_G.*k_X).^2;
subs__2 = PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1;
subs__3 = (u(:,2).*alpha_R)./Rac_Square;
subs__4 = subs__3 + 1;
subs__5 = k_X.^2;
subs__6 = k_G.^2;
subs__7 = PIX.^2;
subs__8 = GIT.^2;
subs__9 = subs__0.^2;
subs__10 = GIT.*PAKtot.*PIX.*k_C.*k_G.*k_X.*subs__0.*subs__4;
subs__11 = PAKtot.*alpha_PAK.*subs__0.*subs__2;
affinity__1=subs__11;
affinity__2=subs__10;
J_gamma_1_2=subs__11 - PAKtot.*alpha_PAK.*subs__2.^2.*subs__3.*subs__9 + 1;
J_gamma_2_2=PAKtot.*u(:,6).*Pax_Square.^2.*Rac_Square.*alpha_R.*k_C.*subs__1.*subs__5.*subs__6.*subs__7.*subs__8;
J_gamma_1_6=PAKtot.*Paxtot.*Pax_Square.*u(:,2).*Rac_Square.^2.*alpha_PAK.*k_C.*subs__1.*subs__5.*subs__6.*subs__7.*subs__8;
J_gamma_2_6=subs__10 - (PAKtot.*u(:,6).*Paxtot.*k_C.^2.*subs__4.^2.*subs__5.*subs__6.*subs__7.*subs__8.*subs__9)./Pax_Square + 1;
subs2__0 = 1./(J_gamma_1_2 + J_gamma_2_6 + J_gamma_1_2.*J_gamma_2_6 - J_gamma_1_6.*J_gamma_2_2 + 1);
subs2__1 = J_gamma_2_2.*f_Raci;
subs2__2 = J_gamma_1_6.*f_Paxi;
subs2__3 = -subs2__2;
subs2__4 = -subs2__1;
subs2__5 = J_gamma_1_2.*f_Raci;
subs2__6 = J_gamma_1_2.*f_Paxi;

Rx = [f_Raci,...
-subs2__0.*(f_Raci + subs2__3 + J_gamma_2_6.*f_Raci),...
f_Rhoi,...
-f_Rhoi,...
f_Paxi,...
-subs2__0.*(f_Paxi + subs2__4 + subs2__6),...
-subs2__0.*(subs2__2 + subs2__5 + J_gamma_1_6.*subs2__4 + J_gamma_2_6.*subs2__5),...
-subs2__0.*(subs2__1 + J_gamma_2_6.*f_Paxi + J_gamma_2_2.*subs2__3 + J_gamma_2_6.*subs2__6)];
adj_empty = ~cell_mask(up) | ~cell_mask(down) | ~cell_mask(left) | ~cell_mask(right); 
adj_full = cell_mask(up) | cell_mask(down) | cell_mask(left) | cell_mask(right); 

bndry_cell = cell_mask & adj_empty;
bndry_empty = ~cell_mask & adj_full;
bndry_mask = bndry_cell | bndry_empty;
bndry = find( bndry_mask );

bndry_up=cell_mask  & ~cell_mask(up);
bndry_down=cell_mask  & ~cell_mask(down);
bndry_l=cell_mask  & ~cell_mask(left);
bndry_r=cell_mask  & ~cell_mask(right);

bndry_ud= bndry_up | bndry_down;
bndry_lr= bndry_l | bndry_r;

bndrys=[bndry_up(:) bndry_down(:) bndry_l(:) bndry_r(:)];is_discrete = all(mod(x(cell_inds(1:A)),1)==0);

if any(cell_maskp~=cell_mask)
    error('not reseting')
end


rho_eq=mean(RhoRatio(find(cell_mask)));
R_eq=mean(RacRatio(find(cell_mask)));
Ncell_mask=squeeze(sum(sum(x))); 
A0=A;

no_holes=false;

while ~no_holes
    vox_trial = bndry(randi(length(bndry)));
    
    r=randi(size(jump,2));
    vox_ref=jump(sub2ind([sz,4],vox_trial,r));
    cell_maskp(vox_trial) = cell_mask(vox_ref);
    
    Per=perim(cell_maskp); 
    A=nnz(cell_maskp); 
    HA=lam_a*(a-A)^2+lam_p*(per-Per)^2+J*Per; 
    dH=HA-H0;
    no_holes = getfield(bwconncomp(cell_maskp,4),'NumObjects')==1 && getfield(bwconncomp(~cell_maskp,4),'NumObjects')==1 ;
    if ~no_holes
        cell_maskp(vox_trial)=cell_mask(vox_trial);
    end
end


reacted = 0;
if  no_holes
    
    grow= cell_maskp(vox_trial) & ~cell_mask(vox_trial);
    shrink= ~cell_maskp(vox_trial) & cell_mask(vox_trial);
    
    
    if grow
        f=1;
        dH_chem=B_rho*(RhoRatio(vox_ref)-rho_eq)-B_R*(RacRatio(vox_ref)-R_eq);
        
    elseif shrink
        f=-1;
        dH_chem=-B_rho*(RhoRatio(vox_trial)-rho_eq)+B_R*(RacRatio(vox_trial)-R_eq);
       
    end
    
    
    if (grow || shrink) && rand<exp(-(dH+dH_chem+Hb)/T)
        reacted=1;
            cm0=cell_mask;
        cell_mask=cell_maskp; 
        
        if shrink
            bndry_up=cell_mask  & ~cell_mask(up);
            bndry_down=cell_mask  & ~cell_mask(down);
            bndry_l=cell_mask  & ~cell_mask(left);
            bndry_r=cell_mask  & ~cell_mask(right);

            bndry_ud= bndry_up | bndry_down;
            bndry_lr= bndry_l | bndry_r;

        end
        
            
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
                dist=max(abs(i0(vox_ref)-i0(inds)),abs(j0(vox_ref)-j0(inds)));
            else
                dist=max(abs(i0(vox_trial)-i0(inds)),abs(j0(vox_trial)-j0(inds)));
            end

            min_dist=5600;
            transport_mask=((D~=0).*D/min(D(D~=0))+(D==0).*prod(shape))*min_dist>dist;

            transport_mask(find(vox_trial==inds),:)=false;
            
        x0=x;
            inds2=inds+i_chem_0;
            i_trial=vox_trial+i_chem_0;
            
            
            if grow
                us=x(vox_ref+i_chem_0);
                Ts=sum(x(inds2));
                f=Ts./(Ts+us);
                x(i_trial)=us;
                inds2=[inds2; i_trial];
            else
                ut=x(i_trial);
                f=1+(ut./sum(x(inds2)));
                x(i_trial)=0;
            end

            if is_discrete
                x(inds2)=floor(x(inds2).*f)+[zeros(1,N_species); diff(floor(cumsum(rem(x(inds2).*f,1.0))+1e-5))]; 
            else
                x(inds2)=x(inds2).*f;
            end
            
        
        
        I=[vox_trial vox_ref]; 
        H0=HA; 
        
        
        Per=perim(cell_mask);
        A=nnz(cell_mask);
        cell_inds(1:A)=find(cell_mask);

        
        if grow
            vox=cell_inds(1:A);
        else
            vox=[cell_inds(1:A); vox_trial];
        end
        
       if length(vox)>1
[tmp,tmp2]=meshgrid(ir0,vox);
    I_rx=tmp+tmp2;
else
    I_rx=vox+ir0;
end
a_c_0=alpha_chem(I_rx);


RacRatio0(vox) = x(vox+1*sz) / Rac_Square;
RacRatio(vox) = (x(vox+1*sz) + x(vox+6*sz)) / Rac_Square;
RhoRatio(vox) = x(vox+3*sz) / Rho_Square;
PaxRatio(vox) = x(vox+5*sz) / Pax_Square;
K_is(vox)=1./((1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(vox)).*(1+alpha_R*RacRatio0(vox))+k_G*k_X*GIT*PIX);
K(vox)=alpha_R*RacRatio0(vox).*K_is(vox).*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio(vox));
I_Ks(vox)=I_K*(1-K_is(vox).*(1+alpha_R*RacRatio0(vox)));
Q_R(vox) = (I_R+I_Ks(vox)).*(L_rho^m./(L_rho^m+RhoRatio(vox).^m));
Q_rho(vox) = I_rho*(L_R^m./(L_R^m +(RacRatio(vox)).^m));
Q_P(vox) = B_1*(K(vox).^n./(L_K^n+K(vox).^n));
alpha_chem(vox+0*sz)=(Q_R(vox)).*x(vox+0*sz);
alpha_chem(vox+1*sz)=(delta_R).*x(vox+1*sz);
alpha_chem(vox+2*sz)=(Q_rho(vox)).*x(vox+2*sz);
alpha_chem(vox+3*sz)=(delta_rho).*x(vox+3*sz);
alpha_chem(vox+4*sz)=(Q_P(vox)).*x(vox+4*sz);
alpha_chem(vox+5*sz)=(delta_P).*x(vox+5*sz);

alpha_rx=alpha_rx+sum(alpha_chem(I_rx)-a_c_0,1);alpha_rx=sum(alpha_chem(ir0 + cell_inds(1:A)));
        if grow
            grow_count=grow_count+1;
        else
            shrink_count=shrink_count+1;
        end
        
    end
end

if ~reacted
    
    cell_maskp=cell_mask;
    Per=perim(cell_mask); 
    A=nnz(cell_mask); 
    cell_inds(1:A)=find(cell_mask);
else
end



Ncell_maskp=squeeze(sum(sum(x)));
if (is_discrete & any(Ncell_mask~=Ncell_maskp)) | (~is_discrete & any(abs(Ncell_mask-Ncell_maskp)>1e-5))
    error('molecule loss')
end

if min(cell_mask(:))<0
    error('Oh no! D: (negtive numbers)')
end
end
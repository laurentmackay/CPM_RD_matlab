function [A,B_1,D,GIT,I_K,I_R,I_rho,L_K,L_R,L_rho,N_dim,N_slow,N_species,PAKtot,...
PIX,Pax_Square,Paxtot,Rac_Square,Rho_Square,alpha_PAK,alpha_R,bndrys,cell_inds,...
cell_mask,delta_P,delta_R,delta_rho,dt,h,jump,k_C,k_G,k_X,m,n,shape,sz,x] = FVM_integrator_2_fun(...
A,B_1,D,GIT,I_K,I_R,I_rho,L_K,L_R,L_rho,N_dim,N_slow,N_species,PAKtot,PIX,Pax_Square,...
Paxtot,Rac_Square,Rho_Square,alpha_PAK,alpha_R,bndrys,cell_inds,cell_mask,delta_P,...
delta_R,delta_rho,dt,h,jump,k_C,k_G,k_X,m,n,shape,sz,x)
T_final=1e4;
dt=0.003;
t_plot=100;


figure(3);clf();

% u = reshape(x,[sz,size(x,3)]);
% u = u()
u = x(cell_inds(1:A) + ((1:N_species)-1)*sz);

% u(i0(:)>40 & u(:,2)>0,2) = u(i0(:)>40 & u(:,2)>0,2)/2;
%
% u(i0(:)<40 & u(:,2)>0,2) = u(i0(:)<40 & u(:,2)>0,2) + sum( u(i0(:)>40 & u(:,2)>0,2))/nnz(u(i0(:)<40 & u(:,2)>0));

interior=~bndrys&cell_mask(:);

vox=repmat((1:sz)',1,N_dim);
vox=vox(interior)';

% ind_0(:) = [ind_0(interior(:,1:2)) ind_0(interior(:,3:4))

%%

dir=1:N_dim;

row = find(interior);
N_ind = length(row);
i = [row'; row'];
j=[ vox(:)';  jump(row)';];

Delta=repmat([-1/h; +1/h],N_ind,1);
u_x = sparse(i,j,Delta,numel(interior),sz);

Delta2 = repelem(repmat([-1/h -1/h],1,N_dim/2),sum(interior));
v=nonzeros(Delta2.*u_x(row,:)');
i2 = repelem(vox,2);
j2 = mod(find(u_x(row,:)')-1,sz)+1;
u_xx = sparse(i2,j2,v,sz,sz);

u_xx=u_xx(cell_inds(1:A),cell_inds(1:A));

%%

%     u_xx = (u_x(jump(:,1),:)-u_x(jump(:,2),:)+...
%         u_x(jump(:,3),:)-u_x(jump(:,4),:))/h;

t=0;
t_last=0;
im=nan(shape);
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
subs__3 = (u(:,2).*alpha_R)./Rac_Square + 1;
subs__4 = k_X.^2;
subs__5 = k_G.^2;
subs__6 = PIX.^2;
subs__7 = GIT.^2;
subs__8 = subs__0.^2;
affinity__1=PAKtot.*alpha_PAK.*subs__0.*subs__2;
affinity__2=GIT.*PAKtot.*PIX.*k_C.*k_G.*k_X.*subs__0.*subs__3;
J_gamma_1_2=1 - (PAKtot.*alpha_R.*alpha_PAK.*subs__2.^2.*subs__8)./Rac_Square;
J_gamma_2_2=PAKtot.*Pax_Square.^2.*Rac_Square.*alpha_R.*k_C.*subs__1.*subs__4.*subs__5.*subs__6.*subs__7;
J_gamma_1_6=PAKtot.*Paxtot.*Pax_Square.*Rac_Square.^2.*alpha_PAK.*k_C.*subs__1.*subs__4.*subs__5.*subs__6.*subs__7;
J_gamma_2_6=1 - (PAKtot.*Paxtot.*k_C.^2.*subs__3.^2.*subs__4.*subs__5.*subs__6.*subs__7.*subs__8)./Pax_Square;
sub2__0 = 1./(J_gamma_1_2 + J_gamma_2_6 + J_gamma_1_2.*J_gamma_2_6 - J_gamma_1_6.*J_gamma_2_2 + 1);

Rx = [f_Raci,...
-sub2__0.*(f_Raci - J_gamma_1_6.*f_Paxi + J_gamma_2_6.*f_Raci),...
f_Rhoi,...
-f_Rhoi,...
f_Paxi,...
-sub2__0.*(f_Paxi + J_gamma_1_2.*f_Paxi - J_gamma_2_2.*f_Raci)];
%  dt * u_xx*(D(1:N_slow).*u(:,1:N_slow))
tic;

while t<=T_final

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
    subs__3 = (u(:,2).*alpha_R)./Rac_Square + 1;
    subs__4 = k_X.^2;
    subs__5 = k_G.^2;
    subs__6 = PIX.^2;
    subs__7 = GIT.^2;
    subs__8 = subs__0.^2;
    affinity__1=PAKtot.*alpha_PAK.*subs__0.*subs__2;
    affinity__2=GIT.*PAKtot.*PIX.*k_C.*k_G.*k_X.*subs__0.*subs__3;
    J_gamma_1_2=1 - (PAKtot.*alpha_R.*alpha_PAK.*subs__2.^2.*subs__8)./Rac_Square;
    J_gamma_2_2=PAKtot.*Pax_Square.^2.*Rac_Square.*alpha_R.*k_C.*subs__1.*subs__4.*subs__5.*subs__6.*subs__7;
    J_gamma_1_6=PAKtot.*Paxtot.*Pax_Square.*Rac_Square.^2.*alpha_PAK.*k_C.*subs__1.*subs__4.*subs__5.*subs__6.*subs__7;
    J_gamma_2_6=1 - (PAKtot.*Paxtot.*k_C.^2.*subs__3.^2.*subs__4.*subs__5.*subs__6.*subs__7.*subs__8)./Pax_Square;
    sub2__0 = 1./(J_gamma_1_2 + J_gamma_2_6 + J_gamma_1_2.*J_gamma_2_6 - J_gamma_1_6.*J_gamma_2_2 + 1);
    
    Rx = [f_Raci,...
    -sub2__0.*(f_Raci - J_gamma_1_6.*f_Paxi + J_gamma_2_6.*f_Raci),...
    f_Rhoi,...
    -f_Rhoi,...
    f_Paxi,...
    -sub2__0.*(f_Paxi + J_gamma_1_2.*f_Paxi - J_gamma_2_2.*f_Raci)];

    u(:,1:N_slow) = u(:,1:N_slow) - dt * u_xx*(D(1:N_slow).*u(:,1:N_slow)) + 3*dt/2*Rx - dt/2*Rx_prev; %two-step adams bashforth

    utot = u(:,2)+u(:,7);
    u(:,2) = utot./(1+affinity__1);
    u(:,7) = utot - u(:,2);
    utot = u(:,6)+u(:,8);
    u(:,6) = utot./(1+affinity__2);
    u(:,8) = utot - u(:,6);


    t=t+dt;
    %
    if t-t_last>=t_plot || t>=T_final
        t_plot/toc
        im(cell_inds(1:A))=u(:,2)/Rac_Square;
        imagesc(im); colorbar();
        title(['time = ' num2str(t)])
        drawnow;
        tic
        t_last=t;
    end
end
end

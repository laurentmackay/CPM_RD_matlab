function [B_1,D,GIT,I_K,I_R,I_rho,L_K,L_R,L_rho,N_dim,N_slow,PAKtot,PIX,Pax_Square,...
Paxtot,RacRatio0,Rac_Square,Rho_Square,Rx,T_final,alpha_PAK,alpha_R,bndry_mask,...
bndrys,cell_mask,delta_P,delta_R,delta_rho,dt,h,jump,k_C,k_G,k_X,m,n,sz,x] = FVM_integrator_fun(...
B_1,D,GIT,I_K,I_R,I_rho,L_K,L_R,L_rho,N_dim,N_slow,PAKtot,PIX,Pax_Square,Paxtot,RacRatio0,...
Rac_Square,Rho_Square,Rx,T_final,alpha_PAK,alpha_R,bndry_mask,bndrys,cell_mask,delta_P,...
delta_R,delta_rho,dt,h,jump,k_C,k_G,k_X,m,n,sz,x)

u = reshape(x,[sz,size(x,3)]);

interior=~bndrys&cell_mask(:);

inds=repmat((1:sz)',N_dim);
inds=inds(interior)';

dir=1:N_dim;
i = [inds; inds];
j=[ jump(interior)'; inds];

u_x = sparse(i,j,repmat([+1; -1],size(j,2),1),sz,sz)/h;


u_xx = (u_x(:,jump(:,1))-u_x(:,jump(:,2))+...
        u_x(:,jump(:,3))-u_x(:,jump(:,4)))/h;

%     u_xx = (u_x(jump(:,1),:)-u_x(jump(:,2),:)+...
%         u_x(jump(:,3),:)-u_x(jump(:,4),:))/h;

    t=0;

RacRatio0(:) = u(:,2) ./ Rac_Square;
RacRatio = (u(:,2) + u(:,7)) ./ Rac_Square;
RhoRatio = u(:,4) ./ Rho_Square;
PaxRatio = u(:,6) ./ Pax_Square;
K_is=1./((1+k_X.*PIX+k_G.*k_X.*k_C.*GIT.*PIX.*Paxtot.*PaxRatio(:)).*(1+alpha_R.*RacRatio0(:))+k_G.*k_X.*GIT.*PIX);
K=alpha_R.*RacRatio0(:).*K_is(:).*(1+k_X.*PIX+k_G.*k_X.*k_C.*Paxtot.*GIT.*PIX.*PaxRatio(:));
I_Ks=I_K.*(1-K_is(:).*(1+alpha_R.*RacRatio0(:)));
Q_R = (I_R+I_Ks(:)).*(L_rho.^m./(L_rho.^m+RhoRatio(:).^m));
Q_rho = I_rho.*(L_R.^m./(L_R.^m +(RacRatio(:)).^m));
Q_P = B_1.*bndry_mask(:).*(K(:).^n./(L_K.^n+K(:).^n));
f_Raci = -(Q_R(:).*u(:,1))+(delta_R.*u(:,2));
f_Rhoi = -(Q_rho(:).*u(:,3))+(delta_rho.*u(:,4));
f_Rho = (Q_rho(:).*u(:,3))-(delta_rho.*u(:,4));
f_Paxi = -(Q_P(:).*u(:,5))+(delta_P.*u(:,6));
J_gamma_plus_ell1_1=1;
J_gamma_plus_ell1_2=1 - (PAKtot.*alpha_R.*alpha_PAK.*(PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1).^2)./(Rac_Square.*(((u(:,2).*alpha_R)./Rac_Square + 1).*(PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1) + GIT.*PIX.*k_G.*k_X).^2);
J_gamma_plus_ell2_2=(GIT.^2.*PAKtot.*PIX.^2.*Pax_Square.^2.*Rac_Square.*alpha_R.*k_C.*k_G.^2.*k_X.^2)./(Pax_Square.*Rac_Square + Pax_Square.*u(:,2).*alpha_R + PIX.*Pax_Square.*Rac_Square.*k_X + PIX.*Pax_Square.*u(:,2).*alpha_R.*k_X + GIT.*PIX.*Pax_Square.*Rac_Square.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*Rac_Square.*k_C.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*u(:,2).*alpha_R.*k_C.*k_G.*k_X).^2;
J_gamma_plus_ell2_5=1;
J_gamma_plus_ell1_6=(GIT.^2.*PAKtot.*PIX.^2.*Paxtot.*Pax_Square.*Rac_Square.^2.*alpha_PAK.*k_C.*k_G.^2.*k_X.^2)./(Pax_Square.*Rac_Square + Pax_Square.*u(:,2).*alpha_R + PIX.*Pax_Square.*Rac_Square.*k_X + PIX.*Pax_Square.*u(:,2).*alpha_R.*k_X + GIT.*PIX.*Pax_Square.*Rac_Square.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*Rac_Square.*k_C.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*u(:,2).*alpha_R.*k_C.*k_G.*k_X).^2;
J_gamma_plus_ell2_6=1 - (GIT.^2.*PAKtot.*PIX.^2.*Paxtot.*k_C.^2.*k_G.^2.*k_X.^2.*((u(:,2).*alpha_R)./Rac_Square + 1).^2)./(Pax_Square.*(((u(:,2).*alpha_R)./Rac_Square + 1).*(PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1) + GIT.*PIX.*k_G.*k_X).^2);

Rx(:) = [f_Raci,...
(J_gamma_plus_ell1_6.*J_gamma_plus_ell2_5.*f_Paxi - J_gamma_plus_ell1_1.*J_gamma_plus_ell2_6.*f_Raci)./(J_gamma_plus_ell1_2.*J_gamma_plus_ell2_6 - J_gamma_plus_ell1_6.*J_gamma_plus_ell2_2),...
f_Rhoi,...
f_Rho,...
f_Paxi,...
-(J_gamma_plus_ell1_2.*J_gamma_plus_ell2_5.*f_Paxi - J_gamma_plus_ell1_1.*J_gamma_plus_ell2_2.*f_Raci)./(J_gamma_plus_ell1_2.*J_gamma_plus_ell2_6 - J_gamma_plus_ell1_6.*J_gamma_plus_ell2_2)];
%  dt * u_xx*(D(1:N_slow).*u(:,1:N_slow))

while t<=T_final

    Rx_prev=Rx;
    RacRatio0(:) = u(:,2) ./ Rac_Square;
    RacRatio = (u(:,2) + u(:,7)) ./ Rac_Square;
    RhoRatio = u(:,4) ./ Rho_Square;
    PaxRatio = u(:,6) ./ Pax_Square;
    K_is=1./((1+k_X.*PIX+k_G.*k_X.*k_C.*GIT.*PIX.*Paxtot.*PaxRatio(:)).*(1+alpha_R.*RacRatio0(:))+k_G.*k_X.*GIT.*PIX);
    K=alpha_R.*RacRatio0(:).*K_is(:).*(1+k_X.*PIX+k_G.*k_X.*k_C.*Paxtot.*GIT.*PIX.*PaxRatio(:));
    I_Ks=I_K.*(1-K_is(:).*(1+alpha_R.*RacRatio0(:)));
    Q_R = (I_R+I_Ks(:)).*(L_rho.^m./(L_rho.^m+RhoRatio(:).^m));
    Q_rho = I_rho.*(L_R.^m./(L_R.^m +(RacRatio(:)).^m));
    Q_P = B_1.*bndry_mask(:).*(K(:).^n./(L_K.^n+K(:).^n));
    f_Raci = -(Q_R(:).*u(:,1))+(delta_R.*u(:,2));
    f_Rhoi = -(Q_rho(:).*u(:,3))+(delta_rho.*u(:,4));
    f_Rho = (Q_rho(:).*u(:,3))-(delta_rho.*u(:,4));
    f_Paxi = -(Q_P(:).*u(:,5))+(delta_P.*u(:,6));
    J_gamma_plus_ell1_1=1;
    J_gamma_plus_ell1_2=1 - (PAKtot.*alpha_R.*alpha_PAK.*(PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1).^2)./(Rac_Square.*(((u(:,2).*alpha_R)./Rac_Square + 1).*(PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1) + GIT.*PIX.*k_G.*k_X).^2);
    J_gamma_plus_ell2_2=(GIT.^2.*PAKtot.*PIX.^2.*Pax_Square.^2.*Rac_Square.*alpha_R.*k_C.*k_G.^2.*k_X.^2)./(Pax_Square.*Rac_Square + Pax_Square.*u(:,2).*alpha_R + PIX.*Pax_Square.*Rac_Square.*k_X + PIX.*Pax_Square.*u(:,2).*alpha_R.*k_X + GIT.*PIX.*Pax_Square.*Rac_Square.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*Rac_Square.*k_C.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*u(:,2).*alpha_R.*k_C.*k_G.*k_X).^2;
    J_gamma_plus_ell2_5=1;
    J_gamma_plus_ell1_6=(GIT.^2.*PAKtot.*PIX.^2.*Paxtot.*Pax_Square.*Rac_Square.^2.*alpha_PAK.*k_C.*k_G.^2.*k_X.^2)./(Pax_Square.*Rac_Square + Pax_Square.*u(:,2).*alpha_R + PIX.*Pax_Square.*Rac_Square.*k_X + PIX.*Pax_Square.*u(:,2).*alpha_R.*k_X + GIT.*PIX.*Pax_Square.*Rac_Square.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*Rac_Square.*k_C.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*u(:,2).*alpha_R.*k_C.*k_G.*k_X).^2;
    J_gamma_plus_ell2_6=1 - (GIT.^2.*PAKtot.*PIX.^2.*Paxtot.*k_C.^2.*k_G.^2.*k_X.^2.*((u(:,2).*alpha_R)./Rac_Square + 1).^2)./(Pax_Square.*(((u(:,2).*alpha_R)./Rac_Square + 1).*(PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1) + GIT.*PIX.*k_G.*k_X).^2);
    
    Rx(:) = [f_Raci,...
    (J_gamma_plus_ell1_6.*J_gamma_plus_ell2_5.*f_Paxi - J_gamma_plus_ell1_1.*J_gamma_plus_ell2_6.*f_Raci)./(J_gamma_plus_ell1_2.*J_gamma_plus_ell2_6 - J_gamma_plus_ell1_6.*J_gamma_plus_ell2_2),...
    f_Rhoi,...
    f_Rho,...
    f_Paxi,...
    -(J_gamma_plus_ell1_2.*J_gamma_plus_ell2_5.*f_Paxi - J_gamma_plus_ell1_1.*J_gamma_plus_ell2_2.*f_Raci)./(J_gamma_plus_ell1_2.*J_gamma_plus_ell2_6 - J_gamma_plus_ell1_6.*J_gamma_plus_ell2_2)];

    u(:,1:N_slow) = u(:,1:N_slow) + dt * u_xx*(D(1:N_slow).*u(:,1:N_slow)) + 3*dt/2*Rx - dt/2*Rx_prev; %two-step adams bashforth

    utot = sum(u(:,2)+u(:,7),2);
    u(:,2) = utot./(1+(1+k_X.*PIX+k_G.*k_X.*k_C.*GIT.*PIX.*Paxtot.*PaxRatio(:)).*alpha_PAK.*PAKtot.*K_is(:));
    u(:,7) = utot - u(:,2);
    utot = sum(u(:,6)+u(:,7),2);
    u(:,6) = utot./(1+(k_G.*k_X.*k_C.*GIT.*PIX.*K_is(:).*PAKtot.*(1+alpha_R.*RacRatio0(:))));
    u(:,8) = utot - u(:,6);


    t=t+dt;

end
end

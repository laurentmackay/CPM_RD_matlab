RacRatio0 = u(:,2) ./ Rac_Square;
RhoRatio = u(:,4) ./ Rho_Square;
PaxRatio = u(:,6) ./ Pax_Square;
K_is=1.0./((1.0+k_X.*PIX+k_G.*k_X.*k_C.*GIT.*PIX.*Paxtot.*PaxRatio).*(1+alpha_R.*RacRatio0)+k_G.*k_X.*GIT.*PIX);
K=alpha_R.*RacRatio0.*K_is.*(1+k_X.*PIX+k_G.*k_X.*k_C.*Paxtot.*GIT.*PIX.*PaxRatio);
I_Ks=I_K.*(1.0-K_is.*(1+alpha_R.*RacRatio0));
P_i=1.0-PaxRatio.*(1+k_G.*k_X.*k_C.*GIT.*PIX.*PAKtot.*K_is.*(1+alpha_R.*RacRatio0));
RacRatio = (u(:,2) + alpha_PAK.*K) ./ Rac_Square;
Q_R = ((1-RacRatio)./(1-RacRatio0)).*(I_R+I_Ks).*(L_rho.^m./(L_rho.^m+RhoRatio.^m));
Q_rho = I_rho.*(L_R.^m./(L_R.^m +(RacRatio).^m));
Q_P = (P_i./(1-PaxRatio)).*B.*(K.^n./(L_K.^n+K.^n));
f_Raci = -(Q_R.*u(:,1))+ (delta_R.*u(:,2));
f_Rhoi = -(Q_rho.*u(:,3))+ (delta_rho.*u(:,4));
f_Paxi = -(Q_P.*u(:,5))+ (delta_P.*u(:,6));

Rx = [f_Raci,...
-f_Raci,...
f_Rhoi,...
-f_Rhoi,...
f_Paxi,...
-f_Paxi];
RacRatio0 = u(:,1) ./ Rac_Square;
PaxRatio = u(:,3) ./ Pax_Square;
K_is=1./((1+k_X.*PIX+k_G.*k_X.*k_C.*GIT.*PIX.*Paxtot.*PaxRatio).*(1+alpha_R.*RacRatio0)+k_G.*k_X.*GIT.*PIX);
RhoRatio = u(:,2) ./ Rho_Square;
K=alpha_R.*RacRatio0.*K_is.*(1+k_X.*PIX+k_G.*k_X.*k_C.*Paxtot.*GIT.*PIX.*PaxRatio);
RacRatio = (u(:,1) + alpha_PAK.*K) ./ Rac_Square;
I_Ks=I_K.*(1-K_is.*(1+alpha_R.*RacRatio0));
Q_R = (I_R+I_Ks).*(L_rho.^m./(L_rho.^m+RhoRatio.^m));
Q_rho = I_rho.*(L_R.^m./(L_R.^m +(RacRatio).^m));
Q_P = B_1.*(K.^n./(L_K.^n+K.^n));
f_Raci = -(Q_R.*(cnsrv_1 - u(:,1)))+ (delta_R.*u(:,1));
f_Rac =  (Q_R.*(cnsrv_1 - u(:,1)))-(delta_R.*u(:,1));
f_Rhoi = -(Q_rho.*(cnsrv_2 - u(:,2)))+ (delta_rho.*u(:,2));
f_Rho =  (Q_rho.*(cnsrv_2 - u(:,2)))-(delta_rho.*u(:,2));
f_Paxi = -(Q_P.*(cnsrv_3 - u(:,3)))+ (delta_P.*u(:,3));
f_Pax =  (Q_P.*(cnsrv_3 - u(:,3)))-(delta_P.*u(:,3));

Rx = [- u(:,1).*delta_R - Q_R.*(u(:,1) - cnsrv_1),...
- u(:,2).*delta_rho - Q_rho.*(u(:,2) - cnsrv_2),...
- u(:,3).*delta_P - Q_P.*(u(:,3) - cnsrv_3)];
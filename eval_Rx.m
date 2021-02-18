R = u(:,2) ./ Rac_Square;
RhoRatio = u(:,4) ./ Rho_Square;
PaxRatio = u(:,6) ./ Pax_Square;
K_is=1.0./((1.0+k_X.*PIX+k_G.*k_X.*k_C.*GIT.*PIX.*Paxtot.*PaxRatio).*(1+alpha_R.*R)+k_G.*k_X.*GIT.*PIX);
K0=alpha_R.*K_is.*(1+k_X.*PIX+k_G.*k_X.*k_C.*Paxtot.*GIT.*PIX.*PaxRatio);
K=R.*K0;
I_Ks=I_K.*(1.0-K_is.*(1+alpha_R.*R));
P_i=1.0-PaxRatio.*(1+k_G.*k_X.*k_C.*GIT.*PIX.*PAKtot.*K_is.*(1+alpha_R.*R));
Rbar = (u(:,2) + ((K0.*u(:,2).*alpha_PAK)./Rac_Square)) ./ Rac_Square;
Q_R = (I_R+I_Ks).*(L_rho.^m./(L_rho.^m+RhoRatio.^m));
Q_rho = I_rho.*(L_R.^m./(L_R.^m +(Rbar).^m));
Q_P = B.*(K.^n./(L_K.^n+K.^n));
f_Raci = -(Q_R.*u(:,1))+ (delta_R.*u(:,2));
f_Rhoi = -(Q_rho.*u(:,3))+ (delta_rho.*u(:,4));
f_Paxi = -(Q_P.*u(:,5))+ (delta_P.*u(:,6));
subs__0 = 1./(((u(:,2).*alpha_R)./Rac_Square + 1).*(PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1) + GIT.*PIX.*k_G.*k_X);
subs__1 = 1./(Pax_Square.*Rac_Square + Pax_Square.*u(:,2).*alpha_R + PIX.*Pax_Square.*Rac_Square.*k_X + PIX.*Pax_Square.*u(:,2).*alpha_R.*k_X + GIT.*PIX.*Pax_Square.*Rac_Square.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*Rac_Square.*k_C.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*u(:,2).*alpha_R.*k_C.*k_G.*k_X).^2;
subs__2 = PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1;
subs__3 = (u(:,2).*alpha_R)./Rac_Square;
subs__4 = 1./Pax_Square;
subs__5 = subs__3 + 1;
subs__6 = k_X.^2;
subs__7 = k_G.^2;
subs__8 = PIX.^2;
subs__9 = GIT.^2;
subs__10 = 1./Rac_Square;
subs__11 = subs__0.^2;
subs__12 = GIT.*PAKtot.*PIX.*k_C.*k_G.*k_X.*subs__0.*subs__4.*subs__5;
J_gamma_1_2=(alpha_R.*alpha_PAK.*(PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1))./(Rac_Square.*(((u(:,2).*alpha_R)./Rac_Square + 1).*(PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1) + GIT.*PIX.*k_G.*k_X)) - (u(:,2).*alpha_R.^2.*alpha_PAK.*(PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1).^2)./(Rac_Square.^2.*(((u(:,2).*alpha_R)./Rac_Square + 1).*(PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1) + GIT.*PIX.*k_G.*k_X).^2);
J_gamma_2_2=(GIT.^2.*PAKtot.*PIX.^2.*u(:,6).*Pax_Square.*Rac_Square.*alpha_R.*k_C.*k_G.^2.*k_X.^2)./(Pax_Square.*Rac_Square + Pax_Square.*u(:,2).*alpha_R + PIX.*Pax_Square.*Rac_Square.*k_X + PIX.*Pax_Square.*u(:,2).*alpha_R.*k_X + GIT.*PIX.*Pax_Square.*Rac_Square.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*Rac_Square.*k_C.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*u(:,2).*alpha_R.*k_C.*k_G.*k_X).^2;
J_gamma_1_6=(GIT.^2.*PIX.^2.*Paxtot.*Pax_Square.*u(:,2).*Rac_Square.*alpha_R.*alpha_PAK.*k_C.*k_G.^2.*k_X.^2)./(Pax_Square.*Rac_Square + Pax_Square.*u(:,2).*alpha_R + PIX.*Pax_Square.*Rac_Square.*k_X + PIX.*Pax_Square.*u(:,2).*alpha_R.*k_X + GIT.*PIX.*Pax_Square.*Rac_Square.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*Rac_Square.*k_C.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*u(:,2).*alpha_R.*k_C.*k_G.*k_X).^2;
J_gamma_2_6=(GIT.*PAKtot.*PIX.*k_C.*k_G.*k_X.*((u(:,2).*alpha_R)./Rac_Square + 1))./(Pax_Square.*(((u(:,2).*alpha_R)./Rac_Square + 1).*(PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1) + GIT.*PIX.*k_G.*k_X)) - (GIT.^2.*PAKtot.*PIX.^2.*u(:,6).*Paxtot.*k_C.^2.*k_G.^2.*k_X.^2.*((u(:,2).*alpha_R)./Rac_Square + 1).^2)./(Pax_Square.^2.*(((u(:,2).*alpha_R)./Rac_Square + 1).*(PIX.*k_X + (GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)./Pax_Square + 1) + GIT.*PIX.*k_G.*k_X).^2);
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
RacRatio0 = u(:,2) ./ Rac_Square;
RacRatio = (u(:,2) + ((u(:,2).*(Pax_Square.*Rac_Square + Pax_Square.*u(:,2).*alpha_R + PIX.*Pax_Square.*Rac_Square.*k_X + PIX.*Pax_Square.*u(:,2).*alpha_R.*k_X + GIT.*PIX.*Pax_Square.*Rac_Square.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*Rac_Square.*k_C.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*u(:,2).*alpha_R.*k_C.*k_G.*k_X))./(PAKtot.*Rac_Square.*alpha_PAK.*(Pax_Square + PIX.*Pax_Square.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X)))) ./ Rac_Square;
RhoRatio = u(:,4) ./ Rho_Square;
PaxRatio = u(:,6) ./ Pax_Square;
K_is=1./((1+k_X.*PIX+k_G.*k_X.*k_C.*GIT.*PIX.*Paxtot.*PaxRatio).*(1+alpha_R.*RacRatio0)+k_G.*k_X.*GIT.*PIX);
K=alpha_R.*RacRatio0.*K_is.*(1+k_X.*PIX+k_G.*k_X.*k_C.*Paxtot.*GIT.*PIX.*PaxRatio);
I_Ks=I_K.*(1-K_is.*(1+alpha_R.*RacRatio0));
Q_R = (I_R+I_Ks).*(L_rho.^m./(L_rho.^m+RhoRatio.^m));
Q_rho = I_rho.*(L_R.^m./(L_R.^m +(RacRatio).^m));
Q_P = B_1.*(K.^n./(L_K.^n+K.^n));
f_Raci = -(Q_R.*(cnsrv_1 - u(:,2) - (u(:,2).*(Pax_Square.*Rac_Square + Pax_Square.*u(:,2).*alpha_R + PIX.*Pax_Square.*Rac_Square.*k_X + PIX.*Pax_Square.*u(:,2).*alpha_R.*k_X + GIT.*PIX.*Pax_Square.*Rac_Square.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*Rac_Square.*k_C.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*u(:,2).*alpha_R.*k_C.*k_G.*k_X))./(PAKtot.*Rac_Square.*alpha_PAK.*(Pax_Square + PIX.*Pax_Square.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X))))+(delta_R.*u(:,2));
f_Rhoi = -(Q_rho.*u(:,3))+(delta_rho.*u(:,4));
f_Paxi = -(Q_P.*u(:,5))+(delta_P.*u(:,6));
subs__0 = GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X;
subs__1 = Pax_Square.*Rac_Square + Rac_Square.*subs__0 + Pax_Square.*u(:,2).*alpha_R + u(:,2).*alpha_R.*subs__0 + PIX.*Pax_Square.*Rac_Square.*k_X + PIX.*Pax_Square.*u(:,2).*alpha_R.*k_X + GIT.*PIX.*Pax_Square.*Rac_Square.*k_G.*k_X;
subs__2 = PIX.*Pax_Square.*k_X;
subs__3 = 1./(Rac_Square + u(:,2).*alpha_R);
subs__4 = 1./PAKtot;
subs__5 = 1./(Pax_Square + subs__0 + subs__2);
subs__6 = 1./k_C;
subs__7 = 1./alpha_PAK;
subs__8 = 1./k_G;
subs__9 = 1./GIT;
subs__10 = 1./subs__2;
subs__11 = 1./Rac_Square;
subs__12 = u(:,2).*alpha_R;
subs__13 = PIX.*Rac_Square.*k_X;
subs__14 = subs__0.*subs__12;
J_gamma_1_2=(Pax_Square.*Rac_Square + 2.*Pax_Square.*u(:,2).*alpha_R + PIX.*Pax_Square.*Rac_Square.*k_X + 2.*PIX.*Pax_Square.*u(:,2).*alpha_R.*k_X + GIT.*PIX.*Pax_Square.*Rac_Square.*k_G.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*Rac_Square.*k_C.*k_G.*k_X + 2.*GIT.*PIX.*u(:,6).*Paxtot.*u(:,2).*alpha_R.*k_C.*k_G.*k_X)./(PAKtot.*Rac_Square.*alpha_PAK.*(Pax_Square + PIX.*Pax_Square.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X));
J_gamma_2_2=-(u(:,6).*Rac_Square.*alpha_R)./(PAKtot.*k_C.*(Rac_Square + u(:,2).*alpha_R).^2);
J_gamma_1_6=-(GIT.^2.*PIX.^2.*Paxtot.*Pax_Square.*u(:,2).*k_C.*k_G.^2.*k_X.^2)./(PAKtot.*alpha_PAK.*(Pax_Square + PIX.*Pax_Square.*k_X + GIT.*PIX.*u(:,6).*Paxtot.*k_C.*k_G.*k_X).^2);
J_gamma_2_6=(Rac_Square + u(:,2).*alpha_R + PIX.*Rac_Square.*k_X + PIX.*u(:,2).*alpha_R.*k_X + GIT.*PIX.*Rac_Square.*k_G.*k_X)./(GIT.*PAKtot.*PIX.*k_C.*k_G.*k_X.*(Rac_Square + u(:,2).*alpha_R)) + (2.*GIT.*PIX.*u(:,6).*Paxtot.*Rac_Square.*k_C.*k_G.*k_X + 2.*GIT.*PIX.*u(:,6).*Paxtot.*u(:,2).*alpha_R.*k_C.*k_G.*k_X)./(GIT.*PAKtot.*PIX.*Pax_Square.*k_C.*k_G.*k_X.*(Rac_Square + u(:,2).*alpha_R));
subs2__0 = 1./(J_gamma_1_2 + J_gamma_2_6 + J_gamma_1_2.*J_gamma_2_6 - J_gamma_1_6.*J_gamma_2_2 + 1);
subs2__1 = J_gamma_2_2.*f_Raci;
subs2__2 = J_gamma_1_6.*f_Paxi;
subs2__3 = -subs2__2;
subs2__4 = -subs2__1;
subs2__5 = J_gamma_1_2.*f_Raci;
subs2__6 = J_gamma_1_2.*f_Paxi;

Rx = [f_Raci,...
-(f_Raci - J_gamma_1_6.*f_Paxi + J_gamma_2_6.*f_Raci)./(J_gamma_1_2 + J_gamma_2_6 + J_gamma_1_2.*J_gamma_2_6 - J_gamma_1_6.*J_gamma_2_2 + 1),...
f_Rhoi,...
-f_Rhoi,...
f_Paxi,...
-(f_Paxi + J_gamma_1_2.*f_Paxi - J_gamma_2_2.*f_Raci)./(J_gamma_1_2 + J_gamma_2_6 + J_gamma_1_2.*J_gamma_2_6 - J_gamma_1_6.*J_gamma_2_2 + 1)];
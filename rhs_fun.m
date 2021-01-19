function Rx = rhs_fun(t,u)
u=transpose(u);
B_1=0.5;
I_rho=0.016;
L_rho=0.34;
delta_rho=0.016;
L_R=0.415;
I_R=0.003;
delta_R=0.025;
alpha_R=15;
Rtot=7.5;
delta_P=0.004;
I_K=0.009;
L_K=5.77;
k_X=41.7;
k_G=5.71;
k_C=5;
GIT=0.11;
PIX=0.069;
Paxtot=2.3;
n=4;
m=4;
alpha_PAK=0.3;
PAKtot= Rtot*alpha_PAK;
Rho_Square= 1.6141e+03;
Rac_Square= 1.6141e+03;
Pax_Square= 494.9785;
RacRatio0 = u(:,2) ./ Rac_Square;
RhoRatio = u(:,4) ./ Rho_Square;
PaxRatio = u(:,6) ./ Pax_Square;
K_is=1./((1+k_X.*PIX+k_G.*k_X.*k_C.*GIT.*PIX.*Paxtot.*PaxRatio).*(1+alpha_R.*RacRatio0)+k_G.*k_X.*GIT.*PIX);
K=alpha_R.*RacRatio0.*K_is.*(1+k_X.*PIX+k_G.*k_X.*k_C.*Paxtot.*GIT.*PIX.*PaxRatio);
I_Ks=I_K.*(1-K_is.*(1+alpha_R.*RacRatio0));
RacRatio = (u(:,2) + alpha_PAK.*K) ./ Rac_Square;
Q_R = (I_R+I_Ks).*(L_rho.^m./(L_rho.^m+RhoRatio.^m));
Q_rho = I_rho.*(L_R.^m./(L_R.^m +(RacRatio + alpha_PAK.*K).^m));
Q_P = B_1.*(K.^n./(L_K.^n+K.^n));
f_Raci = -(Q_R.*u(:,1))+(delta_R.*u(:,2));
f_Rhoi = -(Q_rho.*u(:,3))+(delta_rho.*u(:,4));
f_Paxi = -(Q_P.*u(:,5))+(delta_P.*u(:,6));

Rx = [f_Raci,...
-f_Raci,...
f_Rhoi,...
-f_Rhoi,...
f_Paxi,...
-f_Paxi];
Rx=transpose(Rx);
end
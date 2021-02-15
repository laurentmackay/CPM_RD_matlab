function Rx = rhs_fun(t,u)
u=transpose(u);
B_1=0.500000000;
I_rho=0.016000000;
I_R=0.003000000;
I_K=0.009000000;
L_rho=0.340000000;
L_R=0.340000000;
delta_rho=0.016000000;
delta_R=0.025000000;
delta_P=0.004000000;
alpha_R=15.000000000;
L_K=5.770000000;
k_X=41.700000000;
k_G=5.710000000;
k_C=5.000000000;
GIT=0.110000000;
PIX=0.069000000;
Paxtot=2.300000000;
n=4.000000000;
m=4.000000000;
alpha_PAK=0.300000000;
Rho_Square=1.000000000;
Rac_Square=1.000000000;
Pax_Square=1.000000000;
Rtot=7.5;
PAKtot= 15*0.3;
cnsrv_1=Rac_Square;
cnsrv_2=Rac_Square;
cnsrv_3=Pax_Square;
RacRatio0 = u(:,2) ./ Rac_Square;
PaxRatio = u(:,6) ./ Pax_Square;
K_is=1./((1+k_X.*PIX+k_G.*k_X.*k_C.*GIT.*PIX.*Paxtot.*PaxRatio).*(1+alpha_R.*RacRatio0)+k_G.*k_X.*GIT.*PIX);
RhoRatio = u(:,4) ./ Rho_Square;
K=alpha_R.*RacRatio0.*K_is.*(1+k_X.*PIX+k_G.*k_X.*k_C.*Paxtot.*GIT.*PIX.*PaxRatio);
RacRatio = (u(:,2) + alpha_PAK.*K) ./ Rac_Square;
I_Ks=I_K.*(1-K_is.*(1+alpha_R.*RacRatio0));
Q_R = (I_R+I_Ks).*(L_rho.^m./(L_rho.^m+RhoRatio.^m));
Q_rho = I_rho.*(L_R.^m./(L_R.^m +(RacRatio).^m));
Q_P = B_1.*(K.^n./(L_K.^n+K.^n));
f_Raci = -(Q_R.*u(:,1))+ (delta_R.*u(:,2));
f_Rac =  (Q_R.*u(:,1))-(delta_R.*u(:,2));
f_Rhoi = -(Q_rho.*u(:,3))+ (delta_rho.*u(:,4));
f_Rho =  (Q_rho.*u(:,3))-(delta_rho.*u(:,4));
f_Paxi = -(Q_P.*u(:,5))+ (delta_P.*u(:,6));
f_Pax =  (Q_P.*u(:,5))-(delta_P.*u(:,6));

Rx = [f_Raci,...
-f_Raci,...
f_Rhoi,...
-f_Rhoi,...
f_Paxi,...
-f_Paxi];
Rx=transpose(Rx);
end
if length(vox)>1
[tmp,tmp2]=meshgrid(ir0,vox);
    I_rx=tmp+tmp2;
else
    I_rx=vox+ir0;
end
a_c_0=alpha_chem(I_rx);


RacRatio0(vox) = x(vox+1*sz) / Rac_Square;
RhoRatio(vox) = x(vox+3*sz) / Rho_Square;
PaxRatio(vox) = x(vox+5*sz) / Pax_Square;
K_is(vox)=1.0/((1.0+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(vox))*(1+alpha_R*RacRatio0(vox))+k_G*k_X*GIT*PIX);
K(vox)=alpha_R*RacRatio0(vox)*K_is(vox)*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio(vox));
I_Ks(vox)=I_K*(1.0-K_is(vox)*(1+alpha_R*RacRatio0(vox)));
P_i(vox)=1.0-PaxRatio(vox)*(1+k_G*k_X*k_C*GIT*PIX*PAKtot*K_is(vox)*(1+alpha_R*RacRatio0(vox)));
RacRatio(vox) = (x(vox+1*sz) + alpha_PAK*K(vox)) / Rac_Square;
Q_R(vox) = ((1-RacRatio(vox))/(1-RacRatio0(vox)))*(I_R+I_Ks(vox))*(L_rho^m/(L_rho^m+RhoRatio(vox)^m));
Q_rho(vox) = I_rho*(L_R^m/(L_R^m +(RacRatio(vox))^m));
Q_P(vox) = (P_i(vox)/(1-PaxRatio(vox)))*B*(K(vox)^n/(L_K^n+K(vox)^n));
alpha_chem(vox+0*sz)=(Q_R(vox)).*x(vox+0*sz);
alpha_chem(vox+1*sz)=(delta_R).*x(vox+1*sz);
alpha_chem(vox+2*sz)=(Q_rho(vox)).*x(vox+2*sz);
alpha_chem(vox+3*sz)=(delta_rho).*x(vox+3*sz);
alpha_chem(vox+4*sz)=(Q_P(vox)).*x(vox+4*sz);
alpha_chem(vox+5*sz)=(delta_P).*x(vox+5*sz);

alpha_rx=alpha_rx+sum(alpha_chem(I_rx)-a_c_0,1);
if length(vox)>1
[tmp,tmp2]=meshgrid(ir0,vox);
    I_rx=tmp+tmp2;
else
    I_rx=vox+ir0;
end
a_c_0=alpha_chem(I_rx);


R(vox) = x(vox+1*sz) / Rac_Square;
RhoRatio(vox) = x(vox+3*sz) / Rho_Square;
PaxRatio(vox) = x(vox+5*sz) / Pax_Square;
K_is(vox)=1.0/((1.0+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(vox))*(1+alpha_R*R(vox))+k_G*k_X*GIT*PIX);
K0(vox)=alpha_R*K_is(vox)*(1+k_X*PIX+k_G*k_X*k_C*Paxtot*GIT*PIX*PaxRatio(vox));
K(vox)=R(vox)*K0(vox);
I_Ks(vox)=I_K*(1.0-K_is(vox)*(1+alpha_R*R(vox)));
P_i(vox)=1.0-PaxRatio(vox)*(1+k_G*k_X*k_C*GIT*PIX*PAKtot*K_is(vox)*(1+alpha_R*R(vox)));
Rbar(vox) = (x(vox+1*sz) + x(vox+10*sz)) / Rac_Square;
Q_R(vox) = (I_R+I_Ks(vox)+(x(vox+6*sz)*0.3+0.2)*0.002)*(L_rho^m/(L_rho^m+RhoRatio(vox)^m));
Q_rho(vox) = I_rho*(L_R^m/(L_R^m +(Rbar(vox))^m));
Q_Ps(vox) = B*(K(vox)^n/(L_K^n+K(vox)^n));
Q_Pt(vox) = k1*(x(vox+7*sz)^n/(L_F^n+x(vox+7*sz)^n));
alpha_chem(vox+0*sz)=(Q_R(vox)).*x(vox+0*sz);
alpha_chem(vox+1*sz)=(delta_R).*x(vox+1*sz);
alpha_chem(vox+2*sz)=(Q_rho(vox)).*x(vox+2*sz);
alpha_chem(vox+3*sz)=(delta_rho).*x(vox+3*sz);
alpha_chem(vox+4*sz)=(Q_Ps(vox)).*x(vox+4*sz);
alpha_chem(vox+5*sz)=(delta_P).*x(vox+5*sz);
alpha_chem(vox+6*sz)=(Q_Pt(vox)).*x(vox+4*sz);
alpha_chem(vox+7*sz)=(k2).*x(vox+6*sz);
alpha_chem(vox+8*sz)=(k6).*x(vox+6*sz).*x(vox+7*sz);
alpha_chem(vox+9*sz)=(k7).*x(vox+8*sz);
alpha_chem(vox+10*sz)=(k3).*x(vox+9*sz);
alpha_chem(vox+11*sz)=(k4).*x(vox+7*sz);
alpha_chem(vox+12*sz)=(k8).*x(vox+13*sz);
alpha_chem(vox+13*sz)=(k9).*x(vox+12*sz);

alpha_rx=alpha_rx+sum(alpha_chem(I_rx)-a_c_0,1);
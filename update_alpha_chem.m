if length(vox)>1
[tmp,tmp2]=meshgrid(ir0,vox);
    I_rx=tmp+tmp2;
else
    I_rx=vox+ir0;
end
a_c_0=alpha_chem(I_rx);


RacRatio(vox)=x(vox+1*sz)/Rac_Square;
RhoRatio(vox)=x(vox+3*sz)/Rho_Square;
PaxRatio=0;
alpha_chem(vox+0*sz)=(I_R*(L_rho^m./(L_rho^m+RhoRatio(vox).^m))).*x(vox+0*sz);
alpha_chem(vox+1*sz)=(delta_R).*x(vox+1*sz);
alpha_chem(vox+2*sz)=(I_rho*(L_R^m./(L_R^m +RacRatio(vox).^m))).*x(vox+2*sz);
alpha_chem(vox+3*sz)=(delta_rho).*x(vox+3*sz);

alpha_rx=alpha_rx+sum(alpha_chem(I_rx)-a_c_0,1);
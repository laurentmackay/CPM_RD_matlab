if length(vox)>1
[tmp,tmp2]=meshgrid(ir0,vox);
    I_rx=tmp+tmp2;
else
    I_rx=vox+ir0;
end
a_c_0=alpha_chem(I_rx);


alpha_chem(vox+0*sz)=(a).*x(vox+0*sz).*x(vox+1*sz);
alpha_chem(vox+1*sz)=(d).*x(vox+2*sz);
alpha_chem(vox+2*sz)=(k).*x(vox+2*sz);

alpha_rx=alpha_rx+sum(alpha_chem(I_rx)-a_c_0,1);
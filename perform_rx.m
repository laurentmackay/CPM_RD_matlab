if rx==1
    x(vox+[0  1]*sz)=x(vox+[0  1]*sz)+[-1  1];
elseif rx==2
    x(vox+[0  1]*sz)=x(vox+[0  1]*sz)+[1 -1];
elseif rx==3
    x(vox+[2  3]*sz)=x(vox+[2  3]*sz)+[-1  1];
elseif rx==4
    x(vox+[2  3]*sz)=x(vox+[2  3]*sz)+[1 -1];
elseif rx==5
    x(vox+[4  5]*sz)=x(vox+[4  5]*sz)+[-1  1];
elseif rx==6
    x(vox+[4  5]*sz)=x(vox+[4  5]*sz)+[1 -1];
end

if rx==1||rx==2
    xtot=sum(x(vox+[1  6]*sz));
    x(vox+1*sz)=round(xtot/(1+(1+k_X*PIX+k_G*k_X*k_C*GIT*PIX*Paxtot*PaxRatio(vox))*alpha_PAK*PAKtot*K_is(vox)));
    x(vox+6*sz)=xtot-x(vox+1*sz);
end
if rx==5||rx==6
    xtot=sum(x(vox+[5  7]*sz));
    x(vox+5*sz)=round(xtot/(1+(k_G*k_X*k_C*GIT*PIX*K_is(vox)*PAKtot*(1+alpha_R*RacRatio0(vox)))));
    x(vox+7*sz)=xtot-x(vox+5*sz);
end
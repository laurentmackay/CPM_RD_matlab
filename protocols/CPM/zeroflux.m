function [pl,ql,pr,qr] = zeroflux(xl,ul,xr,ur,t)
D=[ 0.4300    0.0200    0.4300    0.0200    0.0200  0.0200];

ql=ones(size(ul)).*D';
qr=ones(size(ur)).*D';

pl=zeros(size(ul));
pr=zeros(size(ur));
end


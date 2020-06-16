function m = Alg5(M,pT,xi)
l=0;u=M;
fl=betainc(pT,l,M-l+1)-xi;
fu=betainc(pT,u,M-u+1)-xi;


while abs(l-u)>1
r=round((l+u)/2);
fr=betainc(pT,r,M-r+1)-xi;
if sign(fl)==sign(fr)
    l=r;
    fl=fr;
else
    u=r;
    fu=fr;
end

end
m=u;
end


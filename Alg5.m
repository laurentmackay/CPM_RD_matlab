function m = Alg5(M,pT,xi)
f=@(i) betainc(pT,i,M-i+1)-xi;
l=0;u=M;
fl=f(l);
fu=f(u);


while abs(l-u)>1
r=round((l+u)/2);
fr=f(r);
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


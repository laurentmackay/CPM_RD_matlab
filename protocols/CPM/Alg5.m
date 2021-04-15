%determine number of molecules to be transferred
function m = Alg5(M,pT)
l=zeros(size(M),class(M));u=M;
inds=true(size(M));

xi=rand(size(M));

if any(pT>1) || any(pT<0)
    error('wtf mate')
end

fl=betainc(pT,l,M-l+1)-xi;
fu=betainc(pT,u,M-u+1)-xi;

r=round((l+u)/2);
fr=betainc(pT,r,M-r+1)-xi;
samesign=sign(fl)==sign(fr);


while any(inds)
r(inds)=round((l(inds)+u(inds))/2);
fr(inds)=betainc(pT(inds),r(inds),M(inds)-r(inds)+1)-xi(inds);


samesign(inds)= sign(fl(inds))==sign(fr(inds));

inds2=inds&samesign;
l(inds2)=r(inds2);
fl(inds2)=fr(inds2);

inds2=inds&~samesign;
u(inds2)=r(inds2);
fu(inds2)=fr(inds2);

inds=abs(l-u)>1;
end
m=u;
end
function i = Alg23(p,N)
inds=true(1,N);
i=zeros(size(inds));
b=zeros(size(inds));
xi=rand(size(inds));
while any(inds)
    i(inds)=i(inds)+1;
    b(inds)=b(inds)+p(i(inds));
    inds(~inds)=xi(inds)>=b(inds);
end




function mi = sample_mi(m,pi)
mi=zeros(size(pi));
for i=1:m
    j=Alg2(pi);
    mi(j)=mi(j)+1;
end
    
end


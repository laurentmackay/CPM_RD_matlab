function mi = sample_mi(m,pi)
n=size(m,1);

mi=histc(rand(n,m),[zeros(size(pi,1),1) cumsum(pi)],2);
mi(:,end)=[];%delete last bin
end


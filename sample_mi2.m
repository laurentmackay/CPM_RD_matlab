function mi = sample_mi(m,pi)
%SAMPLE_MI: multinomial sampling from bins with probability pi, m samples
%total
%  m = (1xn) vector of integers
%  pi = (kxn) vector of probabilites for each "bin"
%
%  mi = (kxn) vector of integers sampled from pi. sum(mi)=m.

n=size(m,1);
k=size(pi,1);
mi=zeros(size(pi));
for j=1:n
temp=histc(rand(1,m(j)),[0 cumsum(pi(:,j))'],2);
mi(:,j)=temp(1:k);% delete last (empty) bin
end


end


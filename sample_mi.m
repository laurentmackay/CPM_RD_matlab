function mi = sample_mi(m,pi)
%SAMPLE_MI: multinomial sampling from bins with probability pi, m samples
%total
%  m = scalar integer
%  pi = (1xk) vector of probabilites for each "bin"
%
%  mi = (1xk) vector of integers sampled from pi. sum(mi)=m.

mi=histc(rand(1,m),[0 cumsum(pi)],2);
mi=mi(1:end-1);% delete last (empty) bin
end


function mi = sample_mi(m,pi)
%SAMPLE_MI: multinomial sampling from bins with probability pi, m samples
%total
%  m = scalar integer
%  pi = (1xk) vector of probabilites for each "bin"
%
%  mi = (1xk) vector of integers sampled from pi. sum(mi)=m.
f=0.5;
i=1;
mi=zeros(1,length(pi)+1);
m_max=5e2;
mp=sum(mi);
edges=[0 cumsum(pi)];
while mp<m
    m_eff=min(m_max,(m-mp));
    dm=ceil(histc(rand(1,m_eff),edges,2)*(1-f^i)*(m-mp)/m_eff); % hack for speed
    mi=mi+dm; 
    mp=sum(mi);
    i=i+1;
end
mi(:,end)=[];%delete last (empty) bin
end


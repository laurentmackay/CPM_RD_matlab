function t = getHalflife(x,f)
if nargin==1
    f=0.5;
end
N=size(x,2);
t=zeros(1,N);
for i=1:N
    ind=find(x(:,i)<=f,1);
if ~isempty(ind)

    if ind>1
        t(i)=ind-(f-x(ind,i))/(x(ind-1,i)-x(ind,i))-1;
    else
        t(i)=0;
    end

else
    t(i)=inf;
end
end
end


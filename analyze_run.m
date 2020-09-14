%thats (how they do it experimentally)
dx=zeros(size(center,1));
dx(1)=sqrt(sum((center(120/picstep,:)-center(1,:)).^2));
for i=2:length(dx-1)
    dx(i)=sqrt(sum((center(i*120/picstep,:)-center((i-1)*120/picstep,:)).^2));
end
v=dx/120*3600;
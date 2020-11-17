function v = get_instant_velocity(f,iters,n0)

%%load in the results
load(['results/' f]);
iter=min(iter,length(Times));

if nargin==1
    iters=1:iter;
end

center=center(:,iters);
Times=Times(iters);

if nargin<3
    n=1;
end
dr=(center(:,1:end-n0)-center(:,1+n0:end));

v= dr./(Times(1+n0:end)-Times(1:end-n0));

net_vel=sqrt(sum((center(:,1)-center(:,end)).^2,1))/(Times(end)-Times(1))

end


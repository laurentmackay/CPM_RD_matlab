function [ux,uy] = centralGrad(u)
ux=zeros(size(u));
uy=zeros(size(u));
ux(2:end-1,:)=u(3:end,:)-u(1:end-2,:);
ux(:,2:end-1)=u(:,3:end)-u(:,1:end-2);
end


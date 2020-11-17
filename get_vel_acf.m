function [acf, lags] = get_vel_acf(v)
M=max(size(v));
lags=0:M-1;
acf=zeros(size(lags));
for l=lags

acf(l+1)=sum(sum(v(:,l+1:end).*v(:,1:end-l),1))/sum(sum(v(:,l+1:end).^2));

end




end


function [acf, lags] = get_vel_acf(v)
M=max(size(v));
lags=0:M-1;
acf=zeros(size(lags));
norm_factor=sum(sum(v.^2));
for l=lags

acf(l+1)=sum(sum(v(:,l+1:end).*v(:,1:end-l),1))/norm_factor

end
% vcom=(v(1,:)+i*v(2,:))';
% 
% [acf,lags]=xcorr(vcom,conj(vcom),'coeff');


end


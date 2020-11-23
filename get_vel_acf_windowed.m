function [acf, lags] = get_vel_acf_windowed(v,win_len,step)
if nargin==2
    step=1;
end
% type='for'

M=size(v,2);
lags=0:win_len-1;
% if strcmp(type,'for')
max_lag=lags(end)-1;
steps=1:step:M-max_lag;

acf=zeros(size(lags,2),size(steps,2));
j=0;
for i=steps
    norm_factor=sum(sum(v(:,i:i+max_lag).^2));
    j=j+1;
    for l=lags
        acf(l+1,j)=sum(sum(v(:,i:i+max_lag-l).*v(:,l+i:i+max_lag),1))/norm_factor;
    end
end
end


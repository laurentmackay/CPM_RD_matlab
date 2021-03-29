function [acf, lags, steps] = get_vel_acf_windowed(v,win_len,step)
if nargin==2
    step=1;
end
% type='for'

M=size(v,2);
lags=0:win_len-1;
% if strcmp(type,'for')
max_lag=lags(end);
steps=1:step:M;



acf=zeros(size(lags,2),size(steps,2));
j=0;
for i=steps
    %     norm_factor=sum(sum(v(:,i:i+max_lag).^2));
    j=j+1;
    for l=lags
        
        if M>i+l
            
            if i>ceil(max_lag/2)

                    sp=min(i,M-max_lag-l);
            else
                sp=i;
            end
            %         ep=min(sp,M-l);
            inds=sp:sp+max_lag;
            
%             if i>M-max_lag
%                 l
%             end
            
            %             if l==0;
            try
                s1=sqrt(sum(sum(v(:,inds).*v(:,inds),1)));
                s2=sqrt(sum(sum(v(:,inds+l).*v(:,inds+l),1)));
                norm_factor = s1*s2;
                %                end
                acf(l+1,j)=sum(sum(v(:,inds).*v(:,inds+l),1))/norm_factor;
            catch err
                rethrow(err)
            end
        end
        
    end
end
end


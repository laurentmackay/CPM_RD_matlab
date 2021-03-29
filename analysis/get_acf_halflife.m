function [lambda, t] = get_acf_halflife(f, n, win_len)

if nargin<2
    n=3;
end

if nargin<3
    win_len=10;
end


[vtot, t] = get_instant_velocity(f,[],n);
[acf_tot,~,t_acf] = get_vel_acf_windowed(vtot,win_len,1); 
lambda = getHalflife(acf_tot(:,1:end));
t=t(1:length(lambda));
end


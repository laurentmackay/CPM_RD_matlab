function varargout = get_polarization_time(f,thresh)
if nargin<2
    thresh=0.05;
end

[t,y]=timecourse(f,"[min(min(x(inds))), max(max(x(inds)))]",...
        "i_rac = find(strcmp(chems,'Rac')); inds=cell_inds(1:A)+sz*(i_rac-1); ");
delta=abs( y(:,1)-y(:,2));

figure(3);
semilogy(t,delta)

figure(4);
plot(t,y(:,1),t, y(:,2))

i_polarize=find(delta>thresh,1);

if isempty(i_polarize)
    t_polarize=Inf;
else
    t_polarize=t(i_polarize);
end

varargout{1}=t_polarize;

if nargout>1
    varargout{2}=i_polarize;
end
        
end


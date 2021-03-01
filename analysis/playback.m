function playback(f,plot_fun)
addpath(genpath('..'));
if nargin<2 || isempty(plot_uf)

%     addpath(genpath('..'));

    plot_fun = @pic;

%     rmpath(genpath('..'))

end

load(strcat('../results/',f));

plotting=true;
N_steps=size(Results,4);
initialize_pic()

for i=1:iter
   x=Results(:,:,2:end,i);
   time=Times(i);
   eval_model
   plot_fun()
   drawnow()
end

rmpath(genpath('..'))


end
function playback(f,plot_fun)
addpath(genpath('..'));
if nargin<2 || isempty(plot_fun)


    plot_fun = @pic;
    initialize_pic()
    tic();



end

load(strcat('../results/',f));

plotting=true;
N_steps=size(Results,4);

for i=1:iter
   x=Results(:,:,2:end,i);
   time=Times(i);
   eval_model
   plot_fun()
   drawnow()
end

rmpath(genpath('..'))


end
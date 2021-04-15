function playback(f, frame_skip, plot_fun, code)
% addpath(genpath('..'));

load(strcat(results_dir(),f));

if nargin<2
    frame_skip=1;
end

if nargin<3 || isempty(plot_fun)
    
    plot_fun = @pic;
    initialize_pic()
    tic();

end

plotting=true;
N_steps=size(Results,4);

if nargin>3
    if iscell(code) && length(code)>1
        x=Results(:,:,2:end,1);
        cell_mask = Results(:,:,1,1);
        eval(code{1})
        time=Times(1);
        eval_model
        plot_fun()
        drawnow()
        eval(code{1})
        code=code{2}; 
        eval(code);
    end
    
    iplot_0=1+frame_skip;
else
    iplot_0=1;
    
end





for i=iplot_0:frame_skip:iter
   x=Results(:,:,2:end,i);
   cell_mask = Results(:,:,1,i);
   
   time=Times(i);
   eval_model
   plot_fun()
   drawnow()
   eval(code)
   
end

% rmpath(genpath('..'))


end
function playback(f, frame_skip,plot_fun)
% addpath(genpath('..'));
[~, i_polarize] = get_polarization_time(f);
if ~isempty(i_polarize)
    
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
    
    for i=i_polarize:frame_skip:iter
        x=Results(:,:,2:end,i);
        cell_mask = Results(:,:,1,i);
        
        time=Times(i);
        eval_model
        plot_fun()
        drawnow()
    end
    
    % rmpath(genpath('..'))
end

end
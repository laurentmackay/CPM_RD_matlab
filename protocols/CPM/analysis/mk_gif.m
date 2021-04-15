function [outputArg1,outputArg2] = mk_gif(f,frame_skip,plot_fun,frame)
if nargin<2
    frame_skip=1;
end
if nargin<3
    plot_fun=[];
end

if nargin<4
    frame='panel1';
end

nm = regexprep(f,'\.mat$','\.gif');

playback(f,frame_skip, plot_fun, {strcat("delete ", nm ,"; gif('",nm,"','frame',",frame,")"),"gif"})
end


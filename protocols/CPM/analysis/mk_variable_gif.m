function [outputArg1,outputArg2] = mk_variable_gif(file, var,frame_skip)

if nargin<3
    frame_skip=1;
end

load(strcat(results_dir(),file));

i_var = find(strcmp(var, chems));
var_max = max(max(max(Results(:,:,i_var+1,1:iter))));
cell_mask = Results(:,:,1,1);
% var_min = min(min(min(Results(:,:,i_var+1,1:iter))));
figure(1);clf()

subplot(1,1,1);
plot_var
caxis([0 var_max]);
colorbar();

nm = regexprep(file,'\.mat$','\.gif');

playback(file,frame_skip, @plot_var, {strcat("delete ", nm ,"; gif('",nm,"','frame',gcf);var='",var,"';i_var=",int2str(i_var),';'),"gif"})

end


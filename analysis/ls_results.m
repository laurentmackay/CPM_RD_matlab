function varargout=ls_results(base_dir)
global active_model

if nargin<1
    if ~isempty(active_model) 
        base_dir = strcat('../_',active_model);
    else
        base_dir = '..';
    end
end
search_str = strcat(base_dir,'/results/*.mat');
if nargout==1
    varargout{1}=dir(search_str);
else
    dir(search_str)
end
end
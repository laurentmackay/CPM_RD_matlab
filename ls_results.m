function varargout=ls_results(ext)
global active_model

if nargin==0
    ext='mat';
end

if ~isempty(ext)
    search_str = strcat(results_dir(),'*.',ext);
else
    search_str = strcat(results_dir(),'*');
end


if nargout==1
    varargout{1}=dir(search_str);
else
    dir(search_str)
end
end
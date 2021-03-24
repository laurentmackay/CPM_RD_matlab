function dir = work_dir()
global active_model RD_base

if isempty(RD_base)
    p=mfilename('fullpath');
    RD_base = regexprep(p,[ filesep '[^\' filesep ']+$'],filesep);
end

if isempty(active_model)
   error('There is currently no active model. Please activate one using `deploy_model`.') 
end


dir =  strcat(RD_base, '_', active_model, filesep);

end



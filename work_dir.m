function name = work_dir()
global active_model RD_base

if isempty(RD_base)
    p=mfilename('fullpath');
    RD_base = regexprep(p,[ filesep '[^\' filesep ']+$'],filesep);
end

if isempty(active_model)
   error('There is currently no active model. Please activate one using `deploy_model`.') 
end


name =  strcat(RD_base, '_', active_model, filesep);

if isempty(dir(name))
    mkdir(name)
end

end



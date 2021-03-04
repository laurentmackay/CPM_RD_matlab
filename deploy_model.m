function deploy_model(f,force)
global active_model

if nargin<2
    force=false;
end

mk_work_dir = @(f) strcat('./_',f);
[~,fn,~] = fileparts(f);
work_dir = mk_work_dir(fn);

if isempty(active_model)
    active_model=fn;
    addpath(genpath(work_dir));
elseif force || ~strcmp(active_model,fn)
    rmpath(genpath(mk_work_dir(active_model)));
    addpath(genpath(work_dir));
    active_model=fn;
end



if isempty(dir(work_dir))
    mkdir(work_dir)
end

hash_file = strcat(work_dir,filesep,'model.hash');
hash = string2hash(fileread(f));
%check the model-description hash to see if it has changed
if force || isempty(dir(hash_file)) || strcmp(fileread(hash_file), hash)
    
    try
    mk_rxn_files(f, work_dir);
    catch err
        rmpath(genpath(work_dir));
        active_model=[];
        rethrow(err)
    end
    dlmwrite( hash_file, hash, '');
end






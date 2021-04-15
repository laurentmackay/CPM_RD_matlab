function deploy_model(f,force)
global active_model RD_base protocol


if isempty(RD_base)
    get_RD_base();
end

if isempty(protocol) && length(dir(strcat(RD_base,'protocols')))==3 %if there is only one protocol, automatically use that one
    set_protocol(last(dir(strcat(RD_base,'protocols'))).name)
end


if nargin<2
    force=false;
end


[~,fn,~] = fileparts(f);


if isempty(active_model)
    active_model=fn;
    addpath(work_dir());
elseif force || ~strcmp(active_model,fn)
    rmpath(work_dir());
    active_model=fn;
    addpath(work_dir());
    
end





hash_file = strcat(work_dir,filesep,'model.hash');
hash = string2hash(fileread(strcat(RD_base,'models',filesep,f)));
%check the model-description hash to see if it has changed
if force || isempty(dir(hash_file)) || strcmp(fileread(hash_file), hash)
    
    try
    mk_rxn_files(f, work_dir());
    catch err
        rmpath(work_dir());
        active_model=[];
        rethrow(err)
    end
    dlmwrite( hash_file, hash, '');
end






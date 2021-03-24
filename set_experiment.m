function set_experiment(nm)
global experiment
experiment=nm;
if exist(results_dir())~=7
    mkdir(results_dir());
end
end

